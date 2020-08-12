#include <petscsys.h>        /*I  "petscsys.h"  I*/
#include <petsc/private/petscimpl.h>

struct _n_PetscShmComm {
  PetscMPIInt *globranks;       /* global ranks of each rank in the shared memory communicator */
  PetscMPIInt shmsize;          /* size of the shared memory communicator */
  MPI_Comm    globcomm,shmcomm; /* global communicator and shared memory communicator (a sub-communicator of the former) */
};

/*
   Private routine to delete internal shared memory communicator when a communicator is freed.

   This is called by MPI, not by users. This is called by MPI_Comm_free() when the communicator that has this  data as an attribute is freed.

   Note: this is declared extern "C" because it is passed to MPI_Comm_create_keyval()

*/
PETSC_EXTERN PetscMPIInt MPIAPI Petsc_ShmComm_Attr_Delete_Fn(MPI_Comm comm,PetscMPIInt keyval,void *val,void *extra_state)
{
  PetscErrorCode  ierr;
  PetscShmComm p = (PetscShmComm)val;

  PetscFunctionBegin;
  ierr = PetscInfo1(NULL,"Deleting shared memory subcommunicator in a MPI_Comm %ld\n",(long)comm);CHKERRMPI(ierr);
  ierr = MPI_Comm_free(&p->shmcomm);CHKERRMPI(ierr);
  ierr = PetscFree(p->globranks);CHKERRMPI(ierr);
  ierr = PetscFree(val);CHKERRMPI(ierr);
  PetscFunctionReturn(MPI_SUCCESS);
}

/*@C
    PetscShmCommGet - Given a PETSc communicator returns a communicator of all ranks that share a common memory


    Collective.

    Input Parameter:
.   globcomm - MPI_Comm

    Output Parameter:
.   pshmcomm - the PETSc shared memory communicator object

    Level: developer

    Notes:
    This should be called only with an PetscCommDuplicate() communictor

           When used with MPICH, MPICH must be configured with --download-mpich-device=ch3:nemesis

@*/
PetscErrorCode PetscShmCommGet(MPI_Comm globcomm,PetscShmComm *pshmcomm)
{
#ifdef PETSC_HAVE_MPI_PROCESS_SHARED_MEMORY
  PetscErrorCode   ierr;
  MPI_Group        globgroup,shmgroup;
  PetscMPIInt      *shmranks,i,flg;
  PetscCommCounter *counter;

  PetscFunctionBegin;
  ierr = MPI_Comm_get_attr(globcomm,Petsc_Counter_keyval,&counter,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(globcomm,PETSC_ERR_ARG_CORRUPT,"Bad MPI communicator supplied; must be a PETSc communicator");

  ierr = MPI_Comm_get_attr(globcomm,Petsc_ShmComm_keyval,pshmcomm,&flg);CHKERRQ(ierr);
  if (flg) PetscFunctionReturn(0);

  ierr        = PetscNew(pshmcomm);CHKERRQ(ierr);
  (*pshmcomm)->globcomm = globcomm;

  ierr = MPI_Comm_split_type(globcomm, MPI_COMM_TYPE_SHARED,0, MPI_INFO_NULL,&(*pshmcomm)->shmcomm);CHKERRQ(ierr);

  ierr = MPI_Comm_size((*pshmcomm)->shmcomm,&(*pshmcomm)->shmsize);CHKERRQ(ierr);
  ierr = MPI_Comm_group(globcomm, &globgroup);CHKERRQ(ierr);
  ierr = MPI_Comm_group((*pshmcomm)->shmcomm, &shmgroup);CHKERRQ(ierr);
  ierr = PetscMalloc1((*pshmcomm)->shmsize,&shmranks);CHKERRQ(ierr);
  ierr = PetscMalloc1((*pshmcomm)->shmsize,&(*pshmcomm)->globranks);CHKERRQ(ierr);
  for (i=0; i<(*pshmcomm)->shmsize; i++) shmranks[i] = i;
  ierr = MPI_Group_translate_ranks(shmgroup, (*pshmcomm)->shmsize, shmranks, globgroup, (*pshmcomm)->globranks);CHKERRQ(ierr);
  ierr = PetscFree(shmranks);CHKERRQ(ierr);
  ierr = MPI_Group_free(&globgroup);CHKERRQ(ierr);
  ierr = MPI_Group_free(&shmgroup);CHKERRQ(ierr);

  for (i=0; i<(*pshmcomm)->shmsize; i++) {
    ierr = PetscInfo2(NULL,"Shared memory rank %d global rank %d\n",i,(*pshmcomm)->globranks[i]);CHKERRQ(ierr);
  }
  ierr = MPI_Comm_set_attr(globcomm,Petsc_ShmComm_keyval,*pshmcomm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
#else
  SETERRQ(globcomm, PETSC_ERR_SUP, "Shared memory communicators need MPI-3 package support.\nPlease upgrade your MPI or reconfigure with --download-mpich.");
#endif
}

/*@C
    PetscShmCommGlobalToLocal - Given a global rank returns the local rank in the shared memory communicator

    Input Parameters:
+   pshmcomm - the shared memory communicator object
-   grank    - the global rank

    Output Parameter:
.   lrank - the local rank, or MPI_PROC_NULL if it does not exist

    Level: developer

    Developer Notes:
    Assumes the pshmcomm->globranks[] is sorted

    It may be better to rewrite this to map multiple global ranks to local in the same function call

@*/
PetscErrorCode PetscShmCommGlobalToLocal(PetscShmComm pshmcomm,PetscMPIInt grank,PetscMPIInt *lrank)
{
  PetscMPIInt    low,high,t,i;
  PetscBool      flg = PETSC_FALSE;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  *lrank = MPI_PROC_NULL;
  if (grank < pshmcomm->globranks[0]) PetscFunctionReturn(0);
  if (grank > pshmcomm->globranks[pshmcomm->shmsize-1]) PetscFunctionReturn(0);
  ierr = PetscOptionsGetBool(NULL,NULL,"-noshared",&flg,NULL);CHKERRQ(ierr);
  if (flg) PetscFunctionReturn(0);
  low  = 0;
  high = pshmcomm->shmsize;
  while (high-low > 5) {
    t = (low+high)/2;
    if (pshmcomm->globranks[t] > grank) high = t;
    else low = t;
  }
  for (i=low; i<high; i++) {
    if (pshmcomm->globranks[i] > grank) PetscFunctionReturn(0);
    if (pshmcomm->globranks[i] == grank) {
      *lrank = i;
      PetscFunctionReturn(0);
    }
  }
  PetscFunctionReturn(0);
}

/*@C
    PetscShmCommLocalToGlobal - Given a local rank in the shared memory communicator returns the global rank

    Input Parameters:
+   pshmcomm - the shared memory communicator object
-   lrank    - the local rank in the shared memory communicator

    Output Parameter:
.   grank - the global rank in the global communicator where the shared memory communicator is built

    Level: developer

@*/
PetscErrorCode PetscShmCommLocalToGlobal(PetscShmComm pshmcomm,PetscMPIInt lrank,PetscMPIInt *grank)
{
  PetscFunctionBegin;
#ifdef PETSC_USE_DEBUG
  {
    PetscErrorCode ierr;
    if (lrank < 0 || lrank >= pshmcomm->shmsize) { SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"No rank %D in the shared memory communicator",lrank);CHKERRQ(ierr); }
  }
#endif
  *grank = pshmcomm->globranks[lrank];
  PetscFunctionReturn(0);
}

/*@C
    PetscShmCommGetMpiShmComm - Returns the MPI communicator that represents all processes with common shared memory

    Input Parameter:
.   pshmcomm - PetscShmComm object obtained with PetscShmCommGet()

    Output Parameter:
.   comm     - the MPI communicator

    Level: developer

@*/
PetscErrorCode PetscShmCommGetMpiShmComm(PetscShmComm pshmcomm,MPI_Comm *comm)
{
  PetscFunctionBegin;
  *comm = pshmcomm->shmcomm;
  PetscFunctionReturn(0);
}

#if defined(PETSC_HAVE_OPENMP_SUPPORT)
#include <pthread.h>
#include <hwloc.h>
#include <omp.h>

/* Use mmap() to allocate shared mmeory (for the pthread_barrier_t object) if it is available,
   otherwise use MPI_Win_allocate_shared. They should have the same effect except MPI-3 is much
   simpler to use. However, on a Cori Haswell node with Cray MPI, MPI-3 worsened a test's performance
   by 50%. Until the reason is found out, we use mmap() instead.
*/
#define USE_MMAP_ALLOCATE_SHARED_MEMORY

#if defined(USE_MMAP_ALLOCATE_SHARED_MEMORY) && defined(PETSC_HAVE_MMAP)
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#endif

struct _n_PetscOmpCtrl {
  MPI_Comm          omp_comm;        /* a shared memory communicator to spawn omp threads */
  MPI_Comm          omp_master_comm; /* a communicator to give to third party libraries */
  PetscMPIInt       omp_comm_size;   /* size of omp_comm, a kind of OMP_NUM_THREADS */
  PetscBool         is_omp_master;   /* rank 0's in omp_comm */
  MPI_Win           omp_win;         /* a shared memory window containing a barrier */
  pthread_barrier_t *barrier;        /* pointer to the barrier */
  hwloc_topology_t  topology;
  hwloc_cpuset_t    cpuset;          /* cpu bindings of omp master */
  hwloc_cpuset_t    omp_cpuset;      /* union of cpu bindings of ranks in omp_comm */
};


/* Allocate and initialize a pthread_barrier_t object in memory shared by processes in omp_comm
   contained by the controller.

   PETSc OpenMP controller users do not call this function directly. This function exists
   only because we want to separate shared memory allocation methods from other code.
 */
PETSC_STATIC_INLINE PetscErrorCode PetscOmpCtrlCreateBarrier(PetscOmpCtrl ctrl)
{
  PetscErrorCode        ierr;
  MPI_Aint              size;
  void                  *baseptr;
  pthread_barrierattr_t  attr;

#if defined(USE_MMAP_ALLOCATE_SHARED_MEMORY) && defined(PETSC_HAVE_MMAP)
  PetscInt              fd;
  PetscChar             pathname[PETSC_MAX_PATH_LEN];
#else
  PetscMPIInt           disp_unit;
#endif

  PetscFunctionBegin;
#if defined(USE_MMAP_ALLOCATE_SHARED_MEMORY) && defined(PETSC_HAVE_MMAP)
  size = sizeof(pthread_barrier_t);
  if (ctrl->is_omp_master) {
    /* use PETSC_COMM_SELF in PetscGetTmp, since it is a collective call. Using omp_comm would otherwise bcast the partially populated pathname to slaves */
    ierr    = PetscGetTmp(PETSC_COMM_SELF,pathname,PETSC_MAX_PATH_LEN);CHKERRQ(ierr);
    ierr    = PetscStrlcat(pathname,"/petsc-shm-XXXXXX",PETSC_MAX_PATH_LEN);CHKERRQ(ierr);
    /* mkstemp replaces XXXXXX with a unique file name and opens the file for us */
    fd      = mkstemp(pathname); if(fd == -1) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_LIB,"Could not create tmp file %s with mkstemp\n", pathname);
    ierr    = ftruncate(fd,size);CHKERRQ(ierr);
    baseptr = mmap(NULL,size,PROT_READ | PROT_WRITE, MAP_SHARED,fd,0); if (baseptr == MAP_FAILED) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_LIB,"mmap() failed\n");
    ierr    = close(fd);CHKERRQ(ierr);
    ierr    = MPI_Bcast(pathname,PETSC_MAX_PATH_LEN,MPI_CHAR,0,ctrl->omp_comm);CHKERRQ(ierr);
    /* this MPI_Barrier is to wait slaves to open the file before master unlinks it */
    ierr    = MPI_Barrier(ctrl->omp_comm);CHKERRQ(ierr);
    ierr    = unlink(pathname);CHKERRQ(ierr);
  } else {
    ierr    = MPI_Bcast(pathname,PETSC_MAX_PATH_LEN,MPI_CHAR,0,ctrl->omp_comm);CHKERRQ(ierr);
    fd      = open(pathname,O_RDWR); if(fd == -1) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_LIB,"Could not open tmp file %s\n", pathname);
    baseptr = mmap(NULL,size,PROT_READ | PROT_WRITE, MAP_SHARED,fd,0); if (baseptr == MAP_FAILED) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_LIB,"mmap() failed\n");
    ierr    = close(fd);CHKERRQ(ierr);
    ierr    = MPI_Barrier(ctrl->omp_comm);CHKERRQ(ierr);
  }
#else
  size = ctrl->is_omp_master ? sizeof(pthread_barrier_t) : 0;
  ierr = MPI_Win_allocate_shared(size,1,MPI_INFO_NULL,ctrl->omp_comm,&baseptr,&ctrl->omp_win);CHKERRQ(ierr);
  ierr = MPI_Win_shared_query(ctrl->omp_win,0,&size,&disp_unit,&baseptr);CHKERRQ(ierr);
#endif
  ctrl->barrier = (pthread_barrier_t*)baseptr;

  /* omp master initializes the barrier */
  if (ctrl->is_omp_master) {
    ierr = MPI_Comm_size(ctrl->omp_comm,&ctrl->omp_comm_size);CHKERRQ(ierr);
    ierr = pthread_barrierattr_init(&attr);CHKERRQ(ierr);
    ierr = pthread_barrierattr_setpshared(&attr,PTHREAD_PROCESS_SHARED);CHKERRQ(ierr); /* make the barrier also work for processes */
    ierr = pthread_barrier_init(ctrl->barrier,&attr,(unsigned int)ctrl->omp_comm_size);CHKERRQ(ierr);
    ierr = pthread_barrierattr_destroy(&attr);CHKERRQ(ierr);
  }

  /* this MPI_Barrier is to make sure the omp barrier is initialized before slaves use it */
  ierr = MPI_Barrier(ctrl->omp_comm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* Destroy the pthread barrier in the PETSc OpenMP controller */
PETSC_STATIC_INLINE PetscErrorCode PetscOmpCtrlDestroyBarrier(PetscOmpCtrl ctrl)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  /* this MPI_Barrier is to make sure slaves have finished using the omp barrier before master destroys it */
  ierr = MPI_Barrier(ctrl->omp_comm);CHKERRQ(ierr);
  if (ctrl->is_omp_master) { ierr = pthread_barrier_destroy(ctrl->barrier);CHKERRQ(ierr); }

#if defined(USE_MMAP_ALLOCATE_SHARED_MEMORY) && defined(PETSC_HAVE_MMAP)
  ierr = munmap(ctrl->barrier,sizeof(pthread_barrier_t));CHKERRQ(ierr);
#else
  ierr = MPI_Win_free(&ctrl->omp_win);CHKERRQ(ierr);
#endif
  PetscFunctionReturn(0);
}

/*@C
    PetscOmpCtrlCreate - create a PETSc OpenMP controller, which manages PETSc's interaction with third party libraries using OpenMP

    Input Parameter:
+   petsc_comm - a communicator some PETSc object (for example, a matrix) lives in
-   nthreads   - number of threads per MPI rank to spawn in a library using OpenMP. If nthreads = -1, let PETSc decide a suitable value

    Output Parameter:
.   pctrl      - a PETSc OpenMP controller

    Level: developer

    TODO: Possibly use the variable PetscNumOMPThreads to determine the number for threads to use 

.seealso PetscOmpCtrlDestroy()
@*/
PetscErrorCode PetscOmpCtrlCreate(MPI_Comm petsc_comm,PetscInt nthreads,PetscOmpCtrl *pctrl)
{
  PetscErrorCode        ierr;
  PetscOmpCtrl          ctrl;
  unsigned long         *cpu_ulongs=NULL;
  PetscInt              i,nr_cpu_ulongs;
  PetscShmComm          pshmcomm;
  MPI_Comm              shm_comm;
  PetscMPIInt           shm_rank,shm_comm_size,omp_rank,color;
  PetscInt              num_packages,num_cores;

  PetscFunctionBegin;
  ierr = PetscNew(&ctrl);CHKERRQ(ierr);

  /*=================================================================================
    Init hwloc
   ==================================================================================*/
  ierr = hwloc_topology_init(&ctrl->topology);CHKERRQ(ierr);
#if HWLOC_API_VERSION >= 0x00020000
  /* to filter out unneeded info and have faster hwloc_topology_load */
  ierr = hwloc_topology_set_all_types_filter(ctrl->topology,HWLOC_TYPE_FILTER_KEEP_NONE);CHKERRQ(ierr);
  ierr = hwloc_topology_set_type_filter(ctrl->topology,HWLOC_OBJ_CORE,HWLOC_TYPE_FILTER_KEEP_ALL);CHKERRQ(ierr);
#endif
  ierr = hwloc_topology_load(ctrl->topology);CHKERRQ(ierr);

  /*=================================================================================
    Split petsc_comm into multiple omp_comms. Ranks in an omp_comm have access to
    physically shared memory. Rank 0 of each omp_comm is called an OMP master, and
    others are called slaves. OMP Masters make up a new comm called omp_master_comm,
    which is usually passed to third party libraries.
   ==================================================================================*/

  /* fetch the stored shared memory communicator */
  ierr = PetscShmCommGet(petsc_comm,&pshmcomm);CHKERRQ(ierr);
  ierr = PetscShmCommGetMpiShmComm(pshmcomm,&shm_comm);CHKERRQ(ierr);

  ierr = MPI_Comm_rank(shm_comm,&shm_rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(shm_comm,&shm_comm_size);CHKERRQ(ierr);

  /* PETSc decides nthreads, which is the smaller of shm_comm_size or cores per package(socket) */
  if (nthreads == -1) {
    num_packages = hwloc_get_nbobjs_by_type(ctrl->topology,HWLOC_OBJ_PACKAGE) <= 0 ? 1 : hwloc_get_nbobjs_by_type(ctrl->topology,HWLOC_OBJ_PACKAGE);
    num_cores    = hwloc_get_nbobjs_by_type(ctrl->topology,HWLOC_OBJ_CORE) <= 0 ? 1 :  hwloc_get_nbobjs_by_type(ctrl->topology,HWLOC_OBJ_CORE);
    nthreads     = num_cores/num_packages;
    if (nthreads > shm_comm_size) nthreads = shm_comm_size;
  }

  if (nthreads < 1 || nthreads > shm_comm_size) SETERRQ2(petsc_comm,PETSC_ERR_ARG_OUTOFRANGE,"number of OpenMP threads %D can not be < 1 or > the MPI shared memory communicator size %d\n",nthreads,shm_comm_size);
  if (shm_comm_size % nthreads) { ierr = PetscPrintf(petsc_comm,"Warning: number of OpenMP threads %D is not a factor of the MPI shared memory communicator size %d, which may cause load-imbalance!\n",nthreads,shm_comm_size);CHKERRQ(ierr); }

  /* split shm_comm into a set of omp_comms with each of size nthreads. Ex., if
     shm_comm_size=16, nthreads=8, then ranks 0~7 get color 0 and ranks 8~15 get
     color 1. They are put in two omp_comms. Note that petsc_ranks may or may not
     be consecutive in a shm_comm, but shm_ranks always run from 0 to shm_comm_size-1.
     Use 0 as key so that rank ordering wont change in new comm.
   */
  color = shm_rank / nthreads;
  ierr  = MPI_Comm_split(shm_comm,color,0/*key*/,&ctrl->omp_comm);CHKERRQ(ierr);

  /* put rank 0's in omp_comms (i.e., master ranks) into a new comm - omp_master_comm */
  ierr = MPI_Comm_rank(ctrl->omp_comm,&omp_rank);CHKERRQ(ierr);
  if (!omp_rank) {
    ctrl->is_omp_master = PETSC_TRUE;  /* master */
    color = 0;
  } else {
    ctrl->is_omp_master = PETSC_FALSE; /* slave */
    color = MPI_UNDEFINED; /* to make slaves get omp_master_comm = MPI_COMM_NULL in MPI_Comm_split */
  }
  ierr = MPI_Comm_split(petsc_comm,color,0/*key*/,&ctrl->omp_master_comm);CHKERRQ(ierr); /* rank 0 in omp_master_comm is rank 0 in petsc_comm */

  /*=================================================================================
    Each omp_comm has a pthread_barrier_t in its shared memory, which is used to put
    slave ranks in sleep and idle their CPU, so that the master can fork OMP threads
    and run them on the idle CPUs.
   ==================================================================================*/
  ierr = PetscOmpCtrlCreateBarrier(ctrl);CHKERRQ(ierr);

  /*=================================================================================
    omp master logs its cpu binding (i.e., cpu set) and computes a new binding that
    is the union of the bindings of all ranks in the omp_comm
    =================================================================================*/

  ctrl->cpuset = hwloc_bitmap_alloc(); if (!ctrl->cpuset) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_LIB,"hwloc_bitmap_alloc() failed\n");
  ierr = hwloc_get_cpubind(ctrl->topology,ctrl->cpuset, HWLOC_CPUBIND_PROCESS);CHKERRQ(ierr);

  /* hwloc main developer said they will add new APIs hwloc_bitmap_{nr,to,from}_ulongs in 2.1 to help us simplify the following bitmap pack/unpack code */
  nr_cpu_ulongs = (hwloc_bitmap_last(hwloc_topology_get_topology_cpuset (ctrl->topology))+sizeof(unsigned long)*8)/sizeof(unsigned long)/8;
  ierr = PetscMalloc1(nr_cpu_ulongs,&cpu_ulongs);CHKERRQ(ierr);
  if (nr_cpu_ulongs == 1) {
    cpu_ulongs[0] = hwloc_bitmap_to_ulong(ctrl->cpuset);
  } else {
    for (i=0; i<nr_cpu_ulongs; i++) cpu_ulongs[i] = hwloc_bitmap_to_ith_ulong(ctrl->cpuset,(unsigned)i);
  }

  ierr = MPI_Reduce(ctrl->is_omp_master ? MPI_IN_PLACE : cpu_ulongs, cpu_ulongs,nr_cpu_ulongs, MPI_UNSIGNED_LONG,MPI_BOR,0,ctrl->omp_comm);CHKERRQ(ierr);

  if (ctrl->is_omp_master) {
    ctrl->omp_cpuset = hwloc_bitmap_alloc(); if (!ctrl->omp_cpuset) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_LIB,"hwloc_bitmap_alloc() failed\n");
    if (nr_cpu_ulongs == 1) {
#if HWLOC_API_VERSION >= 0x00020000
      ierr = hwloc_bitmap_from_ulong(ctrl->omp_cpuset,cpu_ulongs[0]);CHKERRQ(ierr);
#else
      hwloc_bitmap_from_ulong(ctrl->omp_cpuset,cpu_ulongs[0]);
#endif
    } else {
      for (i=0; i<nr_cpu_ulongs; i++)  {
#if HWLOC_API_VERSION >= 0x00020000
        ierr = hwloc_bitmap_set_ith_ulong(ctrl->omp_cpuset,(unsigned)i,cpu_ulongs[i]);CHKERRQ(ierr);
#else
        hwloc_bitmap_set_ith_ulong(ctrl->omp_cpuset,(unsigned)i,cpu_ulongs[i]);
#endif
      }
    }
  }

  ierr = PetscFree(cpu_ulongs);CHKERRQ(ierr);
  *pctrl = ctrl;
  PetscFunctionReturn(0);
}

/*@C
    PetscOmpCtrlDestroy - destroy the PETSc OpenMP controller

    Input Parameter:
.   pctrl  - a PETSc OpenMP controller

    Level: developer

.seealso PetscOmpCtrlCreate()
@*/
PetscErrorCode PetscOmpCtrlDestroy(PetscOmpCtrl *pctrl)
{
  PetscErrorCode  ierr;
  PetscOmpCtrl    ctrl = *pctrl;

  PetscFunctionBegin;
  hwloc_bitmap_free(ctrl->cpuset);
  hwloc_topology_destroy(ctrl->topology);
  PetscOmpCtrlDestroyBarrier(ctrl);
  ierr = MPI_Comm_free(&ctrl->omp_comm);CHKERRQ(ierr);
  if (ctrl->is_omp_master) {
    hwloc_bitmap_free(ctrl->omp_cpuset);
    ierr = MPI_Comm_free(&ctrl->omp_master_comm);CHKERRQ(ierr);
  }
  ierr = PetscFree(ctrl);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
    PetscOmpCtrlGetOmpComms - Get MPI communicators from a PETSc OMP controller

    Input Parameter:
.   ctrl - a PETSc OMP controller

    Output Parameter:
+   omp_comm         - a communicator that includes a master rank and slave ranks where master spawns threads
.   omp_master_comm  - on master ranks, return a communicator that include master ranks of each omp_comm;
                       on slave ranks, MPI_COMM_NULL will be return in reality.
-   is_omp_master    - true if the calling process is an OMP master rank.

    Notes: any output parameter can be NULL. The parameter is just ignored.

    Level: developer
@*/
PetscErrorCode PetscOmpCtrlGetOmpComms(PetscOmpCtrl ctrl,MPI_Comm *omp_comm,MPI_Comm *omp_master_comm,PetscBool *is_omp_master)
{
  PetscFunctionBegin;
  if (omp_comm)        *omp_comm        = ctrl->omp_comm;
  if (omp_master_comm) *omp_master_comm = ctrl->omp_master_comm;
  if (is_omp_master)   *is_omp_master   = ctrl->is_omp_master;
  PetscFunctionReturn(0);
}

/*@C
    PetscOmpCtrlBarrier - Do barrier on MPI ranks in omp_comm contained by the PETSc OMP controller (to let slave ranks free their CPU)

    Input Parameter:
.   ctrl - a PETSc OMP controller

    Notes:
    this is a pthread barrier on MPI processes. Using MPI_Barrier instead is conceptually correct. But MPI standard does not
    require processes blocked by MPI_Barrier free their CPUs to let other processes progress. In practice, to minilize latency,
    MPI processes stuck in MPI_Barrier keep polling and do not free CPUs. In contrast, pthread_barrier has this requirement.

    A code using PetscOmpCtrlBarrier() would be like this,

    if (is_omp_master) {
      PetscOmpCtrlOmpRegionOnMasterBegin(ctrl);
      Call the library using OpenMP
      PetscOmpCtrlOmpRegionOnMasterEnd(ctrl);
    }
    PetscOmpCtrlBarrier(ctrl);

    Level: developer

.seealso PetscOmpCtrlOmpRegionOnMasterBegin(), PetscOmpCtrlOmpRegionOnMasterEnd()
@*/
PetscErrorCode PetscOmpCtrlBarrier(PetscOmpCtrl ctrl)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = pthread_barrier_wait(ctrl->barrier);
  if (ierr && ierr != PTHREAD_BARRIER_SERIAL_THREAD) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_LIB,"pthread_barrier_wait failed within PetscOmpCtrlBarrier with return code %D\n", ierr);
  PetscFunctionReturn(0);
}

/*@C
    PetscOmpCtrlOmpRegionOnMasterBegin - Mark the beginning of an OpenMP library call on master ranks

    Input Parameter:
.   ctrl - a PETSc OMP controller

    Notes:
    Only master ranks can call this function. Call PetscOmpCtrlGetOmpComms() to know if this is a master rank.
    This function changes CPU binding of master ranks and nthreads-var of OpenMP runtime

    Level: developer

.seealso: PetscOmpCtrlOmpRegionOnMasterEnd()
@*/
PetscErrorCode PetscOmpCtrlOmpRegionOnMasterBegin(PetscOmpCtrl ctrl)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = hwloc_set_cpubind(ctrl->topology,ctrl->omp_cpuset,HWLOC_CPUBIND_PROCESS);CHKERRQ(ierr);
  omp_set_num_threads(ctrl->omp_comm_size); /* may override the OMP_NUM_THREAD env var */
  PetscFunctionReturn(0);
}

/*@C
   PetscOmpCtrlOmpRegionOnMasterEnd - Mark the end of an OpenMP library call on master ranks

   Input Parameter:
.  ctrl - a PETSc OMP controller

   Notes:
   Only master ranks can call this function. Call PetscOmpCtrlGetOmpComms() to know if this is a master rank.
   This function restores the CPU binding of master ranks and set and nthreads-var of OpenMP runtime to 1.

   Level: developer

.seealso: PetscOmpCtrlOmpRegionOnMasterBegin()
@*/
PetscErrorCode PetscOmpCtrlOmpRegionOnMasterEnd(PetscOmpCtrl ctrl)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = hwloc_set_cpubind(ctrl->topology,ctrl->cpuset,HWLOC_CPUBIND_PROCESS);CHKERRQ(ierr);
  omp_set_num_threads(1);
  PetscFunctionReturn(0);
}

#undef USE_MMAP_ALLOCATE_SHARED_MEMORY
#endif /* defined(PETSC_HAVE_OPENMP_SUPPORT) */
