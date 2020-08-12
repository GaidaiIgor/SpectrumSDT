#include <petsc/private/petscimpl.h>
#include <petscvalgrind.h>
#if defined(PETSC_HAVE_CUDA)
#include <cuda_runtime.h>
#include <petsccublas.h>
#endif

static PetscInt petsc_checkpointer_intensity = 1;

/*@
   PetscCheckPointerSetIntensity - An intense pointer check registers a signal handler and attempts to dereference to
   confirm whether the address is valid.  An intensity of 0 never uses signal handlers, 1 uses them when not in a "hot"
   function, and intensity of 2 always uses a signal handler.

   Not Collective

   Input Arguments:
.  intensity - how much to check pointers for validity

   Options Database:
.  -check_pointer_intensity - intensity (0, 1, or 2)

   Level: advanced

.seealso: PetscCheckPointer(), PetscFunctionBeginHot()
@*/
PetscErrorCode PetscCheckPointerSetIntensity(PetscInt intensity)
{

  PetscFunctionBegin;
  switch (intensity) {
  case 0:
  case 1:
  case 2:
    petsc_checkpointer_intensity = intensity;
    break;
  default: SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Intensity %D not in 0,1,2",intensity);
  }
  PetscFunctionReturn(0);
}

/* ---------------------------------------------------------------------------------------*/

#if defined(PETSC_HAVE_SETJMP_H)
#include <setjmp.h>
static jmp_buf PetscSegvJumpBuf;
static PetscBool PetscSegvJumpBuf_set;

/*@C
   PetscSignalSegvCheckPointerOrMpi - To be called from a signal handler for SIGSEGV.  If the signal was received
   while executing PetscCheckPointer()/PetscCheckMpiGpuAwareness(), this function longjmps back there, otherwise returns
   with no effect. This function is called automatically by PetscSignalHandlerDefault().

   Not Collective

   Level: developer

.seealso: PetscPushSignalHandler()
@*/
void PetscSignalSegvCheckPointerOrMpi(void) {
  if (PetscSegvJumpBuf_set) longjmp(PetscSegvJumpBuf,1);
}

/*@C
     PetscCheckPointer - Returns PETSC_TRUE if a pointer points to accessible data

   Not Collective

   Input Parameters:
+     ptr - the pointer
-     dtype - the type of data the pointer is suppose to point to

   Level: developer

.seealso: PetscCheckPointerSetIntensity()
@*/
PetscBool PetscCheckPointer(const void *ptr,PetscDataType dtype)
{

  if (PETSC_RUNNING_ON_VALGRIND) return PETSC_TRUE;
  if (!ptr) return PETSC_FALSE;
  if (petsc_checkpointer_intensity < 1) return PETSC_TRUE;

  /* Skip the verbose check if we are inside a hot function. */
  if (petscstack && petscstack->hotdepth > 0 && petsc_checkpointer_intensity < 2) return PETSC_TRUE;

  PetscSegvJumpBuf_set = PETSC_TRUE;

  if (setjmp(PetscSegvJumpBuf)) {
    /* A segv was triggered in the code below hence we return with an error code */
    PetscSegvJumpBuf_set = PETSC_FALSE;
    return PETSC_FALSE;
  } else {
    switch (dtype) {
    case PETSC_INT:{
      PETSC_UNUSED PetscInt x = (PetscInt)*(volatile PetscInt*)ptr;
      break;
    }
#if defined(PETSC_USE_COMPLEX)
    case PETSC_SCALAR:{         /* C++ is seriously dysfunctional with volatile std::complex. */
#if defined(PETSC_USE_CXXCOMPLEX)
      PetscReal xreal = ((volatile PetscReal*)ptr)[0],ximag = ((volatile PetscReal*)ptr)[1];
      PETSC_UNUSED volatile PetscScalar x = xreal + PETSC_i*ximag;
#else
      PETSC_UNUSED PetscScalar x = *(volatile PetscScalar*)ptr;
#endif
      break;
    }
#endif
    case PETSC_REAL:{
      PETSC_UNUSED PetscReal x = *(volatile PetscReal*)ptr;
      break;
    }
    case PETSC_BOOL:{
      PETSC_UNUSED PetscBool x = *(volatile PetscBool*)ptr;
      break;
    }
    case PETSC_ENUM:{
      PETSC_UNUSED PetscEnum x = *(volatile PetscEnum*)ptr;
      break;
    }
    case PETSC_CHAR:{
      PETSC_UNUSED char x = *(volatile char*)ptr;
      break;
    }
    case PETSC_OBJECT:{
      PETSC_UNUSED volatile PetscClassId classid = ((PetscObject)ptr)->classid;
      break;
    }
    default:;
    }
  }
  PetscSegvJumpBuf_set = PETSC_FALSE;
  return PETSC_TRUE;
}

#if defined (PETSC_HAVE_CUDA)
PetscBool PetscCheckMpiGpuAwareness(void)
{
  cudaError_t cerr=cudaSuccess;
  int         ierr,hbuf[2]={1,0},*dbuf=NULL;
  PetscBool   awareness=PETSC_FALSE;

  cerr = cudaMalloc((void**)&dbuf,sizeof(int)*2);if (cerr != cudaSuccess) return PETSC_FALSE;
  cerr = cudaMemcpy(dbuf,hbuf,sizeof(int)*2,cudaMemcpyHostToDevice);if (cerr != cudaSuccess) return PETSC_FALSE;

  PetscSegvJumpBuf_set = PETSC_TRUE;
  if (setjmp(PetscSegvJumpBuf)) {
    /* If a segv was triggered in the MPI_Allreduce below, it is very likely due to the MPI is not GPU-aware */
    awareness = PETSC_FALSE;
  } else {
    ierr = MPI_Allreduce(dbuf,dbuf+1,1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);
    if (!ierr) awareness = PETSC_TRUE;
  }
  PetscSegvJumpBuf_set = PETSC_FALSE;
  cerr = cudaFree(dbuf);if (cerr != cudaSuccess) return PETSC_FALSE;
  return awareness;
}
#endif
#else
void PetscSignalSegvCheckPointerOrMpi(void) {
  return;
}

PetscBool PetscCheckPointer(const void *ptr,PETSC_UNUSED PetscDataType dtype)
{
  if (!ptr) return PETSC_FALSE;
  return PETSC_TRUE;
}

PetscBool PetscCheckMpiGpuAwareness(void)
{
  /* If no setjmp (rare), return true and let users code run (and segfault if they should) */
  return PETSC_TRUE;
}
#endif
