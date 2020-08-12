
#include <petscsys.h>        /*I  "petscsys.h"  I*/
#include <petsc/private/petscimpl.h>
/*
    Note that tag of 0 is ok because comm is a private communicator
  generated below just for these routines.
*/

PETSC_INTERN PetscErrorCode PetscSequentialPhaseBegin_Private(MPI_Comm comm,int ng)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank,size,tag = 0;
  MPI_Status     status;

  PetscFunctionBegin;
  ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);
  if (size == 1) PetscFunctionReturn(0);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  if (rank) {
    ierr = MPI_Recv(NULL,0,MPI_INT,rank-1,tag,comm,&status);CHKERRQ(ierr);
  }
  /* Send to the next process in the group unless we are the last process */
  if ((rank % ng) < ng - 1 && rank != size - 1) {
    ierr = MPI_Send(NULL,0,MPI_INT,rank + 1,tag,comm);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PETSC_INTERN PetscErrorCode PetscSequentialPhaseEnd_Private(MPI_Comm comm,int ng)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank,size,tag = 0;
  MPI_Status     status;

  PetscFunctionBegin;
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);
  if (size == 1) PetscFunctionReturn(0);

  /* Send to the first process in the next group */
  if ((rank % ng) == ng - 1 || rank == size - 1) {
    ierr = MPI_Send(NULL,0,MPI_INT,(rank + 1) % size,tag,comm);CHKERRQ(ierr);
  }
  if (!rank) {
    ierr = MPI_Recv(NULL,0,MPI_INT,size-1,tag,comm,&status);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/* ---------------------------------------------------------------------*/
/*
    The variable Petsc_Seq_keyval is used to indicate an MPI attribute that
  is attached to a communicator that manages the sequential phase code below.
*/
PetscMPIInt Petsc_Seq_keyval = MPI_KEYVAL_INVALID;

/*@
   PetscSequentialPhaseBegin - Begins a sequential section of code.

   Collective

   Input Parameters:
+  comm - Communicator to sequentialize.
-  ng   - Number in processor group.  This many processes are allowed to execute
   at the same time (usually 1)

   Level: intermediate

   Notes:
   PetscSequentialPhaseBegin() and PetscSequentialPhaseEnd() provide a
   way to force a section of code to be executed by the processes in
   rank order.  Typically, this is done with
.vb
      PetscSequentialPhaseBegin(comm, 1);
      <code to be executed sequentially>
      PetscSequentialPhaseEnd(comm, 1);
.ve

   Often, the sequential code contains output statements (e.g., printf) to
   be executed.  Note that you may need to flush the I/O buffers before
   calling PetscSequentialPhaseEnd().  Also, note that some systems do
   not propagate I/O in any order to the controling terminal (in other words,
   even if you flush the output, you may not get the data in the order
   that you want).

.seealso: PetscSequentialPhaseEnd()

@*/
PetscErrorCode  PetscSequentialPhaseBegin(MPI_Comm comm,int ng)
{
  PetscErrorCode ierr;
  PetscMPIInt    size;
  MPI_Comm       local_comm,*addr_local_comm;

  PetscFunctionBegin;
  ierr = PetscSysInitializePackage();CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);
  if (size == 1) PetscFunctionReturn(0);

  /* Get the private communicator for the sequential operations */
  if (Petsc_Seq_keyval == MPI_KEYVAL_INVALID) {
    ierr = MPI_Comm_create_keyval(MPI_COMM_NULL_COPY_FN,MPI_COMM_NULL_DELETE_FN,&Petsc_Seq_keyval,NULL);CHKERRQ(ierr);
  }

  ierr = MPI_Comm_dup(comm,&local_comm);CHKERRQ(ierr);
  ierr = PetscMalloc1(1,&addr_local_comm);CHKERRQ(ierr);

  *addr_local_comm = local_comm;

  ierr = MPI_Comm_set_attr(comm,Petsc_Seq_keyval,(void*)addr_local_comm);CHKERRQ(ierr);
  ierr = PetscSequentialPhaseBegin_Private(local_comm,ng);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   PetscSequentialPhaseEnd - Ends a sequential section of code.

   Collective

   Input Parameters:
+  comm - Communicator to sequentialize.
-  ng   - Number in processor group.  This many processes are allowed to execute
   at the same time (usually 1)

   Level: intermediate

   Notes:
   See PetscSequentialPhaseBegin() for more details.

.seealso: PetscSequentialPhaseBegin()

@*/
PetscErrorCode  PetscSequentialPhaseEnd(MPI_Comm comm,int ng)
{
  PetscErrorCode ierr;
  PetscMPIInt    size,flag;
  MPI_Comm       local_comm,*addr_local_comm;

  PetscFunctionBegin;
  ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);
  if (size == 1) PetscFunctionReturn(0);

  ierr = MPI_Comm_get_attr(comm,Petsc_Seq_keyval,(void**)&addr_local_comm,&flag);CHKERRQ(ierr);
  if (!flag) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Wrong MPI communicator; must pass in one used with PetscSequentialPhaseBegin()");
  local_comm = *addr_local_comm;

  ierr = PetscSequentialPhaseEnd_Private(local_comm,ng);CHKERRQ(ierr);

  ierr = PetscFree(addr_local_comm);CHKERRQ(ierr);
  ierr = MPI_Comm_free(&local_comm);CHKERRQ(ierr);
  ierr = MPI_Comm_delete_attr(comm,Petsc_Seq_keyval);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
  PetscGlobalMinMaxInt - Get the global min/max from local min/max input

  Collective

  Input Parameter:
. minMaxVal - An array with the local min and max

  Output Parameter:
. minMaxValGlobal - An array with the global min and max

  Level: beginner

.seealso: PetscSplitOwnership()
@*/
PetscErrorCode PetscGlobalMinMaxInt(MPI_Comm comm, PetscInt minMaxVal[2], PetscInt minMaxValGlobal[2])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  minMaxVal[1] = -minMaxVal[1];
  ierr = MPI_Allreduce(minMaxVal, minMaxValGlobal, 2, MPIU_INT, MPI_MIN, comm);CHKERRQ(ierr);
  minMaxValGlobal[1] = -minMaxValGlobal[1];
  PetscFunctionReturn(0);
}

/*@C
  PetscGlobalMinMaxReal - Get the global min/max from local min/max input

  Collective

  Input Parameter:
. minMaxVal - An array with the local min and max

  Output Parameter:
. minMaxValGlobal - An array with the global min and max

  Level: beginner

.seealso: PetscSplitOwnership()
@*/
PetscErrorCode PetscGlobalMinMaxReal(MPI_Comm comm, PetscReal minMaxVal[2], PetscReal minMaxValGlobal[2])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  minMaxVal[1] = -minMaxVal[1];
  ierr = MPI_Allreduce(minMaxVal, minMaxValGlobal, 2, MPIU_REAL, MPI_MIN, comm);CHKERRQ(ierr);
  minMaxValGlobal[1] = -minMaxValGlobal[1];
  PetscFunctionReturn(0);
}
