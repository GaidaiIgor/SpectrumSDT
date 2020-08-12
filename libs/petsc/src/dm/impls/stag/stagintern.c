/* DMStag dimension-independent internal functions. If added to the public API,
   these would move to stagutils.c */

#include <petsc/private/dmstagimpl.h>

/* Note: this is an internal function but we provide a man page in case it's made public */
/*@C
  DMStagDuplicateWithoutSetup - duplicate a DMStag object without setting it up

  Collective

  Input Parameters:
+ dm - The original DM object
- comm - the MPI communicator for the new DM (MPI_COMM_NULL to use the same communicator as dm)

  Output Parameter:
. newdm  - The new DM object

  Developer Notes:
  Copies over all of the state for a DMStag object, except that which is
  populated during DMSetUp().  This function is used within (all) other
  functions that require an un-setup clone, which is common when duplicating,
  coarsening, refining, or creating compatible DMs with different fields.  For
  this reason it also accepts an MPI communicator as an argument (though note
  that at the time of this writing, implementations of DMCoarsen and DMRefine
  don't usually seem to respect their "comm" arguments). This function could be
  pushed up to the general DM API (and perhaps given a different name).

  Level: developer

  .seealso: DMClone(), DMStagCreateCompatibleDMStag(), DMCoarsen(), DMRefine()
@*/
PetscErrorCode DMStagDuplicateWithoutSetup(DM dm, MPI_Comm comm, DM *newdm)
{
  PetscErrorCode  ierr;
  DM_Stag * const stag  = (DM_Stag*)dm->data;
  DM_Stag         *newstag;
  PetscInt        dim;
  MPI_Comm        newcomm;

  PetscValidHeaderSpecificType(dm,DM_CLASSID,1,DMSTAG);
  newcomm = (comm == MPI_COMM_NULL) ? PetscObjectComm((PetscObject)dm) : comm;
  ierr = DMCreate(newcomm,newdm);CHKERRQ(ierr);
  ierr = DMGetDimension(dm,&dim);CHKERRQ(ierr);
  ierr = DMSetDimension(*newdm,dim);CHKERRQ(ierr);

  /* Call routine to define all data required for setup */
  ierr = DMStagInitialize(stag->boundaryType[0],stag->boundaryType[1],stag->boundaryType[2],stag->N[0],stag->N[1],stag->N[2],stag->nRanks[0],stag->nRanks[1],stag->nRanks[2],stag->dof[0],stag->dof[1],stag->dof[2],stag->dof[3],stag->stencilType,stag->stencilWidth,stag->l[0],stag->l[1],stag->l[2],*newdm);CHKERRQ(ierr);

  /* Copy all data unrelated to setup */
  newstag = (DM_Stag*)(*newdm)->data;
  ierr = PetscStrallocpy(stag->coordinateDMType,(char**)&newstag->coordinateDMType);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* Populate data created after DMCreate_Stag() is called, which is used by DMSetUp_Stag(),
   such as the grid dimensions and dof information. Arguments are ignored for dimensions
   less than three. */
PetscErrorCode DMStagInitialize(DMBoundaryType bndx,DMBoundaryType bndy,DMBoundaryType bndz,PetscInt M,PetscInt N,PetscInt P,PetscInt m,PetscInt n,PetscInt p,PetscInt dof0,PetscInt dof1,PetscInt dof2,PetscInt dof3,DMStagStencilType stencilType,PetscInt stencilWidth,const PetscInt lx[],const PetscInt ly[],const PetscInt lz[],DM dm)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMSetType(dm,DMSTAG);CHKERRQ(ierr);
  ierr = DMStagSetBoundaryTypes(dm,bndx,bndy,bndz);CHKERRQ(ierr);
  ierr = DMStagSetGlobalSizes(dm,M,N,P);CHKERRQ(ierr);
  ierr = DMStagSetNumRanks(dm,m,n,p);CHKERRQ(ierr);
  ierr = DMStagSetStencilType(dm,stencilType);CHKERRQ(ierr);
  ierr = DMStagSetStencilWidth(dm,stencilWidth);CHKERRQ(ierr);
  ierr = DMStagSetDOF(dm,dof0,dof1,dof2,dof3);CHKERRQ(ierr);
  ierr = DMStagSetOwnershipRanges(dm,lx,ly,lz);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
