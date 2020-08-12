#include <petsc/private/dmdaimpl.h>    /*I   "petscdmda.h"   I*/

/*@
  DMDASetSizes - Sets the number of grid points in the three dimensional directions

  Logically Collective on da

  Input Parameters:
+ da - the DMDA
. M - the global X size
. N - the global Y size
- P - the global Z size

  Level: intermediate

  Developer Notes:
  Since the dimension may not yet have been set the code cannot error check for non-positive Y and Z number of grid points

.seealso: PetscSplitOwnership()
@*/
PetscErrorCode  DMDASetSizes(DM da, PetscInt M, PetscInt N, PetscInt P)
{
  DM_DA *dd = (DM_DA*)da->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecificType(da, DM_CLASSID, 1,DMDA);
  PetscValidLogicalCollectiveInt(da,M,2);
  PetscValidLogicalCollectiveInt(da,N,3);
  PetscValidLogicalCollectiveInt(da,P,4);
  if (da->setupcalled) SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_ARG_WRONGSTATE,"This function must be called before DMSetUp()");
  if (M < 1) SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_ARG_SIZ,"Number of grid points in X direction must be positive");
  if (N < 0) SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_ARG_SIZ,"Number of grid points in Y direction must be positive");
  if (P < 0) SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_ARG_SIZ,"Number of grid points in Z direction must be positive");

  dd->M = M;
  dd->N = N;
  dd->P = P;
  PetscFunctionReturn(0);
}

/*@
  DMDASetNumProcs - Sets the number of processes in each dimension

  Logically Collective on da

  Input Parameters:
+ da - the DMDA
. m - the number of X procs (or PETSC_DECIDE)
. n - the number of Y procs (or PETSC_DECIDE)
- p - the number of Z procs (or PETSC_DECIDE)

  Level: intermediate

.seealso: DMDASetSizes(), PetscSplitOwnership()
@*/
PetscErrorCode  DMDASetNumProcs(DM da, PetscInt m, PetscInt n, PetscInt p)
{
  DM_DA          *dd = (DM_DA*)da->data;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecificType(da, DM_CLASSID, 1,DMDA);
  PetscValidLogicalCollectiveInt(da,m,2);
  PetscValidLogicalCollectiveInt(da,n,3);
  PetscValidLogicalCollectiveInt(da,p,4);
  if (da->setupcalled) SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_ARG_WRONGSTATE,"This function must be called before DMSetUp()");
  dd->m = m;
  dd->n = n;
  dd->p = p;
  if (da->dim == 2) {
    PetscMPIInt size;
    ierr = MPI_Comm_size(PetscObjectComm((PetscObject)da),&size);CHKERRQ(ierr);
    if ((dd->m > 0) && (dd->n < 0)) {
      dd->n = size/dd->m;
      if (dd->n*dd->m != size) SETERRQ2(PetscObjectComm((PetscObject)da),PETSC_ERR_ARG_OUTOFRANGE,"%D processes in X direction not divisible into comm size %d",m,size);
    }
    if ((dd->n > 0) && (dd->m < 0)) {
      dd->m = size/dd->n;
      if (dd->n*dd->m != size) SETERRQ2(PetscObjectComm((PetscObject)da),PETSC_ERR_ARG_OUTOFRANGE,"%D processes in Y direction not divisible into comm size %d",n,size);
    }
  }
  PetscFunctionReturn(0);
}

/*@
  DMDASetBoundaryType - Sets the type of ghost nodes on domain boundaries.

  Not collective

  Input Parameter:
+ da    - The DMDA
- bx,by,bz - One of DM_BOUNDARY_NONE, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_PERIODIC

  Level: intermediate

.seealso: DMDACreate(), DMDestroy(), DMDA, DMBoundaryType
@*/
PetscErrorCode  DMDASetBoundaryType(DM da,DMBoundaryType bx,DMBoundaryType by,DMBoundaryType bz)
{
  DM_DA *dd = (DM_DA*)da->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecificType(da,DM_CLASSID,1,DMDA);
  PetscValidLogicalCollectiveEnum(da,bx,2);
  PetscValidLogicalCollectiveEnum(da,by,3);
  PetscValidLogicalCollectiveEnum(da,bz,4);
  if (da->setupcalled) SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_ARG_WRONGSTATE,"This function must be called before DMSetUp()");
  dd->bx = bx;
  dd->by = by;
  dd->bz = bz;
  PetscFunctionReturn(0);
}

/*@
  DMDASetDof - Sets the number of degrees of freedom per vertex

  Not collective

  Input Parameters:
+ da  - The DMDA
- dof - Number of degrees of freedom

  Level: intermediate

.seealso: DMDAGetDof(), DMDACreate(), DMDestroy(), DMDA
@*/
PetscErrorCode  DMDASetDof(DM da, PetscInt dof)
{
  DM_DA *dd = (DM_DA*)da->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecificType(da,DM_CLASSID,1,DMDA);
  PetscValidLogicalCollectiveInt(da,dof,2);
  if (da->setupcalled) SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_ARG_WRONGSTATE,"This function must be called before DMSetUp()");
  dd->w  = dof;
  da->bs = dof;
  PetscFunctionReturn(0);
}

/*@
  DMDAGetDof - Gets the number of degrees of freedom per vertex

  Not collective

  Input Parameter:
. da  - The DMDA

  Output Parameter:
. dof - Number of degrees of freedom

  Level: intermediate

.seealso: DMDASetDof(), DMDACreate(), DMDestroy(), DMDA
@*/
PetscErrorCode DMDAGetDof(DM da, PetscInt *dof)
{
  DM_DA *dd = (DM_DA *) da->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecificType(da,DM_CLASSID,1,DMDA);
  PetscValidPointer(dof,2);
  *dof = dd->w;
  PetscFunctionReturn(0);
}

/*@
  DMDAGetOverlap - Gets the size of the per-processor overlap.

  Not collective

  Input Parameters:
. da  - The DMDA

  Output Parameters:
+ x   - Overlap in the x direction
. y   - Overlap in the y direction
- z   - Overlap in the z direction

  Level: intermediate

.seealso: DMDACreateDomainDecomposition(), DMDASetOverlap(), DMDA
@*/
PetscErrorCode  DMDAGetOverlap(DM da,PetscInt *x,PetscInt *y,PetscInt *z)
{
  DM_DA *dd = (DM_DA*)da->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecificType(da,DM_CLASSID,1,DMDA);
  if (x) *x = dd->xol;
  if (y) *y = dd->yol;
  if (z) *z = dd->zol;
  PetscFunctionReturn(0);
}

/*@
  DMDASetOverlap - Sets the size of the per-processor overlap.

  Not collective

  Input Parameters:
+ da  - The DMDA
. x   - Overlap in the x direction
. y   - Overlap in the y direction
- z   - Overlap in the z direction

  Level: intermediate

.seealso: DMDACreateDomainDecomposition(), DMDAGetOverlap(), DMDA
@*/
PetscErrorCode  DMDASetOverlap(DM da,PetscInt x,PetscInt y,PetscInt z)
{
  DM_DA *dd = (DM_DA*)da->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecificType(da,DM_CLASSID,1,DMDA);
  PetscValidLogicalCollectiveInt(da,x,2);
  PetscValidLogicalCollectiveInt(da,y,3);
  PetscValidLogicalCollectiveInt(da,z,4);
  dd->xol = x;
  dd->yol = y;
  dd->zol = z;
  PetscFunctionReturn(0);
}


/*@
  DMDAGetNumLocalSubDomains - Gets the number of local subdomains created upon decomposition.

  Not collective

  Input Parameters:
. da  - The DMDA

  Output Parameters:
. Nsub   - Number of local subdomains created upon decomposition

  Level: intermediate

.seealso: DMDACreateDomainDecomposition(), DMDASetNumLocalSubDomains(), DMDA
@*/
PetscErrorCode  DMDAGetNumLocalSubDomains(DM da,PetscInt *Nsub)
{
  DM_DA *dd = (DM_DA*)da->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecificType(da,DM_CLASSID,1,DMDA);
  if (Nsub) *Nsub = dd->Nsub;
  PetscFunctionReturn(0);
}

/*@
  DMDASetNumLocalSubDomains - Sets the number of local subdomains created upon decomposition.

  Not collective

  Input Parameters:
+ da  - The DMDA
- Nsub - The number of local subdomains requested

  Level: intermediate

.seealso: DMDACreateDomainDecomposition(), DMDAGetNumLocalSubDomains(), DMDA
@*/
PetscErrorCode  DMDASetNumLocalSubDomains(DM da,PetscInt Nsub)
{
  DM_DA *dd = (DM_DA*)da->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecificType(da,DM_CLASSID,1,DMDA);
  PetscValidLogicalCollectiveInt(da,Nsub,2);
  dd->Nsub = Nsub;
  PetscFunctionReturn(0);
}

/*@
  DMDASetOffset - Sets the index offset of the DA.

  Collective on DA

  Input Parameter:
+ da  - The DMDA
. xo  - The offset in the x direction
. yo  - The offset in the y direction
- zo  - The offset in the z direction

  Level: intermediate

  Notes:
    This is used primarily to overlap a computation on a local DA with that on a global DA without
  changing boundary conditions or subdomain features that depend upon the global offsets.

.seealso: DMDAGetOffset(), DMDAVecGetArray()
@*/
PetscErrorCode  DMDASetOffset(DM da, PetscInt xo, PetscInt yo, PetscInt zo, PetscInt Mo, PetscInt No, PetscInt Po)
{
  PetscErrorCode ierr;
  DM_DA          *dd = (DM_DA*)da->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecificType(da,DM_CLASSID,1,DMDA);
  PetscValidLogicalCollectiveInt(da,xo,2);
  PetscValidLogicalCollectiveInt(da,yo,3);
  PetscValidLogicalCollectiveInt(da,zo,4);
  PetscValidLogicalCollectiveInt(da,Mo,5);
  PetscValidLogicalCollectiveInt(da,No,6);
  PetscValidLogicalCollectiveInt(da,Po,7);
  dd->xo = xo;
  dd->yo = yo;
  dd->zo = zo;
  dd->Mo = Mo;
  dd->No = No;
  dd->Po = Po;

  if (da->coordinateDM) {
    ierr = DMDASetOffset(da->coordinateDM,xo,yo,zo,Mo,No,Po);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
  DMDAGetOffset - Gets the index offset of the DA.

  Not collective

  Input Parameter:
. da  - The DMDA

  Output Parameters:
+ xo  - The offset in the x direction
. yo  - The offset in the y direction
. zo  - The offset in the z direction
. Mo  - The global size in the x direction
. No  - The global size in the y direction
- Po  - The global size in the z direction

  Level: intermediate

.seealso: DMDASetOffset(), DMDAVecGetArray()
@*/
PetscErrorCode  DMDAGetOffset(DM da,PetscInt *xo,PetscInt *yo,PetscInt *zo,PetscInt *Mo,PetscInt *No,PetscInt *Po)
{
  DM_DA *dd = (DM_DA*)da->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecificType(da,DM_CLASSID,1,DMDA);
  if (xo) *xo = dd->xo;
  if (yo) *yo = dd->yo;
  if (zo) *zo = dd->zo;
  if (Mo) *Mo = dd->Mo;
  if (No) *No = dd->No;
  if (Po) *Po = dd->Po;
  PetscFunctionReturn(0);
}

/*@
  DMDAGetNonOverlappingRegion - Gets the indices of the nonoverlapping region of a subdomain DM.

  Not collective

  Input Parameter:
. da  - The DMDA

  Output Parameters:
+ xs  - The start of the region in x
. ys  - The start of the region in y
. zs  - The start of the region in z
. xs  - The size of the region in x
. ys  - The size of the region in y
- zs  - The size of the region in z

  Level: intermediate

.seealso: DMDAGetOffset(), DMDAVecGetArray()
@*/
PetscErrorCode  DMDAGetNonOverlappingRegion(DM da, PetscInt *xs, PetscInt *ys, PetscInt *zs, PetscInt *xm, PetscInt *ym, PetscInt *zm)
{
  DM_DA          *dd = (DM_DA*)da->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecificType(da,DM_CLASSID,1,DMDA);
  if (xs) *xs = dd->nonxs;
  if (ys) *ys = dd->nonys;
  if (zs) *zs = dd->nonzs;
  if (xm) *xm = dd->nonxm;
  if (ym) *ym = dd->nonym;
  if (zm) *zm = dd->nonzm;
  PetscFunctionReturn(0);
}


/*@
  DMDASetNonOverlappingRegion - Sets the indices of the nonoverlapping region of a subdomain DM.

  Collective on DA

  Input Parameter:
+ da  - The DMDA
. xs  - The start of the region in x
. ys  - The start of the region in y
. zs  - The start of the region in z
. xs  - The size of the region in x
. ys  - The size of the region in y
- zs  - The size of the region in z

  Level: intermediate

.seealso: DMDAGetOffset(), DMDAVecGetArray()
@*/
PetscErrorCode  DMDASetNonOverlappingRegion(DM da, PetscInt xs, PetscInt ys, PetscInt zs, PetscInt xm, PetscInt ym, PetscInt zm)
{
  DM_DA          *dd = (DM_DA*)da->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecificType(da,DM_CLASSID,1,DMDA);
  PetscValidLogicalCollectiveInt(da,xs,2);
  PetscValidLogicalCollectiveInt(da,ys,3);
  PetscValidLogicalCollectiveInt(da,zs,4);
  PetscValidLogicalCollectiveInt(da,xm,5);
  PetscValidLogicalCollectiveInt(da,ym,6);
  PetscValidLogicalCollectiveInt(da,zm,7);
  dd->nonxs = xs;
  dd->nonys = ys;
  dd->nonzs = zs;
  dd->nonxm = xm;
  dd->nonym = ym;
  dd->nonzm = zm;

  PetscFunctionReturn(0);
}

/*@
  DMDASetStencilType - Sets the type of the communication stencil

  Logically Collective on da

  Input Parameter:
+ da    - The DMDA
- stype - The stencil type, use either DMDA_STENCIL_BOX or DMDA_STENCIL_STAR.

  Level: intermediate

.seealso: DMDACreate(), DMDestroy(), DMDA
@*/
PetscErrorCode  DMDASetStencilType(DM da, DMDAStencilType stype)
{
  DM_DA *dd = (DM_DA*)da->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecificType(da,DM_CLASSID,1,DMDA);
  PetscValidLogicalCollectiveEnum(da,stype,2);
  if (da->setupcalled) SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_ARG_WRONGSTATE,"This function must be called before DMSetUp()");
  dd->stencil_type = stype;
  PetscFunctionReturn(0);
}

/*@
  DMDAGetStencilType - Gets the type of the communication stencil

  Not collective

  Input Parameter:
. da    - The DMDA

  Output Parameter:
. stype - The stencil type, use either DMDA_STENCIL_BOX or DMDA_STENCIL_STAR.

  Level: intermediate

.seealso: DMDACreate(), DMDestroy(), DMDA
@*/
PetscErrorCode DMDAGetStencilType(DM da, DMDAStencilType *stype)
{
  DM_DA *dd = (DM_DA*)da->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecificType(da,DM_CLASSID,1,DMDA);
  PetscValidPointer(stype,2);
  *stype = dd->stencil_type;
  PetscFunctionReturn(0);
}

/*@
  DMDASetStencilWidth - Sets the width of the communication stencil

  Logically Collective on da

  Input Parameter:
+ da    - The DMDA
- width - The stencil width

  Level: intermediate

.seealso: DMDACreate(), DMDestroy(), DMDA
@*/
PetscErrorCode  DMDASetStencilWidth(DM da, PetscInt width)
{
  DM_DA *dd = (DM_DA*)da->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecificType(da,DM_CLASSID,1,DMDA);
  PetscValidLogicalCollectiveInt(da,width,2);
  if (da->setupcalled) SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_ARG_WRONGSTATE,"This function must be called before DMSetUp()");
  dd->s = width;
  PetscFunctionReturn(0);
}

/*@
  DMDAGetStencilWidth - Gets the width of the communication stencil

  Not collective

  Input Parameter:
. da    - The DMDA

  Output Parameter:
. width - The stencil width

  Level: intermediate

.seealso: DMDACreate(), DMDestroy(), DMDA
@*/
PetscErrorCode DMDAGetStencilWidth(DM da, PetscInt *width)
{
  DM_DA *dd = (DM_DA *) da->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecificType(da,DM_CLASSID,1,DMDA);
  PetscValidPointer(width,2);
  *width = dd->s;
  PetscFunctionReturn(0);
}

static PetscErrorCode DMDACheckOwnershipRanges_Private(DM da,PetscInt M,PetscInt m,const PetscInt lx[])
{
  PetscInt i,sum;

  PetscFunctionBegin;
  if (M < 0) SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_ARG_WRONGSTATE,"Global dimension not set");
  for (i=sum=0; i<m; i++) sum += lx[i];
  if (sum != M) SETERRQ2(PetscObjectComm((PetscObject)da),PETSC_ERR_ARG_INCOMP,"Ownership ranges sum to %D but global dimension is %D",sum,M);
  PetscFunctionReturn(0);
}

/*@
  DMDASetOwnershipRanges - Sets the number of nodes in each direction on each process

  Logically Collective on da

  Input Parameter:
+ da - The DMDA
. lx - array containing number of nodes in the X direction on each process, or NULL. If non-null, must be of length da->m
. ly - array containing number of nodes in the Y direction on each process, or NULL. If non-null, must be of length da->n
- lz - array containing number of nodes in the Z direction on each process, or NULL. If non-null, must be of length da->p.

  Level: intermediate

  Note: these numbers are NOT multiplied by the number of dof per node.

.seealso: DMDACreate(), DMDestroy(), DMDA
@*/
PetscErrorCode  DMDASetOwnershipRanges(DM da, const PetscInt lx[], const PetscInt ly[], const PetscInt lz[])
{
  PetscErrorCode ierr;
  DM_DA          *dd = (DM_DA*)da->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecificType(da,DM_CLASSID,1,DMDA);
  if (da->setupcalled) SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_ARG_WRONGSTATE,"This function must be called before DMSetUp()");
  if (lx) {
    if (dd->m < 0) SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_ARG_WRONGSTATE,"Cannot set ownership ranges before setting number of procs");
    ierr = DMDACheckOwnershipRanges_Private(da,dd->M,dd->m,lx);CHKERRQ(ierr);
    if (!dd->lx) {
      ierr = PetscMalloc1(dd->m, &dd->lx);CHKERRQ(ierr);
    }
    ierr = PetscArraycpy(dd->lx, lx, dd->m);CHKERRQ(ierr);
  }
  if (ly) {
    if (dd->n < 0) SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_ARG_WRONGSTATE,"Cannot set ownership ranges before setting number of procs");
    ierr = DMDACheckOwnershipRanges_Private(da,dd->N,dd->n,ly);CHKERRQ(ierr);
    if (!dd->ly) {
      ierr = PetscMalloc1(dd->n, &dd->ly);CHKERRQ(ierr);
    }
    ierr = PetscArraycpy(dd->ly, ly, dd->n);CHKERRQ(ierr);
  }
  if (lz) {
    if (dd->p < 0) SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_ARG_WRONGSTATE,"Cannot set ownership ranges before setting number of procs");
    ierr = DMDACheckOwnershipRanges_Private(da,dd->P,dd->p,lz);CHKERRQ(ierr);
    if (!dd->lz) {
      ierr = PetscMalloc1(dd->p, &dd->lz);CHKERRQ(ierr);
    }
    ierr = PetscArraycpy(dd->lz, lz, dd->p);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
       DMDASetInterpolationType - Sets the type of interpolation that will be
          returned by DMCreateInterpolation()

   Logically Collective on da

   Input Parameter:
+  da - initial distributed array
-  ctype - DMDA_Q1 and DMDA_Q0 are currently the only supported forms

   Level: intermediate

   Notes:
    you should call this on the coarser of the two DMDAs you pass to DMCreateInterpolation()

.seealso: DMDACreate1d(), DMDACreate2d(), DMDACreate3d(), DMDestroy(), DMDA, DMDAInterpolationType
@*/
PetscErrorCode  DMDASetInterpolationType(DM da,DMDAInterpolationType ctype)
{
  DM_DA *dd = (DM_DA*)da->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecificType(da,DM_CLASSID,1,DMDA);
  PetscValidLogicalCollectiveEnum(da,ctype,2);
  dd->interptype = ctype;
  PetscFunctionReturn(0);
}

/*@
       DMDAGetInterpolationType - Gets the type of interpolation that will be
          used by DMCreateInterpolation()

   Not Collective

   Input Parameter:
.  da - distributed array

   Output Parameter:
.  ctype - interpolation type (DMDA_Q1 and DMDA_Q0 are currently the only supported forms)

   Level: intermediate

.seealso: DMDA, DMDAInterpolationType, DMDASetInterpolationType(), DMCreateInterpolation()
@*/
PetscErrorCode  DMDAGetInterpolationType(DM da,DMDAInterpolationType *ctype)
{
  DM_DA *dd = (DM_DA*)da->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecificType(da,DM_CLASSID,1,DMDA);
  PetscValidPointer(ctype,2);
  *ctype = dd->interptype;
  PetscFunctionReturn(0);
}

/*@C
      DMDAGetNeighbors - Gets an array containing the MPI rank of all the current
        processes neighbors.

    Not Collective

   Input Parameter:
.     da - the DMDA object

   Output Parameters:
.     ranks - the neighbors ranks, stored with the x index increasing most rapidly.
              this process itself is in the list

   Notes:
    In 2d the array is of length 9, in 3d of length 27
          Not supported in 1d
          Do not free the array, it is freed when the DMDA is destroyed.

   Fortran Notes:
    In fortran you must pass in an array of the appropriate length.

   Level: intermediate

@*/
PetscErrorCode  DMDAGetNeighbors(DM da,const PetscMPIInt *ranks[])
{
  DM_DA *dd = (DM_DA*)da->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecificType(da,DM_CLASSID,1,DMDA);
  *ranks = dd->neighbors;
  PetscFunctionReturn(0);
}

/*@C
      DMDAGetOwnershipRanges - Gets the ranges of indices in the x, y and z direction that are owned by each process

    Not Collective

   Input Parameter:
.     da - the DMDA object

   Output Parameter:
+     lx - ownership along x direction (optional)
.     ly - ownership along y direction (optional)
-     lz - ownership along z direction (optional)

   Level: intermediate

    Note: these correspond to the optional final arguments passed to DMDACreate(), DMDACreate2d(), DMDACreate3d()

    In Fortran one must pass in arrays lx, ly, and lz that are long enough to hold the values; the sixth, seventh and
    eighth arguments from DMDAGetInfo()

     In C you should not free these arrays, nor change the values in them. They will only have valid values while the
    DMDA they came from still exists (has not been destroyed).

    These numbers are NOT multiplied by the number of dof per node.

.seealso: DMDAGetCorners(), DMDAGetGhostCorners(), DMDACreate(), DMDACreate1d(), DMDACreate2d(), DMDACreate3d(), VecGetOwnershipRanges()
@*/
PetscErrorCode  DMDAGetOwnershipRanges(DM da,const PetscInt *lx[],const PetscInt *ly[],const PetscInt *lz[])
{
  DM_DA *dd = (DM_DA*)da->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecificType(da,DM_CLASSID,1,DMDA);
  if (lx) *lx = dd->lx;
  if (ly) *ly = dd->ly;
  if (lz) *lz = dd->lz;
  PetscFunctionReturn(0);
}

/*@
     DMDASetRefinementFactor - Set the ratios that the DMDA grid is refined

    Logically Collective on da

  Input Parameters:
+    da - the DMDA object
.    refine_x - ratio of fine grid to coarse in x direction (2 by default)
.    refine_y - ratio of fine grid to coarse in y direction (2 by default)
-    refine_z - ratio of fine grid to coarse in z direction (2 by default)

  Options Database:
+  -da_refine_x - refinement ratio in x direction
.  -da_refine_y - refinement ratio in y direction
-  -da_refine_z - refinement ratio in z direction

  Level: intermediate

    Notes:
    Pass PETSC_IGNORE to leave a value unchanged

.seealso: DMRefine(), DMDAGetRefinementFactor()
@*/
PetscErrorCode  DMDASetRefinementFactor(DM da, PetscInt refine_x, PetscInt refine_y,PetscInt refine_z)
{
  DM_DA *dd = (DM_DA*)da->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecificType(da,DM_CLASSID,1,DMDA);
  PetscValidLogicalCollectiveInt(da,refine_x,2);
  PetscValidLogicalCollectiveInt(da,refine_y,3);
  PetscValidLogicalCollectiveInt(da,refine_z,4);

  if (refine_x > 0) dd->refine_x = refine_x;
  if (refine_y > 0) dd->refine_y = refine_y;
  if (refine_z > 0) dd->refine_z = refine_z;
  PetscFunctionReturn(0);
}

/*@C
     DMDAGetRefinementFactor - Gets the ratios that the DMDA grid is refined

    Not Collective

  Input Parameter:
.    da - the DMDA object

  Output Parameters:
+    refine_x - ratio of fine grid to coarse in x direction (2 by default)
.    refine_y - ratio of fine grid to coarse in y direction (2 by default)
-    refine_z - ratio of fine grid to coarse in z direction (2 by default)

  Level: intermediate

    Notes:
    Pass NULL for values you do not need

.seealso: DMRefine(), DMDASetRefinementFactor()
@*/
PetscErrorCode  DMDAGetRefinementFactor(DM da, PetscInt *refine_x, PetscInt *refine_y,PetscInt *refine_z)
{
  DM_DA *dd = (DM_DA*)da->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecificType(da,DM_CLASSID,1,DMDA);
  if (refine_x) *refine_x = dd->refine_x;
  if (refine_y) *refine_y = dd->refine_y;
  if (refine_z) *refine_z = dd->refine_z;
  PetscFunctionReturn(0);
}

/*@C
     DMDASetGetMatrix - Sets the routine used by the DMDA to allocate a matrix.

    Logically Collective on da

  Input Parameters:
+    da - the DMDA object
-    f - the function that allocates the matrix for that specific DMDA

  Level: developer

   Notes:
    See DMDASetBlockFills() that provides a simple way to provide the nonzero structure for
       the diagonal and off-diagonal blocks of the matrix

   Not supported from Fortran

.seealso: DMCreateMatrix(), DMDASetBlockFills()
@*/
PetscErrorCode  DMDASetGetMatrix(DM da,PetscErrorCode (*f)(DM, Mat*))
{
  PetscFunctionBegin;
  PetscValidHeaderSpecificType(da,DM_CLASSID,1,DMDA);
  da->ops->creatematrix = f;
  PetscFunctionReturn(0);
}

/*
  Creates "balanced" ownership ranges after refinement, constrained by the need for the
  fine grid boundaries to fall within one stencil width of the coarse partition.

  Uses a greedy algorithm to handle non-ideal layouts, could probably do something better.
*/
static PetscErrorCode DMDARefineOwnershipRanges(DM da,PetscBool periodic,PetscInt stencil_width,PetscInt ratio,PetscInt m,const PetscInt lc[],PetscInt lf[])
{
  PetscInt       i,totalc = 0,remaining,startc = 0,startf = 0;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (ratio < 1) SETERRQ1(PetscObjectComm((PetscObject)da),PETSC_ERR_USER,"Requested refinement ratio %D must be at least 1",ratio);
  if (ratio == 1) {
    ierr = PetscArraycpy(lf,lc,m);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  for (i=0; i<m; i++) totalc += lc[i];
  remaining = (!periodic) + ratio * (totalc - (!periodic));
  for (i=0; i<m; i++) {
    PetscInt want = remaining/(m-i) + !!(remaining%(m-i));
    if (i == m-1) lf[i] = want;
    else {
      const PetscInt nextc = startc + lc[i];
      /* Move the first fine node of the next subdomain to the right until the coarse node on its left is within one
       * coarse stencil width of the first coarse node in the next subdomain. */
      while ((startf+want)/ratio < nextc - stencil_width) want++;
      /* Move the last fine node in the current subdomain to the left until the coarse node on its right is within one
       * coarse stencil width of the last coarse node in the current subdomain. */
      while ((startf+want-1+ratio-1)/ratio > nextc-1+stencil_width) want--;
      /* Make sure all constraints are satisfied */
      if (want < 0 || want > remaining || ((startf+want)/ratio < nextc - stencil_width)
          || ((startf+want-1+ratio-1)/ratio > nextc-1+stencil_width)) SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_ARG_SIZ,"Could not find a compatible refined ownership range");
    }
    lf[i]      = want;
    startc    += lc[i];
    startf    += lf[i];
    remaining -= lf[i];
  }
  PetscFunctionReturn(0);
}

/*
  Creates "balanced" ownership ranges after coarsening, constrained by the need for the
  fine grid boundaries to fall within one stencil width of the coarse partition.

  Uses a greedy algorithm to handle non-ideal layouts, could probably do something better.
*/
static PetscErrorCode DMDACoarsenOwnershipRanges(DM da,PetscBool periodic,PetscInt stencil_width,PetscInt ratio,PetscInt m,const PetscInt lf[],PetscInt lc[])
{
  PetscInt       i,totalf,remaining,startc,startf;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (ratio < 1) SETERRQ1(PetscObjectComm((PetscObject)da),PETSC_ERR_USER,"Requested refinement ratio %D must be at least 1",ratio);
  if (ratio == 1) {
    ierr = PetscArraycpy(lc,lf,m);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  for (i=0,totalf=0; i<m; i++) totalf += lf[i];
  remaining = (!periodic) + (totalf - (!periodic)) / ratio;
  for (i=0,startc=0,startf=0; i<m; i++) {
    PetscInt want = remaining/(m-i) + !!(remaining%(m-i));
    if (i == m-1) lc[i] = want;
    else {
      const PetscInt nextf = startf+lf[i];
      /* Slide first coarse node of next subdomain to the left until the coarse node to the left of the first fine
       * node is within one stencil width. */
      while (nextf/ratio < startc+want-stencil_width) want--;
      /* Slide the last coarse node of the current subdomain to the right until the coarse node to the right of the last
       * fine node is within one stencil width. */
      while ((nextf-1+ratio-1)/ratio > startc+want-1+stencil_width) want++;
      if (want < 0 || want > remaining
          || (nextf/ratio < startc+want-stencil_width) || ((nextf-1+ratio-1)/ratio > startc+want-1+stencil_width)) SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_ARG_SIZ,"Could not find a compatible coarsened ownership range");
    }
    lc[i]      = want;
    startc    += lc[i];
    startf    += lf[i];
    remaining -= lc[i];
  }
  PetscFunctionReturn(0);
}

PetscErrorCode  DMRefine_DA(DM da,MPI_Comm comm,DM *daref)
{
  PetscErrorCode ierr;
  PetscInt       M,N,P,i,dim;
  DM             da2;
  DM_DA          *dd = (DM_DA*)da->data,*dd2;

  PetscFunctionBegin;
  PetscValidHeaderSpecificType(da,DM_CLASSID,1,DMDA);
  PetscValidPointer(daref,3);

  ierr = DMGetDimension(da, &dim);CHKERRQ(ierr);
  if (dd->bx == DM_BOUNDARY_PERIODIC || dd->interptype == DMDA_Q0) {
    M = dd->refine_x*dd->M;
  } else {
    M = 1 + dd->refine_x*(dd->M - 1);
  }
  if (dd->by == DM_BOUNDARY_PERIODIC || dd->interptype == DMDA_Q0) {
    if (dim > 1) {
      N = dd->refine_y*dd->N;
    } else {
      N = 1;
    }
  } else {
    N = 1 + dd->refine_y*(dd->N - 1);
  }
  if (dd->bz == DM_BOUNDARY_PERIODIC || dd->interptype == DMDA_Q0) {
    if (dim > 2) {
      P = dd->refine_z*dd->P;
    } else {
      P = 1;
    }
  } else {
    P = 1 + dd->refine_z*(dd->P - 1);
  }
  ierr = DMDACreate(PetscObjectComm((PetscObject)da),&da2);CHKERRQ(ierr);
  ierr = DMSetOptionsPrefix(da2,((PetscObject)da)->prefix);CHKERRQ(ierr);
  ierr = DMSetDimension(da2,dim);CHKERRQ(ierr);
  ierr = DMDASetSizes(da2,M,N,P);CHKERRQ(ierr);
  ierr = DMDASetNumProcs(da2,dd->m,dd->n,dd->p);CHKERRQ(ierr);
  ierr = DMDASetBoundaryType(da2,dd->bx,dd->by,dd->bz);CHKERRQ(ierr);
  ierr = DMDASetDof(da2,dd->w);CHKERRQ(ierr);
  ierr = DMDASetStencilType(da2,dd->stencil_type);CHKERRQ(ierr);
  ierr = DMDASetStencilWidth(da2,dd->s);CHKERRQ(ierr);
  if (dim == 3) {
    PetscInt *lx,*ly,*lz;
    ierr = PetscMalloc3(dd->m,&lx,dd->n,&ly,dd->p,&lz);CHKERRQ(ierr);
    ierr = DMDARefineOwnershipRanges(da,(PetscBool)(dd->bx == DM_BOUNDARY_PERIODIC || dd->interptype == DMDA_Q0),dd->s,dd->refine_x,dd->m,dd->lx,lx);CHKERRQ(ierr);
    ierr = DMDARefineOwnershipRanges(da,(PetscBool)(dd->by == DM_BOUNDARY_PERIODIC || dd->interptype == DMDA_Q0),dd->s,dd->refine_y,dd->n,dd->ly,ly);CHKERRQ(ierr);
    ierr = DMDARefineOwnershipRanges(da,(PetscBool)(dd->bz == DM_BOUNDARY_PERIODIC || dd->interptype == DMDA_Q0),dd->s,dd->refine_z,dd->p,dd->lz,lz);CHKERRQ(ierr);
    ierr = DMDASetOwnershipRanges(da2,lx,ly,lz);CHKERRQ(ierr);
    ierr = PetscFree3(lx,ly,lz);CHKERRQ(ierr);
  } else if (dim == 2) {
    PetscInt *lx,*ly;
    ierr = PetscMalloc2(dd->m,&lx,dd->n,&ly);CHKERRQ(ierr);
    ierr = DMDARefineOwnershipRanges(da,(PetscBool)(dd->bx == DM_BOUNDARY_PERIODIC || dd->interptype == DMDA_Q0),dd->s,dd->refine_x,dd->m,dd->lx,lx);CHKERRQ(ierr);
    ierr = DMDARefineOwnershipRanges(da,(PetscBool)(dd->by == DM_BOUNDARY_PERIODIC || dd->interptype == DMDA_Q0),dd->s,dd->refine_y,dd->n,dd->ly,ly);CHKERRQ(ierr);
    ierr = DMDASetOwnershipRanges(da2,lx,ly,NULL);CHKERRQ(ierr);
    ierr = PetscFree2(lx,ly);CHKERRQ(ierr);
  } else if (dim == 1) {
    PetscInt *lx;
    ierr = PetscMalloc1(dd->m,&lx);CHKERRQ(ierr);
    ierr = DMDARefineOwnershipRanges(da,(PetscBool)(dd->bx == DM_BOUNDARY_PERIODIC || dd->interptype == DMDA_Q0),dd->s,dd->refine_x,dd->m,dd->lx,lx);CHKERRQ(ierr);
    ierr = DMDASetOwnershipRanges(da2,lx,NULL,NULL);CHKERRQ(ierr);
    ierr = PetscFree(lx);CHKERRQ(ierr);
  }
  dd2 = (DM_DA*)da2->data;

  /* allow overloaded (user replaced) operations to be inherited by refinement clones */
  da2->ops->creatematrix = da->ops->creatematrix;
  /* da2->ops->createinterpolation = da->ops->createinterpolation; this causes problem with SNESVI */
  da2->ops->getcoloring = da->ops->getcoloring;
  dd2->interptype       = dd->interptype;

  /* copy fill information if given */
  if (dd->dfill) {
    ierr = PetscMalloc1(dd->dfill[dd->w]+dd->w+1,&dd2->dfill);CHKERRQ(ierr);
    ierr = PetscArraycpy(dd2->dfill,dd->dfill,dd->dfill[dd->w]+dd->w+1);CHKERRQ(ierr);
  }
  if (dd->ofill) {
    ierr = PetscMalloc1(dd->ofill[dd->w]+dd->w+1,&dd2->ofill);CHKERRQ(ierr);
    ierr = PetscArraycpy(dd2->ofill,dd->ofill,dd->ofill[dd->w]+dd->w+1);CHKERRQ(ierr);
  }
  /* copy the refine information */
  dd2->coarsen_x = dd2->refine_x = dd->refine_x;
  dd2->coarsen_y = dd2->refine_y = dd->refine_y;
  dd2->coarsen_z = dd2->refine_z = dd->refine_z;

  if (dd->refine_z_hier) {
    if (da->levelup - da->leveldown + 1 > -1 && da->levelup - da->leveldown + 1 < dd->refine_z_hier_n) {
      dd2->refine_z = dd->refine_z_hier[da->levelup - da->leveldown + 1];
    }
    if (da->levelup - da->leveldown > -1 && da->levelup - da->leveldown < dd->refine_z_hier_n) {
      dd2->coarsen_z = dd->refine_z_hier[da->levelup - da->leveldown];
    }
    dd2->refine_z_hier_n = dd->refine_z_hier_n;
    ierr = PetscMalloc1(dd2->refine_z_hier_n,&dd2->refine_z_hier);CHKERRQ(ierr);
    ierr = PetscArraycpy(dd2->refine_z_hier,dd->refine_z_hier,dd2->refine_z_hier_n);CHKERRQ(ierr);
  }
  if (dd->refine_y_hier) {
    if (da->levelup - da->leveldown + 1 > -1 && da->levelup - da->leveldown + 1 < dd->refine_y_hier_n) {
      dd2->refine_y = dd->refine_y_hier[da->levelup - da->leveldown + 1];
    }
    if (da->levelup - da->leveldown > -1 && da->levelup - da->leveldown < dd->refine_y_hier_n) {
      dd2->coarsen_y = dd->refine_y_hier[da->levelup - da->leveldown];
    }
    dd2->refine_y_hier_n = dd->refine_y_hier_n;
    ierr = PetscMalloc1(dd2->refine_y_hier_n,&dd2->refine_y_hier);CHKERRQ(ierr);
    ierr = PetscArraycpy(dd2->refine_y_hier,dd->refine_y_hier,dd2->refine_y_hier_n);CHKERRQ(ierr);
  }
  if (dd->refine_x_hier) {
    if (da->levelup - da->leveldown + 1 > -1 && da->levelup - da->leveldown + 1 < dd->refine_x_hier_n) {
      dd2->refine_x = dd->refine_x_hier[da->levelup - da->leveldown + 1];
    }
    if (da->levelup - da->leveldown > -1 && da->levelup - da->leveldown < dd->refine_x_hier_n) {
      dd2->coarsen_x = dd->refine_x_hier[da->levelup - da->leveldown];
    }
    dd2->refine_x_hier_n = dd->refine_x_hier_n;
    ierr = PetscMalloc1(dd2->refine_x_hier_n,&dd2->refine_x_hier);CHKERRQ(ierr);
    ierr = PetscArraycpy(dd2->refine_x_hier,dd->refine_x_hier,dd2->refine_x_hier_n);CHKERRQ(ierr);
  }


  /* copy vector type information */
  ierr = DMSetVecType(da2,da->vectype);CHKERRQ(ierr);

  dd2->lf = dd->lf;
  dd2->lj = dd->lj;

  da2->leveldown = da->leveldown;
  da2->levelup   = da->levelup + 1;

  ierr = DMSetUp(da2);CHKERRQ(ierr);

  /* interpolate coordinates if they are set on the coarse grid */
  if (da->coordinates) {
    DM  cdaf,cdac;
    Vec coordsc,coordsf;
    Mat II;

    ierr = DMGetCoordinateDM(da,&cdac);CHKERRQ(ierr);
    ierr = DMGetCoordinates(da,&coordsc);CHKERRQ(ierr);
    ierr = DMGetCoordinateDM(da2,&cdaf);CHKERRQ(ierr);
    /* force creation of the coordinate vector */
    ierr = DMDASetUniformCoordinates(da2,0.0,1.0,0.0,1.0,0.0,1.0);CHKERRQ(ierr);
    ierr = DMGetCoordinates(da2,&coordsf);CHKERRQ(ierr);
    ierr = DMCreateInterpolation(cdac,cdaf,&II,NULL);CHKERRQ(ierr);
    ierr = MatInterpolate(II,coordsc,coordsf);CHKERRQ(ierr);
    ierr = MatDestroy(&II);CHKERRQ(ierr);
  }

  for (i=0; i<da->bs; i++) {
    const char *fieldname;
    ierr = DMDAGetFieldName(da,i,&fieldname);CHKERRQ(ierr);
    ierr = DMDASetFieldName(da2,i,fieldname);CHKERRQ(ierr);
  }

  *daref = da2;
  PetscFunctionReturn(0);
}


PetscErrorCode  DMCoarsen_DA(DM dmf, MPI_Comm comm,DM *dmc)
{
  PetscErrorCode ierr;
  PetscInt       M,N,P,i,dim;
  DM             dmc2;
  DM_DA          *dd = (DM_DA*)dmf->data,*dd2;

  PetscFunctionBegin;
  PetscValidHeaderSpecificType(dmf,DM_CLASSID,1,DMDA);
  PetscValidPointer(dmc,3);

  ierr = DMGetDimension(dmf, &dim);CHKERRQ(ierr);
  if (dd->bx == DM_BOUNDARY_PERIODIC || dd->interptype == DMDA_Q0) {
    M = dd->M / dd->coarsen_x;
  } else {
    M = 1 + (dd->M - 1) / dd->coarsen_x;
  }
  if (dd->by == DM_BOUNDARY_PERIODIC || dd->interptype == DMDA_Q0) {
    if (dim > 1) {
      N = dd->N / dd->coarsen_y;
    } else {
      N = 1;
    }
  } else {
    N = 1 + (dd->N - 1) / dd->coarsen_y;
  }
  if (dd->bz == DM_BOUNDARY_PERIODIC || dd->interptype == DMDA_Q0) {
    if (dim > 2) {
      P = dd->P / dd->coarsen_z;
    } else {
      P = 1;
    }
  } else {
    P = 1 + (dd->P - 1) / dd->coarsen_z;
  }
  ierr = DMDACreate(PetscObjectComm((PetscObject)dmf),&dmc2);CHKERRQ(ierr);
  ierr = DMSetOptionsPrefix(dmc2,((PetscObject)dmf)->prefix);CHKERRQ(ierr);
  ierr = DMSetDimension(dmc2,dim);CHKERRQ(ierr);
  ierr = DMDASetSizes(dmc2,M,N,P);CHKERRQ(ierr);
  ierr = DMDASetNumProcs(dmc2,dd->m,dd->n,dd->p);CHKERRQ(ierr);
  ierr = DMDASetBoundaryType(dmc2,dd->bx,dd->by,dd->bz);CHKERRQ(ierr);
  ierr = DMDASetDof(dmc2,dd->w);CHKERRQ(ierr);
  ierr = DMDASetStencilType(dmc2,dd->stencil_type);CHKERRQ(ierr);
  ierr = DMDASetStencilWidth(dmc2,dd->s);CHKERRQ(ierr);
  if (dim == 3) {
    PetscInt *lx,*ly,*lz;
    ierr = PetscMalloc3(dd->m,&lx,dd->n,&ly,dd->p,&lz);CHKERRQ(ierr);
    ierr = DMDACoarsenOwnershipRanges(dmf,(PetscBool)(dd->bx == DM_BOUNDARY_PERIODIC || dd->interptype == DMDA_Q0),dd->s,dd->coarsen_x,dd->m,dd->lx,lx);CHKERRQ(ierr);
    ierr = DMDACoarsenOwnershipRanges(dmf,(PetscBool)(dd->by == DM_BOUNDARY_PERIODIC || dd->interptype == DMDA_Q0),dd->s,dd->coarsen_y,dd->n,dd->ly,ly);CHKERRQ(ierr);
    ierr = DMDACoarsenOwnershipRanges(dmf,(PetscBool)(dd->bz == DM_BOUNDARY_PERIODIC || dd->interptype == DMDA_Q0),dd->s,dd->coarsen_z,dd->p,dd->lz,lz);CHKERRQ(ierr);
    ierr = DMDASetOwnershipRanges(dmc2,lx,ly,lz);CHKERRQ(ierr);
    ierr = PetscFree3(lx,ly,lz);CHKERRQ(ierr);
  } else if (dim == 2) {
    PetscInt *lx,*ly;
    ierr = PetscMalloc2(dd->m,&lx,dd->n,&ly);CHKERRQ(ierr);
    ierr = DMDACoarsenOwnershipRanges(dmf,(PetscBool)(dd->bx == DM_BOUNDARY_PERIODIC || dd->interptype == DMDA_Q0),dd->s,dd->coarsen_x,dd->m,dd->lx,lx);CHKERRQ(ierr);
    ierr = DMDACoarsenOwnershipRanges(dmf,(PetscBool)(dd->by == DM_BOUNDARY_PERIODIC || dd->interptype == DMDA_Q0),dd->s,dd->coarsen_y,dd->n,dd->ly,ly);CHKERRQ(ierr);
    ierr = DMDASetOwnershipRanges(dmc2,lx,ly,NULL);CHKERRQ(ierr);
    ierr = PetscFree2(lx,ly);CHKERRQ(ierr);
  } else if (dim == 1) {
    PetscInt *lx;
    ierr = PetscMalloc1(dd->m,&lx);CHKERRQ(ierr);
    ierr = DMDACoarsenOwnershipRanges(dmf,(PetscBool)(dd->bx == DM_BOUNDARY_PERIODIC || dd->interptype == DMDA_Q0),dd->s,dd->coarsen_x,dd->m,dd->lx,lx);CHKERRQ(ierr);
    ierr = DMDASetOwnershipRanges(dmc2,lx,NULL,NULL);CHKERRQ(ierr);
    ierr = PetscFree(lx);CHKERRQ(ierr);
  }
  dd2 = (DM_DA*)dmc2->data;

  /* allow overloaded (user replaced) operations to be inherited by refinement clones; why are only some inherited and not all? */
  /* dmc2->ops->createinterpolation = dmf->ops->createinterpolation; copying this one causes trouble for DMSetVI */
  dmc2->ops->creatematrix = dmf->ops->creatematrix;
  dmc2->ops->getcoloring  = dmf->ops->getcoloring;
  dd2->interptype         = dd->interptype;

  /* copy fill information if given */
  if (dd->dfill) {
    ierr = PetscMalloc1(dd->dfill[dd->w]+dd->w+1,&dd2->dfill);CHKERRQ(ierr);
    ierr = PetscArraycpy(dd2->dfill,dd->dfill,dd->dfill[dd->w]+dd->w+1);CHKERRQ(ierr);
  }
  if (dd->ofill) {
    ierr = PetscMalloc1(dd->ofill[dd->w]+dd->w+1,&dd2->ofill);CHKERRQ(ierr);
    ierr = PetscArraycpy(dd2->ofill,dd->ofill,dd->ofill[dd->w]+dd->w+1);CHKERRQ(ierr);
  }
  /* copy the refine information */
  dd2->coarsen_x = dd2->refine_x = dd->coarsen_x;
  dd2->coarsen_y = dd2->refine_y = dd->coarsen_y;
  dd2->coarsen_z = dd2->refine_z = dd->coarsen_z;

  if (dd->refine_z_hier) {
    if (dmf->levelup - dmf->leveldown -1 > -1 && dmf->levelup - dmf->leveldown - 1< dd->refine_z_hier_n) {
      dd2->refine_z = dd->refine_z_hier[dmf->levelup - dmf->leveldown - 1];
    }
    if (dmf->levelup - dmf->leveldown - 2 > -1 && dmf->levelup - dmf->leveldown - 2 < dd->refine_z_hier_n) {
      dd2->coarsen_z = dd->refine_z_hier[dmf->levelup - dmf->leveldown - 2];
    }
    dd2->refine_z_hier_n = dd->refine_z_hier_n;
    ierr = PetscMalloc1(dd2->refine_z_hier_n,&dd2->refine_z_hier);CHKERRQ(ierr);
    ierr = PetscArraycpy(dd2->refine_z_hier,dd->refine_z_hier,dd2->refine_z_hier_n);CHKERRQ(ierr);
  }
  if (dd->refine_y_hier) {
    if (dmf->levelup - dmf->leveldown - 1 > -1 && dmf->levelup - dmf->leveldown - 1< dd->refine_y_hier_n) {
      dd2->refine_y = dd->refine_y_hier[dmf->levelup - dmf->leveldown - 1];
    }
    if (dmf->levelup - dmf->leveldown - 2 > -1 && dmf->levelup - dmf->leveldown - 2 < dd->refine_y_hier_n) {
      dd2->coarsen_y = dd->refine_y_hier[dmf->levelup - dmf->leveldown - 2];
    }
    dd2->refine_y_hier_n = dd->refine_y_hier_n;
    ierr = PetscMalloc1(dd2->refine_y_hier_n,&dd2->refine_y_hier);CHKERRQ(ierr);
    ierr = PetscArraycpy(dd2->refine_y_hier,dd->refine_y_hier,dd2->refine_y_hier_n);CHKERRQ(ierr);
  }
  if (dd->refine_x_hier) {
    if (dmf->levelup - dmf->leveldown - 1 > -1 && dmf->levelup - dmf->leveldown - 1 < dd->refine_x_hier_n) {
      dd2->refine_x = dd->refine_x_hier[dmf->levelup - dmf->leveldown - 1];
    }
    if (dmf->levelup - dmf->leveldown - 2 > -1 && dmf->levelup - dmf->leveldown - 2 < dd->refine_x_hier_n) {
      dd2->coarsen_x = dd->refine_x_hier[dmf->levelup - dmf->leveldown - 2];
    }
    dd2->refine_x_hier_n = dd->refine_x_hier_n;
    ierr = PetscMalloc1(dd2->refine_x_hier_n,&dd2->refine_x_hier);CHKERRQ(ierr);
    ierr = PetscArraycpy(dd2->refine_x_hier,dd->refine_x_hier,dd2->refine_x_hier_n);CHKERRQ(ierr);
  }

  /* copy vector type information */
  ierr = DMSetVecType(dmc2,dmf->vectype);CHKERRQ(ierr);

  dd2->lf = dd->lf;
  dd2->lj = dd->lj;

  dmc2->leveldown = dmf->leveldown + 1;
  dmc2->levelup   = dmf->levelup;

  ierr = DMSetUp(dmc2);CHKERRQ(ierr);

  /* inject coordinates if they are set on the fine grid */
  if (dmf->coordinates) {
    DM         cdaf,cdac;
    Vec        coordsc,coordsf;
    Mat        inject;
    VecScatter vscat;

    ierr = DMGetCoordinateDM(dmf,&cdaf);CHKERRQ(ierr);
    ierr = DMGetCoordinates(dmf,&coordsf);CHKERRQ(ierr);
    ierr = DMGetCoordinateDM(dmc2,&cdac);CHKERRQ(ierr);
    /* force creation of the coordinate vector */
    ierr = DMDASetUniformCoordinates(dmc2,0.0,1.0,0.0,1.0,0.0,1.0);CHKERRQ(ierr);
    ierr = DMGetCoordinates(dmc2,&coordsc);CHKERRQ(ierr);

    ierr = DMCreateInjection(cdac,cdaf,&inject);CHKERRQ(ierr);
    ierr = MatScatterGetVecScatter(inject,&vscat);CHKERRQ(ierr);
    ierr = VecScatterBegin(vscat,coordsf,coordsc,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecScatterEnd(vscat,coordsf,coordsc,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = MatDestroy(&inject);CHKERRQ(ierr);
  }

  for (i=0; i<dmf->bs; i++) {
    const char *fieldname;
    ierr = DMDAGetFieldName(dmf,i,&fieldname);CHKERRQ(ierr);
    ierr = DMDASetFieldName(dmc2,i,fieldname);CHKERRQ(ierr);
  }

  *dmc = dmc2;
  PetscFunctionReturn(0);
}

PetscErrorCode  DMRefineHierarchy_DA(DM da,PetscInt nlevels,DM daf[])
{
  PetscErrorCode ierr;
  PetscInt       i,n,*refx,*refy,*refz;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(da,DM_CLASSID,1);
  if (nlevels < 0) SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_ARG_OUTOFRANGE,"nlevels cannot be negative");
  if (nlevels == 0) PetscFunctionReturn(0);
  PetscValidPointer(daf,3);

  /* Get refinement factors, defaults taken from the coarse DMDA */
  ierr = PetscMalloc3(nlevels,&refx,nlevels,&refy,nlevels,&refz);CHKERRQ(ierr);
  for (i=0; i<nlevels; i++) {
    ierr = DMDAGetRefinementFactor(da,&refx[i],&refy[i],&refz[i]);CHKERRQ(ierr);
  }
  n    = nlevels;
  ierr = PetscOptionsGetIntArray(((PetscObject)da)->options,((PetscObject)da)->prefix,"-da_refine_hierarchy_x",refx,&n,NULL);CHKERRQ(ierr);
  n    = nlevels;
  ierr = PetscOptionsGetIntArray(((PetscObject)da)->options,((PetscObject)da)->prefix,"-da_refine_hierarchy_y",refy,&n,NULL);CHKERRQ(ierr);
  n    = nlevels;
  ierr = PetscOptionsGetIntArray(((PetscObject)da)->options,((PetscObject)da)->prefix,"-da_refine_hierarchy_z",refz,&n,NULL);CHKERRQ(ierr);

  ierr = DMDASetRefinementFactor(da,refx[0],refy[0],refz[0]);CHKERRQ(ierr);
  ierr = DMRefine(da,PetscObjectComm((PetscObject)da),&daf[0]);CHKERRQ(ierr);
  for (i=1; i<nlevels; i++) {
    ierr = DMDASetRefinementFactor(daf[i-1],refx[i],refy[i],refz[i]);CHKERRQ(ierr);
    ierr = DMRefine(daf[i-1],PetscObjectComm((PetscObject)da),&daf[i]);CHKERRQ(ierr);
  }
  ierr = PetscFree3(refx,refy,refz);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode  DMCoarsenHierarchy_DA(DM da,PetscInt nlevels,DM dac[])
{
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(da,DM_CLASSID,1);
  if (nlevels < 0) SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_ARG_OUTOFRANGE,"nlevels cannot be negative");
  if (nlevels == 0) PetscFunctionReturn(0);
  PetscValidPointer(dac,3);
  ierr = DMCoarsen(da,PetscObjectComm((PetscObject)da),&dac[0]);CHKERRQ(ierr);
  for (i=1; i<nlevels; i++) {
    ierr = DMCoarsen(dac[i-1],PetscObjectComm((PetscObject)da),&dac[i]);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DMDASetGLLCoordinates_1d(DM dm,PetscInt n,PetscReal *nodes)
{
  PetscErrorCode ierr;
  PetscInt       i,j,xs,xn,q;
  PetscScalar    *xx;
  PetscReal      h;
  Vec            x;
  DM_DA          *da = (DM_DA*)dm->data;

  PetscFunctionBegin;
  if (da->bx != DM_BOUNDARY_PERIODIC) {
    ierr = DMDAGetInfo(dm,NULL,&q,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    q    = (q-1)/(n-1);  /* number of spectral elements */
    h    = 2.0/q;
    ierr = DMDAGetCorners(dm,&xs,NULL,NULL,&xn,NULL,NULL);CHKERRQ(ierr);
    xs   = xs/(n-1);
    xn   = xn/(n-1);
    ierr = DMDASetUniformCoordinates(dm,-1.,1.,0.,0.,0.,0.);CHKERRQ(ierr);
    ierr = DMGetCoordinates(dm,&x);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(dm,x,&xx);CHKERRQ(ierr);

    /* loop over local spectral elements */
    for (j=xs; j<xs+xn; j++) {
      /*
       Except for the first process, each process starts on the second GLL point of the first element on that process
       */
      for (i= (j == xs && xs > 0)? 1 : 0; i<n; i++) {
        xx[j*(n-1) + i] = -1.0 + h*j + h*(nodes[i]+1.0)/2.;
      }
    }
    ierr = DMDAVecRestoreArray(dm,x,&xx);CHKERRQ(ierr);
  } else SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_SUP,"Not yet implemented for periodic");
  PetscFunctionReturn(0);
}

/*@

     DMDASetGLLCoordinates - Sets the global coordinates from -1 to 1 to the GLL points of as many GLL elements that fit the number of grid points

   Collective on da

   Input Parameters:
+   da - the DMDA object
-   n - the number of GLL nodes
-   nodes - the GLL nodes

   Notes:
    the parallel decomposition of grid points must correspond to the degree of the GLL. That is, the number of grid points
          on each process much be divisible by the number of GLL elements needed per process. This depends on whether the DM is
          periodic or not.

   Level: advanced

.seealso:   DMDACreate(), PetscDTGaussLobattoLegendreQuadrature(), DMGetCoordinates()
@*/
PetscErrorCode DMDASetGLLCoordinates(DM da,PetscInt n,PetscReal *nodes)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (da->dim == 1) {
    ierr = DMDASetGLLCoordinates_1d(da,n,nodes);CHKERRQ(ierr);
  } else SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_SUP,"Not yet implemented for 2 or 3d");
  PetscFunctionReturn(0);
}

PETSC_INTERN PetscErrorCode DMGetCompatibility_DA(DM da1,DM dm2,PetscBool *compatible,PetscBool *set)
{
  PetscErrorCode ierr;
  DM_DA          *dd1 = (DM_DA*)da1->data,*dd2;
  DM             da2;
  DMType         dmtype2;
  PetscBool      isda,compatibleLocal;
  PetscInt       i;

  PetscFunctionBegin;
  if (!da1->setupcalled) SETERRQ(PetscObjectComm((PetscObject)da1),PETSC_ERR_ARG_WRONGSTATE,"DMSetUp() must be called on first DM before DMGetCompatibility()");
  ierr = DMGetType(dm2,&dmtype2);CHKERRQ(ierr);
  ierr = PetscStrcmp(dmtype2,DMDA,&isda);CHKERRQ(ierr);
  if (isda) {
    da2 = dm2;
    dd2 = (DM_DA*)da2->data;
    if (!da2->setupcalled) SETERRQ(PetscObjectComm((PetscObject)da2),PETSC_ERR_ARG_WRONGSTATE,"DMSetUp() must be called on second DM before DMGetCompatibility()");
    compatibleLocal = (PetscBool)(da1->dim == da2->dim);
    if (compatibleLocal) compatibleLocal = (PetscBool)(compatibleLocal && (dd1->s == dd2->s)); /* Stencil width */
    /*                                                                           Global size              ranks               Boundary type */
    if (compatibleLocal)                 compatibleLocal = (PetscBool)(compatibleLocal && (dd1->M == dd2->M) && (dd1->m == dd2->m) && (dd1->bx == dd2->bx));
    if (compatibleLocal && da1->dim > 1) compatibleLocal = (PetscBool)(compatibleLocal && (dd1->N == dd2->N) && (dd1->n == dd2->n) && (dd1->by == dd2->by));
    if (compatibleLocal && da1->dim > 2) compatibleLocal = (PetscBool)(compatibleLocal && (dd1->P == dd2->P) && (dd1->p == dd2->p) && (dd1->bz == dd2->bz));
    if (compatibleLocal) {
      for (i=0; i<dd1->m; ++i) {
        compatibleLocal = (PetscBool)(compatibleLocal && (dd1->lx[i] == dd2->lx[i]));           /* Local size     */
      }
    }
    if (compatibleLocal && da1->dim > 1) {
      for (i=0; i<dd1->n; ++i) {
        compatibleLocal = (PetscBool)(compatibleLocal && (dd1->ly[i] == dd2->ly[i]));
      }
    }
    if (compatibleLocal && da1->dim > 2) {
      for (i=0; i<dd1->p; ++i) {
        compatibleLocal = (PetscBool)(compatibleLocal && (dd1->lz[i] == dd2->lz[i]));
      }
    }
    *compatible = compatibleLocal;
    *set = PETSC_TRUE;
  } else {
    /* Decline to determine compatibility with other DM types */
    *set = PETSC_FALSE;
  }
  PetscFunctionReturn(0);
}
