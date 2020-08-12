/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Basic RG routines
*/

#include <slepc/private/rgimpl.h>      /*I "slepcrg.h" I*/

PetscFunctionList RGList = 0;
PetscBool         RGRegisterAllCalled = PETSC_FALSE;
PetscClassId      RG_CLASSID = 0;
static PetscBool  RGPackageInitialized = PETSC_FALSE;

/*@C
   RGFinalizePackage - This function destroys everything in the Slepc interface
   to the RG package. It is called from SlepcFinalize().

   Level: developer

.seealso: SlepcFinalize()
@*/
PetscErrorCode RGFinalizePackage(void)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFunctionListDestroy(&RGList);CHKERRQ(ierr);
  RGPackageInitialized = PETSC_FALSE;
  RGRegisterAllCalled  = PETSC_FALSE;
  PetscFunctionReturn(0);
}

/*@C
  RGInitializePackage - This function initializes everything in the RG package.
  It is called from PetscDLLibraryRegister() when using dynamic libraries, and
  on the first call to RGCreate() when using static libraries.

  Level: developer

.seealso: SlepcInitialize()
@*/
PetscErrorCode RGInitializePackage(void)
{
  char           logList[256];
  PetscBool      opt,pkg;
  PetscClassId   classids[1];
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (RGPackageInitialized) PetscFunctionReturn(0);
  RGPackageInitialized = PETSC_TRUE;
  /* Register Classes */
  ierr = PetscClassIdRegister("Region",&RG_CLASSID);CHKERRQ(ierr);
  /* Register Constructors */
  ierr = RGRegisterAll();CHKERRQ(ierr);
  /* Process Info */
  classids[0] = RG_CLASSID;
  ierr = PetscInfoProcessClass("rg",1,&classids[0]);CHKERRQ(ierr);
  /* Process summary exclusions */
  ierr = PetscOptionsGetString(NULL,NULL,"-log_exclude",logList,sizeof(logList),&opt);CHKERRQ(ierr);
  if (opt) {
    ierr = PetscStrInList("rg",logList,',',&pkg);CHKERRQ(ierr);
    if (pkg) { ierr = PetscLogEventDeactivateClass(RG_CLASSID);CHKERRQ(ierr); }
  }
  /* Register package finalizer */
  ierr = PetscRegisterFinalize(RGFinalizePackage);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   RGCreate - Creates an RG context.

   Collective

   Input Parameter:
.  comm - MPI communicator

   Output Parameter:
.  newrg - location to put the RG context

   Level: beginner

.seealso: RGDestroy(), RG
@*/
PetscErrorCode RGCreate(MPI_Comm comm,RG *newrg)
{
  RG             rg;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidPointer(newrg,2);
  *newrg = 0;
  ierr = RGInitializePackage();CHKERRQ(ierr);
  ierr = SlepcHeaderCreate(rg,RG_CLASSID,"RG","Region","RG",comm,RGDestroy,RGView);CHKERRQ(ierr);
  rg->complement = PETSC_FALSE;
  rg->sfactor    = 1.0;
  rg->osfactor   = 0.0;
  rg->data       = NULL;

  *newrg = rg;
  PetscFunctionReturn(0);
}

/*@C
   RGSetOptionsPrefix - Sets the prefix used for searching for all
   RG options in the database.

   Logically Collective on rg

   Input Parameters:
+  rg     - the region context
-  prefix - the prefix string to prepend to all RG option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the
   hyphen.

   Level: advanced

.seealso: RGAppendOptionsPrefix()
@*/
PetscErrorCode RGSetOptionsPrefix(RG rg,const char *prefix)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  ierr = PetscObjectSetOptionsPrefix((PetscObject)rg,prefix);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   RGAppendOptionsPrefix - Appends to the prefix used for searching for all
   RG options in the database.

   Logically Collective on rg

   Input Parameters:
+  rg     - the region context
-  prefix - the prefix string to prepend to all RG option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the hyphen.

   Level: advanced

.seealso: RGSetOptionsPrefix()
@*/
PetscErrorCode RGAppendOptionsPrefix(RG rg,const char *prefix)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  ierr = PetscObjectAppendOptionsPrefix((PetscObject)rg,prefix);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   RGGetOptionsPrefix - Gets the prefix used for searching for all
   RG options in the database.

   Not Collective

   Input Parameters:
.  rg - the region context

   Output Parameters:
.  prefix - pointer to the prefix string used is returned

   Note:
   On the Fortran side, the user should pass in a string 'prefix' of
   sufficient length to hold the prefix.

   Level: advanced

.seealso: RGSetOptionsPrefix(), RGAppendOptionsPrefix()
@*/
PetscErrorCode RGGetOptionsPrefix(RG rg,const char *prefix[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscValidPointer(prefix,2);
  ierr = PetscObjectGetOptionsPrefix((PetscObject)rg,prefix);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   RGSetType - Selects the type for the RG object.

   Logically Collective on rg

   Input Parameter:
+  rg   - the region context
-  type - a known type

   Level: intermediate

.seealso: RGGetType()
@*/
PetscErrorCode RGSetType(RG rg,RGType type)
{
  PetscErrorCode ierr,(*r)(RG);
  PetscBool      match;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscValidCharPointer(type,2);

  ierr = PetscObjectTypeCompare((PetscObject)rg,type,&match);CHKERRQ(ierr);
  if (match) PetscFunctionReturn(0);

  ierr =  PetscFunctionListFind(RGList,type,&r);CHKERRQ(ierr);
  if (!r) SETERRQ1(PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_UNKNOWN_TYPE,"Unable to find requested RG type %s",type);

  if (rg->ops->destroy) { ierr = (*rg->ops->destroy)(rg);CHKERRQ(ierr); }
  ierr = PetscMemzero(rg->ops,sizeof(struct _RGOps));CHKERRQ(ierr);

  ierr = PetscObjectChangeTypeName((PetscObject)rg,type);CHKERRQ(ierr);
  ierr = (*r)(rg);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   RGGetType - Gets the RG type name (as a string) from the RG context.

   Not Collective

   Input Parameter:
.  rg - the region context

   Output Parameter:
.  name - name of the region

   Level: intermediate

.seealso: RGSetType()
@*/
PetscErrorCode RGGetType(RG rg,RGType *type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscValidPointer(type,2);
  *type = ((PetscObject)rg)->type_name;
  PetscFunctionReturn(0);
}

/*@
   RGSetFromOptions - Sets RG options from the options database.

   Collective on rg

   Input Parameters:
.  rg - the region context

   Notes:
   To see all options, run your program with the -help option.

   Level: beginner
@*/
PetscErrorCode RGSetFromOptions(RG rg)
{
  PetscErrorCode ierr;
  char           type[256];
  PetscBool      flg;
  PetscReal      sfactor;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  ierr = RGRegisterAll();CHKERRQ(ierr);
  ierr = PetscObjectOptionsBegin((PetscObject)rg);CHKERRQ(ierr);
    ierr = PetscOptionsFList("-rg_type","Region type","RGSetType",RGList,(char*)(((PetscObject)rg)->type_name?((PetscObject)rg)->type_name:RGINTERVAL),type,256,&flg);CHKERRQ(ierr);
    if (flg) {
      ierr = RGSetType(rg,type);CHKERRQ(ierr);
    } else if (!((PetscObject)rg)->type_name) {
      ierr = RGSetType(rg,RGINTERVAL);CHKERRQ(ierr);
    }

    ierr = PetscOptionsBool("-rg_complement","Whether region is complemented or not","RGSetComplement",rg->complement,&rg->complement,NULL);CHKERRQ(ierr);

    ierr = PetscOptionsReal("-rg_scale","Scaling factor","RGSetScale",1.0,&sfactor,&flg);CHKERRQ(ierr);
    if (flg) { ierr = RGSetScale(rg,sfactor);CHKERRQ(ierr); }

    if (rg->ops->setfromoptions) {
      ierr = (*rg->ops->setfromoptions)(PetscOptionsObject,rg);CHKERRQ(ierr);
    }
    ierr = PetscObjectProcessOptionsHandlers(PetscOptionsObject,(PetscObject)rg);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   RGView - Prints the RG data structure.

   Collective on rg

   Input Parameters:
+  rg - the region context
-  viewer - optional visualization context

   Note:
   The available visualization contexts include
+     PETSC_VIEWER_STDOUT_SELF - standard output (default)
-     PETSC_VIEWER_STDOUT_WORLD - synchronized standard
         output where only the first processor opens
         the file.  All other processors send their
         data to the first processor to print.

   The user can open an alternative visualization context with
   PetscViewerASCIIOpen() - output to a specified file.

   Level: beginner
@*/
PetscErrorCode RGView(RG rg,PetscViewer viewer)
{
  PetscBool      isdraw,isascii;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)rg),&viewer);CHKERRQ(ierr);
  }
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(rg,1,viewer,2);
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERDRAW,&isdraw);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    ierr = PetscObjectPrintClassNamePrefixType((PetscObject)rg,viewer);CHKERRQ(ierr);
    if (rg->ops->view) {
      ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
      ierr = (*rg->ops->view)(rg,viewer);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
    }
    if (rg->complement) {
      ierr = PetscViewerASCIIPrintf(viewer,"  selected region is the complement of the specified one\n");CHKERRQ(ierr);
    }
    if (rg->sfactor!=1.0) {
      ierr = PetscViewerASCIIPrintf(viewer,"  scaling factor = %g\n",(double)rg->sfactor);CHKERRQ(ierr);
    }
  } else if (isdraw) {
    if (rg->ops->view) { ierr = (*rg->ops->view)(rg,viewer);CHKERRQ(ierr); }
  }
  PetscFunctionReturn(0);
}

/*@C
   RGViewFromOptions - View from options

   Collective on RG

   Input Parameters:
+  rg   - the region context
.  obj  - optional object
-  name - command line option

   Level: intermediate

.seealso: RGView(), RGCreate()
@*/
PetscErrorCode RGViewFromOptions(RG rg,PetscObject obj,const char name[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  ierr = PetscObjectViewFromOptions((PetscObject)rg,obj,name);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   RGIsTrivial - Whether it is the trivial region (whole complex plane).

   Not Collective

   Input Parameter:
.  rg - the region context

   Output Parameter:
.  trivial - true if the region is equal to the whole complex plane, e.g.,
             an interval region with all four endpoints unbounded or an
             ellipse with infinite radius.

   Level: beginner
@*/
PetscErrorCode RGIsTrivial(RG rg,PetscBool *trivial)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscValidType(rg,1);
  PetscValidBoolPointer(trivial,2);
  if (rg->ops->istrivial) {
    ierr = (*rg->ops->istrivial)(rg,trivial);CHKERRQ(ierr);
  } else *trivial = PETSC_FALSE;
  PetscFunctionReturn(0);
}

/*@
   RGCheckInside - Determines if a set of given points are inside the region or not.

   Not Collective

   Input Parameters:
+  rg - the region context
.  n  - number of points to check
.  ar - array of real parts
-  ai - array of imaginary parts

   Output Parameter:
.  inside - array of results (1=inside, 0=on the contour, -1=outside)

   Note:
   The point a is expressed as a couple of PetscScalar variables ar,ai.
   If built with complex scalars, the point is supposed to be stored in ar,
   otherwise ar,ai contain the real and imaginary parts, respectively.

   If a scaling factor was set, the points are scaled before checking.

   Level: intermediate

.seealso: RGSetScale(), RGSetComplement()
@*/
PetscErrorCode RGCheckInside(RG rg,PetscInt n,PetscScalar *ar,PetscScalar *ai,PetscInt *inside)
{
  PetscErrorCode ierr;
  PetscReal      px,py;
  PetscInt       i;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscValidType(rg,1);
  PetscValidScalarPointer(ar,3);
#if !defined(PETSC_USE_COMPLEX)
  PetscValidScalarPointer(ai,4);
#endif
  PetscValidIntPointer(inside,5);

  for (i=0;i<n;i++) {
#if defined(PETSC_USE_COMPLEX)
    px = PetscRealPart(ar[i]);
    py = PetscImaginaryPart(ar[i]);
#else
    px = ar[i];
    py = ai[i];
#endif
    if (rg->sfactor != 1.0) {
      px /= rg->sfactor;
      py /= rg->sfactor;
    }
    ierr = (*rg->ops->checkinside)(rg,px,py,inside+i);CHKERRQ(ierr);
    if (rg->complement) inside[i] = -inside[i];
  }
  PetscFunctionReturn(0);
}

/*@
   RGComputeContour - Computes the coordinates of several points lying on the
   contour of the region.

   Not Collective

   Input Parameters:
+  rg - the region context
-  n  - number of points to compute

   Output Parameters:
+  cr - location to store real parts
-  ci - location to store imaginary parts

   Level: intermediate
@*/
PetscErrorCode RGComputeContour(RG rg,PetscInt n,PetscScalar cr[],PetscScalar ci[])
{
  PetscInt       i;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscValidType(rg,1);
  PetscValidScalarPointer(cr,3);
#if !defined(PETSC_USE_COMPLEX)
  PetscValidScalarPointer(ci,4);
#endif
  if (rg->complement) SETERRQ(PetscObjectComm((PetscObject)rg),PETSC_ERR_SUP,"Cannot compute contour of region with complement flag set");
  ierr = (*rg->ops->computecontour)(rg,n,cr,ci);CHKERRQ(ierr);
  for (i=0;i<n;i++) {
    cr[i] *= rg->sfactor;
    ci[i] *= rg->sfactor;
  }
  PetscFunctionReturn(0);
}

/*@
   RGComputeBoundingBox - Determines the endpoints of a rectangle in the complex plane that
   contains the region.

   Not Collective

   Input Parameters:
.  rg - the region context

   Output Parameters:
+  a,b - endpoints of the bounding box in the real axis
-  c,d - endpoints of the bounding box in the imaginary axis

   Notes:
   The bounding box is defined as [a,b]x[c,d]. In regions that are not bounded (e.g. an
   open interval) or with the complement flag set, it makes no sense to compute a bounding
   box, so the return values are infinite.

   Level: intermediate

.seealso: RGSetComplement()
@*/
PetscErrorCode RGComputeBoundingBox(RG rg,PetscReal *a,PetscReal *b,PetscReal *c,PetscReal *d)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscValidType(rg,1);

  if (rg->complement) {  /* cannot compute bounding box */
    if (a) *a = -PETSC_MAX_REAL;
    if (b) *b =  PETSC_MAX_REAL;
    if (c) *c = -PETSC_MAX_REAL;
    if (d) *d =  PETSC_MAX_REAL;
  } else {
    ierr = (*rg->ops->computebbox)(rg,a,b,c,d);CHKERRQ(ierr);
    if (a && *a!=-PETSC_MAX_REAL) *a *= rg->sfactor;
    if (b && *b!= PETSC_MAX_REAL) *b *= rg->sfactor;
    if (c && *c!=-PETSC_MAX_REAL) *c *= rg->sfactor;
    if (d && *d!= PETSC_MAX_REAL) *d *= rg->sfactor;
  }
  PetscFunctionReturn(0);
}

/*@
   RGSetComplement - Sets a flag to indicate that the region is the complement
   of the specified one.

   Logically Collective on rg

   Input Parameters:
+  rg  - the region context
-  flg - the boolean flag

   Options Database Key:
.  -rg_complement <bool> - Activate/deactivate the complementation of the region

   Level: intermediate

.seealso: RGGetComplement()
@*/
PetscErrorCode RGSetComplement(RG rg,PetscBool flg)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscValidLogicalCollectiveBool(rg,flg,2);
  rg->complement = flg;
  PetscFunctionReturn(0);
}

/*@
   RGGetComplement - Gets a flag that that indicates whether the region
   is complemented or not.

   Not Collective

   Input Parameter:
.  rg - the region context

   Output Parameter:
.  flg - the flag

   Level: intermediate

.seealso: RGSetComplement()
@*/
PetscErrorCode RGGetComplement(RG rg,PetscBool *flg)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscValidBoolPointer(flg,2);
  *flg = rg->complement;
  PetscFunctionReturn(0);
}

/*@
   RGSetScale - Sets the scaling factor to be used when checking that a
   point is inside the region and when computing the contour.

   Logically Collective on rg

   Input Parameters:
+  rg      - the region context
-  sfactor - the scaling factor

   Options Database Key:
.  -rg_scale <real> - Sets the scaling factor

   Level: advanced

.seealso: RGGetScale(), RGCheckInside()
@*/
PetscErrorCode RGSetScale(RG rg,PetscReal sfactor)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscValidLogicalCollectiveReal(rg,sfactor,2);
  if (sfactor == PETSC_DEFAULT || sfactor == PETSC_DECIDE) rg->sfactor = 1.0;
  else {
    if (sfactor<=0.0) SETERRQ(PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of scaling factor. Must be > 0");
    rg->sfactor = sfactor;
  }
  PetscFunctionReturn(0);
}

/*@
   RGGetScale - Gets the scaling factor.

   Not Collective

   Input Parameter:
.  rg - the region context

   Output Parameter:
.  flg - the flag

   Level: advanced

.seealso: RGSetScale()
@*/
PetscErrorCode RGGetScale(RG rg,PetscReal *sfactor)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscValidRealPointer(sfactor,2);
  *sfactor = rg->sfactor;
  PetscFunctionReturn(0);
}

/*@
   RGPushScale - Sets an additional scaling factor, that will multiply the
   user-defined scaling factor.

   Logically Collective on rg

   Input Parameters:
+  rg      - the region context
-  sfactor - the scaling factor

   Notes:
   The current implementation does not allow pushing several scaling factors.

   This is intended for internal use, for instance in polynomial eigensolvers
   that use parameter scaling.

   Level: developer

.seealso: RGPopScale(), RGSetScale()
@*/
PetscErrorCode RGPushScale(RG rg,PetscReal sfactor)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscValidLogicalCollectiveReal(rg,sfactor,2);
  if (sfactor<=0.0) SETERRQ(PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of scaling factor. Must be > 0");
  if (rg->osfactor) SETERRQ(PetscObjectComm((PetscObject)rg),PETSC_ERR_SUP,"Current implementation does not allow pushing several scaling factors");
  rg->osfactor = rg->sfactor;
  rg->sfactor *= sfactor;
  PetscFunctionReturn(0);
}

/*@
   RGPopScale - Pops the scaling factor set with RGPushScale().

   Not Collective

   Input Parameter:
.  rg - the region context

   Level: developer

.seealso: RGPushScale()
@*/
PetscErrorCode RGPopScale(RG rg)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  if (!rg->osfactor) SETERRQ(PetscObjectComm((PetscObject)rg),PETSC_ERR_ORDER,"Must call RGPushScale first");
  rg->sfactor  = rg->osfactor;
  rg->osfactor = 0.0;
  PetscFunctionReturn(0);
}

/*@
   RGDestroy - Destroys RG context that was created with RGCreate().

   Collective on rg

   Input Parameter:
.  rg - the region context

   Level: beginner

.seealso: RGCreate()
@*/
PetscErrorCode RGDestroy(RG *rg)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!*rg) PetscFunctionReturn(0);
  PetscValidHeaderSpecific(*rg,RG_CLASSID,1);
  if (--((PetscObject)(*rg))->refct > 0) { *rg = 0; PetscFunctionReturn(0); }
  if ((*rg)->ops->destroy) { ierr = (*(*rg)->ops->destroy)(*rg);CHKERRQ(ierr); }
  ierr = PetscHeaderDestroy(rg);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   RGRegister - Adds a region to the RG package.

   Not collective

   Input Parameters:
+  name - name of a new user-defined RG
-  function - routine to create context

   Notes:
   RGRegister() may be called multiple times to add several user-defined regions.

   Level: advanced

.seealso: RGRegisterAll()
@*/
PetscErrorCode RGRegister(const char *name,PetscErrorCode (*function)(RG))
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = RGInitializePackage();CHKERRQ(ierr);
  ierr = PetscFunctionListAdd(&RGList,name,function);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

