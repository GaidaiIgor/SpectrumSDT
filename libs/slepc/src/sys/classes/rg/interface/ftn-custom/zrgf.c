/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <petsc/private/fortranimpl.h>
#include <slepcrg.h>

#if defined(PETSC_HAVE_FORTRAN_CAPS)
#define rgsettype_                RGSETTYPE
#define rggettype_                RGGETTYPE
#define rgsetoptionsprefix_       RGSETOPTIONSPREFIX
#define rgappendoptionsprefix_    RGAPPENDOPTIONSPREFIX
#define rggetoptionsprefix_       RGGETOPTIONSPREFIX
#define rgview_                   RGVIEW
#define rgviewfromoptions_        RGVIEWFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define rgsettype_                rgsettype
#define rggettype_                rggettype
#define rgsetoptionsprefix_       rgsetoptionsprefix
#define rgappendoptionsprefix_    rgappendoptionsprefix
#define rggetoptionsprefix_       rggetoptionsprefix
#define rgview_                   rgview
#define rgviewfromoptions_        rgviewfromoptions
#endif

SLEPC_EXTERN void rgsettype_(RG *rg,char *type,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(type,len,t);
  *ierr = RGSetType(*rg,t);
  FREECHAR(type,t);
}

SLEPC_EXTERN void rggettype_(RG *rg,char *name,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  RGType tname;

  *ierr = RGGetType(*rg,&tname); if (*ierr) return;
  *ierr = PetscStrncpy(name,tname,len);if (*ierr) return;
  FIXRETURNCHAR(PETSC_TRUE,name,len);
}

SLEPC_EXTERN void rgsetoptionsprefix_(RG *rg,char *prefix,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(prefix,len,t);
  *ierr = RGSetOptionsPrefix(*rg,t);if (*ierr) return;
  FREECHAR(prefix,t);
}

SLEPC_EXTERN void rgappendoptionsprefix_(RG *rg,char *prefix,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(prefix,len,t);
  *ierr = RGAppendOptionsPrefix(*rg,t);if (*ierr) return;
  FREECHAR(prefix,t);
}

SLEPC_EXTERN void rggetoptionsprefix_(RG *rg,char *prefix,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  const char *tname;

  *ierr = RGGetOptionsPrefix(*rg,&tname); if (*ierr) return;
  *ierr = PetscStrncpy(prefix,tname,len);if (*ierr) return;
  FIXRETURNCHAR(PETSC_TRUE,prefix,len);
}

SLEPC_EXTERN void rgview_(RG *rg,PetscViewer *viewer,PetscErrorCode *ierr)
{
  PetscViewer v;
  PetscPatchDefaultViewers_Fortran(viewer,v);
  *ierr = RGView(*rg,v);
}

SLEPC_EXTERN void rgviewfromoptions_(RG *rg,PetscObject obj,char* type,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(type,len,t);
  *ierr = RGViewFromOptions(*rg,obj,t);if (*ierr) return;
  FREECHAR(type,t);
}

