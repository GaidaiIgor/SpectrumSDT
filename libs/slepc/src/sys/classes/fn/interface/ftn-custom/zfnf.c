/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <petsc/private/fortranimpl.h>
#include <slepc/private/slepcimpl.h>
#include <slepc/private/fnimpl.h>

#if defined(PETSC_HAVE_FORTRAN_CAPS)
#define fnview_                    FNVIEW
#define fnviewfromoptions_         FNVIEWFROMOPTIONS
#define fnsetoptionsprefix_        FNSETOPTIONSPREFIX
#define fnappendoptionsprefix_     FNAPPENDOPTIONSPREFIX
#define fngetoptionsprefix_        FNGETOPTIONSPREFIX
#define fnsettype_                 FNSETTYPE
#define fngettype_                 FNGETTYPE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define fnview_                    fnview
#define fnviewfromoptions_         fnviewfromoptions
#define fnsetoptionsprefix_        fnsetoptionsprefix
#define fnappendoptionsprefix_     fnappendoptionsprefix
#define fngetoptionsprefix_        fngetoptionsprefix
#define fnsettype_                 fnsettype
#define fngettype_                 fngettype
#endif

SLEPC_EXTERN void fnview_(FN *fn,PetscViewer *viewer,PetscErrorCode *ierr)
{
  PetscViewer v;
  PetscPatchDefaultViewers_Fortran(viewer,v);
  *ierr = FNView(*fn,v);
}

SLEPC_EXTERN void fnviewfromoptions_(FN *fn,PetscObject obj,char* type,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(type,len,t);
  *ierr = FNViewFromOptions(*fn,obj,t);if (*ierr) return;
  FREECHAR(type,t);
}

SLEPC_EXTERN void fnsetoptionsprefix_(FN *fn,char *prefix,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(prefix,len,t);
  *ierr = FNSetOptionsPrefix(*fn,t);if (*ierr) return;
  FREECHAR(prefix,t);
}

SLEPC_EXTERN void fnappendoptionsprefix_(FN *fn,char *prefix,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(prefix,len,t);
  *ierr = FNAppendOptionsPrefix(*fn,t);if (*ierr) return;
  FREECHAR(prefix,t);
}

SLEPC_EXTERN void fngetoptionsprefix_(FN *fn,char *prefix,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  const char *tname;

  *ierr = FNGetOptionsPrefix(*fn,&tname); if (*ierr) return;
  *ierr = PetscStrncpy(prefix,tname,len);if (*ierr) return;
  FIXRETURNCHAR(PETSC_TRUE,prefix,len);
}

SLEPC_EXTERN void fnsettype_(FN *fn,char *type,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(type,len,t);
  *ierr = FNSetType(*fn,t);if (*ierr) return;
  FREECHAR(type,t);
}

SLEPC_EXTERN void fngettype_(FN *fn,char *name,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  FNType tname;

  *ierr = FNGetType(*fn,&tname);if (*ierr) return;
  *ierr = PetscStrncpy(name,tname,len);if (*ierr) return;
  FIXRETURNCHAR(PETSC_TRUE,name,len);
}


