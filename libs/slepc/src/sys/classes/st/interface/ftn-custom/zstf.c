/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <petsc/private/fortranimpl.h>
#include <slepcst.h>

#if defined(PETSC_HAVE_FORTRAN_CAPS)
#define stsettype_                STSETTYPE
#define stgettype_                STGETTYPE
#define stsetoptionsprefix_       STSETOPTIONSPREFIX
#define stappendoptionsprefix_    STAPPENDOPTIONSPREFIX
#define stgetoptionsprefix_       STGETOPTIONSPREFIX
#define stview_                   STVIEW
#define stviewfromoptions_        STVIEWFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define stsettype_                stsettype
#define stgettype_                stgettype
#define stsetoptionsprefix_       stsetoptionsprefix
#define stappendoptionsprefix_    stappendoptionsprefix
#define stgetoptionsprefix_       stgetoptionsprefix
#define stview_                   stview
#define stviewfromoptions_        stviewfromoptions
#endif

SLEPC_EXTERN void stsettype_(ST *st,char *type,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(type,len,t);
  *ierr = STSetType(*st,t);if (*ierr) return;
  FREECHAR(type,t);
}

SLEPC_EXTERN void stgettype_(ST *st,char *name,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  STType tname;

  *ierr = STGetType(*st,&tname); if (*ierr) return;
  *ierr = PetscStrncpy(name,tname,len);if (*ierr) return;
  FIXRETURNCHAR(PETSC_TRUE,name,len);
}

SLEPC_EXTERN void stsetoptionsprefix_(ST *st,char *prefix,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(prefix,len,t);
  *ierr = STSetOptionsPrefix(*st,t);if (*ierr) return;
  FREECHAR(prefix,t);
}

SLEPC_EXTERN void stappendoptionsprefix_(ST *st,char *prefix,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(prefix,len,t);
  *ierr = STAppendOptionsPrefix(*st,t);if (*ierr) return;
  FREECHAR(prefix,t);
}

SLEPC_EXTERN void stgetoptionsprefix_(ST *st,char *prefix,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  const char *tname;

  *ierr = STGetOptionsPrefix(*st,&tname); if (*ierr) return;
  *ierr = PetscStrncpy(prefix,tname,len);if (*ierr) return;
  FIXRETURNCHAR(PETSC_TRUE,prefix,len);
}

SLEPC_EXTERN void stview_(ST *st,PetscViewer *viewer,PetscErrorCode *ierr)
{
  PetscViewer v;
  PetscPatchDefaultViewers_Fortran(viewer,v);
  *ierr = STView(*st,v);
}

SLEPC_EXTERN void stviewfromoptions_(ST *st,PetscObject obj,char* type,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(type,len,t);
  *ierr = STViewFromOptions(*st,obj,t);if (*ierr) return;
  FREECHAR(type,t);
}

