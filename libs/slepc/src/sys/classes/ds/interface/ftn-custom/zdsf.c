/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <petsc/private/fortranimpl.h>
#include <slepcds.h>

#if defined(PETSC_HAVE_FORTRAN_CAPS)
#define dssettype_                DSSETTYPE
#define dsgettype_                DSGETTYPE
#define dssetoptionsprefix_       DSSETOPTIONSPREFIX
#define dsappendoptionsprefix_    DSAPPENDOPTIONSPREFIX
#define dsgetoptionsprefix_       DSGETOPTIONSPREFIX
#define dsview_                   DSVIEW
#define dsviewfromoptions_        DSVIEWFROMOPTIONS
#define dsviewmat_                DSVIEWMAT
#define dsvectors_                DSVECTORS
#define dssort_                   DSSORT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define dssettype_                dssettype
#define dsgettype_                dsgettype
#define dssetoptionsprefix_       dssetoptionsprefix
#define dsappendoptionsprefix_    dsappendoptionsprefix
#define dsgetoptionsprefix_       dsgetoptionsprefix
#define dsview_                   dsview
#define dsviewfromoptions_        dsviewfromoptions
#define dsviewmat_                dsviewmat
#define dsvectors_                dsvectors
#define dssort_                   dssort
#endif

SLEPC_EXTERN void dssettype_(DS *ds,char *type,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(type,len,t);
  *ierr = DSSetType(*ds,t);if (*ierr) return;
  FREECHAR(type,t);
}

SLEPC_EXTERN void dsgettype_(DS *ds,char *name,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  DSType tname;

  *ierr = DSGetType(*ds,&tname);if (*ierr) return;
  *ierr = PetscStrncpy(name,tname,len);if (*ierr) return;
  FIXRETURNCHAR(PETSC_TRUE,name,len);
}

SLEPC_EXTERN void dssetoptionsprefix_(DS *ds,char *prefix,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(prefix,len,t);
  *ierr = DSSetOptionsPrefix(*ds,t);if (*ierr) return;
  FREECHAR(prefix,t);
}

SLEPC_EXTERN void dsappendoptionsprefix_(DS *ds,char *prefix,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(prefix,len,t);
  *ierr = DSAppendOptionsPrefix(*ds,t);if (*ierr) return;
  FREECHAR(prefix,t);
}

SLEPC_EXTERN void dsgetoptionsprefix_(DS *ds,char *prefix,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  const char *tname;

  *ierr = DSGetOptionsPrefix(*ds,&tname); if (*ierr) return;
  *ierr = PetscStrncpy(prefix,tname,len);if (*ierr) return;
  FIXRETURNCHAR(PETSC_TRUE,prefix,len);
}

SLEPC_EXTERN void dsview_(DS *ds,PetscViewer *viewer,PetscErrorCode *ierr)
{
  PetscViewer v;
  PetscPatchDefaultViewers_Fortran(viewer,v);
  *ierr = DSView(*ds,v);
}

SLEPC_EXTERN void dsviewfromoptions_(DS *ds,PetscObject obj,char* type,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(type,len,t);
  *ierr = DSViewFromOptions(*ds,obj,t);if (*ierr) return;
  FREECHAR(type,t);
}

SLEPC_EXTERN void dsviewmat_(DS *ds,PetscViewer *viewer,DSMatType *m,PetscErrorCode *ierr)
{
  PetscViewer v;
  PetscPatchDefaultViewers_Fortran(viewer,v);
  *ierr = DSViewMat(*ds,v,*m);
}

SLEPC_EXTERN void dsvectors_(DS *ds,DSMatType *mat,PetscInt *j,PetscReal *rnorm,PetscErrorCode *ierr)
{
  CHKFORTRANNULLINTEGER(j);
  CHKFORTRANNULLREAL(rnorm);
  *ierr = DSVectors(*ds,*mat,j,rnorm);
}

SLEPC_EXTERN void dssort_(DS *ds,PetscScalar *eigr,PetscScalar *eigi,PetscScalar *rr,PetscScalar *ri,PetscInt *k,PetscErrorCode *ierr)
{
  CHKFORTRANNULLSCALAR(eigr);
  CHKFORTRANNULLSCALAR(eigi);
  CHKFORTRANNULLSCALAR(rr);
  CHKFORTRANNULLSCALAR(ri);
  CHKFORTRANNULLINTEGER(k);
  *ierr = DSSort(*ds,eigr,eigi,rr,ri,k);
}

