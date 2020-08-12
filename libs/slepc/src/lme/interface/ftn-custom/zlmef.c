/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <petsc/private/fortranimpl.h>
#include <slepc/private/lmeimpl.h>

#if defined(PETSC_HAVE_FORTRAN_CAPS)
#define lmeview_                          LMEVIEW
#define lmeviewfromoptions_               LMEVIEWFROMOPTIONS
#define lmereasonview_                    LMEREASONVIEW
#define lmesetoptionsprefix_              LMESETOPTIONSPREFIX
#define lmeappendoptionsprefix_           LMEAPPENDOPTIONSPREFIX
#define lmegetoptionsprefix_              LMEGETOPTIONSPREFIX
#define lmesettype_                       LMESETTYPE
#define lmegettype_                       LMEGETTYPE
#define lmemonitordefault_                LMEMONITORDEFAULT
#define lmemonitorlg_                     LMEMONITORLG
#define lmemonitorset_                    LMEMONITORSET
#define lmegettolerances00_               LMEGETTOLERANCES00
#define lmegettolerances10_               LMEGETTOLERANCES10
#define lmegettolerances01_               LMEGETTOLERANCES01
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define lmeview_                          lmeview
#define lmeviewfromoptions_               lmeviewfromoptions
#define lmereasonview_                    lmereasonview
#define lmesetoptionsprefix_              lmesetoptionsprefix
#define lmeappendoptionsprefix_           lmeappendoptionsprefix
#define lmegetoptionsprefix_              lmegetoptionsprefix
#define lmesettype_                       lmesettype
#define lmegettype_                       lmegettype
#define lmemonitordefault_                lmemonitordefault
#define lmemonitorlg_                     lmemonitorlg
#define lmemonitorset_                    lmemonitorset
#define lmegettolerances00_               lmegettolerances00
#define lmegettolerances10_               lmegettolerances10
#define lmegettolerances01_               lmegettolerances01
#endif

/*
   These are not usually called from Fortran but allow Fortran users
   to transparently set these monitors from .F code
*/
SLEPC_EXTERN void lmemonitordefault_(LME *lme,PetscInt *it,PetscReal *errest,PetscViewerAndFormat **ctx,PetscErrorCode *ierr)
{
  *ierr = LMEMonitorDefault(*lme,*it,*errest,*ctx);
}

SLEPC_EXTERN void lmemonitorlg_(LME *lme,PetscInt *it,PetscReal *errest,void *ctx,PetscErrorCode *ierr)
{
  *ierr = LMEMonitorLG(*lme,*it,*errest,ctx);
}

static struct {
  PetscFortranCallbackId monitor;
  PetscFortranCallbackId monitordestroy;
} _cb;

/* These are not extern C because they are passed into non-extern C user level functions */
static PetscErrorCode ourmonitor(LME lme,PetscInt i,PetscReal d,void* ctx)
{
  PetscObjectUseFortranCallback(lme,_cb.monitor,(LME*,PetscInt*,PetscReal*,void*,PetscErrorCode*),(&lme,&i,&d,_ctx,&ierr));
}

static PetscErrorCode ourdestroy(void** ctx)
{
  LME lme = (LME)*ctx;
  PetscObjectUseFortranCallback(lme,_cb.monitordestroy,(void*,PetscErrorCode*),(_ctx,&ierr));
}

SLEPC_EXTERN void lmeview_(LME *lme,PetscViewer *viewer,PetscErrorCode *ierr)
{
  PetscViewer v;
  PetscPatchDefaultViewers_Fortran(viewer,v);
  *ierr = LMEView(*lme,v);
}

SLEPC_EXTERN void lmeviewfromoptions_(LME *lme,PetscObject obj,char* type,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(type,len,t);
  *ierr = LMEViewFromOptions(*lme,obj,t);if (*ierr) return;
  FREECHAR(type,t);
}

SLEPC_EXTERN void lmereasonview_(LME *lme,PetscViewer *viewer,PetscErrorCode *ierr)
{
  PetscViewer v;
  PetscPatchDefaultViewers_Fortran(viewer,v);
  *ierr = LMEReasonView(*lme,v);
}

SLEPC_EXTERN void lmesettype_(LME *lme,char *type,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(type,len,t);
  *ierr = LMESetType(*lme,t);if (*ierr) return;
  FREECHAR(type,t);
}

SLEPC_EXTERN void lmegettype_(LME *lme,char *name,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  LMEType tname;

  *ierr = LMEGetType(*lme,&tname);if (*ierr) return;
  *ierr = PetscStrncpy(name,tname,len);if (*ierr) return;
  FIXRETURNCHAR(PETSC_TRUE,name,len);
}

SLEPC_EXTERN void lmesetoptionsprefix_(LME *lme,char *prefix,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(prefix,len,t);
  *ierr = LMESetOptionsPrefix(*lme,t);if (*ierr) return;
  FREECHAR(prefix,t);
}

SLEPC_EXTERN void lmeappendoptionsprefix_(LME *lme,char *prefix,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(prefix,len,t);
  *ierr = LMEAppendOptionsPrefix(*lme,t);if (*ierr) return;
  FREECHAR(prefix,t);
}

SLEPC_EXTERN void lmegetoptionsprefix_(LME *lme,char *prefix,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  const char *tname;

  *ierr = LMEGetOptionsPrefix(*lme,&tname); if (*ierr) return;
  *ierr = PetscStrncpy(prefix,tname,len);if (*ierr) return;
  FIXRETURNCHAR(PETSC_TRUE,prefix,len);
}

SLEPC_EXTERN void lmemonitorset_(LME *lme,void (*monitor)(LME*,PetscInt*,PetscReal*,void*,PetscErrorCode*),void *mctx,void (*monitordestroy)(void *,PetscErrorCode*),PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(mctx);
  CHKFORTRANNULLFUNCTION(monitordestroy);
  if ((PetscVoidFunction)monitor == (PetscVoidFunction)lmemonitordefault_) {
    *ierr = LMEMonitorSet(*lme,(PetscErrorCode (*)(LME,PetscInt,PetscReal,void*))LMEMonitorDefault,*(PetscViewerAndFormat**)mctx,(PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy);
  } else if ((PetscVoidFunction)monitor == (PetscVoidFunction)lmemonitorlg_) {
    *ierr = LMEMonitorSet(*lme,LMEMonitorLG,0,0);
  } else {
    *ierr = PetscObjectSetFortranCallback((PetscObject)*lme,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.monitor,(PetscVoidFunction)monitor,mctx); if (*ierr) return;
    *ierr = PetscObjectSetFortranCallback((PetscObject)*lme,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.monitordestroy,(PetscVoidFunction)monitordestroy,mctx); if (*ierr) return;
    *ierr = LMEMonitorSet(*lme,ourmonitor,*lme,ourdestroy);
  }
}

SLEPC_EXTERN void lmegettolerances_(LME *lme,PetscReal *tol,PetscInt *maxits,PetscErrorCode *ierr)
{
  CHKFORTRANNULLREAL(tol);
  CHKFORTRANNULLINTEGER(maxits);
  *ierr = LMEGetTolerances(*lme,tol,maxits);
}

SLEPC_EXTERN void lmegettolerances00_(LME *lme,PetscReal *tol,PetscInt *maxits,PetscErrorCode *ierr)
{
  lmegettolerances_(lme,tol,maxits,ierr);
}

SLEPC_EXTERN void lmegettolerances10_(LME *lme,PetscReal *tol,PetscInt *maxits,PetscErrorCode *ierr)
{
  lmegettolerances_(lme,tol,maxits,ierr);
}

SLEPC_EXTERN void lmegettolerances01_(LME *lme,PetscReal *tol,PetscInt *maxits,PetscErrorCode *ierr)
{
  lmegettolerances_(lme,tol,maxits,ierr);
}

