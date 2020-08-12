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
#include <slepc/private/epsimpl.h>

#if defined(PETSC_HAVE_FORTRAN_CAPS)
#define epsview_                          EPSVIEW
#define epsviewfromoptions_               EPSVIEWFROMOPTIONS
#define epserrorview_                     EPSERRORVIEW
#define epsreasonview_                    EPSREASONVIEW
#define epsvaluesview_                    EPSVALUESVIEW
#define epsvectorsview_                   EPSVECTORSVIEW
#define epssetoptionsprefix_              EPSSETOPTIONSPREFIX
#define epsappendoptionsprefix_           EPSAPPENDOPTIONSPREFIX
#define epsgetoptionsprefix_              EPSGETOPTIONSPREFIX
#define epssettype_                       EPSSETTYPE
#define epsgettype_                       EPSGETTYPE
#define epsmonitorall_                    EPSMONITORALL
#define epsmonitorlg_                     EPSMONITORLG
#define epsmonitorlgall_                  EPSMONITORLGALL
#define epsmonitorset_                    EPSMONITORSET
#define epsmonitorconverged_              EPSMONITORCONVERGED
#define epsmonitorfirst_                  EPSMONITORFIRST
#define epsconvergedabsolute_             EPSCONVERGEDABSOLUTE
#define epsconvergedrelative_             EPSCONVERGEDRELATIVE
#define epsconvergednorm_                 EPSCONVERGEDNORM
#define epssetconvergencetestfunction_    EPSSETCONVERGENCETESTFUNCTION
#define epssetstoppingtestfunction_       EPSSETSTOPPINGTESTFUNCTION
#define epsseteigenvaluecomparison_       EPSSETEIGENVALUECOMPARISON
#define epssetarbitraryselection_         EPSSETARBITRARYSELECTION
#define epsgetdimensions000_              EPSGETDIMENSIONS000
#define epsgetdimensions100_              EPSGETDIMENSIONS100
#define epsgetdimensions010_              EPSGETDIMENSIONS010
#define epsgetdimensions001_              EPSGETDIMENSIONS001
#define epsgetdimensions110_              EPSGETDIMENSIONS110
#define epsgetdimensions011_              EPSGETDIMENSIONS011
#define epsgetdimensions101_              EPSGETDIMENSIONS101
#define epsgeteigenpair00_                EPSGETEIGENPAIR00
#define epsgeteigenpair10_                EPSGETEIGENPAIR10
#define epsgeteigenpair01_                EPSGETEIGENPAIR01
#define epsgeteigenpair11_                EPSGETEIGENPAIR11
#define epsgeteigenvalue00_               EPSGETEIGENVALUE00
#define epsgeteigenvalue10_               EPSGETEIGENVALUE10
#define epsgeteigenvalue01_               EPSGETEIGENVALUE01
#define epsgettolerances00_               EPSGETTOLERANCES00
#define epsgettolerances10_               EPSGETTOLERANCES10
#define epsgettolerances01_               EPSGETTOLERANCES01
#define epsgetbalance000_                 EPSGETBALANCE000
#define epsgetbalance100_                 EPSGETBALANCE100
#define epsgetbalance010_                 EPSGETBALANCE010
#define epsgetbalance001_                 EPSGETBALANCE001
#define epsgetbalance110_                 EPSGETBALANCE110
#define epsgetbalance011_                 EPSGETBALANCE011
#define epsgetbalance101_                 EPSGETBALANCE101
#define epssetdeflationspace0_            EPSSETDEFLATIONSPACE0
#define epssetdeflationspace1_            EPSSETDEFLATIONSPACE1
#define epssetinitialspace0_              EPSSETINITIALSPACE0
#define epssetinitialspace1_              EPSSETINITIALSPACE1
#define epssetleftinitialspace0_          EPSSETLEFTINITIALSPACE0
#define epssetleftinitialspace1_          EPSSETLEFTINITIALSPACE1
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define epsview_                          epsview
#define epsviewfromoptions_               epsviewfromoptions
#define epserrorview_                     epserrorview
#define epsreasonview_                    epsreasonview
#define epsvaluesview_                    epsvaluesview
#define epsvectorsview_                   epsvectorsview
#define epssetoptionsprefix_              epssetoptionsprefix
#define epsappendoptionsprefix_           epsappendoptionsprefix
#define epsgetoptionsprefix_              epsgetoptionsprefix
#define epssettype_                       epssettype
#define epsgettype_                       epsgettype
#define epsmonitorall_                    epsmonitorall
#define epsmonitorlg_                     epsmonitorlg
#define epsmonitorlgall_                  epsmonitorlgall
#define epsmonitorset_                    epsmonitorset
#define epsmonitorconverged_              epsmonitorconverged
#define epsmonitorfirst_                  epsmonitorfirst
#define epsconvergedabsolute_             epsconvergedabsolute
#define epsconvergedrelative_             epsconvergedrelative
#define epsconvergednorm_                 epsconvergednorm
#define epssetconvergencetestfunction_    epssetconvergencetestfunction
#define epssetstoppingtestfunction_       epssetstoppingtestfunction
#define epsseteigenvaluecomparison_       epsseteigenvaluecomparison
#define epssetarbitraryselection_         epssetarbitraryselection
#define epsgetdimensions000_              epsgetdimensions000
#define epsgetdimensions100_              epsgetdimensions100
#define epsgetdimensions010_              epsgetdimensions010
#define epsgetdimensions001_              epsgetdimensions001
#define epsgetdimensions110_              epsgetdimensions110
#define epsgetdimensions011_              epsgetdimensions011
#define epsgetdimensions101_              epsgetdimensions101
#define epsgeteigenpair00_                epsgeteigenpair00
#define epsgeteigenpair10_                epsgeteigenpair10
#define epsgeteigenpair01_                epsgeteigenpair01
#define epsgeteigenpair11_                epsgeteigenpair11
#define epsgeteigenvalue00_               epsgeteigenvalue00
#define epsgeteigenvalue10_               epsgeteigenvalue10
#define epsgeteigenvalue01_               epsgeteigenvalue01
#define epsgettolerances00_               epsgettolerances00
#define epsgettolerances10_               epsgettolerances10
#define epsgettolerances01_               epsgettolerances01
#define epsgetbalance000_                 epsgetbalance000
#define epsgetbalance100_                 epsgetbalance100
#define epsgetbalance010_                 epsgetbalance010
#define epsgetbalance001_                 epsgetbalance001
#define epsgetbalance110_                 epsgetbalance110
#define epsgetbalance011_                 epsgetbalance011
#define epsgetbalance101_                 epsgetbalance101
#define epssetdeflationspace0_            epssetdeflationspace0
#define epssetdeflationspace1_            epssetdeflationspace1
#define epssetinitialspace0_              epssetinitialspace0
#define epssetinitialspace1_              epssetinitialspace1
#define epssetleftinitialspace0_          epssetleftinitialspace0
#define epssetleftinitialspace1_          epssetleftinitialspace1
#endif

/*
   These are not usually called from Fortran but allow Fortran users
   to transparently set these monitors from .F code
*/
SLEPC_EXTERN void epsmonitorall_(EPS *eps,PetscInt *it,PetscInt *nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt *nest,PetscViewerAndFormat **ctx,PetscErrorCode *ierr)
{
  *ierr = EPSMonitorAll(*eps,*it,*nconv,eigr,eigi,errest,*nest,*ctx);
}

SLEPC_EXTERN void epsmonitorconverged_(EPS *eps,PetscInt *it,PetscInt *nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt *nest,SlepcConvMonitor *ctx,PetscErrorCode *ierr)
{
  *ierr = EPSMonitorConverged(*eps,*it,*nconv,eigr,eigi,errest,*nest,*ctx);
}

SLEPC_EXTERN void epsmonitorfirst_(EPS *eps,PetscInt *it,PetscInt *nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt *nest,PetscViewerAndFormat **ctx,PetscErrorCode *ierr)
{
  *ierr = EPSMonitorFirst(*eps,*it,*nconv,eigr,eigi,errest,*nest,*ctx);
}

SLEPC_EXTERN void epsmonitorlg_(EPS *eps,PetscInt *it,PetscInt *nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt *nest,void *ctx,PetscErrorCode *ierr)
{
  *ierr = EPSMonitorLG(*eps,*it,*nconv,eigr,eigi,errest,*nest,ctx);
}

SLEPC_EXTERN void epsmonitorlgall_(EPS *eps,PetscInt *it,PetscInt *nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt *nest,void *ctx,PetscErrorCode *ierr)
{
  *ierr = EPSMonitorLGAll(*eps,*it,*nconv,eigr,eigi,errest,*nest,ctx);
}

static struct {
  PetscFortranCallbackId monitor;
  PetscFortranCallbackId monitordestroy;
  PetscFortranCallbackId convergence;
  PetscFortranCallbackId convdestroy;
  PetscFortranCallbackId stopping;
  PetscFortranCallbackId stopdestroy;
  PetscFortranCallbackId comparison;
  PetscFortranCallbackId arbitrary;
} _cb;

/* These are not extern C because they are passed into non-extern C user level functions */
static PetscErrorCode ourmonitor(EPS eps,PetscInt i,PetscInt nc,PetscScalar *er,PetscScalar *ei,PetscReal *d,PetscInt l,void* ctx)
{
  PetscObjectUseFortranCallback(eps,_cb.monitor,(EPS*,PetscInt*,PetscInt*,PetscScalar*,PetscScalar*,PetscReal*,PetscInt*,void*,PetscErrorCode*),(&eps,&i,&nc,er,ei,d,&l,_ctx,&ierr));
}

static PetscErrorCode ourdestroy(void** ctx)
{
  EPS eps = (EPS)*ctx;
  PetscObjectUseFortranCallback(eps,_cb.monitordestroy,(void*,PetscErrorCode*),(_ctx,&ierr));
}

static PetscErrorCode ourconvergence(EPS eps,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscObjectUseFortranCallback(eps,_cb.convergence,(EPS*,PetscScalar*,PetscScalar*,PetscReal*,PetscReal*,void*,PetscErrorCode*),(&eps,&eigr,&eigi,&res,errest,_ctx,&ierr));
}

static PetscErrorCode ourconvdestroy(void *ctx)
{
  EPS eps = (EPS)ctx;
  PetscObjectUseFortranCallback(eps,_cb.convdestroy,(void*,PetscErrorCode*),(_ctx,&ierr));
}

static PetscErrorCode ourstopping(EPS eps,PetscInt its,PetscInt max_it,PetscInt nconv,PetscInt nev,EPSConvergedReason *reason,void *ctx)
{
  PetscObjectUseFortranCallback(eps,_cb.stopping,(EPS*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,EPSConvergedReason*,void*,PetscErrorCode*),(&eps,&its,&max_it,&nconv,&nev,reason,_ctx,&ierr));
}

static PetscErrorCode ourstopdestroy(void *ctx)
{
  EPS eps = (EPS)ctx;
  PetscObjectUseFortranCallback(eps,_cb.stopdestroy,(void*,PetscErrorCode*),(_ctx,&ierr));
}

static PetscErrorCode oureigenvaluecomparison(PetscScalar ar,PetscScalar ai,PetscScalar br,PetscScalar bi,PetscInt *r,void *ctx)
{
  EPS eps = (EPS)ctx;
  PetscObjectUseFortranCallback(eps,_cb.comparison,(PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt*,void*,PetscErrorCode*),(&ar,&ai,&br,&bi,r,_ctx,&ierr));
}

static PetscErrorCode ourarbitraryfunc(PetscScalar er,PetscScalar ei,Vec xr,Vec xi,PetscScalar *rr,PetscScalar *ri,void *ctx)
{
  EPS eps = (EPS)ctx;
  PetscObjectUseFortranCallback(eps,_cb.arbitrary,(PetscScalar*,PetscScalar*,Vec*,Vec*,PetscScalar*,PetscScalar*,void*,PetscErrorCode*),(&er,&ei,&xr,&xi,rr,ri,_ctx,&ierr));
}

SLEPC_EXTERN void epsview_(EPS *eps,PetscViewer *viewer,PetscErrorCode *ierr)
{
  PetscViewer v;
  PetscPatchDefaultViewers_Fortran(viewer,v);
  *ierr = EPSView(*eps,v);
}

SLEPC_EXTERN void epsviewfromoptions_(EPS *eps,PetscObject obj,char* type,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(type,len,t);
  *ierr = EPSViewFromOptions(*eps,obj,t);if (*ierr) return;
  FREECHAR(type,t);
}

SLEPC_EXTERN void epsreasonview_(EPS *eps,PetscViewer *viewer,PetscErrorCode *ierr)
{
  PetscViewer v;
  PetscPatchDefaultViewers_Fortran(viewer,v);
  *ierr = EPSReasonView(*eps,v);
}

SLEPC_EXTERN void epserrorview_(EPS *eps,EPSErrorType *etype,PetscViewer *viewer,PetscErrorCode *ierr)
{
  PetscViewer v;
  PetscPatchDefaultViewers_Fortran(viewer,v);
  *ierr = EPSErrorView(*eps,*etype,v);
}

SLEPC_EXTERN void epsvaluesview_(EPS *eps,PetscViewer *viewer,PetscErrorCode *ierr)
{
  PetscViewer v;
  PetscPatchDefaultViewers_Fortran(viewer,v);
  *ierr = EPSValuesView(*eps,v);
}

SLEPC_EXTERN void epsvectorsview_(EPS *eps,PetscViewer *viewer,PetscErrorCode *ierr)
{
  PetscViewer v;
  PetscPatchDefaultViewers_Fortran(viewer,v);
  *ierr = EPSVectorsView(*eps,v);
}

SLEPC_EXTERN void epssettype_(EPS *eps,char *type,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(type,len,t);
  *ierr = EPSSetType(*eps,t);if (*ierr) return;
  FREECHAR(type,t);
}

SLEPC_EXTERN void epsgettype_(EPS *eps,char *name,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  EPSType tname;

  *ierr = EPSGetType(*eps,&tname);if (*ierr) return;
  *ierr = PetscStrncpy(name,tname,len);if (*ierr) return;
  FIXRETURNCHAR(PETSC_TRUE,name,len);
}

SLEPC_EXTERN void epssetoptionsprefix_(EPS *eps,char *prefix,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(prefix,len,t);
  *ierr = EPSSetOptionsPrefix(*eps,t);if (*ierr) return;
  FREECHAR(prefix,t);
}

SLEPC_EXTERN void epsappendoptionsprefix_(EPS *eps,char *prefix,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(prefix,len,t);
  *ierr = EPSAppendOptionsPrefix(*eps,t);if (*ierr) return;
  FREECHAR(prefix,t);
}

SLEPC_EXTERN void epsgetoptionsprefix_(EPS *eps,char *prefix,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  const char *tname;

  *ierr = EPSGetOptionsPrefix(*eps,&tname); if (*ierr) return;
  *ierr = PetscStrncpy(prefix,tname,len);if (*ierr) return;
  FIXRETURNCHAR(PETSC_TRUE,prefix,len);
}

SLEPC_EXTERN void epsmonitorset_(EPS *eps,void (*monitor)(EPS*,PetscInt*,PetscInt*,PetscScalar*,PetscScalar*,PetscReal*,PetscInt*,void*,PetscErrorCode*),void *mctx,void (*monitordestroy)(void *,PetscErrorCode*),PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(mctx);
  CHKFORTRANNULLFUNCTION(monitordestroy);
  if ((PetscVoidFunction)monitor == (PetscVoidFunction)epsmonitorall_) {
    *ierr = EPSMonitorSet(*eps,(PetscErrorCode (*)(EPS,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*))EPSMonitorAll,*(PetscViewerAndFormat**)mctx,(PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy);
  } else if ((PetscVoidFunction)monitor == (PetscVoidFunction)epsmonitorconverged_) {
    *ierr = EPSMonitorSet(*eps,(PetscErrorCode (*)(EPS,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*))EPSMonitorConverged,*(SlepcConvMonitor*)mctx,(PetscErrorCode (*)(void**))SlepcConvMonitorDestroy);
  } else if ((PetscVoidFunction)monitor == (PetscVoidFunction)epsmonitorfirst_) {
    *ierr = EPSMonitorSet(*eps,(PetscErrorCode (*)(EPS,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*))EPSMonitorFirst,*(PetscViewerAndFormat**)mctx,(PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy);
  } else if ((PetscVoidFunction)monitor == (PetscVoidFunction)epsmonitorlg_) {
    *ierr = EPSMonitorSet(*eps,EPSMonitorLG,0,0);
  } else if ((PetscVoidFunction)monitor == (PetscVoidFunction)epsmonitorlgall_) {
    *ierr = EPSMonitorSet(*eps,EPSMonitorLGAll,0,0);
  } else {
    *ierr = PetscObjectSetFortranCallback((PetscObject)*eps,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.monitor,(PetscVoidFunction)monitor,mctx); if (*ierr) return;
    *ierr = PetscObjectSetFortranCallback((PetscObject)*eps,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.monitordestroy,(PetscVoidFunction)monitordestroy,mctx); if (*ierr) return;
    *ierr = EPSMonitorSet(*eps,ourmonitor,*eps,ourdestroy);
  }
}

SLEPC_EXTERN void epsconvergedabsolute_(EPS *eps,PetscScalar *eigr,PetscScalar *eigi,PetscReal *res,PetscReal *errest,void *ctx,PetscErrorCode *ierr)
{
  *ierr = EPSConvergedAbsolute(*eps,*eigr,*eigi,*res,errest,ctx);
}

SLEPC_EXTERN void epsconvergedrelative_(EPS *eps,PetscScalar *eigr,PetscScalar *eigi,PetscReal *res,PetscReal *errest,void *ctx,PetscErrorCode *ierr)
{
  *ierr = EPSConvergedRelative(*eps,*eigr,*eigi,*res,errest,ctx);
}

SLEPC_EXTERN void epsconvergednorm_(EPS *eps,PetscScalar *eigr,PetscScalar *eigi,PetscReal *res,PetscReal *errest,void *ctx,PetscErrorCode *ierr)
{
  *ierr = EPSConvergedNorm(*eps,*eigr,*eigi,*res,errest,ctx);
}

SLEPC_EXTERN void epssetconvergencetestfunction_(EPS *eps,void (*func)(EPS*,PetscScalar*,PetscScalar*,PetscReal*,PetscReal*,void*,PetscErrorCode*),void* ctx,void (*destroy)(void*,PetscErrorCode*),PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(ctx);
  CHKFORTRANNULLFUNCTION(destroy);
  if ((PetscVoidFunction)func == (PetscVoidFunction)epsconvergedabsolute_) {
    *ierr = EPSSetConvergenceTest(*eps,EPS_CONV_ABS);
  } else if ((PetscVoidFunction)func == (PetscVoidFunction)epsconvergedrelative_) {
    *ierr = EPSSetConvergenceTest(*eps,EPS_CONV_REL);
  } else if ((PetscVoidFunction)func == (PetscVoidFunction)epsconvergednorm_) {
    *ierr = EPSSetConvergenceTest(*eps,EPS_CONV_NORM);
  } else {
    *ierr = PetscObjectSetFortranCallback((PetscObject)*eps,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.convergence,(PetscVoidFunction)func,ctx); if (*ierr) return;
    *ierr = PetscObjectSetFortranCallback((PetscObject)*eps,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.convdestroy,(PetscVoidFunction)destroy,ctx); if (*ierr) return;
    *ierr = EPSSetConvergenceTestFunction(*eps,ourconvergence,*eps,ourconvdestroy);
  }
}

SLEPC_EXTERN void epsstoppingbasic_(EPS *eps,PetscInt *its,PetscInt *max_it,PetscInt *nconv,PetscInt *nev,EPSConvergedReason *reason,void *ctx,PetscErrorCode *ierr)
{
  *ierr = EPSStoppingBasic(*eps,*its,*max_it,*nconv,*nev,reason,ctx);
}

SLEPC_EXTERN void epssetstoppingtestfunction_(EPS *eps,void (*func)(EPS*,PetscInt,PetscInt,PetscInt,PetscInt,EPSConvergedReason*,void*,PetscErrorCode*),void* ctx,void (*destroy)(void*,PetscErrorCode*),PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(ctx);
  CHKFORTRANNULLFUNCTION(destroy);
  if ((PetscVoidFunction)func == (PetscVoidFunction)epsstoppingbasic_) {
    *ierr = EPSSetStoppingTest(*eps,EPS_STOP_BASIC);
  } else {
    *ierr = PetscObjectSetFortranCallback((PetscObject)*eps,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.stopping,(PetscVoidFunction)func,ctx); if (*ierr) return;
    *ierr = PetscObjectSetFortranCallback((PetscObject)*eps,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.stopdestroy,(PetscVoidFunction)destroy,ctx); if (*ierr) return;
    *ierr = EPSSetStoppingTestFunction(*eps,ourstopping,*eps,ourstopdestroy);
  }
}

SLEPC_EXTERN void epsseteigenvaluecomparison_(EPS *eps,void (*func)(PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt*,void*),void* ctx,PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(ctx);
  *ierr = PetscObjectSetFortranCallback((PetscObject)*eps,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.comparison,(PetscVoidFunction)func,ctx); if (*ierr) return;
  *ierr = EPSSetEigenvalueComparison(*eps,oureigenvaluecomparison,*eps);
}

SLEPC_EXTERN void epssetarbitraryselection_(EPS *eps,void (*func)(PetscScalar*,PetscScalar*,Vec*,Vec*,PetscScalar*,PetscScalar*,void*,PetscErrorCode*),void *ctx,PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(ctx);
  *ierr = PetscObjectSetFortranCallback((PetscObject)*eps,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.arbitrary,(PetscVoidFunction)func,ctx); if (*ierr) return;
  *ierr = EPSSetArbitrarySelection(*eps,ourarbitraryfunc,*eps);
}

SLEPC_EXTERN void epsgetdimensions_(EPS *eps,PetscInt *nev,PetscInt *ncv,PetscInt *mpd,PetscErrorCode *ierr)
{
  CHKFORTRANNULLINTEGER(nev);
  CHKFORTRANNULLINTEGER(ncv);
  CHKFORTRANNULLINTEGER(mpd);
  *ierr = EPSGetDimensions(*eps,nev,ncv,mpd);
}

SLEPC_EXTERN void epsgetdimensions000_(EPS *eps,PetscInt *nev,PetscInt *ncv,PetscInt *mpd,PetscErrorCode *ierr)
{
  epsgetdimensions_(eps,nev,ncv,mpd,ierr);
}

SLEPC_EXTERN void epsgetdimensions100_(EPS *eps,PetscInt *nev,PetscInt *ncv,PetscInt *mpd,PetscErrorCode *ierr)
{
  epsgetdimensions_(eps,nev,ncv,mpd,ierr);
}

SLEPC_EXTERN void epsgetdimensions010_(EPS *eps,PetscInt *nev,PetscInt *ncv,PetscInt *mpd,PetscErrorCode *ierr)
{
  epsgetdimensions_(eps,nev,ncv,mpd,ierr);
}

SLEPC_EXTERN void epsgetdimensions001_(EPS *eps,PetscInt *nev,PetscInt *ncv,PetscInt *mpd,PetscErrorCode *ierr)
{
  epsgetdimensions_(eps,nev,ncv,mpd,ierr);
}

SLEPC_EXTERN void epsgetdimensions110_(EPS *eps,PetscInt *nev,PetscInt *ncv,PetscInt *mpd,PetscErrorCode *ierr)
{
  epsgetdimensions_(eps,nev,ncv,mpd,ierr);
}

SLEPC_EXTERN void epsgetdimensions011_(EPS *eps,PetscInt *nev,PetscInt *ncv,PetscInt *mpd,PetscErrorCode *ierr)
{
  epsgetdimensions_(eps,nev,ncv,mpd,ierr);
}

SLEPC_EXTERN void epsgetdimensions101_(EPS *eps,PetscInt *nev,PetscInt *ncv,PetscInt *mpd,PetscErrorCode *ierr)
{
  epsgetdimensions_(eps,nev,ncv,mpd,ierr);
}

SLEPC_EXTERN void epsgeteigenpair_(EPS *eps,PetscInt *i,PetscScalar *eigr,PetscScalar *eigi,Vec *Vr,Vec *Vi,PetscErrorCode *ierr)
{
  CHKFORTRANNULLSCALAR(eigr);
  CHKFORTRANNULLSCALAR(eigi);
  *ierr = EPSGetEigenpair(*eps,*i,eigr,eigi,*Vr,*Vi);
}

SLEPC_EXTERN void epsgeteigenpair00_(EPS *eps,PetscInt *i,PetscScalar *eigr,PetscScalar *eigi,Vec *Vr,Vec *Vi,PetscErrorCode *ierr)
{
  epsgeteigenpair_(eps,i,eigr,eigi,Vr,Vi,ierr);
}

SLEPC_EXTERN void epsgeteigenpair10_(EPS *eps,PetscInt *i,PetscScalar *eigr,PetscScalar *eigi,Vec *Vr,Vec *Vi,PetscErrorCode *ierr)
{
  epsgeteigenpair_(eps,i,eigr,eigi,Vr,Vi,ierr);
}

SLEPC_EXTERN void epsgeteigenpair01_(EPS *eps,PetscInt *i,PetscScalar *eigr,PetscScalar *eigi,Vec *Vr,Vec *Vi,PetscErrorCode *ierr)
{
  epsgeteigenpair_(eps,i,eigr,eigi,Vr,Vi,ierr);
}

SLEPC_EXTERN void epsgeteigenpair11_(EPS *eps,PetscInt *i,PetscScalar *eigr,PetscScalar *eigi,Vec *Vr,Vec *Vi,PetscErrorCode *ierr)
{
  epsgeteigenpair_(eps,i,eigr,eigi,Vr,Vi,ierr);
}

SLEPC_EXTERN void epsgeteigenvalue_(EPS *eps,PetscInt *i,PetscScalar *eigr,PetscScalar *eigi,PetscErrorCode *ierr)
{
  CHKFORTRANNULLSCALAR(eigr);
  CHKFORTRANNULLSCALAR(eigi);
  *ierr = EPSGetEigenvalue(*eps,*i,eigr,eigi);
}

SLEPC_EXTERN void epsgeteigenvalue00_(EPS *eps,PetscInt *i,PetscScalar *eigr,PetscScalar *eigi,PetscErrorCode *ierr)
{
  epsgeteigenvalue_(eps,i,eigr,eigi,ierr);
}

SLEPC_EXTERN void epsgeteigenvalue10_(EPS *eps,PetscInt *i,PetscScalar *eigr,PetscScalar *eigi,PetscErrorCode *ierr)
{
  epsgeteigenvalue_(eps,i,eigr,eigi,ierr);
}

SLEPC_EXTERN void epsgeteigenvalue01_(EPS *eps,PetscInt *i,PetscScalar *eigr,PetscScalar *eigi,PetscErrorCode *ierr)
{
  epsgeteigenvalue_(eps,i,eigr,eigi,ierr);
}

SLEPC_EXTERN void epsgettolerances_(EPS *eps,PetscReal *tol,PetscInt *maxits,PetscErrorCode *ierr)
{
  CHKFORTRANNULLREAL(tol);
  CHKFORTRANNULLINTEGER(maxits);
  *ierr = EPSGetTolerances(*eps,tol,maxits);
}

SLEPC_EXTERN void epsgettolerances00_(EPS *eps,PetscReal *tol,PetscInt *maxits,PetscErrorCode *ierr)
{
  epsgettolerances_(eps,tol,maxits,ierr);
}

SLEPC_EXTERN void epsgettolerances10_(EPS *eps,PetscReal *tol,PetscInt *maxits,PetscErrorCode *ierr)
{
  epsgettolerances_(eps,tol,maxits,ierr);
}

SLEPC_EXTERN void epsgettolerances01_(EPS *eps,PetscReal *tol,PetscInt *maxits,PetscErrorCode *ierr)
{
  epsgettolerances_(eps,tol,maxits,ierr);
}

SLEPC_EXTERN void epsgetbalance_(EPS *eps,EPSBalance *bal,PetscInt *its,PetscReal *cutoff,PetscErrorCode *ierr)
{
  CHKFORTRANNULLINTEGER(bal);
  CHKFORTRANNULLINTEGER(its);
  CHKFORTRANNULLREAL(cutoff);
  *ierr = EPSGetBalance(*eps,bal,its,cutoff);
}

SLEPC_EXTERN void epsgetbalance000_(EPS *eps,EPSBalance *bal,PetscInt *its,PetscReal *cutoff,PetscErrorCode *ierr)
{
  epsgetbalance_(eps,bal,its,cutoff,ierr);
}

SLEPC_EXTERN void epsgetbalance100_(EPS *eps,EPSBalance *bal,PetscInt *its,PetscReal *cutoff,PetscErrorCode *ierr)
{
  epsgetbalance_(eps,bal,its,cutoff,ierr);
}

SLEPC_EXTERN void epsgetbalance010_(EPS *eps,EPSBalance *bal,PetscInt *its,PetscReal *cutoff,PetscErrorCode *ierr)
{
  epsgetbalance_(eps,bal,its,cutoff,ierr);
}

SLEPC_EXTERN void epsgetbalance001_(EPS *eps,EPSBalance *bal,PetscInt *its,PetscReal *cutoff,PetscErrorCode *ierr)
{
  epsgetbalance_(eps,bal,its,cutoff,ierr);
}

SLEPC_EXTERN void epsgetbalance110_(EPS *eps,EPSBalance *bal,PetscInt *its,PetscReal *cutoff,PetscErrorCode *ierr)
{
  epsgetbalance_(eps,bal,its,cutoff,ierr);
}

SLEPC_EXTERN void epsgetbalance011_(EPS *eps,EPSBalance *bal,PetscInt *its,PetscReal *cutoff,PetscErrorCode *ierr)
{
  epsgetbalance_(eps,bal,its,cutoff,ierr);
}

SLEPC_EXTERN void epsgetbalance101_(EPS *eps,EPSBalance *bal,PetscInt *its,PetscReal *cutoff,PetscErrorCode *ierr)
{
  epsgetbalance_(eps,bal,its,cutoff,ierr);
}

SLEPC_EXTERN void epssetdeflationspace0_(EPS *eps,PetscInt *n,Vec *ds,PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(ds);
  *ierr = EPSSetDeflationSpace(*eps,*n,ds);
}

SLEPC_EXTERN void epssetdeflationspace1_(EPS *eps,PetscInt *n,Vec *ds,PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(ds);
  *ierr = EPSSetDeflationSpace(*eps,*n,ds);
}

SLEPC_EXTERN void epssetinitialspace0_(EPS *eps,PetscInt *n,Vec *is,PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(is);
  *ierr = EPSSetInitialSpace(*eps,*n,is);
}

SLEPC_EXTERN void epssetinitialspace1_(EPS *eps,PetscInt *n,Vec *is,PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(is);
  *ierr = EPSSetInitialSpace(*eps,*n,is);
}

SLEPC_EXTERN void epssetleftinitialspace0_(EPS *eps,PetscInt *n,Vec *isl,PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(isl);
  *ierr = EPSSetLeftInitialSpace(*eps,*n,isl);
}

SLEPC_EXTERN void epssetleftinitialspace1_(EPS *eps,PetscInt *n,Vec *isl,PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(isl);
  *ierr = EPSSetLeftInitialSpace(*eps,*n,isl);
}

