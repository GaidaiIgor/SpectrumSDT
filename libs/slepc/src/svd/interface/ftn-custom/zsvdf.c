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
#include <slepcsvd.h>

#if defined(PETSC_HAVE_FORTRAN_CAPS)
#define svdmonitorall_                    SVDMONITORALL
#define svdmonitorlg_                     SVDMONITORLG
#define svdmonitorlgall_                  SVDMONITORLGALL
#define svdmonitorconverged_              SVDMONITORCONVERGED
#define svdmonitorfirst_                  SVDMONITORFIRST
#define svdview_                          SVDVIEW
#define svderrorview_                     SVDERRORVIEW
#define svdreasonview_                    SVDREASONVIEW
#define svdvaluesview_                    SVDVALUESVIEW
#define svdvectorsview_                   SVDVECTORSVIEW
#define svdsettype_                       SVDSETTYPE
#define svdgettype_                       SVDGETTYPE
#define svdmonitorset_                    SVDMONITORSET
#define svdsetoptionsprefix_              SVDSETOPTIONSPREFIX
#define svdappendoptionsprefix_           SVDAPPENDOPTIONSPREFIX
#define svdgetoptionsprefix_              SVDGETOPTIONSPREFIX
#define svdconvergedabsolute_             SVDCONVERGEDABSOLUTE
#define svdconvergedrelative_             SVDCONVERGEDRELATIVE
#define svdsetconvergencetestfunction_    SVDSETCONVERGENCETESTFUNCTION
#define svdsetstoppingtestfunction_       SVDSETSTOPPINGTESTFUNCTION
#define svdgetdimensions000_              SVDGETDIMENSIONS000
#define svdgetdimensions100_              SVDGETDIMENSIONS100
#define svdgetdimensions010_              SVDGETDIMENSIONS010
#define svdgetdimensions001_              SVDGETDIMENSIONS001
#define svdgetdimensions110_              SVDGETDIMENSIONS110
#define svdgetdimensions011_              SVDGETDIMENSIONS011
#define svdgetdimensions101_              SVDGETDIMENSIONS101
#define svdgetsingulartriplet0_           SVDGETSINGULARTRIPLET0
#define svdgetsingulartriplet1_           SVDGETSINGULARTRIPLET1
#define svdgettolerances00_               SVDGETTOLERANCES00
#define svdgettolerances10_               SVDGETTOLERANCES10
#define svdgettolerances01_               SVDGETTOLERANCES01
#define svdsetinitialspaces00_            SVDSETINITIALSPACES00
#define svdsetinitialspaces01_            SVDSETINITIALSPACES01
#define svdsetinitialspaces10_            SVDSETINITIALSPACES10
#define svdsetinitialspaces11_            SVDSETINITIALSPACES11
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define svdmonitorall_                    svdmonitorall
#define svdmonitorlg_                     svdmonitorlg
#define svdmonitorlgall_                  svdmonitorlgall
#define svdmonitorconverged_              svdmonitorconverged
#define svdmonitorfirst_                  svdmonitorfirst
#define svdview_                          svdview
#define svderrorview_                     svderrorview
#define svdreasonview_                    svdreasonview
#define svdvaluesview_                    svdvaluesview
#define svdvectorsview_                   svdvectorsview
#define svdsettype_                       svdsettype
#define svdgettype_                       svdgettype
#define svdmonitorset_                    svdmonitorset
#define svdsetoptionsprefix_              svdsetoptionsprefix
#define svdappendoptionsprefix_           svdappendoptionsprefix
#define svdgetoptionsprefix_              svdgetoptionsprefix
#define svdconvergedabsolute_             svdconvergedabsolute
#define svdconvergedrelative_             svdconvergedrelative
#define svdsetconvergencetestfunction_    svdsetconvergencetestfunction
#define svdsetstoppingtestfunction_       svdsetstoppingtestfunction
#define svdgetdimensions000_              svdgetdimensions000
#define svdgetdimensions100_              svdgetdimensions100
#define svdgetdimensions010_              svdgetdimensions010
#define svdgetdimensions001_              svdgetdimensions001
#define svdgetdimensions110_              svdgetdimensions110
#define svdgetdimensions011_              svdgetdimensions011
#define svdgetdimensions101_              svdgetdimensions101
#define svdgetsingulartriplet0_           svdgetsingulartriplet0
#define svdgetsingulartriplet1_           svdgetsingulartriplet1
#define svdgettolerances00_               svdgettolerances00
#define svdgettolerances10_               svdgettolerances10
#define svdgettolerances01_               svdgettolerances01
#define svdsetinitialspaces00_            svdsetinitialspaces00
#define svdsetinitialspaces01_            svdsetinitialspaces01
#define svdsetinitialspaces10_            svdsetinitialspaces10
#define svdsetinitialspaces11_            svdsetinitialspaces11
#endif

/*
   These are not usually called from Fortran but allow Fortran users
   to transparently set these monitors from .F code
*/
SLEPC_EXTERN void svdmonitorall_(SVD *svd,PetscInt *it,PetscInt *nconv,PetscReal *sigma,PetscReal *errest,PetscInt *nest,PetscViewerAndFormat **ctx,PetscErrorCode *ierr)
{
  *ierr = SVDMonitorAll(*svd,*it,*nconv,sigma,errest,*nest,*ctx);
}

SLEPC_EXTERN void svdmonitorconverged_(SVD *svd,PetscInt *it,PetscInt *nconv,PetscReal *sigma,PetscReal *errest,PetscInt *nest,SlepcConvMonitor *ctx,PetscErrorCode *ierr)
{
  *ierr = SVDMonitorConverged(*svd,*it,*nconv,sigma,errest,*nest,*ctx);
}

SLEPC_EXTERN void svdmonitorfirst_(SVD *svd,PetscInt *it,PetscInt *nconv,PetscReal *sigma,PetscReal *errest,PetscInt *nest,PetscViewerAndFormat **ctx,PetscErrorCode *ierr)
{
  *ierr = SVDMonitorFirst(*svd,*it,*nconv,sigma,errest,*nest,*ctx);
}

SLEPC_EXTERN void svdmonitorlg_(SVD *svd,PetscInt *it,PetscInt *nconv,PetscReal *sigma,PetscReal *errest,PetscInt *nest,void *ctx,PetscErrorCode *ierr)
{
  *ierr = SVDMonitorLG(*svd,*it,*nconv,sigma,errest,*nest,ctx);
}

SLEPC_EXTERN void svdmonitorlgall_(SVD *svd,PetscInt *it,PetscInt *nconv,PetscReal *sigma,PetscReal *errest,PetscInt *nest,void *ctx,PetscErrorCode *ierr)
{
  *ierr = SVDMonitorLGAll(*svd,*it,*nconv,sigma,errest,*nest,ctx);
}

static struct {
  PetscFortranCallbackId monitor;
  PetscFortranCallbackId monitordestroy;
  PetscFortranCallbackId convergence;
  PetscFortranCallbackId convdestroy;
  PetscFortranCallbackId stopping;
  PetscFortranCallbackId stopdestroy;
} _cb;

/* These are not extern C because they are passed into non-extern C user level functions */
static PetscErrorCode ourmonitor(SVD svd,PetscInt i,PetscInt nc,PetscReal *sigma,PetscReal *d,PetscInt l,void* ctx)
{
  PetscObjectUseFortranCallback(svd,_cb.monitor,(SVD*,PetscInt*,PetscInt*,PetscReal*,PetscReal*,PetscInt*,void*,PetscErrorCode*),(&svd,&i,&nc,sigma,d,&l,_ctx,&ierr));
}

static PetscErrorCode ourdestroy(void** ctx)
{
  SVD svd = (SVD)*ctx;
  PetscObjectUseFortranCallback(svd,_cb.monitordestroy,(void*,PetscErrorCode*),(_ctx,&ierr));
}

static PetscErrorCode ourconvergence(SVD svd,PetscReal sigma,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscObjectUseFortranCallback(svd,_cb.convergence,(SVD*,PetscReal*,PetscReal*,PetscReal*,void*,PetscErrorCode*),(&svd,&sigma,&res,errest,_ctx,&ierr));
}

static PetscErrorCode ourconvdestroy(void *ctx)
{
  SVD svd = (SVD)ctx;
  PetscObjectUseFortranCallback(svd,_cb.convdestroy,(void*,PetscErrorCode*),(_ctx,&ierr));
}

static PetscErrorCode ourstopping(SVD svd,PetscInt its,PetscInt max_it,PetscInt nconv,PetscInt nsv,SVDConvergedReason *reason,void *ctx)
{
  PetscObjectUseFortranCallback(svd,_cb.stopping,(SVD*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,SVDConvergedReason*,void*,PetscErrorCode*),(&svd,&its,&max_it,&nconv,&nsv,reason,_ctx,&ierr));
}

static PetscErrorCode ourstopdestroy(void *ctx)
{
  SVD svd = (SVD)ctx;
  PetscObjectUseFortranCallback(svd,_cb.stopdestroy,(void*,PetscErrorCode*),(_ctx,&ierr));
}

SLEPC_EXTERN void svdview_(SVD *svd,PetscViewer *viewer,PetscErrorCode *ierr)
{
  PetscViewer v;
  PetscPatchDefaultViewers_Fortran(viewer,v);
  *ierr = SVDView(*svd,v);
}

SLEPC_EXTERN void svdreasonview_(SVD *svd,PetscViewer *viewer,PetscErrorCode *ierr)
{
  PetscViewer v;
  PetscPatchDefaultViewers_Fortran(viewer,v);
  *ierr = SVDReasonView(*svd,v);
}

SLEPC_EXTERN void svderrorview_(SVD *svd,SVDErrorType *etype,PetscViewer *viewer,PetscErrorCode *ierr)
{
  PetscViewer v;
  PetscPatchDefaultViewers_Fortran(viewer,v);
  *ierr = SVDErrorView(*svd,*etype,v);
}

SLEPC_EXTERN void svdvaluesview_(SVD *svd,PetscViewer *viewer,PetscErrorCode *ierr)
{
  PetscViewer v;
  PetscPatchDefaultViewers_Fortran(viewer,v);
  *ierr = SVDValuesView(*svd,v);
}

SLEPC_EXTERN void svdvectorsview_(SVD *svd,PetscViewer *viewer,PetscErrorCode *ierr)
{
  PetscViewer v;
  PetscPatchDefaultViewers_Fortran(viewer,v);
  *ierr = SVDVectorsView(*svd,v);
}

SLEPC_EXTERN void svdsettype_(SVD *svd,char *type,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(type,len,t);
  *ierr = SVDSetType(*svd,t);if (*ierr) return;
  FREECHAR(type,t);
}

SLEPC_EXTERN void svdgettype_(SVD *svd,char *name,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  SVDType tname;

  *ierr = SVDGetType(*svd,&tname);if (*ierr) return;
  *ierr = PetscStrncpy(name,tname,len);if (*ierr) return;
  FIXRETURNCHAR(PETSC_TRUE,name,len);
}

SLEPC_EXTERN void svdmonitorset_(SVD *svd,void (*monitor)(SVD*,PetscInt*,PetscInt*,PetscReal*,PetscReal*,PetscInt*,void*,PetscErrorCode*),void *mctx,void (*monitordestroy)(void *,PetscErrorCode*),PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(mctx);
  CHKFORTRANNULLFUNCTION(monitordestroy);
  if ((PetscVoidFunction)monitor == (PetscVoidFunction)svdmonitorall_) {
    *ierr = SVDMonitorSet(*svd,(PetscErrorCode (*)(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,void*))SVDMonitorAll,*(PetscViewerAndFormat**)mctx,(PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy);
  } else if ((PetscVoidFunction)monitor == (PetscVoidFunction)svdmonitorconverged_) {
    *ierr = SVDMonitorSet(*svd,(PetscErrorCode (*)(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,void*))SVDMonitorConverged,*(SlepcConvMonitor*)mctx,(PetscErrorCode (*)(void**))SlepcConvMonitorDestroy);
  } else if ((PetscVoidFunction)monitor == (PetscVoidFunction)svdmonitorfirst_) {
    *ierr = SVDMonitorSet(*svd,(PetscErrorCode (*)(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,void*))SVDMonitorFirst,*(PetscViewerAndFormat**)mctx,(PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy);
  } else if ((PetscVoidFunction)monitor == (PetscVoidFunction)svdmonitorlg_) {
    *ierr = SVDMonitorSet(*svd,SVDMonitorLG,0,0);
  } else if ((PetscVoidFunction)monitor == (PetscVoidFunction)svdmonitorlgall_) {
    *ierr = SVDMonitorSet(*svd,SVDMonitorLGAll,0,0);
  } else {
    *ierr = PetscObjectSetFortranCallback((PetscObject)*svd,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.monitor,(PetscVoidFunction)monitor,mctx); if (*ierr) return;
    *ierr = PetscObjectSetFortranCallback((PetscObject)*svd,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.monitordestroy,(PetscVoidFunction)monitordestroy,mctx); if (*ierr) return;
    *ierr = SVDMonitorSet(*svd,ourmonitor,*svd,ourdestroy);
  }
}

SLEPC_EXTERN void svdsetoptionsprefix_(SVD *svd,char *prefix,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(prefix,len,t);
  *ierr = SVDSetOptionsPrefix(*svd,t);if (*ierr) return;
  FREECHAR(prefix,t);
}

SLEPC_EXTERN void svdappendoptionsprefix_(SVD *svd,char *prefix,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(prefix,len,t);
  *ierr = SVDAppendOptionsPrefix(*svd,t);if (*ierr) return;
  FREECHAR(prefix,t);
}

SLEPC_EXTERN void svdgetoptionsprefix_(SVD *svd,char *prefix,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  const char *tname;

  *ierr = SVDGetOptionsPrefix(*svd,&tname); if (*ierr) return;
  *ierr = PetscStrncpy(prefix,tname,len);if (*ierr) return;
  FIXRETURNCHAR(PETSC_TRUE,prefix,len);
}

SLEPC_EXTERN void svdconvergedabsolute_(SVD *svd,PetscReal *sigma,PetscReal *res,PetscReal *errest,void *ctx,PetscErrorCode *ierr)
{
  *ierr = SVDConvergedAbsolute(*svd,*sigma,*res,errest,ctx);
}

SLEPC_EXTERN void svdconvergedrelative_(SVD *svd,PetscReal *sigma,PetscReal *res,PetscReal *errest,void *ctx,PetscErrorCode *ierr)
{
  *ierr = SVDConvergedRelative(*svd,*sigma,*res,errest,ctx);
}

SLEPC_EXTERN void svdsetconvergencetestfunction_(SVD *svd,void (*func)(SVD*,PetscReal*,PetscReal*,PetscReal*,void*,PetscErrorCode*),void* ctx,void (*destroy)(void*,PetscErrorCode*),PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(ctx);
  CHKFORTRANNULLFUNCTION(destroy);
  if ((PetscVoidFunction)func == (PetscVoidFunction)svdconvergedabsolute_) {
    *ierr = SVDSetConvergenceTest(*svd,SVD_CONV_ABS);
  } else if ((PetscVoidFunction)func == (PetscVoidFunction)svdconvergedrelative_) {
    *ierr = SVDSetConvergenceTest(*svd,SVD_CONV_REL);
  } else {
    *ierr = PetscObjectSetFortranCallback((PetscObject)*svd,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.convergence,(PetscVoidFunction)func,ctx); if (*ierr) return;
    *ierr = PetscObjectSetFortranCallback((PetscObject)*svd,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.convdestroy,(PetscVoidFunction)destroy,ctx); if (*ierr) return;
    *ierr = SVDSetConvergenceTestFunction(*svd,ourconvergence,*svd,ourconvdestroy);
  }
}

SLEPC_EXTERN void svdstoppingbasic_(SVD *svd,PetscInt *its,PetscInt *max_it,PetscInt *nconv,PetscInt *nsv,SVDConvergedReason *reason,void *ctx,PetscErrorCode *ierr)
{
  *ierr = SVDStoppingBasic(*svd,*its,*max_it,*nconv,*nsv,reason,ctx);
}

SLEPC_EXTERN void svdsetstoppingtestfunction_(SVD *svd,void (*func)(SVD*,PetscInt,PetscInt,PetscInt,PetscInt,SVDConvergedReason*,void*,PetscErrorCode*),void* ctx,void (*destroy)(void*,PetscErrorCode*),PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(ctx);
  CHKFORTRANNULLFUNCTION(destroy);
  if ((PetscVoidFunction)func == (PetscVoidFunction)svdstoppingbasic_) {
    *ierr = SVDSetStoppingTest(*svd,SVD_STOP_BASIC);
  } else {
    *ierr = PetscObjectSetFortranCallback((PetscObject)*svd,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.stopping,(PetscVoidFunction)func,ctx); if (*ierr) return;
    *ierr = PetscObjectSetFortranCallback((PetscObject)*svd,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.stopdestroy,(PetscVoidFunction)destroy,ctx); if (*ierr) return;
    *ierr = SVDSetStoppingTestFunction(*svd,ourstopping,*svd,ourstopdestroy);
  }
}

SLEPC_EXTERN void svdgetdimensions_(SVD *svd,PetscInt *nev,PetscInt *ncv,PetscInt *mpd,PetscErrorCode *ierr)
{
  CHKFORTRANNULLINTEGER(nev);
  CHKFORTRANNULLINTEGER(ncv);
  CHKFORTRANNULLINTEGER(mpd);
  *ierr = SVDGetDimensions(*svd,nev,ncv,mpd);
}

SLEPC_EXTERN void svdgetdimensions000_(SVD *svd,PetscInt *nev,PetscInt *ncv,PetscInt *mpd,PetscErrorCode *ierr)
{
  svdgetdimensions_(svd,nev,ncv,mpd,ierr);
}

SLEPC_EXTERN void svdgetdimensions100_(SVD *svd,PetscInt *nev,PetscInt *ncv,PetscInt *mpd,PetscErrorCode *ierr)
{
  svdgetdimensions_(svd,nev,ncv,mpd,ierr);
}

SLEPC_EXTERN void svdgetdimensions010_(SVD *svd,PetscInt *nev,PetscInt *ncv,PetscInt *mpd,PetscErrorCode *ierr)
{
  svdgetdimensions_(svd,nev,ncv,mpd,ierr);
}

SLEPC_EXTERN void svdgetdimensions001_(SVD *svd,PetscInt *nev,PetscInt *ncv,PetscInt *mpd,PetscErrorCode *ierr)
{
  svdgetdimensions_(svd,nev,ncv,mpd,ierr);
}

SLEPC_EXTERN void svdgetdimensions110_(SVD *svd,PetscInt *nev,PetscInt *ncv,PetscInt *mpd,PetscErrorCode *ierr)
{
  svdgetdimensions_(svd,nev,ncv,mpd,ierr);
}

SLEPC_EXTERN void svdgetdimensions011_(SVD *svd,PetscInt *nev,PetscInt *ncv,PetscInt *mpd,PetscErrorCode *ierr)
{
  svdgetdimensions_(svd,nev,ncv,mpd,ierr);
}

SLEPC_EXTERN void svdgetdimensions101_(SVD *svd,PetscInt *nev,PetscInt *ncv,PetscInt *mpd,PetscErrorCode *ierr)
{
  svdgetdimensions_(svd,nev,ncv,mpd,ierr);
}

SLEPC_EXTERN void svdgetsingulartriplet_(SVD *svd,PetscInt *i,PetscReal *sigma,Vec *u,Vec *v,PetscErrorCode *ierr)
{
  CHKFORTRANNULLREAL(sigma);
  *ierr = SVDGetSingularTriplet(*svd,*i,sigma,*u,*v);
}

SLEPC_EXTERN void svdgetsingulartriplet0_(SVD *svd,PetscInt *i,PetscReal *sigma,Vec *u,Vec *v, PetscErrorCode *ierr)
{
  svdgetsingulartriplet_(svd,i,sigma,u,v,ierr);
}

SLEPC_EXTERN void svdgetsingulartriplet1_(SVD *svd,PetscInt *i,PetscReal *sigma,Vec *u,Vec *v, PetscErrorCode *ierr)
{
  svdgetsingulartriplet_(svd,i,sigma,u,v,ierr);
}

SLEPC_EXTERN void svdgettolerances_(SVD *svd,PetscReal *tol,PetscInt *maxits,PetscErrorCode *ierr)
{
  CHKFORTRANNULLREAL(tol);
  CHKFORTRANNULLINTEGER(maxits);
  *ierr = SVDGetTolerances(*svd,tol,maxits);
}

SLEPC_EXTERN void svdgettolerances00_(SVD *svd,PetscReal *tol,PetscInt *maxits,PetscErrorCode *ierr)
{
  svdgettolerances_(svd,tol,maxits,ierr);
}

SLEPC_EXTERN void svdgettolerances10_(SVD *svd,PetscReal *tol,PetscInt *maxits,PetscErrorCode *ierr)
{
  svdgettolerances_(svd,tol,maxits,ierr);
}

SLEPC_EXTERN void svdgettolerances01_(SVD *svd,PetscReal *tol,PetscInt *maxits,PetscErrorCode *ierr)
{
  svdgettolerances_(svd,tol,maxits,ierr);
}

SLEPC_EXTERN void svdsetinitialspaces_(SVD *svd,PetscInt *nr,Vec *isr,PetscInt *nl,Vec *isl,PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(isr);
  CHKFORTRANNULLOBJECT(isl);
  *ierr = SVDSetInitialSpaces(*svd,*nr,isr,*nl,isl);
}

SLEPC_EXTERN void svdsetinitialspaces00_(SVD *svd,PetscInt *nr,Vec *isr,PetscInt *nl,Vec *isl,PetscErrorCode *ierr)
{
  svdsetinitialspaces_(svd,nr,isr,nl,isl,ierr);
}

SLEPC_EXTERN void svdsetinitialspaces01_(SVD *svd,PetscInt *nr,Vec *isr,PetscInt *nl,Vec *isl,PetscErrorCode *ierr)
{
  svdsetinitialspaces_(svd,nr,isr,nl,isl,ierr);
}

SLEPC_EXTERN void svdsetinitialspaces10_(SVD *svd,PetscInt *nr,Vec *isr,PetscInt *nl,Vec *isl,PetscErrorCode *ierr)
{
  svdsetinitialspaces_(svd,nr,isr,nl,isl,ierr);
}

SLEPC_EXTERN void svdsetinitialspaces11_(SVD *svd,PetscInt *nr,Vec *isr,PetscInt *nl,Vec *isl,PetscErrorCode *ierr)
{
  svdsetinitialspaces_(svd,nr,isr,nl,isl,ierr);
}

