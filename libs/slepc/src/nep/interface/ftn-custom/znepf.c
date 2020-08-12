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
#include <slepc/private/nepimpl.h>
#include <petsc/private/f90impl.h>

#if defined(PETSC_HAVE_FORTRAN_CAPS)
#define nepview_                          NEPVIEW
#define nepviewfromoptions_               NEPVIEWFROMOPTIONS
#define neperrorview_                     NEPERRORVIEW
#define nepreasonview_                    NEPREASONVIEW
#define nepvaluesview_                    NEPVALUESVIEW
#define nepvectorsview_                   NEPVECTORSVIEW
#define nepsetoptionsprefix_              NEPSETOPTIONSPREFIX
#define nepappendoptionsprefix_           NEPAPPENDOPTIONSPREFIX
#define nepgetoptionsprefix_              NEPGETOPTIONSPREFIX
#define nepsettype_                       NEPSETTYPE
#define nepgettype_                       NEPGETTYPE
#define nepmonitorall_                    NEPMONITORALL
#define nepmonitorlg_                     NEPMONITORLG
#define nepmonitorlgall_                  NEPMONITORLGALL
#define nepmonitorset_                    NEPMONITORSET
#define nepmonitorconverged_              NEPMONITORCONVERGED
#define nepmonitorfirst_                  NEPMONITORFIRST
#define nepconvergedabsolute_             NEPCONVERGEDABSOLUTE
#define nepconvergedrelative_             NEPCONVERGEDRELATIVE
#define nepsetconvergencetestfunction_    NEPSETCONVERGENCETESTFUNCTION
#define nepsetstoppingtestfunction_       NEPSETSTOPPINGTESTFUNCTION
#define nepseteigenvaluecomparison_       NEPSETEIGENVALUECOMPARISON
#define nepsetfunction_                   NEPSETFUNCTION
#define nepgetfunction_                   NEPGETFUNCTION
#define nepsetjacobian_                   NEPSETJACOBIAN
#define nepgetjacobian_                   NEPGETJACOBIAN
#define nepgetdimensions000_              NEPGETDIMENSIONS000
#define nepgetdimensions100_              NEPGETDIMENSIONS100
#define nepgetdimensions010_              NEPGETDIMENSIONS010
#define nepgetdimensions001_              NEPGETDIMENSIONS001
#define nepgetdimensions110_              NEPGETDIMENSIONS110
#define nepgetdimensions011_              NEPGETDIMENSIONS011
#define nepgetdimensions101_              NEPGETDIMENSIONS101
#define nepgeteigenpair00_                NEPGETEIGENPAIR00
#define nepgeteigenpair10_                NEPGETEIGENPAIR10
#define nepgeteigenpair01_                NEPGETEIGENPAIR01
#define nepgeteigenpair11_                NEPGETEIGENPAIR11
#define nepgettolerances00_               NEPGETTOLERANCES00
#define nepgettolerances10_               NEPGETTOLERANCES10
#define nepgettolerances01_               NEPGETTOLERANCES01
#define nepgetrefine000_                  NEPGETREFINE000
#define nepgetrefine100_                  NEPGETREFINE100
#define nepgetrefine010_                  NEPGETREFINE010
#define nepgetrefine001_                  NEPGETREFINE001
#define nepgetrefine110_                  NEPGETREFINE110
#define nepgetrefine011_                  NEPGETREFINE011
#define nepgetrefine101_                  NEPGETREFINE101
#define nepgetrefine111_                  NEPGETREFINE111
#define nepsetinitialspace0_              NEPSETINITIALSPACE0
#define nepsetinitialspace1_              NEPSETINITIALSPACE1
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define nepview_                          nepview
#define nepviewfromoptions_               nepviewfromoptions
#define neperrorview_                     neperrorview
#define nepreasonview_                    nepreasonview
#define nepvaluesview_                    nepvaluesview
#define nepvectorsview_                   nepvectorsview
#define nepsetoptionsprefix_              nepsetoptionsprefix
#define nepappendoptionsprefix_           nepappendoptionsprefix
#define nepgetoptionsprefix_              nepgetoptionsprefix
#define nepsettype_                       nepsettype
#define nepgettype_                       nepgettype
#define nepmonitorall_                    nepmonitorall
#define nepmonitorlg_                     nepmonitorlg
#define nepmonitorlgall_                  nepmonitorlgall
#define nepmonitorset_                    nepmonitorset
#define nepmonitorconverged_              nepmonitorconverged
#define nepmonitorfirst_                  nepmonitorfirst
#define nepconvergedabsolute_             nepconvergedabsolute
#define nepconvergedrelative_             nepconvergedrelative
#define nepsetconvergencetestfunction_    nepsetconvergencetestfunction
#define nepsetstoppingtestfunction_       nepsetstoppingtestfunction
#define nepseteigenvaluecomparison_       nepseteigenvaluecomparison
#define nepsetfunction_                   nepsetfunction
#define nepgetfunction_                   nepgetfunction
#define nepsetjacobian_                   nepsetjacobian
#define nepgetjacobian_                   nepgetjacobian
#define nepgetdimensions000_              nepgetdimensions000
#define nepgetdimensions100_              nepgetdimensions100
#define nepgetdimensions010_              nepgetdimensions010
#define nepgetdimensions001_              nepgetdimensions001
#define nepgetdimensions110_              nepgetdimensions110
#define nepgetdimensions011_              nepgetdimensions011
#define nepgetdimensions101_              nepgetdimensions101
#define nepgeteigenpair00_                nepgeteigenpair00
#define nepgeteigenpair10_                nepgeteigenpair10
#define nepgeteigenpair01_                nepgeteigenpair01
#define nepgeteigenpair11_                nepgeteigenpair11
#define nepgettolerances00_               nepgettolerances00
#define nepgettolerances10_               nepgettolerances10
#define nepgettolerances01_               nepgettolerances01
#define nepgetrefine000_                  nepgetrefine000
#define nepgetrefine100_                  nepgetrefine100
#define nepgetrefine010_                  nepgetrefine010
#define nepgetrefine001_                  nepgetrefine001
#define nepgetrefine110_                  nepgetrefine110
#define nepgetrefine011_                  nepgetrefine011
#define nepgetrefine101_                  nepgetrefine101
#define nepgetrefine111_                  nepgetrefine111
#define nepsetinitialspace0_              nepsetinitialspace0
#define nepsetinitialspace1_              nepsetinitialspace1
#endif

/*
   These are not usually called from Fortran but allow Fortran users
   to transparently set these monitors from .F code
*/
SLEPC_EXTERN void nepmonitorall_(NEP *nep,PetscInt *it,PetscInt *nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt *nest,PetscViewerAndFormat **ctx,PetscErrorCode *ierr)
{
  *ierr = NEPMonitorAll(*nep,*it,*nconv,eigr,eigi,errest,*nest,*ctx);
}

SLEPC_EXTERN void nepmonitorconverged_(NEP *nep,PetscInt *it,PetscInt *nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt *nest,SlepcConvMonitor *ctx,PetscErrorCode *ierr)
{
  *ierr = NEPMonitorConverged(*nep,*it,*nconv,eigr,eigi,errest,*nest,*ctx);
}

SLEPC_EXTERN void nepmonitorfirst_(NEP *nep,PetscInt *it,PetscInt *nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt *nest,PetscViewerAndFormat **ctx,PetscErrorCode *ierr)
{
  *ierr = NEPMonitorFirst(*nep,*it,*nconv,eigr,eigi,errest,*nest,*ctx);
}

SLEPC_EXTERN void nepmonitorlg_(NEP *nep,PetscInt *it,PetscInt *nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt *nest,void *ctx,PetscErrorCode *ierr)
{
  *ierr = NEPMonitorLG(*nep,*it,*nconv,eigr,eigi,errest,*nest,ctx);
}

SLEPC_EXTERN void nepmonitorlgall_(NEP *nep,PetscInt *it,PetscInt *nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt *nest,void *ctx,PetscErrorCode *ierr)
{
  *ierr = NEPMonitorLGAll(*nep,*it,*nconv,eigr,eigi,errest,*nest,ctx);
}

static struct {
  PetscFortranCallbackId monitor;
  PetscFortranCallbackId monitordestroy;
  PetscFortranCallbackId convergence;
  PetscFortranCallbackId convdestroy;
  PetscFortranCallbackId stopping;
  PetscFortranCallbackId stopdestroy;
  PetscFortranCallbackId comparison;
  PetscFortranCallbackId function;
  PetscFortranCallbackId jacobian;
#if defined(PETSC_HAVE_F90_2PTR_ARG)
  PetscFortranCallbackId function_pgiptr;
  PetscFortranCallbackId jacobian_pgiptr;
#endif
} _cb;

/* These are not extern C because they are passed into non-extern C user level functions */
static PetscErrorCode ourmonitor(NEP nep,PetscInt i,PetscInt nc,PetscScalar *er,PetscScalar *ei,PetscReal *d,PetscInt l,void* ctx)
{
  PetscObjectUseFortranCallback(nep,_cb.monitor,(NEP*,PetscInt*,PetscInt*,PetscScalar*,PetscScalar*,PetscReal*,PetscInt*,void*,PetscErrorCode*),(&nep,&i,&nc,er,ei,d,&l,_ctx,&ierr));
}

static PetscErrorCode ourdestroy(void** ctx)
{
  NEP nep = (NEP)*ctx;
  PetscObjectUseFortranCallback(nep,_cb.monitordestroy,(void*,PetscErrorCode*),(_ctx,&ierr));
}

static PetscErrorCode ourconvergence(NEP nep,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscObjectUseFortranCallback(nep,_cb.convergence,(NEP*,PetscScalar*,PetscScalar*,PetscReal*,PetscReal*,void*,PetscErrorCode*),(&nep,&eigr,&eigi,&res,errest,_ctx,&ierr));
}

static PetscErrorCode ourconvdestroy(void *ctx)
{
  NEP nep = (NEP)ctx;
  PetscObjectUseFortranCallback(nep,_cb.convdestroy,(void*,PetscErrorCode*),(_ctx,&ierr));
}

static PetscErrorCode ourstopping(NEP nep,PetscInt its,PetscInt max_it,PetscInt nconv,PetscInt nev,NEPConvergedReason *reason,void *ctx)
{
  PetscObjectUseFortranCallback(nep,_cb.stopping,(NEP*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,NEPConvergedReason*,void*,PetscErrorCode*),(&nep,&its,&max_it,&nconv,&nev,reason,_ctx,&ierr));
}

static PetscErrorCode ourstopdestroy(void *ctx)
{
  NEP nep = (NEP)ctx;
  PetscObjectUseFortranCallback(nep,_cb.stopdestroy,(void*,PetscErrorCode*),(_ctx,&ierr));
}

static PetscErrorCode oureigenvaluecomparison(PetscScalar ar,PetscScalar ai,PetscScalar br,PetscScalar bi,PetscInt *r,void *ctx)
{
  NEP eps = (NEP)ctx;
  PetscObjectUseFortranCallback(eps,_cb.comparison,(PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt*,void*,PetscErrorCode*),(&ar,&ai,&br,&bi,r,_ctx,&ierr));
}

static PetscErrorCode ournepfunction(NEP nep,PetscScalar lambda,Mat T,Mat P,void *ctx)
{
#if defined(PETSC_HAVE_F90_2PTR_ARG)
  void* ptr;
  PetscObjectGetFortranCallback((PetscObject)nep,PETSC_FORTRAN_CALLBACK_CLASS,_cb.function_pgiptr,NULL,&ptr);
#endif
  PetscObjectUseFortranCallback(nep,_cb.function,(NEP*,PetscScalar*,Mat*,Mat*,void*,PetscErrorCode* PETSC_F90_2PTR_PROTO_NOVAR),(&nep,&lambda,&T,&P,_ctx,&ierr PETSC_F90_2PTR_PARAM(ptr)));
}

static PetscErrorCode ournepjacobian(NEP nep,PetscScalar lambda,Mat J,void *ctx)
{
#if defined(PETSC_HAVE_F90_2PTR_ARG)
  void* ptr;
  PetscObjectGetFortranCallback((PetscObject)nep,PETSC_FORTRAN_CALLBACK_CLASS,_cb.jacobian_pgiptr,NULL,&ptr);
#endif
  PetscObjectUseFortranCallback(nep,_cb.jacobian,(NEP*,PetscScalar*,Mat*,void*,PetscErrorCode* PETSC_F90_2PTR_PROTO_NOVAR),(&nep,&lambda,&J,_ctx,&ierr PETSC_F90_2PTR_PARAM(ptr)));
}

SLEPC_EXTERN void nepview_(NEP *nep,PetscViewer *viewer,PetscErrorCode *ierr)
{
  PetscViewer v;
  PetscPatchDefaultViewers_Fortran(viewer,v);
  *ierr = NEPView(*nep,v);
}

SLEPC_EXTERN void nepviewfromoptions_(NEP *nep,PetscObject obj,char* type,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(type,len,t);
  *ierr = NEPViewFromOptions(*nep,obj,t);if (*ierr) return;
  FREECHAR(type,t);
}

SLEPC_EXTERN void nepreasonview_(NEP *nep,PetscViewer *viewer,PetscErrorCode *ierr)
{
  PetscViewer v;
  PetscPatchDefaultViewers_Fortran(viewer,v);
  *ierr = NEPReasonView(*nep,v);
}

SLEPC_EXTERN void neperrorview_(NEP *nep,NEPErrorType *etype,PetscViewer *viewer,PetscErrorCode *ierr)
{
  PetscViewer v;
  PetscPatchDefaultViewers_Fortran(viewer,v);
  *ierr = NEPErrorView(*nep,*etype,v);
}

SLEPC_EXTERN void nepvaluesview_(NEP *nep,PetscViewer *viewer,PetscErrorCode *ierr)
{
  PetscViewer v;
  PetscPatchDefaultViewers_Fortran(viewer,v);
  *ierr = NEPValuesView(*nep,v);
}

SLEPC_EXTERN void nepvectorsview_(NEP *nep,PetscViewer *viewer,PetscErrorCode *ierr)
{
  PetscViewer v;
  PetscPatchDefaultViewers_Fortran(viewer,v);
  *ierr = NEPVectorsView(*nep,v);
}

SLEPC_EXTERN void nepsettype_(NEP *nep,char *type,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(type,len,t);
  *ierr = NEPSetType(*nep,t);if (*ierr) return;
  FREECHAR(type,t);
}

SLEPC_EXTERN void nepgettype_(NEP *nep,char *name,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  NEPType tname;

  *ierr = NEPGetType(*nep,&tname);if (*ierr) return;
  *ierr = PetscStrncpy(name,tname,len);if (*ierr) return;
  FIXRETURNCHAR(PETSC_TRUE,name,len);
}

SLEPC_EXTERN void nepsetoptionsprefix_(NEP *nep,char *prefix,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(prefix,len,t);
  *ierr = NEPSetOptionsPrefix(*nep,t);if (*ierr) return;
  FREECHAR(prefix,t);
}

SLEPC_EXTERN void nepappendoptionsprefix_(NEP *nep,char *prefix,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(prefix,len,t);
  *ierr = NEPAppendOptionsPrefix(*nep,t);if (*ierr) return;
  FREECHAR(prefix,t);
}

SLEPC_EXTERN void nepgetoptionsprefix_(NEP *nep,char *prefix,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  const char *tname;

  *ierr = NEPGetOptionsPrefix(*nep,&tname); if (*ierr) return;
  *ierr = PetscStrncpy(prefix,tname,len);if (*ierr) return;
  FIXRETURNCHAR(PETSC_TRUE,prefix,len);
}

SLEPC_EXTERN void nepmonitorset_(NEP *nep,void (*monitor)(NEP*,PetscInt*,PetscInt*,PetscScalar*,PetscScalar*,PetscReal*,PetscInt*,void*,PetscErrorCode*),void *mctx,void (*monitordestroy)(void *,PetscErrorCode*),PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(mctx);
  CHKFORTRANNULLFUNCTION(monitordestroy);
  if ((PetscVoidFunction)monitor == (PetscVoidFunction)nepmonitorall_) {
    *ierr = NEPMonitorSet(*nep,(PetscErrorCode (*)(NEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*))NEPMonitorAll,*(PetscViewerAndFormat**)mctx,(PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy);
  } else if ((PetscVoidFunction)monitor == (PetscVoidFunction)nepmonitorconverged_) {
    *ierr = NEPMonitorSet(*nep,(PetscErrorCode (*)(NEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*))NEPMonitorConverged,*(SlepcConvMonitor*)mctx,(PetscErrorCode (*)(void**))SlepcConvMonitorDestroy);
  } else if ((PetscVoidFunction)monitor == (PetscVoidFunction)nepmonitorfirst_) {
    *ierr = NEPMonitorSet(*nep,(PetscErrorCode (*)(NEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*))NEPMonitorFirst,*(PetscViewerAndFormat**)mctx,(PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy);
  } else if ((PetscVoidFunction)monitor == (PetscVoidFunction)nepmonitorlg_) {
    *ierr = NEPMonitorSet(*nep,NEPMonitorLG,0,0);
  } else if ((PetscVoidFunction)monitor == (PetscVoidFunction)nepmonitorlgall_) {
    *ierr = NEPMonitorSet(*nep,NEPMonitorLGAll,0,0);
  } else {
    *ierr = PetscObjectSetFortranCallback((PetscObject)*nep,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.monitor,(PetscVoidFunction)monitor,mctx); if (*ierr) return;
    *ierr = PetscObjectSetFortranCallback((PetscObject)*nep,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.monitordestroy,(PetscVoidFunction)monitordestroy,mctx); if (*ierr) return;
    *ierr = NEPMonitorSet(*nep,ourmonitor,*nep,ourdestroy);
  }
}

SLEPC_EXTERN void nepconvergedabsolute_(NEP *nep,PetscScalar *eigr,PetscScalar *eigi,PetscReal *res,PetscReal *errest,void *ctx,PetscErrorCode *ierr)
{
  *ierr = NEPConvergedAbsolute(*nep,*eigr,*eigi,*res,errest,ctx);
}

SLEPC_EXTERN void nepconvergedrelative_(NEP *nep,PetscScalar *eigr,PetscScalar *eigi,PetscReal *res,PetscReal *errest,void *ctx,PetscErrorCode *ierr)
{
  *ierr = NEPConvergedRelative(*nep,*eigr,*eigi,*res,errest,ctx);
}

SLEPC_EXTERN void nepsetconvergencetestfunction_(NEP *nep,void (*func)(NEP*,PetscScalar*,PetscScalar*,PetscReal*,PetscReal*,void*,PetscErrorCode*),void* ctx,void (*destroy)(void*,PetscErrorCode*),PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(ctx);
  CHKFORTRANNULLFUNCTION(destroy);
  if ((PetscVoidFunction)func == (PetscVoidFunction)nepconvergedabsolute_) {
    *ierr = NEPSetConvergenceTest(*nep,NEP_CONV_ABS);
  } else if ((PetscVoidFunction)func == (PetscVoidFunction)nepconvergedrelative_) {
    *ierr = NEPSetConvergenceTest(*nep,NEP_CONV_REL);
  } else {
    *ierr = PetscObjectSetFortranCallback((PetscObject)*nep,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.convergence,(PetscVoidFunction)func,ctx); if (*ierr) return;
    *ierr = PetscObjectSetFortranCallback((PetscObject)*nep,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.convdestroy,(PetscVoidFunction)destroy,ctx); if (*ierr) return;
    *ierr = NEPSetConvergenceTestFunction(*nep,ourconvergence,*nep,ourconvdestroy);
  }
}

SLEPC_EXTERN void nepstoppingbasic_(NEP *nep,PetscInt *its,PetscInt *max_it,PetscInt *nconv,PetscInt *nev,NEPConvergedReason *reason,void *ctx,PetscErrorCode *ierr)
{
  *ierr = NEPStoppingBasic(*nep,*its,*max_it,*nconv,*nev,reason,ctx);
}

SLEPC_EXTERN void nepsetstoppingtestfunction_(NEP *nep,void (*func)(NEP*,PetscInt,PetscInt,PetscInt,PetscInt,NEPConvergedReason*,void*,PetscErrorCode*),void* ctx,void (*destroy)(void*,PetscErrorCode*),PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(ctx);
  CHKFORTRANNULLFUNCTION(destroy);
  if ((PetscVoidFunction)func == (PetscVoidFunction)nepstoppingbasic_) {
    *ierr = NEPSetStoppingTest(*nep,NEP_STOP_BASIC);
  } else {
    *ierr = PetscObjectSetFortranCallback((PetscObject)*nep,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.stopping,(PetscVoidFunction)func,ctx); if (*ierr) return;
    *ierr = PetscObjectSetFortranCallback((PetscObject)*nep,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.stopdestroy,(PetscVoidFunction)destroy,ctx); if (*ierr) return;
    *ierr = NEPSetStoppingTestFunction(*nep,ourstopping,*nep,ourstopdestroy);
  }
}

SLEPC_EXTERN void nepseteigenvaluecomparison_(NEP *nep,void (*func)(PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt*,void*),void* ctx,PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(ctx);
  *ierr = PetscObjectSetFortranCallback((PetscObject)*nep,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.comparison,(PetscVoidFunction)func,ctx); if (*ierr) return;
  *ierr = NEPSetEigenvalueComparison(*nep,oureigenvaluecomparison,*nep);
}

SLEPC_EXTERN void nepsetfunction_(NEP *nep,Mat *A,Mat *B,void (*func)(NEP*,PetscScalar*,Mat*,Mat*,void*,PetscErrorCode*),void *ctx,PetscErrorCode *ierr PETSC_F90_2PTR_PROTO(ptr))
{
  *ierr = PetscObjectSetFortranCallback((PetscObject)*nep,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.function,(PetscVoidFunction)func,ctx);if (*ierr) return;
#if defined(PETSC_HAVE_F90_2PTR_ARG)
  *ierr = PetscObjectSetFortranCallback((PetscObject)*nep,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.function_pgiptr,NULL,ptr);if (*ierr) return;
#endif
  *ierr = NEPSetFunction(*nep,*A,*B,ournepfunction,NULL);
}

/* func is currently ignored from Fortran */
SLEPC_EXTERN void nepgetfunction_(NEP *nep,Mat *A,Mat *B,void *func,void **ctx,PetscErrorCode *ierr)
{
  CHKFORTRANNULLINTEGER(ctx);
  CHKFORTRANNULLOBJECT(A);
  CHKFORTRANNULLOBJECT(B);
  *ierr = NEPGetFunction(*nep,A,B,NULL,NULL); if (*ierr) return;
  *ierr = PetscObjectGetFortranCallback((PetscObject)*nep,PETSC_FORTRAN_CALLBACK_CLASS,_cb.function,NULL,ctx);
}

SLEPC_EXTERN void nepsetjacobian_(NEP *nep,Mat *J,void (*func)(NEP*,PetscScalar*,Mat*,void*,PetscErrorCode*),void *ctx,PetscErrorCode *ierr PETSC_F90_2PTR_PROTO(ptr))
{
  *ierr = PetscObjectSetFortranCallback((PetscObject)*nep,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.jacobian,(PetscVoidFunction)func,ctx);if (*ierr) return;
#if defined(PETSC_HAVE_F90_2PTR_ARG)
  *ierr = PetscObjectSetFortranCallback((PetscObject)*nep,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.jacobian_pgiptr,NULL,ptr);if (*ierr) return;
#endif
  *ierr = NEPSetJacobian(*nep,*J,ournepjacobian,NULL);
}

/* func is currently ignored from Fortran */
SLEPC_EXTERN void nepgetjacobian_(NEP *nep,Mat *J,void *func,void **ctx,PetscErrorCode *ierr)
{
  CHKFORTRANNULLINTEGER(ctx);
  CHKFORTRANNULLOBJECT(J);
  *ierr = NEPGetJacobian(*nep,J,NULL,NULL); if (*ierr) return;
  *ierr = PetscObjectGetFortranCallback((PetscObject)*nep,PETSC_FORTRAN_CALLBACK_CLASS,_cb.jacobian,NULL,ctx);
}

SLEPC_EXTERN void nepgetdimensions_(NEP *nep,PetscInt *nev,PetscInt *ncv,PetscInt *mpd,PetscErrorCode *ierr)
{
  CHKFORTRANNULLINTEGER(nev);
  CHKFORTRANNULLINTEGER(ncv);
  CHKFORTRANNULLINTEGER(mpd);
  *ierr = NEPGetDimensions(*nep,nev,ncv,mpd);
}

SLEPC_EXTERN void nepgetdimensions000_(NEP *nep,PetscInt *nev,PetscInt *ncv,PetscInt *mpd,PetscErrorCode *ierr)
{
  nepgetdimensions_(nep,nev,ncv,mpd,ierr);
}

SLEPC_EXTERN void nepgetdimensions100_(NEP *nep,PetscInt *nev,PetscInt *ncv,PetscInt *mpd,PetscErrorCode *ierr)
{
  nepgetdimensions_(nep,nev,ncv,mpd,ierr);
}

SLEPC_EXTERN void nepgetdimensions010_(NEP *nep,PetscInt *nev,PetscInt *ncv,PetscInt *mpd,PetscErrorCode *ierr)
{
  nepgetdimensions_(nep,nev,ncv,mpd,ierr);
}

SLEPC_EXTERN void nepgetdimensions001_(NEP *nep,PetscInt *nev,PetscInt *ncv,PetscInt *mpd,PetscErrorCode *ierr)
{
  nepgetdimensions_(nep,nev,ncv,mpd,ierr);
}

SLEPC_EXTERN void nepgetdimensions110_(NEP *nep,PetscInt *nev,PetscInt *ncv,PetscInt *mpd,PetscErrorCode *ierr)
{
  nepgetdimensions_(nep,nev,ncv,mpd,ierr);
}

SLEPC_EXTERN void nepgetdimensions011_(NEP *nep,PetscInt *nev,PetscInt *ncv,PetscInt *mpd,PetscErrorCode *ierr)
{
  nepgetdimensions_(nep,nev,ncv,mpd,ierr);
}

SLEPC_EXTERN void nepgeteigenpair_(NEP *nep,PetscInt *i,PetscScalar *eigr,PetscScalar *eigi,Vec *Vr,Vec *Vi,PetscErrorCode *ierr)
{
  CHKFORTRANNULLSCALAR(eigr);
  CHKFORTRANNULLSCALAR(eigi);
  *ierr = NEPGetEigenpair(*nep,*i,eigr,eigi,*Vr,*Vi);
}

SLEPC_EXTERN void nepgeteigenpair00_(NEP *nep,PetscInt *i,PetscScalar *eigr,PetscScalar *eigi,Vec *Vr,Vec *Vi,PetscErrorCode *ierr)
{
  nepgeteigenpair_(nep,i,eigr,eigi,Vr,Vi,ierr);
}

SLEPC_EXTERN void nepgeteigenpair10_(NEP *nep,PetscInt *i,PetscScalar *eigr,PetscScalar *eigi,Vec *Vr,Vec *Vi,PetscErrorCode *ierr)
{
  nepgeteigenpair_(nep,i,eigr,eigi,Vr,Vi,ierr);
}

SLEPC_EXTERN void nepgeteigenpair01_(NEP *nep,PetscInt *i,PetscScalar *eigr,PetscScalar *eigi,Vec *Vr,Vec *Vi,PetscErrorCode *ierr)
{
  nepgeteigenpair_(nep,i,eigr,eigi,Vr,Vi,ierr);
}

SLEPC_EXTERN void nepgeteigenpair11_(NEP *nep,PetscInt *i,PetscScalar *eigr,PetscScalar *eigi,Vec *Vr,Vec *Vi,PetscErrorCode *ierr)
{
  nepgeteigenpair_(nep,i,eigr,eigi,Vr,Vi,ierr);
}

SLEPC_EXTERN void nepgettolerances_(NEP *nep,PetscReal *tol,PetscInt *maxits,PetscErrorCode *ierr)
{
  CHKFORTRANNULLREAL(tol);
  CHKFORTRANNULLINTEGER(maxits);
  *ierr = NEPGetTolerances(*nep,tol,maxits);
}

SLEPC_EXTERN void nepgettolerances00_(NEP *nep,PetscReal *tol,PetscInt *maxits,PetscErrorCode *ierr)
{
  nepgettolerances_(nep,tol,maxits,ierr);
}

SLEPC_EXTERN void nepgettolerances10_(NEP *nep,PetscReal *tol,PetscInt *maxits,PetscErrorCode *ierr)
{
  nepgettolerances_(nep,tol,maxits,ierr);
}

SLEPC_EXTERN void nepgettolerances01_(NEP *nep,PetscReal *tol,PetscInt *maxits,PetscErrorCode *ierr)
{
  nepgettolerances_(nep,tol,maxits,ierr);
}

SLEPC_EXTERN void nepgetrefine_(NEP *nep,NEPRefine *refine,PetscInt *npart,PetscReal *tol,PetscInt *its,NEPRefineScheme *scheme,PetscErrorCode *ierr)
{
  CHKFORTRANNULLINTEGER(npart);
  CHKFORTRANNULLREAL(tol);
  CHKFORTRANNULLINTEGER(its);
  *ierr = NEPGetRefine(*nep,refine,npart,tol,its,scheme);
}

SLEPC_EXTERN void nepgetrefine000_(NEP *nep,NEPRefine *refine,PetscInt *npart,PetscReal *tol,PetscInt *its,NEPRefineScheme *scheme,PetscErrorCode *ierr)
{
  nepgetrefine_(nep,refine,npart,tol,its,scheme,ierr);
}

SLEPC_EXTERN void nepgetrefine100_(NEP *nep,NEPRefine *refine,PetscInt *npart,PetscReal *tol,PetscInt *its,NEPRefineScheme *scheme,PetscErrorCode *ierr)
{
  nepgetrefine_(nep,refine,npart,tol,its,scheme,ierr);
}

SLEPC_EXTERN void nepgetrefine010_(NEP *nep,NEPRefine *refine,PetscInt *npart,PetscReal *tol,PetscInt *its,NEPRefineScheme *scheme,PetscErrorCode *ierr)
{
  nepgetrefine_(nep,refine,npart,tol,its,scheme,ierr);
}

SLEPC_EXTERN void nepgetrefine001_(NEP *nep,NEPRefine *refine,PetscInt *npart,PetscReal *tol,PetscInt *its,NEPRefineScheme *scheme,PetscErrorCode *ierr)
{
  nepgetrefine_(nep,refine,npart,tol,its,scheme,ierr);
}

SLEPC_EXTERN void nepgetrefine110_(NEP *nep,NEPRefine *refine,PetscInt *npart,PetscReal *tol,PetscInt *its,NEPRefineScheme *scheme,PetscErrorCode *ierr)
{
  nepgetrefine_(nep,refine,npart,tol,its,scheme,ierr);
}

SLEPC_EXTERN void nepgetrefine011_(NEP *nep,NEPRefine *refine,PetscInt *npart,PetscReal *tol,PetscInt *its,NEPRefineScheme *scheme,PetscErrorCode *ierr)
{
  nepgetrefine_(nep,refine,npart,tol,its,scheme,ierr);
}

SLEPC_EXTERN void nepgetrefine101_(NEP *nep,NEPRefine *refine,PetscInt *npart,PetscReal *tol,PetscInt *its,NEPRefineScheme *scheme,PetscErrorCode *ierr)
{
  nepgetrefine_(nep,refine,npart,tol,its,scheme,ierr);
}

SLEPC_EXTERN void nepgetrefine111_(NEP *nep,NEPRefine *refine,PetscInt *npart,PetscReal *tol,PetscInt *its,NEPRefineScheme *scheme,PetscErrorCode *ierr)
{
  nepgetrefine_(nep,refine,npart,tol,its,scheme,ierr);
}

SLEPC_EXTERN void nepsetinitialspace0_(NEP *nep,PetscInt *n,Vec *is,PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(is);
  *ierr = NEPSetInitialSpace(*nep,*n,is);
}

SLEPC_EXTERN void nepsetinitialspace1_(NEP *nep,PetscInt *n,Vec *is,PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(is);
  *ierr = NEPSetInitialSpace(*nep,*n,is);
}

