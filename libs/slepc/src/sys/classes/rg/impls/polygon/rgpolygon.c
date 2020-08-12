/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Polygonal region defined by a set of vertices
*/

#include <slepc/private/rgimpl.h>      /*I "slepcrg.h" I*/
#include <petscdraw.h>

#define VERTMAX 30

typedef struct {
  PetscInt    n;         /* number of vertices */
  PetscScalar *vr,*vi;   /* array of vertices (vi not used in complex scalars) */
} RG_POLYGON;

PetscErrorCode RGComputeBoundingBox_Polygon(RG,PetscReal*,PetscReal*,PetscReal*,PetscReal*);

#if !defined(PETSC_USE_COMPLEX)
static PetscBool CheckSymmetry(PetscInt n,PetscScalar *vr,PetscScalar *vi)
{
  PetscInt i,j,k;
  /* find change of sign in imaginary part */
  j = vi[0]!=0.0? 0: 1;
  for (k=j+1;k<n;k++) {
    if (vi[k]!=0.0) {
      if (vi[k]*vi[j]<0.0) break;
      j++;
    }
  }
  if (k==n) return (j==1)? PETSC_TRUE: PETSC_FALSE;
  /* check pairing vertices */
  for (i=0;i<n/2;i++) {
    if (vr[k]!=vr[j] || vi[k]!=-vi[j]) return PETSC_FALSE;
    k = (k+1)%n;
    j = (j-1+n)%n;
  }
  return PETSC_TRUE;
}
#endif

static PetscErrorCode RGPolygonSetVertices_Polygon(RG rg,PetscInt n,PetscScalar *vr,PetscScalar *vi)
{
  PetscErrorCode ierr;
  PetscInt       i;
  RG_POLYGON     *ctx = (RG_POLYGON*)rg->data;

  PetscFunctionBegin;
  if (n<3) SETERRQ1(PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_OUTOFRANGE,"At least 3 vertices required, you provided %s",n);
  if (n>VERTMAX) SETERRQ1(PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_OUTOFRANGE,"Too many points, maximum allowed is %d",VERTMAX);
#if !defined(PETSC_USE_COMPLEX)
  if (!CheckSymmetry(n,vr,vi)) SETERRQ(PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_WRONG,"In real scalars the region must be symmetric wrt real axis");
#endif
  if (ctx->n) {
    ierr = PetscFree(ctx->vr);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
    ierr = PetscFree(ctx->vi);CHKERRQ(ierr);
#endif
  }
  ctx->n = n;
  ierr = PetscMalloc1(n,&ctx->vr);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
  ierr = PetscMalloc1(n,&ctx->vi);CHKERRQ(ierr);
#endif
  for (i=0;i<n;i++) {
    ctx->vr[i] = vr[i];
#if !defined(PETSC_USE_COMPLEX)
    ctx->vi[i] = vi[i];
#endif
  }
  PetscFunctionReturn(0);
}

/*@
   RGPolygonSetVertices - Sets the vertices that define the polygon region.

   Logically Collective on rg

   Input Parameters:
+  rg - the region context
.  n  - number of vertices
.  vr - array of vertices
-  vi - array of vertices (imaginary part)

   Options Database Keys:
+  -rg_polygon_vertices - Sets the vertices
-  -rg_polygon_verticesi - Sets the vertices (imaginary part)

   Notes:
   In the case of complex scalars, only argument vr is used, containing
   the complex vertices; the list of vertices can be provided in the
   command line with a comma-separated list of complex values
   [+/-][realnumber][+/-]realnumberi with no spaces.

   When PETSc is built with real scalars, the real and imaginary parts of
   the vertices must be provided in two separate arrays (or two lists in
   the command line). In this case, the region must be symmetric with
   respect to the real axis.

   Level: advanced

.seealso: RGPolygonGetVertices()
@*/
PetscErrorCode RGPolygonSetVertices(RG rg,PetscInt n,PetscScalar vr[],PetscScalar vi[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscValidLogicalCollectiveInt(rg,n,2);
  PetscValidScalarPointer(vr,3);
#if !defined(PETSC_USE_COMPLEX)
  PetscValidScalarPointer(vi,4);
#endif
  ierr = PetscTryMethod(rg,"RGPolygonSetVertices_C",(RG,PetscInt,PetscScalar*,PetscScalar*),(rg,n,vr,vi));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode RGPolygonGetVertices_Polygon(RG rg,PetscInt *n,PetscScalar **vr,PetscScalar **vi)
{
  PetscErrorCode ierr;
  RG_POLYGON     *ctx = (RG_POLYGON*)rg->data;
  PetscInt       i;

  PetscFunctionBegin;
  if (n) *n  = ctx->n;
  if (vr) {
    if (!ctx->n) *vr = NULL;
    else {
      ierr = PetscMalloc1(ctx->n,vr);CHKERRQ(ierr);
      for (i=0;i<ctx->n;i++) (*vr)[i] = ctx->vr[i];
    }
  }
#if !defined(PETSC_USE_COMPLEX)
  if (vi) {
    if (!ctx->n) *vi = NULL;
    else {
      ierr = PetscMalloc1(ctx->n,vi);CHKERRQ(ierr);
      for (i=0;i<ctx->n;i++) (*vi)[i] = ctx->vi[i];
    }
  }
#endif
  PetscFunctionReturn(0);
}

/*@C
   RGPolygonGetVertices - Gets the vertices that define the polygon region.

   Not Collective

   Input Parameter:
.  rg     - the region context

   Output Parameters:
+  n  - number of vertices
.  vr - array of vertices
-  vi - array of vertices (imaginary part)

   Notes:
   The values passed by user with RGPolygonSetVertices() are returned (or null
   pointers otherwise).
   The returned arrays should be freed by the user when no longer needed.

   Level: advanced

.seealso: RGPolygonSetVertices()
@*/
PetscErrorCode RGPolygonGetVertices(RG rg,PetscInt *n,PetscScalar **vr,PetscScalar **vi)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  ierr = PetscUseMethod(rg,"RGPolygonGetVertices_C",(RG,PetscInt*,PetscScalar**,PetscScalar**),(rg,n,vr,vi));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode RGView_Polygon(RG rg,PetscViewer viewer)
{
  PetscErrorCode ierr;
  RG_POLYGON     *ctx = (RG_POLYGON*)rg->data;
  PetscBool      isdraw,isascii;
  int            winw,winh;
  PetscDraw      draw;
  PetscDrawAxis  axis;
  PetscReal      a,b,c,d,ab,cd,lx,ly,w,x0,y0,x1,y1,scale=1.2;
  PetscInt       i;
  char           str[50];

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERDRAW,&isdraw);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"  vertices: ");CHKERRQ(ierr);
    ierr = PetscViewerASCIIUseTabs(viewer,PETSC_FALSE);CHKERRQ(ierr);
    for (i=0;i<ctx->n;i++) {
#if defined(PETSC_USE_COMPLEX)
      ierr = SlepcSNPrintfScalar(str,50,ctx->vr[i],PETSC_FALSE);CHKERRQ(ierr);
#else
      if (ctx->vi[i]!=0.0) {
        ierr = PetscSNPrintf(str,50,"%g%+gi",(double)ctx->vr[i],(double)ctx->vi[i]);CHKERRQ(ierr);
      } else {
        ierr = PetscSNPrintf(str,50,"%g",(double)ctx->vr[i]);CHKERRQ(ierr);
      }
#endif
      ierr = PetscViewerASCIIPrintf(viewer,"%s%s",str,(i<ctx->n-1)?", ":"");CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
    ierr = PetscViewerASCIIUseTabs(viewer,PETSC_TRUE);CHKERRQ(ierr);
  } else if (isdraw) {
    ierr = PetscViewerDrawGetDraw(viewer,0,&draw);CHKERRQ(ierr);
    ierr = PetscDrawCheckResizedWindow(draw);CHKERRQ(ierr);
    ierr = PetscDrawGetWindowSize(draw,&winw,&winh);CHKERRQ(ierr);
    winw = PetscMax(winw,1); winh = PetscMax(winh,1);
    ierr = PetscDrawClear(draw);CHKERRQ(ierr);
    ierr = PetscDrawSetTitle(draw,"Polygonal region");CHKERRQ(ierr);
    ierr = PetscDrawAxisCreate(draw,&axis);CHKERRQ(ierr);
    ierr = RGComputeBoundingBox_Polygon(rg,&a,&b,&c,&d);CHKERRQ(ierr);
    a *= rg->sfactor;
    b *= rg->sfactor;
    c *= rg->sfactor;
    d *= rg->sfactor;
    lx = b-a;
    ly = d-c;
    ab = (a+b)/2;
    cd = (c+d)/2;
    w  = scale*PetscMax(lx/winw,ly/winh)/2;
    ierr = PetscDrawAxisSetLimits(axis,ab-w*winw,ab+w*winw,cd-w*winh,cd+w*winh);CHKERRQ(ierr);
    ierr = PetscDrawAxisDraw(axis);CHKERRQ(ierr);
    ierr = PetscDrawAxisDestroy(&axis);CHKERRQ(ierr);
    for (i=0;i<ctx->n;i++) {
#if defined(PETSC_USE_COMPLEX)
      x0 = PetscRealPart(ctx->vr[i]); y0 = PetscImaginaryPart(ctx->vr[i]);
      if (i<ctx->n-1) {
        x1 = PetscRealPart(ctx->vr[i+1]); y1 = PetscImaginaryPart(ctx->vr[i+1]);
      } else {
        x1 = PetscRealPart(ctx->vr[0]); y1 = PetscImaginaryPart(ctx->vr[0]);
      }
#else
      x0 = ctx->vr[i]; y0 = ctx->vi[i];
      if (i<ctx->n-1) {
        x1 = ctx->vr[i+1]; y1 = ctx->vi[i+1];
      } else {
        x1 = ctx->vr[0]; y1 = ctx->vi[0];
      }
#endif
      ierr = PetscDrawLine(draw,x0*rg->sfactor,y0*rg->sfactor,x1*rg->sfactor,y1*rg->sfactor,PETSC_DRAW_MAGENTA);CHKERRQ(ierr);
    }
    ierr = PetscDrawFlush(draw);CHKERRQ(ierr);
    ierr = PetscDrawSave(draw);CHKERRQ(ierr);
    ierr = PetscDrawPause(draw);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode RGIsTrivial_Polygon(RG rg,PetscBool *trivial)
{
  RG_POLYGON *ctx = (RG_POLYGON*)rg->data;

  PetscFunctionBegin;
  *trivial = PetscNot(ctx->n);
  PetscFunctionReturn(0);
}

PetscErrorCode RGComputeContour_Polygon(RG rg,PetscInt n,PetscScalar *cr,PetscScalar *ci)
{
  RG_POLYGON  *ctx = (RG_POLYGON*)rg->data;
  PetscReal   length,h,d,rem=0.0;
  PetscInt    k=1,idx=ctx->n-1,i;
  PetscBool   ini=PETSC_FALSE;
  PetscScalar incr;
#if !defined(PETSC_USE_COMPLEX)
  PetscScalar inci;
#endif

  PetscFunctionBegin;
  if (!ctx->n) SETERRQ(PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_WRONGSTATE,"No vertices have been set yet");
  length = SlepcAbsEigenvalue(ctx->vr[0]-ctx->vr[ctx->n-1],ctx->vi[0]-ctx->vi[ctx->n-1]);
  for (i=0;i<ctx->n-1;i++) length += SlepcAbsEigenvalue(ctx->vr[i]-ctx->vr[i+1],ctx->vi[i]-ctx->vi[i+1]);
  h = length/n;
  cr[0] = ctx->vr[0];
#if !defined(PETSC_USE_COMPLEX)
  ci[0] = ctx->vi[0];
#endif
  incr = ctx->vr[ctx->n-1]-ctx->vr[0];
#if !defined(PETSC_USE_COMPLEX)
  inci = ctx->vi[ctx->n-1]-ctx->vi[0];
#endif
  d = SlepcAbsEigenvalue(incr,inci);
  incr /= d;
#if !defined(PETSC_USE_COMPLEX)
  inci /= d;
#endif
  while (k<n) {
    if (ini) {
      incr = ctx->vr[idx]-ctx->vr[idx+1];
#if !defined(PETSC_USE_COMPLEX)
      inci = ctx->vi[idx]-ctx->vi[idx+1];
#endif
      d = SlepcAbsEigenvalue(incr,inci);
      incr /= d;
#if !defined(PETSC_USE_COMPLEX)
      inci /= d;
#endif
      if (rem+d>h) {
        cr[k] = ctx->vr[idx+1]+incr*(h-rem);
#if !defined(PETSC_USE_COMPLEX)
        ci[k] = ctx->vi[idx+1]+inci*(h-rem);
#endif
        k++;
        ini = PETSC_FALSE;
      } else {rem += d; idx--;}
    } else {
#if !defined(PETSC_USE_COMPLEX)
      rem = SlepcAbsEigenvalue(ctx->vr[idx]-cr[k-1],ctx->vi[idx]-ci[k-1]);
#else
      rem = PetscAbsScalar(ctx->vr[idx]-cr[k-1]);
#endif
      if (rem>h) {
        cr[k] = cr[k-1]+incr*h;
#if !defined(PETSC_USE_COMPLEX)
        ci[k] = ci[k-1]+inci*h;
#endif
        k++;
      } else {ini = PETSC_TRUE; idx--;}
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode RGComputeBoundingBox_Polygon(RG rg,PetscReal *a,PetscReal *b,PetscReal *c,PetscReal *d)
{
  RG_POLYGON *ctx = (RG_POLYGON*)rg->data;
  PetscInt   i;

  PetscFunctionBegin;
  *a =  PETSC_MAX_REAL;
  *b = -PETSC_MAX_REAL;
  *c =  PETSC_MAX_REAL;
  *d = -PETSC_MAX_REAL;
  for (i=0;i<ctx->n;i++) {
#if defined(PETSC_USE_COMPLEX)
    if (a) *a = PetscMin(*a,PetscRealPart(ctx->vr[i]));
    if (b) *b = PetscMax(*b,PetscRealPart(ctx->vr[i]));
    if (c) *c = PetscMin(*c,PetscImaginaryPart(ctx->vr[i]));
    if (d) *d = PetscMax(*d,PetscImaginaryPart(ctx->vr[i]));
#else
    if (a) *a = PetscMin(*a,ctx->vr[i]);
    if (b) *b = PetscMax(*b,ctx->vr[i]);
    if (c) *c = PetscMin(*c,ctx->vi[i]);
    if (d) *d = PetscMax(*d,ctx->vi[i]);
#endif
  }
  PetscFunctionReturn(0);
}

PetscErrorCode RGCheckInside_Polygon(RG rg,PetscReal px,PetscReal py,PetscInt *inout)
{
  RG_POLYGON *ctx = (RG_POLYGON*)rg->data;
  PetscReal  val,x[VERTMAX],y[VERTMAX];
  PetscBool  mx,my,nx,ny;
  PetscInt   i,j;

  PetscFunctionBegin;
  for (i=0;i<ctx->n;i++) {
#if defined(PETSC_USE_COMPLEX)
    x[i] = PetscRealPart(ctx->vr[i])-px;
    y[i] = PetscImaginaryPart(ctx->vr[i])-py;
#else
    x[i] = ctx->vr[i]-px;
    y[i] = ctx->vi[i]-py;
#endif
  }
  *inout = -1;
  for (i=0;i<ctx->n;i++) {
    j = (i+1)%ctx->n;
    mx = PetscNot(x[i]<0.0);
    nx = PetscNot(x[j]<0.0);
    my = PetscNot(y[i]<0.0);
    ny = PetscNot(y[j]<0.0);
    if (!((my||ny) && (mx||nx)) || (mx&&nx)) continue;
    if (((my && ny && (mx||nx)) && (!(mx&&nx)))) {
      *inout = -*inout;
      continue;
    }
    val = (y[i]*x[j]-x[i]*y[j])/(x[j]-x[i]);
    if (PetscAbs(val)<10*PETSC_MACHINE_EPSILON) {
      *inout = 0;
      PetscFunctionReturn(0);
    } else if (val>0.0) *inout = -*inout;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode RGSetFromOptions_Polygon(PetscOptionItems *PetscOptionsObject,RG rg)
{
  PetscErrorCode ierr;
  PetscScalar    array[VERTMAX];
  PetscInt       i,k;
  PetscBool      flg,flgi=PETSC_FALSE;
#if !defined(PETSC_USE_COMPLEX)
  PetscScalar    arrayi[VERTMAX];
  PetscInt       ki;
#else
  PetscScalar    *arrayi=NULL;
#endif

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"RG Polygon Options");CHKERRQ(ierr);

    k = VERTMAX;
    for (i=0;i<k;i++) array[i] = 0;
    ierr = PetscOptionsScalarArray("-rg_polygon_vertices","Vertices of polygon","RGPolygonSetVertices",array,&k,&flg);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
    ki = VERTMAX;
    for (i=0;i<ki;i++) arrayi[i] = 0;
    ierr = PetscOptionsScalarArray("-rg_polygon_verticesi","Vertices of polygon (imaginary part)","RGPolygonSetVertices",arrayi,&ki,&flgi);CHKERRQ(ierr);
    if (ki!=k) SETERRQ2(PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_SIZ,"The number of real %D and imaginary %D parts do not match",k,ki);
#endif
    if (flg || flgi) { ierr = RGPolygonSetVertices(rg,k,array,arrayi);CHKERRQ(ierr); }

  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode RGDestroy_Polygon(RG rg)
{
  PetscErrorCode ierr;
  RG_POLYGON     *ctx = (RG_POLYGON*)rg->data;

  PetscFunctionBegin;
  if (ctx->n) {
    ierr = PetscFree(ctx->vr);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
    ierr = PetscFree(ctx->vi);CHKERRQ(ierr);
#endif
  }
  ierr = PetscFree(rg->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)rg,"RGPolygonSetVertices_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)rg,"RGPolygonGetVertices_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode RGCreate_Polygon(RG rg)
{
  RG_POLYGON     *polygon;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNewLog(rg,&polygon);CHKERRQ(ierr);
  rg->data = (void*)polygon;

  rg->ops->istrivial      = RGIsTrivial_Polygon;
  rg->ops->computecontour = RGComputeContour_Polygon;
  rg->ops->computebbox    = RGComputeBoundingBox_Polygon;
  rg->ops->checkinside    = RGCheckInside_Polygon;
  rg->ops->setfromoptions = RGSetFromOptions_Polygon;
  rg->ops->view           = RGView_Polygon;
  rg->ops->destroy        = RGDestroy_Polygon;
  ierr = PetscObjectComposeFunction((PetscObject)rg,"RGPolygonSetVertices_C",RGPolygonSetVertices_Polygon);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)rg,"RGPolygonGetVertices_C",RGPolygonGetVertices_Polygon);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

