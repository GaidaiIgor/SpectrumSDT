/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Ring region, similar to the ellipse but with a start and end angle,
   together with the width
*/

#include <slepc/private/rgimpl.h>      /*I "slepcrg.h" I*/
#include <petscdraw.h>

typedef struct {
  PetscScalar center;     /* center of the ellipse */
  PetscReal   radius;     /* radius of the ellipse */
  PetscReal   vscale;     /* vertical scale of the ellipse */
  PetscReal   start_ang;  /* start angle */
  PetscReal   end_ang;    /* end angle */
  PetscReal   width;      /* ring width */
} RG_RING;

static PetscErrorCode RGRingSetParameters_Ring(RG rg,PetscScalar center,PetscReal radius,PetscReal vscale,PetscReal start_ang,PetscReal end_ang,PetscReal width)
{
  RG_RING *ctx = (RG_RING*)rg->data;

  PetscFunctionBegin;
  ctx->center = center;
  if (radius == PETSC_DEFAULT) {
    ctx->radius = 1.0;
  } else {
    if (radius<=0.0) SETERRQ(PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_OUTOFRANGE,"The radius argument must be > 0.0");
    ctx->radius = radius;
  }
  if (vscale<=0.0) SETERRQ(PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_OUTOFRANGE,"The vscale argument must be > 0.0");
  ctx->vscale = vscale;
  if (start_ang == PETSC_DEFAULT) {
    ctx->start_ang = 0.0;
  } else {
    if (start_ang<0.0) SETERRQ(PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_OUTOFRANGE,"The right-hand side angle argument must be >= 0.0");
    if (start_ang>1.0) SETERRQ(PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_OUTOFRANGE,"The right-hand side angle argument must be <= 1.0");
    ctx->start_ang = start_ang;
  }
  if (end_ang == PETSC_DEFAULT) {
    ctx->end_ang = 1.0;
  } else {
    if (end_ang<0.0) SETERRQ(PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_OUTOFRANGE,"The left-hand side angle argument must be >= 0.0");
    if (end_ang>1.0) SETERRQ(PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_OUTOFRANGE,"The left-hand side angle argument must be <= 1.0");
    ctx->end_ang = end_ang;
  }
#if !defined(PETSC_USE_COMPLEX)
  if (ctx->start_ang+ctx->end_ang!=1.0) SETERRQ(PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_WRONG,"In real scalars the region must be symmetric wrt real axis");
#endif
  if (width == PETSC_DEFAULT) {
    ctx->width = 0.1;
  } else {
    if (width<=0.0) SETERRQ(PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_OUTOFRANGE,"The width argument must be > 0.0");
    ctx->width = width;
  }
  PetscFunctionReturn(0);
}

/*@
   RGRingSetParameters - Sets the parameters defining the ring region.

   Logically Collective on rg

   Input Parameters:
+  rg        - the region context
.  center    - center of the ellipse
.  radius    - radius of the ellipse
.  vscale    - vertical scale of the ellipse
.  start_ang - the right-hand side angle
.  end_ang   - the left-hand side angle
-  width     - width of the ring

   Options Database Keys:
+  -rg_ring_center     - Sets the center
.  -rg_ring_radius     - Sets the radius
.  -rg_ring_vscale     - Sets the vertical scale
.  -rg_ring_startangle - Sets the right-hand side angle
.  -rg_ring_endangle   - Sets the left-hand side angle
-  -rg_ring_width      - Sets the width of the ring

   Notes:
   The values of center, radius and vscale have the same meaning as in the
   ellipse region. The startangle and endangle define the span of the ring
   (by default it is the whole ring), while the width is the separation
   between the two concentric ellipses (above and below the radius by
   width/2).

   The start and end angles are expressed as a fraction of the circumference:
   the allowed range is [0..1], with 0 corresponding to 0 radians, 0.25 to
   pi/2 radians, and so on. It is allowed to have startangle>endangle, in
   which case the ring region crosses over the zero angle.

   In the case of complex scalars, a complex center can be provided in the
   command line with [+/-][realnumber][+/-]realnumberi with no spaces, e.g.
   -rg_ring_center 1.0+2.0i

   When PETSc is built with real scalars, the center is restricted to a real value,
   and the start and end angles must be such that the region is symmetric with
   respect to the real axis.

   Level: advanced

.seealso: RGRingGetParameters()
@*/
PetscErrorCode RGRingSetParameters(RG rg,PetscScalar center,PetscReal radius,PetscReal vscale,PetscReal start_ang,PetscReal end_ang,PetscReal width)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscValidLogicalCollectiveScalar(rg,center,2);
  PetscValidLogicalCollectiveReal(rg,radius,3);
  PetscValidLogicalCollectiveReal(rg,vscale,4);
  PetscValidLogicalCollectiveReal(rg,start_ang,5);
  PetscValidLogicalCollectiveReal(rg,end_ang,6);
  PetscValidLogicalCollectiveReal(rg,width,7);
  ierr = PetscTryMethod(rg,"RGRingSetParameters_C",(RG,PetscScalar,PetscReal,PetscReal,PetscReal,PetscReal,PetscReal),(rg,center,radius,vscale,start_ang,end_ang,width));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode RGRingGetParameters_Ring(RG rg,PetscScalar *center,PetscReal *radius,PetscReal *vscale,PetscReal *start_ang,PetscReal *end_ang,PetscReal *width)
{
  RG_RING *ctx = (RG_RING*)rg->data;

  PetscFunctionBegin;
  if (center)    *center    = ctx->center;
  if (radius)    *radius    = ctx->radius;
  if (vscale)    *vscale    = ctx->vscale;
  if (start_ang) *start_ang = ctx->start_ang;
  if (end_ang)   *end_ang   = ctx->end_ang;
  if (width)     *width     = ctx->width;
  PetscFunctionReturn(0);
}

/*@
   RGRingGetParameters - Gets the parameters that define the ring region.

   Not Collective

   Input Parameter:
.  rg     - the region context

   Output Parameters:
+  center    - center of the region
.  radius    - radius of the region
.  vscale    - vertical scale of the region
.  start_ang - the right-hand side angle
.  end_ang   - the left-hand side angle
-  width     - width of the ring

   Level: advanced

.seealso: RGRingSetParameters()
@*/
PetscErrorCode RGRingGetParameters(RG rg,PetscScalar *center,PetscReal *radius,PetscReal *vscale,PetscReal *start_ang,PetscReal *end_ang,PetscReal *width)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  ierr = PetscUseMethod(rg,"RGRingGetParameters_C",(RG,PetscScalar*,PetscReal*,PetscReal*,PetscReal*,PetscReal*,PetscReal*),(rg,center,radius,vscale,start_ang,end_ang,width));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode RGView_Ring(RG rg,PetscViewer viewer)
{
  PetscErrorCode ierr;
  RG_RING        *ctx = (RG_RING*)rg->data;
  int            winw,winh;
  PetscBool      isdraw,isascii;
  PetscDraw      draw;
  PetscDrawAxis  axis;
  PetscReal      cx,cy,radius,width,ab,cd,lx,ly,w,end_ang,x1,y1,x2,y2,r,theta,scale=1.2;
  char           str[50];

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERDRAW,&isdraw);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    ierr = SlepcSNPrintfScalar(str,50,ctx->center,PETSC_FALSE);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  center: %s, radius: %g, vscale: %g, start angle: %g, end angle: %g, ring width: %g\n",str,RGShowReal(ctx->radius),RGShowReal(ctx->vscale),(double)ctx->start_ang,(double)ctx->end_ang,(double)ctx->width);CHKERRQ(ierr);
  } else if (isdraw) {
    ierr = PetscViewerDrawGetDraw(viewer,0,&draw);CHKERRQ(ierr);
    ierr = PetscDrawCheckResizedWindow(draw);CHKERRQ(ierr);
    ierr = PetscDrawGetWindowSize(draw,&winw,&winh);CHKERRQ(ierr);
    winw = PetscMax(winw,1); winh = PetscMax(winh,1);
    ierr = PetscDrawClear(draw);CHKERRQ(ierr);
    ierr = PetscDrawSetTitle(draw,"Ring region");CHKERRQ(ierr);
    ierr = PetscDrawAxisCreate(draw,&axis);CHKERRQ(ierr);
    cx = PetscRealPart(ctx->center)*rg->sfactor;
    cy = PetscImaginaryPart(ctx->center)*rg->sfactor;
    radius = ctx->radius*rg->sfactor;
    width  = ctx->width*rg->sfactor;
    lx = 2*(radius+width);
    ly = 2*(radius+width)*ctx->vscale;
    ab = cx;
    cd = cy;
    w  = scale*PetscMax(lx/winw,ly/winh)/2;
    ierr = PetscDrawAxisSetLimits(axis,ab-w*winw,ab+w*winw,cd-w*winh,cd+w*winh);CHKERRQ(ierr);
    ierr = PetscDrawAxisDraw(axis);CHKERRQ(ierr);
    ierr = PetscDrawAxisDestroy(&axis);CHKERRQ(ierr);
    /* draw outer ellipse */
    ierr = PetscDrawEllipse(draw,cx,cy,2*(radius+width),2*(radius+width)*ctx->vscale,PETSC_DRAW_ORANGE);CHKERRQ(ierr);
    /* remove inner part */
    ierr = PetscDrawEllipse(draw,cx,cy,2*(radius-width),2*(radius-width)*ctx->vscale,PETSC_DRAW_WHITE);CHKERRQ(ierr);
    if (ctx->start_ang!=ctx->end_ang) {
      /* remove section from end_ang to start_ang */
      end_ang = (ctx->start_ang<ctx->end_ang)? ctx->end_ang-1: ctx->end_ang;
      theta = end_ang;
      r = scale*(radius+width);
      if (ctx->vscale>1) r *= ctx->vscale;
      x1 = PetscMin(PetscMax(ab+r*PetscCosReal(2.0*PETSC_PI*theta),ab-w*winw),ab+w*winw);
      y1 = PetscMin(PetscMax(cd+r*PetscSinReal(2.0*PETSC_PI*theta),cd-w*winh),cd+w*winh);
      do {
        theta = PetscMin(PetscFloorReal(8*theta+1)/8,ctx->start_ang);
        x2 = PetscMin(PetscMax(ab+r*PetscCosReal(2.0*PETSC_PI*theta),ab-w*winw),ab+w*winw);
        y2 = PetscMin(PetscMax(cd+r*PetscSinReal(2.0*PETSC_PI*theta),cd-w*winh),cd+w*winh);
        ierr = PetscDrawTriangle(draw,cx,cy,x1,y1,x2,y2,PETSC_DRAW_WHITE,PETSC_DRAW_WHITE,PETSC_DRAW_WHITE);CHKERRQ(ierr);
        x1 = x2; y1 = y2;
      } while (theta<ctx->start_ang);
    }
    ierr = PetscDrawFlush(draw);CHKERRQ(ierr);
    ierr = PetscDrawSave(draw);CHKERRQ(ierr);
    ierr = PetscDrawPause(draw);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode RGIsTrivial_Ring(RG rg,PetscBool *trivial)
{
  RG_RING *ctx = (RG_RING*)rg->data;

  PetscFunctionBegin;
  if (rg->complement) *trivial = PetscNot(ctx->radius);
  else *trivial = PetscNot(ctx->radius<PETSC_MAX_REAL);
  PetscFunctionReturn(0);
}

PetscErrorCode RGComputeContour_Ring(RG rg,PetscInt n,PetscScalar *cr,PetscScalar *ci)
{
  RG_RING   *ctx = (RG_RING*)rg->data;
  PetscReal theta,start_ang;
  PetscInt  i,n2=n/2;

  PetscFunctionBegin;
  start_ang = (ctx->start_ang>ctx->end_ang)? ctx->start_ang-1: ctx->start_ang;
  for (i=0;i<n;i++) {
    if (i < n2) {
      theta = ((ctx->end_ang-start_ang)*i/n2 + start_ang)*2.0*PETSC_PI;
#if defined(PETSC_USE_COMPLEX)
      cr[i] = ctx->center + (ctx->radius+ctx->width/2.0)*PetscCMPLX(PetscCosReal(theta),ctx->vscale*PetscSinReal(theta));
#else
      cr[i] = ctx->center + (ctx->radius+ctx->width/2.0)*PetscCosReal(theta);
      ci[i] = (ctx->radius+ctx->width/2.0)*ctx->vscale*PetscSinReal(theta);
#endif
    } else {
      theta = ((ctx->end_ang-start_ang)*(n-i)/n2 + start_ang)*2.0*PETSC_PI;
#if defined(PETSC_USE_COMPLEX)
      cr[i] = ctx->center + (ctx->radius-ctx->width/2.0)*PetscCMPLX(PetscCosReal(theta),ctx->vscale*PetscSinReal(theta));
#else
      cr[i] = ctx->center + (ctx->radius-ctx->width/2.0)*PetscCosReal(theta);
      ci[i] = (ctx->radius-ctx->width/2.0)*ctx->vscale*PetscSinReal(theta);
#endif
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode RGComputeBoundingBox_Ring(RG rg,PetscReal *a,PetscReal *b,PetscReal *c,PetscReal *d)
{
  RG_RING *ctx = (RG_RING*)rg->data;

  PetscFunctionBegin;
  /* current implementation does not return a tight bounding box */
  if (a) *a = PetscRealPart(ctx->center) - (ctx->radius+ctx->width/2.0);
  if (b) *b = PetscRealPart(ctx->center) + (ctx->radius+ctx->width/2.0);
  if (c) *c = PetscImaginaryPart(ctx->center) - (ctx->radius+ctx->width/2.0)*ctx->vscale;
  if (d) *d = PetscImaginaryPart(ctx->center) + (ctx->radius+ctx->width/2.0)*ctx->vscale;
  PetscFunctionReturn(0);
}

PetscErrorCode RGCheckInside_Ring(RG rg,PetscReal px,PetscReal py,PetscInt *inside)
{
  RG_RING   *ctx = (RG_RING*)rg->data;
  PetscReal dx,dy,r;

  PetscFunctionBegin;
  /* outer ellipse */
#if defined(PETSC_USE_COMPLEX)
  dx = (px-PetscRealPart(ctx->center))/(ctx->radius+ctx->width/2.0);
  dy = (py-PetscImaginaryPart(ctx->center))/(ctx->radius+ctx->width/2.0);
#else
  dx = (px-ctx->center)/(ctx->radius+ctx->width/2.0);
  dy = py/(ctx->radius+ctx->width/2.0);
#endif
  r = 1.0-dx*dx-(dy*dy)/(ctx->vscale*ctx->vscale);
  *inside = PetscSign(r);
  /* inner ellipse */
#if defined(PETSC_USE_COMPLEX)
  dx = (px-PetscRealPart(ctx->center))/(ctx->radius-ctx->width/2.0);
  dy = (py-PetscImaginaryPart(ctx->center))/(ctx->radius-ctx->width/2.0);
#else
  dx = (px-ctx->center)/(ctx->radius-ctx->width/2.0);
  dy = py/(ctx->radius-ctx->width/2.0);
#endif
  r = -1.0+dx*dx+(dy*dy)/(ctx->vscale*ctx->vscale);
  *inside *= PetscSign(r);
  if (*inside == 1) {  /* check angles */
#if defined(PETSC_USE_COMPLEX)
    dx = (px-PetscRealPart(ctx->center));
    dy = (py-PetscImaginaryPart(ctx->center));
#else
    dx = px-ctx->center;
    dy = py;
#endif
    if (dx == 0) {
      if (dy == 0) r = -1;
      else if (dy > 0) r = 0.25;
      else r = 0.75;
    } else if (dx > 0) {
      r = PetscAtanReal((dy/ctx->vscale)/dx);
      if (dy >= 0) r /= 2*PETSC_PI;
      else r = r/(2*PETSC_PI)+1;
    } else r = PetscAtanReal((dy/ctx->vscale)/dx)/(2*PETSC_PI)+0.5;
    if (ctx->start_ang>ctx->end_ang) {
      if (r>ctx->end_ang && r<ctx->start_ang) *inside = -1;
    } else {
      if (r<ctx->start_ang || r>ctx->end_ang) *inside = -1;
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode RGSetFromOptions_Ring(PetscOptionItems *PetscOptionsObject,RG rg)
{
  PetscErrorCode ierr;
  PetscScalar    s;
  PetscReal      r1,r2,r3,r4,r5;
  PetscBool      flg1,flg2,flg3,flg4,flg5,flg6;

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"RG Ring Options");CHKERRQ(ierr);

    ierr = RGRingGetParameters(rg,&s,&r1,&r2,&r3,&r4,&r5);CHKERRQ(ierr);
    ierr = PetscOptionsScalar("-rg_ring_center","Center of ellipse","RGRingSetParameters",s,&s,&flg1);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-rg_ring_radius","Radius of ellipse","RGRingSetParameters",r1,&r1,&flg2);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-rg_ring_vscale","Vertical scale of ellipse","RGRingSetParameters",r2,&r2,&flg3);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-rg_ring_startangle","Right-hand side angle","RGRingSetParameters",r3,&r3,&flg4);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-rg_ring_endangle","Left-hand side angle","RGRingSetParameters",r4,&r4,&flg5);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-rg_ring_width","Width of ring","RGRingSetParameters",r5,&r5,&flg6);CHKERRQ(ierr);
    if (flg1 || flg2 || flg3 || flg4 || flg5 || flg6) { ierr = RGRingSetParameters(rg,s,r1,r2,r3,r4,r5);CHKERRQ(ierr); }

  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode RGDestroy_Ring(RG rg)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFree(rg->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)rg,"RGRingSetParameters_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)rg,"RGRingGetParameters_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode RGCreate_Ring(RG rg)
{
  RG_RING        *ring;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNewLog(rg,&ring);CHKERRQ(ierr);
  ring->center    = 0.0;
  ring->radius    = PETSC_MAX_REAL;
  ring->vscale    = 1.0;
  ring->start_ang = 0.0;
  ring->end_ang   = 1.0;
  ring->width     = 0.1;
  rg->data = (void*)ring;

  rg->ops->istrivial      = RGIsTrivial_Ring;
  rg->ops->computecontour = RGComputeContour_Ring;
  rg->ops->computebbox    = RGComputeBoundingBox_Ring;
  rg->ops->checkinside    = RGCheckInside_Ring;
  rg->ops->setfromoptions = RGSetFromOptions_Ring;
  rg->ops->view           = RGView_Ring;
  rg->ops->destroy        = RGDestroy_Ring;
  ierr = PetscObjectComposeFunction((PetscObject)rg,"RGRingSetParameters_C",RGRingSetParameters_Ring);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)rg,"RGRingGetParameters_C",RGRingGetParameters_Ring);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

