/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test the ring region.\n\n";

#include <slepcrg.h>

#define NPOINTS 11

PetscErrorCode CheckPoint(RG rg,PetscReal re,PetscReal im)
{
  PetscErrorCode ierr;
  PetscInt       inside;
  PetscScalar    ar,ai;

  PetscFunctionBeginUser;
#if defined(PETSC_USE_COMPLEX)
  ar = PetscCMPLX(re,im);
#else
  ar = re; ai = im;
#endif
  ierr = RGCheckInside(rg,1,&ar,&ai,&inside);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Point (%g,%g) is %s the region\n",(double)re,(double)im,(inside>=0)?"inside":"outside");
  PetscFunctionReturn(0);
}

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  RG             rg;
  RGType         rtype;
  PetscInt       i;
  PetscBool      triv;
  PetscReal      re,im,radius,vscale,start_ang,end_ang,width;
  PetscScalar    center,cr[NPOINTS],ci[NPOINTS];

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = RGCreate(PETSC_COMM_WORLD,&rg);CHKERRQ(ierr);

  ierr = RGSetType(rg,RGRING);CHKERRQ(ierr);
  ierr = RGIsTrivial(rg,&triv);CHKERRQ(ierr);
  if (!triv) SETERRQ(PETSC_COMM_WORLD,1,"Region should be trivial before setting parameters");
  ierr = RGRingSetParameters(rg,2,PETSC_DEFAULT,0.5,0.25,0.75,0.1);CHKERRQ(ierr);
  ierr = RGSetFromOptions(rg);CHKERRQ(ierr);
  ierr = RGIsTrivial(rg,&triv);CHKERRQ(ierr);
  if (triv) SETERRQ(PETSC_COMM_WORLD,1,"Region should be non-trivial after setting parameters");
  ierr = RGView(rg,NULL);CHKERRQ(ierr);
  ierr = RGViewFromOptions(rg,NULL,"-rg_view");CHKERRQ(ierr);

  ierr = RGGetType(rg,&rtype);CHKERRQ(ierr);
  ierr = RGRingGetParameters(rg,&center,&radius,&vscale,&start_ang,&end_ang,&width);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%s region: \n  center=%g, radius=%g, vscale=%g\n  start angle=%g, end angle=%g, width=%g\n\n",rtype,(double)PetscRealPart(center),(double)radius,(double)vscale,(double)start_ang,(double)end_ang,(double)width);

  ierr = CheckPoint(rg,3.0,0.3);CHKERRQ(ierr);
  ierr = CheckPoint(rg,1.1747,0.28253);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nContour points: ");
  ierr = RGComputeContour(rg,NPOINTS,cr,ci);CHKERRQ(ierr);
  for (i=0;i<NPOINTS;i++) {
#if defined(PETSC_USE_COMPLEX)
    re = PetscRealPart(cr[i]);
    im = PetscImaginaryPart(cr[i]);
#else
    re = cr[i];
    im = ci[i];
#endif
    ierr = PetscPrintf(PETSC_COMM_WORLD,"(%.3g,%.3g) ",(double)re,(double)im);
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");

  ierr = RGDestroy(&rg);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      args: -rg_ring_width 0.015

   test:
      suffix: 2
      args: -rg_ring_width 0.015 -rg_scale 1.5

   test:
      suffix: 3
      args: -rg_view draw:tikz:test2_3_ring.tikz
      filter: cat - test2_3_ring.tikz
      requires: !single

TEST*/
