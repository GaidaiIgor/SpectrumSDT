/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test VecComp.\n\n";

#include <slepcsys.h>

int main(int argc,char **argv)
{
  Vec            v,w,x,y,vc,wc,xc,yc,vparent,vchild[2],vecs[2];
  const Vec      *varray;
  PetscMPIInt    size,rank;
  PetscInt       i,n,k,Nx[2];
  PetscReal      norm,normc,norm12[2],norm12c[2],vmax,vmin;
  PetscScalar    dot[2],dotc[2];
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  if (size > 2) SETERRQ(PETSC_COMM_WORLD,1,"This test needs one or two processes");
  ierr = PetscPrintf(PETSC_COMM_WORLD,"VecComp test\n");CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create standard vectors
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = VecCreate(PETSC_COMM_WORLD,&v);CHKERRQ(ierr);
  ierr = VecSetSizes(v,8/size,8);CHKERRQ(ierr);
  ierr = VecSetFromOptions(v);CHKERRQ(ierr);

  if (!rank) {
    ierr = VecSetValue(v,0,2.0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(v,1,-1.0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(v,2,3.0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(v,3,3.5,INSERT_VALUES);CHKERRQ(ierr);
  }
  if ((!rank && size==1) || (rank && size==2)) {
    ierr = VecSetValue(v,4,1.2,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(v,5,1.8,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(v,6,-2.2,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(v,7,2.0,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = VecAssemblyBegin(v);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(v);CHKERRQ(ierr);
  ierr = VecDuplicate(v,&w);CHKERRQ(ierr);
  ierr = VecSet(w,1.0);CHKERRQ(ierr);
  ierr = VecDuplicate(v,&x);CHKERRQ(ierr);
  ierr = VecDuplicate(v,&y);CHKERRQ(ierr);
  if (!rank) {
    ierr = VecSetValue(y,0,1.0,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = VecAssemblyBegin(y);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(y);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create veccomp vectors
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = VecCreate(PETSC_COMM_WORLD,&vparent);CHKERRQ(ierr);
  ierr = VecSetSizes(vparent,4/size,4);CHKERRQ(ierr);
  ierr = VecSetFromOptions(vparent);CHKERRQ(ierr);

  /* create a veccomp vector with two subvectors */
  ierr = VecDuplicate(vparent,&vchild[0]);CHKERRQ(ierr);
  ierr = VecDuplicate(vparent,&vchild[1]);CHKERRQ(ierr);
  if (!rank) {
    ierr = VecSetValue(vchild[0],0,2.0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vchild[0],1,-1.0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vchild[1],0,1.2,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vchild[1],1,1.8,INSERT_VALUES);CHKERRQ(ierr);
  }
  if ((!rank && size==1) || (rank && size==2)) {
    ierr = VecSetValue(vchild[0],2,3.0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vchild[0],3,3.5,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vchild[1],2,-2.2,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vchild[1],3,2.0,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = VecAssemblyBegin(vchild[0]);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(vchild[1]);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(vchild[0]);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(vchild[1]);CHKERRQ(ierr);
  ierr = VecCreateCompWithVecs(vchild,2,vparent,&vc);CHKERRQ(ierr);
  ierr = VecDestroy(&vchild[0]);CHKERRQ(ierr);
  ierr = VecDestroy(&vchild[1]);CHKERRQ(ierr);
  ierr = VecView(vc,NULL);CHKERRQ(ierr);

  ierr = VecGetSize(vc,&k);CHKERRQ(ierr);
  if (k!=8) SETERRQ(PETSC_COMM_WORLD,1,"Vector global length should be 8");

  /* create an empty veccomp vector with two subvectors */
  Nx[0] = 4;
  Nx[1] = 4;
  ierr = VecCreateComp(PETSC_COMM_WORLD,Nx,2,VECSTANDARD,vparent,&wc);CHKERRQ(ierr);
  ierr = VecCompGetSubVecs(wc,&n,&varray);CHKERRQ(ierr);
  if (n!=2) SETERRQ(PETSC_COMM_WORLD,1,"n should be 2");
  for (i=0;i<2;i++) {
    ierr = VecSet(varray[i],1.0);CHKERRQ(ierr);
  }

  ierr = VecGetSize(wc,&k);CHKERRQ(ierr);
  if (k!=8) SETERRQ(PETSC_COMM_WORLD,1,"Vector global length should be 8");

  /* duplicate a veccomp */
  ierr = VecDuplicate(vc,&xc);CHKERRQ(ierr);

  /* create a veccomp via VecSetType */
  ierr = VecCreate(PETSC_COMM_WORLD,&yc);CHKERRQ(ierr);
  ierr = VecSetType(yc,VECCOMP);CHKERRQ(ierr);
  ierr = VecSetSizes(yc,8/size,8);CHKERRQ(ierr);
  ierr = VecCompSetSubVecs(yc,2,NULL);CHKERRQ(ierr);

  ierr = VecCompGetSubVecs(yc,&n,&varray);CHKERRQ(ierr);
  if (n!=2) SETERRQ(PETSC_COMM_WORLD,1,"n should be 2");
  if (!rank) {
    ierr = VecSetValue(varray[0],0,1.0,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = VecAssemblyBegin(varray[0]);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(varray[0]);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Operate with vectors
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = VecCopy(w,x);CHKERRQ(ierr);
  ierr = VecAXPBY(x,1.0,-2.0,v);CHKERRQ(ierr);
  ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
  ierr = VecCopy(wc,xc);CHKERRQ(ierr);
  ierr = VecAXPBY(xc,1.0,-2.0,vc);CHKERRQ(ierr);
  ierr = VecNorm(xc,NORM_2,&normc);CHKERRQ(ierr);
  if (PetscAbsReal(norm-normc)>10*PETSC_MACHINE_EPSILON) SETERRQ(PETSC_COMM_WORLD,1,"Norms are different");

  ierr = VecCopy(w,x);CHKERRQ(ierr);
  ierr = VecWAXPY(x,-2.0,w,v);CHKERRQ(ierr);
  ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
  ierr = VecCopy(wc,xc);CHKERRQ(ierr);
  ierr = VecWAXPY(xc,-2.0,wc,vc);CHKERRQ(ierr);
  ierr = VecNorm(xc,NORM_2,&normc);CHKERRQ(ierr);
  if (PetscAbsReal(norm-normc)>10*PETSC_MACHINE_EPSILON) SETERRQ(PETSC_COMM_WORLD,1,"Norms are different");

  ierr = VecAXPBYPCZ(y,3.0,-1.0,1.0,w,v);CHKERRQ(ierr);
  ierr = VecNorm(y,NORM_2,&norm);CHKERRQ(ierr);
  ierr = VecAXPBYPCZ(yc,3.0,-1.0,1.0,wc,vc);CHKERRQ(ierr);
  ierr = VecNorm(yc,NORM_2,&normc);CHKERRQ(ierr);
  if (PetscAbsReal(norm-normc)>10*PETSC_MACHINE_EPSILON) SETERRQ(PETSC_COMM_WORLD,1,"Norms are different");

  ierr = VecMax(xc,NULL,&vmax);CHKERRQ(ierr);
  ierr = VecMin(xc,NULL,&vmin);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"xc has max value %g min value %g\n",(double)vmax,(double)vmin);CHKERRQ(ierr);

  ierr = VecMaxPointwiseDivide(wc,xc,&vmax);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"wc/xc has max value %g\n",(double)vmax);CHKERRQ(ierr);

  ierr = VecDot(x,y,&dot[0]);CHKERRQ(ierr);
  ierr = VecDot(xc,yc,&dotc[0]);CHKERRQ(ierr);
  if (PetscAbsScalar(dot[0]-dotc[0])>10*PETSC_MACHINE_EPSILON) SETERRQ(PETSC_COMM_WORLD,1,"Dots are different");
  ierr = VecTDot(x,y,&dot[0]);CHKERRQ(ierr);
  ierr = VecTDot(xc,yc,&dotc[0]);CHKERRQ(ierr);
  if (PetscAbsScalar(dot[0]-dotc[0])>10*PETSC_MACHINE_EPSILON) SETERRQ(PETSC_COMM_WORLD,1,"Dots are different");

  vecs[0] = w; vecs[1] = y;
  ierr = VecMDot(x,2,vecs,dot);CHKERRQ(ierr);
  vecs[0] = wc; vecs[1] = yc;
  ierr = VecMDot(xc,2,vecs,dotc);CHKERRQ(ierr);
  if (PetscAbsScalar(dot[0]-dotc[0])>10*PETSC_MACHINE_EPSILON || PetscAbsScalar(dot[1]-dotc[1])>10*PETSC_MACHINE_EPSILON) SETERRQ(PETSC_COMM_WORLD,1,"Dots are different");
  vecs[0] = w; vecs[1] = y;
  ierr = VecMTDot(x,2,vecs,dot);CHKERRQ(ierr);
  vecs[0] = wc; vecs[1] = yc;
  ierr = VecMTDot(xc,2,vecs,dotc);CHKERRQ(ierr);
  if (PetscAbsScalar(dot[0]-dotc[0])>10*PETSC_MACHINE_EPSILON || PetscAbsScalar(dot[1]-dotc[1])>10*PETSC_MACHINE_EPSILON) SETERRQ(PETSC_COMM_WORLD,1,"Dots are different");

  ierr = VecDotNorm2(x,y,&dot[0],&norm);CHKERRQ(ierr);
  ierr = VecDotNorm2(xc,yc,&dotc[0],&normc);CHKERRQ(ierr);
  if (PetscAbsScalar(dot[0]-dotc[0])>10*PETSC_MACHINE_EPSILON) SETERRQ(PETSC_COMM_WORLD,1,"Dots are different");
  if (PetscAbsReal(norm-normc)>100*PETSC_MACHINE_EPSILON) SETERRQ(PETSC_COMM_WORLD,1,"Norms are different");

  ierr = VecAbs(w);CHKERRQ(ierr);
  ierr = VecAbs(wc);CHKERRQ(ierr);
  ierr = VecConjugate(x);CHKERRQ(ierr);
  ierr = VecConjugate(xc);CHKERRQ(ierr);
  ierr = VecShift(y,0.5);CHKERRQ(ierr);
  ierr = VecShift(yc,0.5);CHKERRQ(ierr);
  ierr = VecReciprocal(y);CHKERRQ(ierr);
  ierr = VecReciprocal(yc);CHKERRQ(ierr);
  ierr = VecExp(y);CHKERRQ(ierr);
  ierr = VecExp(yc);CHKERRQ(ierr);
  ierr = VecLog(y);CHKERRQ(ierr);
  ierr = VecLog(yc);CHKERRQ(ierr);
  ierr = VecNorm(y,NORM_1,&norm);CHKERRQ(ierr);
  ierr = VecNorm(yc,NORM_1,&normc);CHKERRQ(ierr);
  if (PetscAbsReal(norm-normc)>10*PETSC_MACHINE_EPSILON) SETERRQ(PETSC_COMM_WORLD,1,"Norms are different");

  ierr = VecPointwiseMult(w,x,y);CHKERRQ(ierr);
  ierr = VecPointwiseMult(wc,xc,yc);CHKERRQ(ierr);
  ierr = VecNorm(w,NORM_INFINITY,&norm);CHKERRQ(ierr);
  ierr = VecNorm(wc,NORM_INFINITY,&normc);CHKERRQ(ierr);
  if (PetscAbsReal(norm-normc)>10*PETSC_MACHINE_EPSILON) SETERRQ(PETSC_COMM_WORLD,1,"Norms are different");

  ierr = VecPointwiseMax(w,x,y);CHKERRQ(ierr);
  ierr = VecPointwiseMax(wc,xc,yc);CHKERRQ(ierr);
  ierr = VecNorm(w,NORM_INFINITY,&norm);CHKERRQ(ierr);
  ierr = VecNorm(wc,NORM_INFINITY,&normc);CHKERRQ(ierr);
  if (PetscAbsReal(norm-normc)>10*PETSC_MACHINE_EPSILON) SETERRQ(PETSC_COMM_WORLD,1,"Norms are different");

  ierr = VecSwap(x,y);CHKERRQ(ierr);
  ierr = VecSwap(xc,yc);CHKERRQ(ierr);
  ierr = VecPointwiseDivide(w,x,y);CHKERRQ(ierr);
  ierr = VecPointwiseDivide(wc,xc,yc);CHKERRQ(ierr);
  ierr = VecScale(w,0.3);CHKERRQ(ierr);
  ierr = VecScale(wc,0.3);CHKERRQ(ierr);
  ierr = VecSqrtAbs(w);CHKERRQ(ierr);
  ierr = VecSqrtAbs(wc);CHKERRQ(ierr);
  ierr = VecNorm(w,NORM_1_AND_2,norm12);CHKERRQ(ierr);
  ierr = VecNorm(wc,NORM_1_AND_2,norm12c);CHKERRQ(ierr);
  if (PetscAbsReal(norm12[0]-norm12c[0])>10*PETSC_MACHINE_EPSILON || PetscAbsReal(norm12[1]-norm12c[1])>10*PETSC_MACHINE_EPSILON) SETERRQ(PETSC_COMM_WORLD,1,"Norms are different");

  ierr = VecPointwiseMin(w,x,y);CHKERRQ(ierr);
  ierr = VecPointwiseMin(wc,xc,yc);CHKERRQ(ierr);
  ierr = VecPointwiseMaxAbs(x,y,w);CHKERRQ(ierr);
  ierr = VecPointwiseMaxAbs(xc,yc,wc);CHKERRQ(ierr);
  ierr = VecNorm(x,NORM_INFINITY,&norm);CHKERRQ(ierr);
  ierr = VecNorm(xc,NORM_INFINITY,&normc);CHKERRQ(ierr);
  if (PetscAbsReal(norm-normc)>10*PETSC_MACHINE_EPSILON) SETERRQ(PETSC_COMM_WORLD,1,"Norms are different");

  ierr = VecSetRandom(wc,NULL);CHKERRQ(ierr);

  ierr = VecDestroy(&v);CHKERRQ(ierr);
  ierr = VecDestroy(&w);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&y);CHKERRQ(ierr);
  ierr = VecDestroy(&vparent);CHKERRQ(ierr);
  ierr = VecDestroy(&vc);CHKERRQ(ierr);
  ierr = VecDestroy(&wc);CHKERRQ(ierr);
  ierr = VecDestroy(&xc);CHKERRQ(ierr);
  ierr = VecDestroy(&yc);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1

   test:
      suffix: 2
      nsize: 2

TEST*/
