/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test exponential function.\n\n";

#include <slepcfn.h>

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  FN             fn,fncopy;
  PetscScalar    x,y,yp,tau,eta,alpha,beta;
  char           strx[50],str[50];

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = FNCreate(PETSC_COMM_WORLD,&fn);CHKERRQ(ierr);
  ierr = FNSetFromOptions(fn);CHKERRQ(ierr);

  /* plain exponential exp(x) */
  ierr = FNSetType(fn,FNEXP);CHKERRQ(ierr);
  ierr = FNView(fn,NULL);CHKERRQ(ierr);
  x = 2.2;
  ierr = SlepcSNPrintfScalar(strx,50,x,PETSC_FALSE);CHKERRQ(ierr);
  ierr = FNEvaluateFunction(fn,x,&y);CHKERRQ(ierr);
  ierr = FNEvaluateDerivative(fn,x,&yp);CHKERRQ(ierr);
  ierr = SlepcSNPrintfScalar(str,50,y,PETSC_FALSE);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"  f(%s)=%s\n",strx,str);CHKERRQ(ierr);
  ierr = SlepcSNPrintfScalar(str,50,yp,PETSC_FALSE);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"  f'(%s)=%s\n",strx,str);CHKERRQ(ierr);

  /* exponential with scaling factors eta*exp(tau*x) */
  ierr = FNSetType(fn,FNEXP);CHKERRQ(ierr);
  tau = -0.2;
  eta = 1.3;
  ierr = FNSetScale(fn,tau,eta);CHKERRQ(ierr);
  ierr = FNView(fn,NULL);CHKERRQ(ierr);
  x = 2.2;
  ierr = SlepcSNPrintfScalar(strx,50,x,PETSC_FALSE);CHKERRQ(ierr);
  ierr = FNEvaluateFunction(fn,x,&y);CHKERRQ(ierr);
  ierr = FNEvaluateDerivative(fn,x,&yp);CHKERRQ(ierr);
  ierr = SlepcSNPrintfScalar(str,50,y,PETSC_FALSE);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"  f(%s)=%s\n",strx,str);CHKERRQ(ierr);
  ierr = SlepcSNPrintfScalar(str,50,yp,PETSC_FALSE);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"  f'(%s)=%s\n",strx,str);CHKERRQ(ierr);

  /* test FNDuplicate */
  ierr = FNDuplicate(fn,PetscObjectComm((PetscObject)fn),&fncopy);CHKERRQ(ierr);

  /* test FNGetScale */
  ierr = FNGetScale(fncopy,&alpha,&beta);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Parameters:\n - alpha: ");CHKERRQ(ierr);
  ierr = SlepcSNPrintfScalar(str,50,alpha,PETSC_FALSE);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%s ",str);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n - beta: ");CHKERRQ(ierr);
  ierr = SlepcSNPrintfScalar(str,50,beta,PETSC_FALSE);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%s ",str);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);

  ierr = FNDestroy(&fn);CHKERRQ(ierr);
  ierr = FNDestroy(&fncopy);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      nsize: 1
      filter: grep -v "computing matrix functions"

TEST*/
