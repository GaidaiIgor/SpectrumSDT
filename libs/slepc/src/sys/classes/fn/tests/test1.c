/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test rational function.\n\n";

#include <slepcfn.h>

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  FN             fn;
  PetscInt       i,na,nb;
  PetscScalar    x,y,yp,p[10],q[10],five=5.0,*pp,*qq;
  char           strx[50],str[50];

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = FNCreate(PETSC_COMM_WORLD,&fn);CHKERRQ(ierr);

  /* polynomial p(x) */
  na = 5;
  p[0] = -3.1; p[1] = 1.1; p[2] = 1.0; p[3] = -2.0; p[4] = 3.5;
  ierr = FNSetType(fn,FNRATIONAL);CHKERRQ(ierr);
  ierr = FNRationalSetNumerator(fn,na,p);CHKERRQ(ierr);
  ierr = FNView(fn,NULL);CHKERRQ(ierr);
  x = 2.2;
  ierr = SlepcSNPrintfScalar(strx,50,x,PETSC_FALSE);CHKERRQ(ierr);
  ierr = FNEvaluateFunction(fn,x,&y);CHKERRQ(ierr);
  ierr = FNEvaluateDerivative(fn,x,&yp);CHKERRQ(ierr);
  ierr = SlepcSNPrintfScalar(str,50,y,PETSC_FALSE);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"  f(%s)=%s\n",strx,str);CHKERRQ(ierr);
  ierr = SlepcSNPrintfScalar(str,50,yp,PETSC_FALSE);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"  f'(%s)=%s\n",strx,str);CHKERRQ(ierr);

  /* inverse of polynomial 1/q(x) */
  nb = 3;
  q[0] = -3.1; q[1] = 1.1; q[2] = 1.0;
  ierr = FNSetType(fn,FNRATIONAL);CHKERRQ(ierr);
  ierr = FNRationalSetNumerator(fn,0,NULL);CHKERRQ(ierr);  /* reset previous values */
  ierr = FNRationalSetDenominator(fn,nb,q);CHKERRQ(ierr);
  ierr = FNView(fn,NULL);CHKERRQ(ierr);
  x = 2.2;
  ierr = SlepcSNPrintfScalar(strx,50,x,PETSC_FALSE);CHKERRQ(ierr);
  ierr = FNEvaluateFunction(fn,x,&y);CHKERRQ(ierr);
  ierr = FNEvaluateDerivative(fn,x,&yp);CHKERRQ(ierr);
  ierr = SlepcSNPrintfScalar(str,50,y,PETSC_FALSE);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"  f(%s)=%s\n",strx,str);CHKERRQ(ierr);
  ierr = SlepcSNPrintfScalar(str,50,yp,PETSC_FALSE);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"  f'(%s)=%s\n",strx,str);CHKERRQ(ierr);

  /* rational p(x)/q(x) */
  na = 2; nb = 3;
  p[0] = 1.1; p[1] = 1.1;
  q[0] = 1.0; q[1] = -2.0; q[2] = 3.5;
  ierr = FNSetType(fn,FNRATIONAL);CHKERRQ(ierr);
  ierr = FNRationalSetNumerator(fn,na,p);CHKERRQ(ierr);
  ierr = FNRationalSetDenominator(fn,nb,q);CHKERRQ(ierr);
  ierr = FNSetScale(fn,1.2,0.5);CHKERRQ(ierr);
  ierr = FNView(fn,NULL);CHKERRQ(ierr);
  x = 2.2;
  ierr = SlepcSNPrintfScalar(strx,50,x,PETSC_FALSE);CHKERRQ(ierr);
  ierr = FNEvaluateFunction(fn,x,&y);CHKERRQ(ierr);
  ierr = FNEvaluateDerivative(fn,x,&yp);CHKERRQ(ierr);
  ierr = SlepcSNPrintfScalar(str,50,y,PETSC_FALSE);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"  f(%s)=%s\n",strx,str);CHKERRQ(ierr);
  ierr = SlepcSNPrintfScalar(str,50,yp,PETSC_FALSE);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"  f'(%s)=%s\n",strx,str);CHKERRQ(ierr);

  ierr = FNRationalGetNumerator(fn,&na,&pp);CHKERRQ(ierr);
  ierr = FNRationalGetDenominator(fn,&nb,&qq);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Coefficients:\n  Numerator: ");CHKERRQ(ierr);
  for (i=0;i<na;i++) {
    ierr = SlepcSNPrintfScalar(str,50,pp[i],PETSC_FALSE);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"%s ",str);CHKERRQ(ierr);
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n  Denominator: ");CHKERRQ(ierr);
  for (i=0;i<nb;i++) {
    ierr = SlepcSNPrintfScalar(str,50,qq[i],PETSC_FALSE);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"%s ",str);CHKERRQ(ierr);
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
  ierr = PetscFree(pp);CHKERRQ(ierr);
  ierr = PetscFree(qq);CHKERRQ(ierr);

  /* constant */
  ierr = FNSetType(fn,FNRATIONAL);CHKERRQ(ierr);
  ierr = FNRationalSetNumerator(fn,1,&five);CHKERRQ(ierr);
  ierr = FNRationalSetDenominator(fn,0,NULL);CHKERRQ(ierr);  /* reset previous values */
  ierr = FNView(fn,NULL);CHKERRQ(ierr);
  x = 2.2;
  ierr = SlepcSNPrintfScalar(strx,50,x,PETSC_FALSE);CHKERRQ(ierr);
  ierr = FNEvaluateFunction(fn,x,&y);CHKERRQ(ierr);
  ierr = FNEvaluateDerivative(fn,x,&yp);CHKERRQ(ierr);
  ierr = SlepcSNPrintfScalar(str,50,y,PETSC_FALSE);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"  f(%s)=%s\n",strx,str);CHKERRQ(ierr);
  ierr = SlepcSNPrintfScalar(str,50,yp,PETSC_FALSE);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"  f'(%s)=%s\n",strx,str);CHKERRQ(ierr);

  ierr = FNDestroy(&fn);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      nsize: 1

TEST*/
