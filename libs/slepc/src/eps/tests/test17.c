/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test interface functions of spectrum-slicing Krylov-Schur.\n\n"
  "This is based on ex12.c. The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions in x dimension.\n"
  "  -m <m>, where <m> = number of grid subdivisions in y dimension.\n\n";

#include <slepceps.h>

int main(int argc,char **argv)
{
  Mat            A,B;         /* matrices */
  Mat            As,Bs;       /* matrices distributed in subcommunicators */
  Mat            Au;          /* matrix used to modify A on subcommunicators */
  EPS            eps;         /* eigenproblem solver context */
  ST             st;          /* spectral transformation context */
  KSP            ksp;
  PC             pc;
  Vec            v;
  PetscMPIInt    size,rank;
  PetscInt       N,n=35,m,Istart,Iend,II,nev,ncv,mpd,i,j,k,*inertias,npart,nval,nloc,nlocs,mlocs;
  PetscBool      flag,showinertia=PETSC_TRUE,lock,detect;
  PetscReal      int0,int1,*shifts,keep,*subint,*evals;
  PetscScalar    lambda;
  char           vlist[4000];
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

  ierr = PetscOptionsGetBool(NULL,NULL,"-showinertia",&showinertia,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-m",&m,&flag);CHKERRQ(ierr);
  if (!flag) m=n;
  N = n*m;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nSpectrum-slicing test, N=%D (%Dx%D grid)\n\n",N,n,m);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the matrices that define the eigensystem, Ax=kBx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD,&B);CHKERRQ(ierr);
  ierr = MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(B);CHKERRQ(ierr);
  ierr = MatSetUp(B);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
  for (II=Istart;II<Iend;II++) {
    i = II/n; j = II-i*n;
    if (i>0) { ierr = MatSetValue(A,II,II-n,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (i<m-1) { ierr = MatSetValue(A,II,II+n,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (j>0) { ierr = MatSetValue(A,II,II-1,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (j<n-1) { ierr = MatSetValue(A,II,II+1,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    ierr = MatSetValue(A,II,II,4.0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatSetValue(B,II,II,2.0,INSERT_VALUES);CHKERRQ(ierr);
  }
  if (Istart==0) {
    ierr = MatSetValue(B,0,0,6.0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatSetValue(B,0,1,-1.0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatSetValue(B,1,0,-1.0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatSetValue(B,1,1,1.0,INSERT_VALUES);CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
  ierr = EPSSetOperators(eps,A,B);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_GHEP);CHKERRQ(ierr);
  ierr = EPSSetType(eps,EPSKRYLOVSCHUR);CHKERRQ(ierr);

  /*
     Set interval and other settings for spectrum slicing
  */
  ierr = EPSSetWhichEigenpairs(eps,EPS_ALL);CHKERRQ(ierr);
  int0 = 1.1; int1 = 1.3;
  ierr = EPSSetInterval(eps,int0,int1);CHKERRQ(ierr);
  ierr = EPSGetST(eps,&st);CHKERRQ(ierr);
  ierr = STSetType(st,STSINVERT);CHKERRQ(ierr);
  ierr = STGetKSP(st,&ksp);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPPREONLY);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCCHOLESKY);CHKERRQ(ierr);

  /*
     Test interface functions of Krylov-Schur solver
  */
  ierr = EPSKrylovSchurGetRestart(eps,&keep);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Restart parameter before changing = %g",(double)keep);CHKERRQ(ierr);
  ierr = EPSKrylovSchurSetRestart(eps,0.4);CHKERRQ(ierr);
  ierr = EPSKrylovSchurGetRestart(eps,&keep);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," ... changed to %g\n",(double)keep);CHKERRQ(ierr);

  ierr = EPSKrylovSchurGetDetectZeros(eps,&detect);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Detect zeros before changing = %d",(int)detect);CHKERRQ(ierr);
  ierr = EPSKrylovSchurSetDetectZeros(eps,PETSC_TRUE);CHKERRQ(ierr);
  ierr = EPSKrylovSchurGetDetectZeros(eps,&detect);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," ... changed to %d\n",(int)detect);CHKERRQ(ierr);

  ierr = EPSKrylovSchurGetLocking(eps,&lock);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Locking flag before changing = %d",(int)lock);CHKERRQ(ierr);
  ierr = EPSKrylovSchurSetLocking(eps,PETSC_FALSE);CHKERRQ(ierr);
  ierr = EPSKrylovSchurGetLocking(eps,&lock);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," ... changed to %d\n",(int)lock);CHKERRQ(ierr);

  ierr = EPSKrylovSchurGetDimensions(eps,&nev,&ncv,&mpd);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Sub-solve dimensions before changing = [%D,%D,%D]",nev,ncv,mpd);CHKERRQ(ierr);
  ierr = EPSKrylovSchurSetDimensions(eps,30,60,60);CHKERRQ(ierr);
  ierr = EPSKrylovSchurGetDimensions(eps,&nev,&ncv,&mpd);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," ... changed to [%D,%D,%D]\n",nev,ncv,mpd);CHKERRQ(ierr);

  if (size>1) {
    ierr = EPSKrylovSchurSetPartitions(eps,size);CHKERRQ(ierr);
    ierr = EPSKrylovSchurGetPartitions(eps,&npart);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Using %D partitions\n",npart);CHKERRQ(ierr);

    ierr = PetscMalloc1(npart+1,&subint);CHKERRQ(ierr);
    subint[0] = int0;
    subint[npart] = int1;
    for (i=1;i<npart;i++) subint[i] = int0+i*(int1-int0)/npart;
    ierr = EPSKrylovSchurSetSubintervals(eps,subint);CHKERRQ(ierr);
    ierr = PetscFree(subint);CHKERRQ(ierr);
    ierr = EPSKrylovSchurGetSubintervals(eps,&subint);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Using sub-interval separations = ");CHKERRQ(ierr);
    for (i=1;i<npart;i++) { ierr = PetscPrintf(PETSC_COMM_WORLD," %g",(double)subint[i]);CHKERRQ(ierr); }
    ierr = PetscFree(subint);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
  }

  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           Compute all eigenvalues in interval and display info
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSSetUp(eps);CHKERRQ(ierr);
  ierr = EPSKrylovSchurGetInertias(eps,&k,&shifts,&inertias);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Inertias after EPSSetUp:\n");CHKERRQ(ierr);
  for (i=0;i<k;i++) {
    ierr = PetscPrintf(PETSC_COMM_WORLD," .. %g (%D)\n",(double)shifts[i],inertias[i]);CHKERRQ(ierr);
  }
  ierr = PetscFree(shifts);CHKERRQ(ierr);
  ierr = PetscFree(inertias);CHKERRQ(ierr);

  ierr = EPSSolve(eps);CHKERRQ(ierr);
  ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
  ierr = EPSGetInterval(eps,&int0,&int1);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Found %D eigenvalues in interval [%g,%g]\n",nev,(double)int0,(double)int1);CHKERRQ(ierr);

  if (showinertia) {
    ierr = EPSKrylovSchurGetInertias(eps,&k,&shifts,&inertias);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Used %D shifts (inertia):\n",k);CHKERRQ(ierr);
    for (i=0;i<k;i++) {
      ierr = PetscPrintf(PETSC_COMM_WORLD," .. %g (%D)\n",(double)shifts[i],inertias[i]);CHKERRQ(ierr);
    }
    ierr = PetscFree(shifts);CHKERRQ(ierr);
    ierr = PetscFree(inertias);CHKERRQ(ierr);
  }

  ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL);CHKERRQ(ierr);

  if (size>1) {
    ierr = EPSKrylovSchurGetSubcommInfo(eps,&k,&nval,&v);CHKERRQ(ierr);
    ierr = PetscMalloc1(nval,&evals);CHKERRQ(ierr);
    for (i=0;i<nval;i++) {
      ierr = EPSKrylovSchurGetSubcommPairs(eps,i,&lambda,v);CHKERRQ(ierr);
      evals[i] = PetscRealPart(lambda);
    }
    ierr = PetscFormatRealArray(vlist,sizeof(vlist),"%f",nval,evals);CHKERRQ(ierr);
    ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD," Process %d has worked in sub-interval %D, containing %D eigenvalues: %s\n",(int)rank,k,nval,vlist);CHKERRQ(ierr);
    ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);CHKERRQ(ierr);
    ierr = VecDestroy(&v);CHKERRQ(ierr);
    ierr = PetscFree(evals);CHKERRQ(ierr);

    ierr = EPSKrylovSchurGetSubcommMats(eps,&As,&Bs);CHKERRQ(ierr);
    ierr = MatGetLocalSize(A,&nloc,NULL);CHKERRQ(ierr);
    ierr = MatGetLocalSize(As,&nlocs,&mlocs);CHKERRQ(ierr);
    ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD," Process %d owns %D rows of the global matrices, and %D rows in the subcommunicator\n",(int)rank,nloc,nlocs);CHKERRQ(ierr);
    ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);CHKERRQ(ierr);

    /* modify A on subcommunicators */
    ierr = MatCreate(PetscObjectComm((PetscObject)As),&Au);CHKERRQ(ierr);
    ierr = MatSetSizes(Au,nlocs,mlocs,N,N);CHKERRQ(ierr);
    ierr = MatSetFromOptions(Au);CHKERRQ(ierr);
    ierr = MatSetUp(Au);CHKERRQ(ierr);
    ierr = MatGetOwnershipRange(Au,&Istart,&Iend);CHKERRQ(ierr);
    for (II=Istart;II<Iend;II++) {
      ierr = MatSetValue(Au,II,II,0.5,INSERT_VALUES);CHKERRQ(ierr);
    }
    ierr = MatAssemblyBegin(Au,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Au,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = EPSKrylovSchurUpdateSubcommMats(eps,1.0,-1.0,Au,0.0,0.0,NULL,DIFFERENT_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
    ierr = MatDestroy(&Au);CHKERRQ(ierr);
  }

  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = MatDestroy(&B);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      nsize: 2
      args: -showinertia 0 -log_exclude eps,st,rg,bv,ds
      requires: !single

   test:
      suffix: 2
      nsize: 1
      args: -showinertia 0 -log_exclude eps,st,rg,bv,ds
      requires: !single

TEST*/
