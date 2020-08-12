/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   HZ iteration for generalized symmetric-indefinite eigenproblem.
   Based on Matlab code from David Watkins.

   References:

       [1] D.S. Watkins, The Matrix Eigenvalue Problem: GR and Krylov Subspace
           Methods, SIAM, 2007.

       [2] M.A. Brebner, J. Grad, "Eigenvalues of Ax = lambda Bx for real
           symmetric matrices A and B computed by reduction to pseudosymmetric
           form and the HR process", Linear Alg. Appl. 43:99-118, 1982.
*/

#include <slepc/private/dsimpl.h>
#include <slepcblaslapack.h>

/*
   Sets up a 2-by-2 matrix to eliminate y in the vector [x y]'.
   Transformation is rotator if sygn = 1 and hyperbolic if sygn = -1.
*/
static PetscErrorCode UnifiedRotation(PetscReal x,PetscReal y,PetscReal sygn,PetscReal *rot,PetscReal *rcond,PetscBool *swap)
{
  PetscReal nrm,c,s;

  PetscFunctionBegin;
  *swap = PETSC_FALSE;
  if (y == 0) {
    rot[0] = 1.0; rot[1] = 0.0; rot[2] = 0.0; rot[3] = 1.0;
    *rcond = 1.0;
  } else {
    nrm = PetscMax(PetscAbs(x),PetscAbs(y));
    c = x/nrm; s = y/nrm;
    if (sygn == 1.0) {  /* set up a rotator */
      nrm = PetscSqrtReal(c*c+s*s);
      c = c/nrm; s = s/nrm;
      /* rot = [c s; -s c]; */
      rot[0] = c; rot[1] = -s; rot[2] = s; rot[3] = c;
      *rcond = 1.0;
    } else if (sygn == -1) {  /* set up a hyperbolic transformation */
      nrm = c*c-s*s;
      if (nrm > 0) nrm = PetscSqrtReal(nrm);
      else if (nrm < 0) {
        nrm = PetscSqrtReal(-nrm);
        *swap = PETSC_TRUE;
      } else SETERRQ(PETSC_COMM_SELF,1,"Breakdown in construction of hyperbolic transformation");
      c = c/nrm; s = s/nrm;
      /* rot = [c -s; -s c]; */
      rot[0] = c; rot[1] = -s; rot[2] = -s; rot[3] = c;
      *rcond = PetscAbs(PetscAbs(s)-PetscAbs(c))/(PetscAbs(s)+PetscAbs(c));
    } else SETERRQ(PETSC_COMM_SELF,1,"Value of sygn sent to transetup must be 1 or -1");
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode HZStep(PetscBLASInt ntop,PetscBLASInt nn,PetscReal tr,PetscReal dt,PetscReal *aa,PetscReal *bb,PetscReal *dd,PetscScalar *uu,PetscInt n,PetscInt ld,PetscBool *flag)
{
  PetscErrorCode ierr;
  PetscBLASInt   one=1;
  PetscInt       k,jj,ii;
  PetscBLASInt   n_;
  PetscReal      bulge10,bulge20,bulge30,bulge31,bulge41,bulge42;
  PetscReal      sygn,rcond=1.0,worstcond,rot[4],buf[2],t;
  PetscScalar    rtmp;
  PetscBool      swap;

  PetscFunctionBegin;
  *flag = PETSC_FALSE;
  worstcond = 1.0;
  ierr = PetscBLASIntCast(n,&n_);CHKERRQ(ierr);

  /* Build initial bulge that sets step in motion */
  bulge10 = dd[ntop+1]*(aa[ntop]*(aa[ntop] - dd[ntop]*tr) + dt*dd[ntop]*dd[ntop]) + dd[ntop]*bb[ntop]*bb[ntop];
  bulge20 = bb[ntop]*(dd[ntop+1]*aa[ntop] + dd[ntop]*aa[ntop+1] - dd[ntop]*dd[ntop+1]*tr);
  bulge30 = bb[ntop]*bb[ntop+1]*dd[ntop];
  bulge31 = 0.0;
  bulge41 = 0.0;
  bulge42 = 0.0;

  /* Chase the bulge */
  for (jj=ntop;jj<nn-1;jj++) {

    /* Check for trivial bulge */
    if (jj>ntop && PetscMax(PetscMax(PetscAbs(bulge10),PetscAbs(bulge20)),PetscAbs(bulge30))<PETSC_MACHINE_EPSILON*(PetscAbs(aa[jj]) + PetscAbs(aa[jj+1]))) {
      bb[jj-1] = 0.0;  /* deflate and move on */

    } else { /* carry out the step */

      /* Annihilate tip entry bulge30 */
      if (bulge30 != 0.0) {

        /* Make an interchange if necessary to ensure that the
           first transformation is othogonal, not hyperbolic.  */
        if (dd[jj+1] != dd[jj+2]) { /* make an interchange */
          if (dd[jj] != dd[jj+1]) {  /* interchange 1st and 2nd */
            buf[0] = bulge20; bulge20 = bulge10; bulge10 = buf[0];
            buf[0] = aa[jj]; aa[jj] = aa[jj+1]; aa[jj+1] = buf[0];
            buf[0] = bb[jj+1]; bb[jj+1] = bulge31; bulge31 = buf[0];
            buf[0] = dd[jj]; dd[jj] = dd[jj+1]; dd[jj+1] = buf[0];
            for (k=0;k<n;k++) {
              rtmp = uu[k+jj*ld]; uu[k+jj*ld] = uu[k+(jj+1)*ld]; uu[k+(jj+1)*ld] = rtmp;
            }
          } else {  /* interchange 1st and 3rd */
            buf[0] = bulge30; bulge30 = bulge10; bulge10 = buf[0];
            buf[0] = aa[jj]; aa[jj] = aa[jj+2]; aa[jj+2] = buf[0];
            buf[0] = bb[jj]; bb[jj] = bb[jj+1]; bb[jj+1] = buf[0];
            buf[0] = dd[jj]; dd[jj] = dd[jj+2]; dd[jj+2] = buf[0];
            if (jj + 2 < nn-1) {
              bulge41 = bb[jj+2];
              bb[jj+2] = 0;
            }
            for (k=0;k<n;k++) {
              rtmp = uu[k+jj*ld]; uu[k+jj*ld] = uu[k+(jj+2)*ld]; uu[k+(jj+2)*ld] = rtmp;
            }
          }
        }

        /* Set up transforming matrix rot. */
        ierr = UnifiedRotation(bulge20,bulge30,1,rot,&rcond,&swap);CHKERRQ(ierr);

        /* Apply transforming matrix rot to T. */
        bulge20 = rot[0]*bulge20 + rot[2]*bulge30;
        buf[0] = rot[0]*bb[jj] + rot[2]*bulge31;
        buf[1] = rot[1]*bb[jj] + rot[3]*bulge31;
        bb[jj] = buf[0];
        bulge31 = buf[1];
        buf[0] = rot[0]*rot[0]*aa[jj+1] + 2.0*rot[0]*rot[2]*bb[jj+1] + rot[2]*rot[2]*aa[jj+2];
        buf[1] = rot[1]*rot[1]*aa[jj+1] + 2.0*rot[3]*rot[1]*bb[jj+1] + rot[3]*rot[3]*aa[jj+2];
        bb[jj+1] = rot[1]*rot[0]*aa[jj+1] + rot[3]*rot[2]*aa[jj+2] + (rot[3]*rot[0] + rot[1]*rot[2])*bb[jj+1];
        aa[jj+1] = buf[0];
        aa[jj+2] = buf[1];
        if (jj + 2 < nn-1) {
          bulge42 = bb[jj+2]*rot[2];
          bb[jj+2] = bb[jj+2]*rot[3];
        }

        /* Accumulate transforming matrix */
        PetscStackCallBLAS("BLASrot",BLASMIXEDrot_(&n_,uu+(jj+1)*ld,&one,uu+(jj+2)*ld,&one,&rot[0],&rot[2]));
      }

      /* Annihilate inner entry bulge20 */
      if (bulge20 != 0.0) {

        /* Begin by setting up transforming matrix rot */
        sygn = dd[jj]*dd[jj+1];
        ierr = UnifiedRotation(bulge10,bulge20,sygn,rot,&rcond,&swap);CHKERRQ(ierr);
        if (rcond<PETSC_MACHINE_EPSILON) {
          *flag = PETSC_TRUE;
          PetscFunctionReturn(0);
        }
        if (rcond < worstcond) worstcond = rcond;

        /* Apply transforming matrix rot to T */
        if (jj > ntop) bb[jj-1] = rot[0]*bulge10 + rot[2]*bulge20;
        buf[0] = rot[0]*rot[0]*aa[jj] + 2*rot[0]*rot[2]*bb[jj] + rot[2]*rot[2]*aa[jj+1];
        buf[1] = rot[1]*rot[1]*aa[jj] + 2*rot[3]*rot[1]*bb[jj] + rot[3]*rot[3]*aa[jj+1];
        bb[jj] = rot[1]*rot[0]*aa[jj] + rot[3]*rot[2]*aa[jj+1] + (rot[3]*rot[0] + rot[1]*rot[2])*bb[jj];
        aa[jj] = buf[0];
        aa[jj+1] = buf[1];
        if (jj + 1 < nn-1) {
          /* buf = [ bulge31 bb(jj+1) ] * rot' */
          buf[0] = rot[0]*bulge31 + rot[2]*bb[jj+1];
          buf[1] = rot[1]*bulge31 + rot[3]*bb[jj+1];
          bulge31 = buf[0];
          bb[jj+1] = buf[1];
        }
        if (jj + 2 < nn-1) {
          /* buf = [bulge41 bulge42] * rot' */
          buf[0] = rot[0]*bulge41 + rot[2]*bulge42;
          buf[1] = rot[1]*bulge41 + rot[3]*bulge42;
          bulge41 = buf[0];
          bulge42 = buf[1];
        }

        /* Apply transforming matrix rot to D */
        if (swap == 1) {
          buf[0] = dd[jj]; dd[jj] = dd[jj+1]; dd[jj+1] = buf[0];
        }

        /* Accumulate transforming matrix, uu(jj:jj+1,:) = rot*uu(jj:jj+1,:) */
        if (sygn==1) {
          PetscStackCallBLAS("BLASrot",BLASMIXEDrot_(&n_,uu+jj*ld,&one,uu+(jj+1)*ld,&one,&rot[0],&rot[2]));
        } else {
          if (PetscAbsReal(rot[0])>PetscAbsReal(rot[1])) { /* Type I */
            t = rot[1]/rot[0];
            for (ii=0;ii<n;ii++) {
              uu[jj*ld+ii] = rot[0]*uu[jj*ld+ii] + rot[1]*uu[(jj+1)*ld+ii];
              uu[(jj+1)*ld+ii] = t*uu[jj*ld+ii] + uu[(jj+1)*ld+ii]/rot[0];
            }
          } else { /* Type II */
            t = rot[0]/rot[1];
            for (ii=0;ii<n;ii++) {
              rtmp = uu[jj*ld+ii];
              uu[jj*ld+ii] = rot[0]*uu[jj*ld+ii] + rot[1]*uu[(jj+1)*ld+ii];
              uu[(jj+1)*ld+ii] = t*uu[jj*ld+ii] + rtmp/rot[1];
            }
          }
        }
      }
    }

    /* Adjust bulge for next step */
    bulge10 = bb[jj];
    bulge20 = bulge31;
    bulge30 = bulge41;
    bulge31 = bulge42;
    bulge41 = 0.0;
    bulge42 = 0.0;
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode HZIteration(PetscBLASInt nn,PetscBLASInt cgd,PetscReal *aa,PetscReal *bb,PetscReal *dd,PetscScalar *uu,PetscBLASInt ld)
{
  PetscErrorCode ierr;
  PetscBLASInt   j2,one=1;
  PetscInt       its,nits,nstop,jj,ntop,nbot,ntry;
  PetscReal      htr,det,dis,dif,tn,kt,c,s,tr,dt;
  PetscBool      flag=PETSC_FALSE;

  PetscFunctionBegin;
  its = 0;
  nbot = nn-1;
  nits = 0;
  nstop = 40*(nn - cgd);

  while (nbot >= cgd && nits < nstop) {

    /* Check for zeros on the subdiagonal */
    jj = nbot - 1;
    while (jj>=cgd && PetscAbs(bb[jj])>PETSC_MACHINE_EPSILON*(PetscAbs(aa[jj])+PetscAbs(aa[jj+1]))) jj = jj-1;
    if (jj>=cgd) bb[jj]=0;
    ntop = jj + 1;  /* starting point for step */
    if (ntop == nbot) {  /* isolate single eigenvalue */
      nbot = ntop - 1;
      its = 0;
    } else if (ntop+1 == nbot) {  /* isolate pair of eigenvalues */
      htr = 0.5*(aa[ntop]*dd[ntop] + aa[nbot]*dd[nbot]);
      det = dd[ntop]*dd[nbot]*(aa[ntop]*aa[nbot]-bb[ntop]*bb[ntop]);
      dis = htr*htr - det;
      if (dis > 0) {  /* distinct real eigenvalues */
        if (dd[ntop] == dd[nbot]) {  /* separate the eigenvalues by a Jacobi rotator */
          dif = aa[ntop]-aa[nbot];
          if (2.0*PetscAbs(bb[ntop])<=dif) {
            tn = 2*bb[ntop]/dif;
            tn = tn/(1.0 + PetscSqrtScalar(1.0+tn*tn));
          } else {
            kt = dif/(2.0*bb[ntop]);
            tn = PetscSign(kt)/(PetscAbs(kt)+PetscSqrtScalar(1.0+kt*kt));
          }
          c = 1.0/PetscSqrtScalar(1.0 + tn*tn);
          s = c*tn;
          aa[ntop] = aa[ntop] + tn*bb[ntop];
          aa[nbot] = aa[nbot] - tn*bb[ntop];
          bb[ntop] = 0;
          j2 = nn-cgd;
          PetscStackCallBLAS("BLASrot",BLASMIXEDrot_(&j2,uu+ntop*ld+cgd,&one,uu+nbot*ld+cgd,&one,&c,&s));
        }
      }
      nbot = ntop - 1;
    } else {  /* Do an HZ iteration */
      its = its + 1;
      nits = nits + 1;
      tr = aa[nbot-1]*dd[nbot-1] + aa[nbot]*dd[nbot];
      dt = dd[nbot-1]*dd[nbot]*(aa[nbot-1]*aa[nbot]-bb[nbot-1]*bb[nbot-1]);
      for (ntry=1;ntry<=6;ntry++) {
        ierr = HZStep(ntop,nbot+1,tr,dt,aa,bb,dd,uu,nn,ld,&flag);CHKERRQ(ierr);
        if (!flag) break;
        else if (ntry == 6) SETERRQ(PETSC_COMM_SELF,1,"Unable to complete hz step on six tries");
        else {
          tr = 0.9*tr; dt = 0.81*dt;
        }
      }
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DSSolve_GHIEP_HZ(DS ds,PetscScalar *wr,PetscScalar *wi)
{
  PetscErrorCode ierr;
  PetscInt       i,off;
  PetscBLASInt   n1,ld;
  PetscScalar    *A,*B,*Q;
  PetscReal      *d,*e,*s;

  PetscFunctionBegin;
#if !defined(PETSC_USE_COMPLEX)
  PetscValidScalarPointer(wi,3);
#endif
  ierr = PetscBLASIntCast(ds->ld,&ld);CHKERRQ(ierr);
  n1  = ds->n - ds->l;
  off = ds->l + ds->l*ld;
  A   = ds->mat[DS_MAT_A];
  B   = ds->mat[DS_MAT_B];
  Q   = ds->mat[DS_MAT_Q];
  d   = ds->rmat[DS_MAT_T];
  e   = ds->rmat[DS_MAT_T] + ld;
  s   = ds->rmat[DS_MAT_D];
  /* Quick return */
  if (n1 == 1) {
    for (i=0;i<=ds->l;i++) Q[i+i*ld] = 1.0;
    ierr = DSGHIEPComplexEigs(ds,0,ds->l,wr,wi);CHKERRQ(ierr);
    if (ds->compact) {
      wr[ds->l] = d[ds->l]/s[ds->l]; wi[ds->l] = 0.0;
    } else {
      d[ds->l] = PetscRealPart(A[off]); s[ds->l] = PetscRealPart(B[off]);
      wr[ds->l] = d[ds->l]/s[ds->l]; wi[ds->l] = 0.0;
    }
    PetscFunctionReturn(0);
  }
  /* Reduce to pseudotriadiagonal form */
  ierr = DSIntermediate_GHIEP(ds);CHKERRQ(ierr);
  ierr = HZIteration(ds->n,ds->l,d,e,s,Q,ld);CHKERRQ(ierr);
  if (!ds->compact) {
    ierr = DSSwitchFormat_GHIEP(ds,PETSC_FALSE);CHKERRQ(ierr);
  }
  /* Undo from diagonal the blocks whith real eigenvalues*/
  ierr = DSGHIEPRealBlocks(ds);CHKERRQ(ierr);

  /* Recover eigenvalues from diagonal */
  ierr = DSGHIEPComplexEigs(ds,0,ds->n,wr,wi);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
  if (wi) {
    for (i=ds->l;i<ds->n;i++) wi[i] = 0.0;
  }
#endif
  ds->t = ds->n;
  PetscFunctionReturn(0);
}

