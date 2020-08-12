/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Power grid small signal stability analysis of WECC 9 bus system.\n\
This example is based on the 9-bus (node) example given in the book Power\n\
Systems Dynamics and Stability (Chapter 8) by P. Sauer and M. A. Pai.\n\
The power grid in this example consists of 9 buses (nodes), 3 generators,\n\
3 loads, and 9 transmission lines. The network equations are written\n\
in current balance form using rectangular coordinates. It uses the SLEPc\n\
package to calculate the eigenvalues for small signal stability analysis\n\n";

/*
   This example is based on PETSc's ex9bus example (under TS).

   The equations for the stability analysis are described by the DAE

   \dot{x} = f(x,y,t)
     0     = g(x,y,t)

   where the generators are described by differential equations, while the algebraic
   constraints define the network equations.

   The generators are modeled with a 4th order differential equation describing the electrical
   and mechanical dynamics. Each generator also has an exciter system modeled by 3rd order
   diff. eqns. describing the exciter, voltage regulator, and the feedback stabilizer
   mechanism.

   The network equations are described by nodal current balance equations.
    I(x,y) - Y*V = 0

   where:
    I(x,y) is the current injected from generators and loads.
      Y    is the admittance matrix, and
      V    is the voltage vector

   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   The linearized equations for the eigenvalue analysis are

     \dot{\delta{x}} = f_x*\delta{x} + f_y*\delta{y}
             0       = g_x*\delta{x} + g_y*\delta{y}

   This gives the linearized sensitivity matrix
     A = | f_x  f_y |
         | g_x  g_y |

   We are interested in the eigenvalues of the Schur complement of A
     \hat{A} = f_x - g_x*inv(g_y)*f_y


   Example contributed by: Shrirang Abhyankar
*/

#include <petscdm.h>
#include <petscdmda.h>
#include <petscdmcomposite.h>
#include <slepceps.h>

#define freq 60
#define w_s (2*PETSC_PI*freq)

/* Sizes and indices */
const PetscInt nbus    = 9; /* Number of network buses */
const PetscInt ngen    = 3; /* Number of generators */
const PetscInt nload   = 3; /* Number of loads */
const PetscInt gbus[3] = {0,1,2}; /* Buses at which generators are incident */
const PetscInt lbus[3] = {4,5,7}; /* Buses at which loads are incident */

/* Generator real and reactive powers (found via loadflow) */
const PetscScalar PG[3] = {0.716786142395021,1.630000000000000,0.850000000000000};
const PetscScalar QG[3] = {0.270702180178785,0.066120127797275,-0.108402221791588};
/* Generator constants */
const PetscScalar H[3]    = {23.64,6.4,3.01};   /* Inertia constant */
const PetscScalar Rs[3]   = {0.0,0.0,0.0}; /* Stator Resistance */
const PetscScalar Xd[3]   = {0.146,0.8958,1.3125};  /* d-axis reactance */
const PetscScalar Xdp[3]  = {0.0608,0.1198,0.1813}; /* d-axis transient reactance */
const PetscScalar Xq[3]   = {0.0969,0.8645,1.2578}; /* q-axis reactance Xq(1) set to 0.4360, value given in text 0.0969 */
const PetscScalar Xqp[3]  = {0.0969,0.1969,0.25}; /* q-axis transient reactance */
const PetscScalar Td0p[3] = {8.96,6.0,5.89}; /* d-axis open circuit time constant */
const PetscScalar Tq0p[3] = {0.31,0.535,0.6}; /* q-axis open circuit time constant */
PetscScalar M[3]; /* M = 2*H/w_s */
PetscScalar D[3]; /* D = 0.1*M */

PetscScalar TM[3]; /* Mechanical Torque */
/* Exciter system constants */
const PetscScalar KA[3] = {20.0,20.0,20.0};  /* Voltage regulartor gain constant */
const PetscScalar TA[3] = {0.2,0.2,0.2};     /* Voltage regulator time constant */
const PetscScalar KE[3] = {1.0,1.0,1.0};     /* Exciter gain constant */
const PetscScalar TE[3] = {0.314,0.314,0.314}; /* Exciter time constant */
const PetscScalar KF[3] = {0.063,0.063,0.063};  /* Feedback stabilizer gain constant */
const PetscScalar TF[3] = {0.35,0.35,0.35};    /* Feedback stabilizer time constant */
const PetscScalar k1[3] = {0.0039,0.0039,0.0039};
const PetscScalar k2[3] = {1.555,1.555,1.555};  /* k1 and k2 for calculating the saturation function SE = k1*exp(k2*Efd) */

PetscScalar Vref[3];
/* Load constants
  We use a composite load model that describes the load and reactive powers at each time instant as follows
  P(t) = \sum\limits_{i=0}^ld_nsegsp \ld_alphap_i*P_D0(\frac{V_m(t)}{V_m0})^\ld_betap_i
  Q(t) = \sum\limits_{i=0}^ld_nsegsq \ld_alphaq_i*Q_D0(\frac{V_m(t)}{V_m0})^\ld_betaq_i
  where
    ld_nsegsp,ld_nsegsq - Number of individual load models for real and reactive power loads
    ld_alphap,ld_alphap - Percentage contribution (weights) or loads
    P_D0                - Real power load
    Q_D0                - Reactive power load
    V_m(t)              - Voltage magnitude at time t
    V_m0                - Voltage magnitude at t = 0
    ld_betap, ld_betaq  - exponents describing the load model for real and reactive part

    Note: All loads have the same characteristic currently.
*/
const PetscScalar PD0[3] = {1.25,0.9,1.0};
const PetscScalar QD0[3] = {0.5,0.3,0.35};
const PetscInt    ld_nsegsp[3] = {3,3,3};
const PetscScalar ld_alphap[3] = {0.0,0.0,1.0};
const PetscScalar ld_betap[3]  = {2.0,1.0,0.0};
const PetscInt    ld_nsegsq[3] = {3,3,3};
const PetscScalar ld_alphaq[3] = {0.0,0.0,1.0};
const PetscScalar ld_betaq[3]  = {2.0,1.0,0.0};

typedef struct {
  DM       dmgen, dmnet; /* DMs to manage generator and network subsystem */
  DM       dmpgrid;      /* Composite DM to manage the entire power grid */
  Mat      Ybus;         /* Network admittance matrix */
  Vec      V0;           /* Initial voltage vector (Power flow solution) */
  PetscInt neqs_gen,neqs_net,neqs_pgrid;
  IS       is_diff;      /* indices for differential equations */
  IS       is_alg;       /* indices for algebraic equations */
} Userctx;

/* Converts from machine frame (dq) to network (phase a real,imag) reference frame */
PetscErrorCode dq2ri(PetscScalar Fd,PetscScalar Fq,PetscScalar delta,PetscScalar *Fr,PetscScalar *Fi)
{
  PetscFunctionBegin;
  *Fr =  Fd*PetscSinScalar(delta) + Fq*PetscCosScalar(delta);
  *Fi = -Fd*PetscCosScalar(delta) + Fq*PetscSinScalar(delta);
  PetscFunctionReturn(0);
}

/* Converts from network frame ([phase a real,imag) to machine (dq) reference frame */
PetscErrorCode ri2dq(PetscScalar Fr,PetscScalar Fi,PetscScalar delta,PetscScalar *Fd,PetscScalar *Fq)
{
  PetscFunctionBegin;
  *Fd =  Fr*PetscSinScalar(delta) - Fi*PetscCosScalar(delta);
  *Fq =  Fr*PetscCosScalar(delta) + Fi*PetscSinScalar(delta);
  PetscFunctionReturn(0);
}

PetscErrorCode SetInitialGuess(Vec X,Userctx *user)
{
  PetscErrorCode ierr;
  Vec            Xgen,Xnet;
  PetscScalar    *xgen,*xnet;
  PetscInt       i,idx=0;
  PetscScalar    Vr,Vi,IGr,IGi,Vm,Vm2;
  PetscScalar    Eqp,Edp,delta;
  PetscScalar    Efd,RF,VR; /* Exciter variables */
  PetscScalar    Id,Iq;  /* Generator dq axis currents */
  PetscScalar    theta,Vd,Vq,SE;

  PetscFunctionBegin;
  M[0] = 2*H[0]/w_s; M[1] = 2*H[1]/w_s; M[2] = 2*H[2]/w_s;
      /*      D[0] = 0.1*M[0]; D[1] = 0.1*M[1]; D[2] = 0.1*M[2];
       */
  D[0] = D[1] = D[2] = 0.0;
  ierr = DMCompositeGetLocalVectors(user->dmpgrid,&Xgen,&Xnet);CHKERRQ(ierr);

  /* Network subsystem initialization */
  ierr = VecCopy(user->V0,Xnet);CHKERRQ(ierr);

  /* Generator subsystem initialization */
  ierr = VecGetArray(Xgen,&xgen);CHKERRQ(ierr);
  ierr = VecGetArray(Xnet,&xnet);CHKERRQ(ierr);

  for (i=0; i < ngen; i++) {
    Vr  = xnet[2*gbus[i]]; /* Real part of generator terminal voltage */
    Vi  = xnet[2*gbus[i]+1]; /* Imaginary part of the generator terminal voltage */
    Vm  = PetscSqrtScalar(Vr*Vr + Vi*Vi); Vm2 = Vm*Vm;
    IGr = (Vr*PG[i] + Vi*QG[i])/Vm2;
    IGi = (Vi*PG[i] - Vr*QG[i])/Vm2;

    delta = atan2(Vi+Xq[i]*IGr,Vr-Xq[i]*IGi); /* Machine angle */

    theta = PETSC_PI/2.0 - delta;

    Id = IGr*PetscCosScalar(theta) - IGi*PetscSinScalar(theta); /* d-axis stator current */
    Iq = IGr*PetscSinScalar(theta) + IGi*PetscCosScalar(theta); /* q-axis stator current */

    Vd = Vr*PetscCosScalar(theta) - Vi*PetscSinScalar(theta);
    Vq = Vr*PetscSinScalar(theta) + Vi*PetscCosScalar(theta);

    Edp = Vd + Rs[i]*Id - Xqp[i]*Iq; /* d-axis transient EMF */
    Eqp = Vq + Rs[i]*Iq + Xdp[i]*Id; /* q-axis transient EMF */

    TM[i] = PG[i];

    /* The generator variables are ordered as [Eqp,Edp,delta,w,Id,Iq] */
    xgen[idx]   = Eqp;
    xgen[idx+1] = Edp;
    xgen[idx+2] = delta;
    xgen[idx+3] = w_s;

    idx = idx + 4;

    xgen[idx]   = Id;
    xgen[idx+1] = Iq;

    idx = idx + 2;

    /* Exciter */
    Efd = Eqp + (Xd[i] - Xdp[i])*Id;
    SE  = k1[i]*PetscExpScalar(k2[i]*Efd);
    VR  =  KE[i]*Efd + SE;
    RF  =  KF[i]*Efd/TF[i];

    xgen[idx]   = Efd;
    xgen[idx+1] = RF;
    xgen[idx+2] = VR;

    Vref[i] = Vm + (VR/KA[i]);

    idx = idx + 3;
  }

  ierr = VecRestoreArray(Xgen,&xgen);CHKERRQ(ierr);
  ierr = VecRestoreArray(Xnet,&xnet);CHKERRQ(ierr);

  /* ierr = VecView(Xgen,0);CHKERRQ(ierr); */
  ierr = DMCompositeGather(user->dmpgrid,INSERT_VALUES,X,Xgen,Xnet);CHKERRQ(ierr);
  ierr = DMCompositeRestoreLocalVectors(user->dmpgrid,&Xgen,&Xnet);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode PreallocateJacobian(Mat J,Userctx *user)
{
  PetscErrorCode ierr;
  PetscInt       *d_nnz;
  PetscInt       i,idx=0,start=0;
  PetscInt       ncols;

  PetscFunctionBegin;
  ierr = PetscMalloc1(user->neqs_pgrid,&d_nnz);CHKERRQ(ierr);
  for (i=0; i<user->neqs_pgrid; i++) d_nnz[i] = 0;
  /* Generator subsystem */
  for (i=0; i < ngen; i++) {

    d_nnz[idx]   += 3;
    d_nnz[idx+1] += 2;
    d_nnz[idx+2] += 2;
    d_nnz[idx+3] += 5;
    d_nnz[idx+4] += 6;
    d_nnz[idx+5] += 6;

    d_nnz[user->neqs_gen+2*gbus[i]]   += 3;
    d_nnz[user->neqs_gen+2*gbus[i]+1] += 3;

    d_nnz[idx+6] += 2;
    d_nnz[idx+7] += 2;
    d_nnz[idx+8] += 5;

    idx = idx + 9;
  }

  start = user->neqs_gen;

  for (i=0; i < nbus; i++) {
    ierr = MatGetRow(user->Ybus,2*i,&ncols,NULL,NULL);CHKERRQ(ierr);
    d_nnz[start+2*i]   += ncols;
    d_nnz[start+2*i+1] += ncols;
    ierr = MatRestoreRow(user->Ybus,2*i,&ncols,NULL,NULL);CHKERRQ(ierr);
  }

  ierr = MatSeqAIJSetPreallocation(J,0,d_nnz);CHKERRQ(ierr);

  ierr = PetscFree(d_nnz);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   J = [-df_dx, -df_dy
        dg_dx, dg_dy]
*/
PetscErrorCode ResidualJacobian(Vec X,Mat J,void *ctx)
{
  PetscErrorCode ierr;
  Userctx        *user=(Userctx*)ctx;
  Vec            Xgen,Xnet;
  PetscScalar    *xgen,*xnet;
  PetscInt       i,idx=0;
  PetscScalar    Vr,Vi,Vm,Vm2;
  PetscScalar    Eqp,Edp,delta; /* Generator variables */
  PetscScalar    Efd;
  PetscScalar    Id,Iq;  /* Generator dq axis currents */
  PetscScalar    Vd,Vq;
  PetscScalar    val[10];
  PetscInt       row[2],col[10];
  PetscInt       net_start=user->neqs_gen;
  PetscScalar    Zdq_inv[4],det;
  PetscScalar    dVd_dVr,dVd_dVi,dVq_dVr,dVq_dVi,dVd_ddelta,dVq_ddelta;
  PetscScalar    dIGr_ddelta,dIGi_ddelta,dIGr_dId,dIGr_dIq,dIGi_dId,dIGi_dIq;
  PetscScalar    dSE_dEfd;
  PetscScalar    dVm_dVd,dVm_dVq,dVm_dVr,dVm_dVi;
  PetscInt          ncols;
  const PetscInt    *cols;
  const PetscScalar *yvals;
  PetscInt          k;
  PetscScalar PD,QD,Vm0,*v0,Vm4;
  PetscScalar dPD_dVr,dPD_dVi,dQD_dVr,dQD_dVi;
  PetscScalar dIDr_dVr,dIDr_dVi,dIDi_dVr,dIDi_dVi;


  PetscFunctionBegin;
  ierr = MatZeroEntries(J);CHKERRQ(ierr);
  ierr = DMCompositeGetLocalVectors(user->dmpgrid,&Xgen,&Xnet);CHKERRQ(ierr);
  ierr = DMCompositeScatter(user->dmpgrid,X,Xgen,Xnet);CHKERRQ(ierr);

  ierr = VecGetArray(Xgen,&xgen);CHKERRQ(ierr);
  ierr = VecGetArray(Xnet,&xnet);CHKERRQ(ierr);

  /* Generator subsystem */
  for (i=0; i < ngen; i++) {
    Eqp   = xgen[idx];
    Edp   = xgen[idx+1];
    delta = xgen[idx+2];
    Id    = xgen[idx+4];
    Iq    = xgen[idx+5];
    Efd   = xgen[idx+6];

    /*    fgen[idx]   = (Eqp + (Xd[i] - Xdp[i])*Id - Efd)/Td0p[i]; */
    row[0] = idx;
    col[0] = idx;           col[1] = idx+4;          col[2] = idx+6;
    val[0] = 1/ Td0p[i]; val[1] = (Xd[i] - Xdp[i])/ Td0p[i]; val[2] = -1/Td0p[i];

    ierr = MatSetValues(J,1,row,3,col,val,INSERT_VALUES);CHKERRQ(ierr);

    /*    fgen[idx+1] = (Edp - (Xq[i] - Xqp[i])*Iq)/Tq0p[i]; */
    row[0] = idx + 1;
    col[0] = idx + 1;       col[1] = idx+5;
    val[0] = 1/Tq0p[i]; val[1] = -(Xq[i] - Xqp[i])/Tq0p[i];
    ierr   = MatSetValues(J,1,row,2,col,val,INSERT_VALUES);CHKERRQ(ierr);

    /*    fgen[idx+2] = - w + w_s; */
    row[0] = idx + 2;
    col[0] = idx + 2; col[1] = idx + 3;
    val[0] = 0;       val[1] = -1;
    ierr   = MatSetValues(J,1,row,2,col,val,INSERT_VALUES);CHKERRQ(ierr);

    /*    fgen[idx+3] = (-TM[i] + Edp*Id + Eqp*Iq + (Xqp[i] - Xdp[i])*Id*Iq + D[i]*(w - w_s))/M[i]; */
    row[0] = idx + 3;
    col[0] = idx; col[1] = idx + 1; col[2] = idx + 3;       col[3] = idx + 4;                  col[4] = idx + 5;
    val[0] = Iq/M[i];  val[1] = Id/M[i];      val[2] = D[i]/M[i]; val[3] = (Edp + (Xqp[i]-Xdp[i])*Iq)/M[i]; val[4] = (Eqp + (Xqp[i] - Xdp[i])*Id)/M[i];
    ierr   = MatSetValues(J,1,row,5,col,val,INSERT_VALUES);CHKERRQ(ierr);

    Vr   = xnet[2*gbus[i]]; /* Real part of generator terminal voltage */
    Vi   = xnet[2*gbus[i]+1]; /* Imaginary part of the generator terminal voltage */
    ierr = ri2dq(Vr,Vi,delta,&Vd,&Vq);CHKERRQ(ierr);

    det = Rs[i]*Rs[i] + Xdp[i]*Xqp[i];

    Zdq_inv[0] = Rs[i]/det;
    Zdq_inv[1] = Xqp[i]/det;
    Zdq_inv[2] = -Xdp[i]/det;
    Zdq_inv[3] = Rs[i]/det;

    dVd_dVr    = PetscSinScalar(delta); dVd_dVi = -PetscCosScalar(delta);
    dVq_dVr    = PetscCosScalar(delta); dVq_dVi = PetscSinScalar(delta);
    dVd_ddelta = Vr*PetscCosScalar(delta) + Vi*PetscSinScalar(delta);
    dVq_ddelta = -Vr*PetscSinScalar(delta) + Vi*PetscCosScalar(delta);

    /*    fgen[idx+4] = Zdq_inv[0]*(-Edp + Vd) + Zdq_inv[1]*(-Eqp + Vq) + Id; */
    row[0] = idx+4;
    col[0] = idx;         col[1] = idx+1;        col[2] = idx + 2;
    val[0] = -Zdq_inv[1]; val[1] = -Zdq_inv[0];  val[2] = Zdq_inv[0]*dVd_ddelta + Zdq_inv[1]*dVq_ddelta;
    col[3] = idx + 4; col[4] = net_start+2*gbus[i];                     col[5] = net_start + 2*gbus[i]+1;
    val[3] = 1;       val[4] = Zdq_inv[0]*dVd_dVr + Zdq_inv[1]*dVq_dVr; val[5] = Zdq_inv[0]*dVd_dVi + Zdq_inv[1]*dVq_dVi;
    ierr   = MatSetValues(J,1,row,6,col,val,INSERT_VALUES);CHKERRQ(ierr);

    /*  fgen[idx+5] = Zdq_inv[2]*(-Edp + Vd) + Zdq_inv[3]*(-Eqp + Vq) + Iq; */
    row[0] = idx+5;
    col[0] = idx;         col[1] = idx+1;        col[2] = idx + 2;
    val[0] = -Zdq_inv[3]; val[1] = -Zdq_inv[2];  val[2] = Zdq_inv[2]*dVd_ddelta + Zdq_inv[3]*dVq_ddelta;
    col[3] = idx + 5; col[4] = net_start+2*gbus[i];                     col[5] = net_start + 2*gbus[i]+1;
    val[3] = 1;       val[4] = Zdq_inv[2]*dVd_dVr + Zdq_inv[3]*dVq_dVr; val[5] = Zdq_inv[2]*dVd_dVi + Zdq_inv[3]*dVq_dVi;
    ierr   = MatSetValues(J,1,row,6,col,val,INSERT_VALUES);CHKERRQ(ierr);

    dIGr_ddelta = Id*PetscCosScalar(delta) - Iq*PetscSinScalar(delta);
    dIGi_ddelta = Id*PetscSinScalar(delta) + Iq*PetscCosScalar(delta);
    dIGr_dId    = PetscSinScalar(delta);  dIGr_dIq = PetscCosScalar(delta);
    dIGi_dId    = -PetscCosScalar(delta); dIGi_dIq = PetscSinScalar(delta);

    /* fnet[2*gbus[i]]   -= IGi; */
    row[0] = net_start + 2*gbus[i];
    col[0] = idx+2;        col[1] = idx + 4;   col[2] = idx + 5;
    val[0] = -dIGi_ddelta; val[1] = -dIGi_dId; val[2] = -dIGi_dIq;
    ierr = MatSetValues(J,1,row,3,col,val,INSERT_VALUES);CHKERRQ(ierr);

    /* fnet[2*gbus[i]+1]   -= IGr; */
    row[0] = net_start + 2*gbus[i]+1;
    col[0] = idx+2;        col[1] = idx + 4;   col[2] = idx + 5;
    val[0] = -dIGr_ddelta; val[1] = -dIGr_dId; val[2] = -dIGr_dIq;
    ierr   = MatSetValues(J,1,row,3,col,val,INSERT_VALUES);CHKERRQ(ierr);

    Vm = PetscSqrtScalar(Vd*Vd + Vq*Vq); Vm2 = Vm*Vm;

    /*    fgen[idx+6] = (KE[i]*Efd + SE - VR)/TE[i]; */
    /*    SE  = k1[i]*PetscExpScalar(k2[i]*Efd); */

    dSE_dEfd = k1[i]*k2[i]*PetscExpScalar(k2[i]*Efd);

    row[0] = idx + 6;
    col[0] = idx + 6;                     col[1] = idx + 8;
    val[0] = (KE[i] + dSE_dEfd)/TE[i];  val[1] = -1/TE[i];
    ierr   = MatSetValues(J,1,row,2,col,val,INSERT_VALUES);CHKERRQ(ierr);

    /* Exciter differential equations */

    /*    fgen[idx+7] = (RF - KF[i]*Efd/TF[i])/TF[i]; */
    row[0] = idx + 7;
    col[0] = idx + 6;       col[1] = idx + 7;
    val[0] = (-KF[i]/TF[i])/TF[i];  val[1] = 1/TF[i];
    ierr   = MatSetValues(J,1,row,2,col,val,INSERT_VALUES);CHKERRQ(ierr);

    /*    fgen[idx+8] = (VR - KA[i]*RF + KA[i]*KF[i]*Efd/TF[i] - KA[i]*(Vref[i] - Vm))/TA[i]; */
    /* Vm = (Vd^2 + Vq^2)^0.5; */

    dVm_dVd = Vd/Vm; dVm_dVq = Vq/Vm;
    dVm_dVr = dVm_dVd*dVd_dVr + dVm_dVq*dVq_dVr;
    dVm_dVi = dVm_dVd*dVd_dVi + dVm_dVq*dVq_dVi;
    row[0]  = idx + 8;
    col[0]  = idx + 6;           col[1] = idx + 7; col[2] = idx + 8;
    val[0]  = (KA[i]*KF[i]/TF[i])/TA[i]; val[1] = -KA[i]/TA[i];  val[2] = 1/TA[i];
    col[3]  = net_start + 2*gbus[i]; col[4] = net_start + 2*gbus[i]+1;
    val[3]  = KA[i]*dVm_dVr/TA[i];         val[4] = KA[i]*dVm_dVi/TA[i];
    ierr    = MatSetValues(J,1,row,5,col,val,INSERT_VALUES);CHKERRQ(ierr);
    idx     = idx + 9;
  }

  for (i=0; i<nbus; i++) {
    ierr   = MatGetRow(user->Ybus,2*i,&ncols,&cols,&yvals);CHKERRQ(ierr);
    row[0] = net_start + 2*i;
    for (k=0; k<ncols; k++) {
      col[k] = net_start + cols[k];
      val[k] = yvals[k];
    }
    ierr = MatSetValues(J,1,row,ncols,col,val,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatRestoreRow(user->Ybus,2*i,&ncols,&cols,&yvals);CHKERRQ(ierr);

    ierr   = MatGetRow(user->Ybus,2*i+1,&ncols,&cols,&yvals);CHKERRQ(ierr);
    row[0] = net_start + 2*i+1;
    for (k=0; k<ncols; k++) {
      col[k] = net_start + cols[k];
      val[k] = yvals[k];
    }
    ierr = MatSetValues(J,1,row,ncols,col,val,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatRestoreRow(user->Ybus,2*i+1,&ncols,&cols,&yvals);CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(J,MAT_FLUSH_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(J,MAT_FLUSH_ASSEMBLY);CHKERRQ(ierr);

  ierr = VecGetArray(user->V0,&v0);CHKERRQ(ierr);
  for (i=0; i < nload; i++) {
    Vr      = xnet[2*lbus[i]]; /* Real part of load bus voltage */
    Vi      = xnet[2*lbus[i]+1]; /* Imaginary part of the load bus voltage */
    Vm      = PetscSqrtScalar(Vr*Vr + Vi*Vi); Vm2 = Vm*Vm; Vm4 = Vm2*Vm2;
    Vm0     = PetscSqrtScalar(v0[2*lbus[i]]*v0[2*lbus[i]] + v0[2*lbus[i]+1]*v0[2*lbus[i]+1]);
    PD      = QD = 0.0;
    dPD_dVr = dPD_dVi = dQD_dVr = dQD_dVi = 0.0;
    for (k=0; k < ld_nsegsp[i]; k++) {
      PD      += ld_alphap[k]*PD0[i]*PetscPowScalar((Vm/Vm0),ld_betap[k]);
      dPD_dVr += ld_alphap[k]*ld_betap[k]*PD0[i]*PetscPowScalar((1/Vm0),ld_betap[k])*Vr*PetscPowScalar(Vm,(ld_betap[k]-2));
      dPD_dVi += ld_alphap[k]*ld_betap[k]*PD0[i]*PetscPowScalar((1/Vm0),ld_betap[k])*Vi*PetscPowScalar(Vm,(ld_betap[k]-2));
    }
    for (k=0; k < ld_nsegsq[i]; k++) {
      QD      += ld_alphaq[k]*QD0[i]*PetscPowScalar((Vm/Vm0),ld_betaq[k]);
      dQD_dVr += ld_alphaq[k]*ld_betaq[k]*QD0[i]*PetscPowScalar((1/Vm0),ld_betaq[k])*Vr*PetscPowScalar(Vm,(ld_betaq[k]-2));
      dQD_dVi += ld_alphaq[k]*ld_betaq[k]*QD0[i]*PetscPowScalar((1/Vm0),ld_betaq[k])*Vi*PetscPowScalar(Vm,(ld_betaq[k]-2));
    }

    /*    IDr = (PD*Vr + QD*Vi)/Vm2; */
    /*    IDi = (-QD*Vr + PD*Vi)/Vm2; */

    dIDr_dVr = (dPD_dVr*Vr + dQD_dVr*Vi + PD)/Vm2 - ((PD*Vr + QD*Vi)*2*Vr)/Vm4;
    dIDr_dVi = (dPD_dVi*Vr + dQD_dVi*Vi + QD)/Vm2 - ((PD*Vr + QD*Vi)*2*Vi)/Vm4;

    dIDi_dVr = (-dQD_dVr*Vr + dPD_dVr*Vi - QD)/Vm2 - ((-QD*Vr + PD*Vi)*2*Vr)/Vm4;
    dIDi_dVi = (-dQD_dVi*Vr + dPD_dVi*Vi + PD)/Vm2 - ((-QD*Vr + PD*Vi)*2*Vi)/Vm4;


    /*    fnet[2*lbus[i]]   += IDi; */
    row[0] = net_start + 2*lbus[i];
    col[0] = net_start + 2*lbus[i];  col[1] = net_start + 2*lbus[i]+1;
    val[0] = dIDi_dVr;               val[1] = dIDi_dVi;
    ierr   = MatSetValues(J,1,row,2,col,val,ADD_VALUES);CHKERRQ(ierr);
    /*    fnet[2*lbus[i]+1] += IDr; */
    row[0] = net_start + 2*lbus[i]+1;
    col[0] = net_start + 2*lbus[i];  col[1] = net_start + 2*lbus[i]+1;
    val[0] = dIDr_dVr;               val[1] = dIDr_dVi;
    ierr   = MatSetValues(J,1,row,2,col,val,ADD_VALUES);CHKERRQ(ierr);
  }
  ierr = VecRestoreArray(user->V0,&v0);CHKERRQ(ierr);

  ierr = VecRestoreArray(Xgen,&xgen);CHKERRQ(ierr);
  ierr = VecRestoreArray(Xnet,&xnet);CHKERRQ(ierr);

  ierr = DMCompositeRestoreLocalVectors(user->dmpgrid,&Xgen,&Xnet);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

int main(int argc,char **argv)
{
  EPS            eps;
  EPSType        type;
  PetscErrorCode ierr;
  PetscMPIInt    size;
  Userctx        user;
  PetscViewer    Xview,Ybusview;
  Vec            X,Xr,Xi;
  Mat            J,Jred=NULL;
  IS             is0,is1;
  PetscInt       i,*idx2,its,nev,nconv;
  PetscReal      error,re,im;
  PetscScalar    kr,ki;
  PetscBool      terse;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  if (size > 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Only for sequential runs");
  /* show detailed info unless -terse option is given by user */
  ierr = PetscOptionsHasName(NULL,NULL,"-terse",&terse);CHKERRQ(ierr);

  user.neqs_gen   = 9*ngen; /* # eqs. for generator subsystem */
  user.neqs_net   = 2*nbus; /* # eqs. for network subsystem   */
  user.neqs_pgrid = user.neqs_gen + user.neqs_net;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nStability analysis in a network with %D buses and %D generators\n\n",nbus,ngen);CHKERRQ(ierr);

  /* Create indices for differential and algebraic equations */
  ierr = PetscMalloc1(7*ngen,&idx2);CHKERRQ(ierr);
  for (i=0; i<ngen; i++) {
    idx2[7*i]   = 9*i;   idx2[7*i+1] = 9*i+1; idx2[7*i+2] = 9*i+2; idx2[7*i+3] = 9*i+3;
    idx2[7*i+4] = 9*i+6; idx2[7*i+5] = 9*i+7; idx2[7*i+6] = 9*i+8;
  }
  ierr = ISCreateGeneral(PETSC_COMM_WORLD,7*ngen,idx2,PETSC_COPY_VALUES,&user.is_diff);CHKERRQ(ierr);
  ierr = ISComplement(user.is_diff,0,user.neqs_pgrid,&user.is_alg);CHKERRQ(ierr);
  ierr = PetscFree(idx2);CHKERRQ(ierr);

  /* Read initial voltage vector and Ybus */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"X.bin",FILE_MODE_READ,&Xview);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Ybus.bin",FILE_MODE_READ,&Ybusview);CHKERRQ(ierr);

  ierr = VecCreate(PETSC_COMM_WORLD,&user.V0);CHKERRQ(ierr);
  ierr = VecSetSizes(user.V0,PETSC_DECIDE,user.neqs_net);CHKERRQ(ierr);
  ierr = VecLoad(user.V0,Xview);CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD,&user.Ybus);CHKERRQ(ierr);
  ierr = MatSetSizes(user.Ybus,PETSC_DECIDE,PETSC_DECIDE,user.neqs_net,user.neqs_net);CHKERRQ(ierr);
  ierr = MatSetType(user.Ybus,MATBAIJ);CHKERRQ(ierr);
  /*  ierr = MatSetBlockSize(user.Ybus,2);CHKERRQ(ierr); */
  ierr = MatLoad(user.Ybus,Ybusview);CHKERRQ(ierr);

  ierr = PetscViewerDestroy(&Xview);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&Ybusview);CHKERRQ(ierr);

  /* Create DMs for generator and network subsystems */
  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,user.neqs_gen,1,1,NULL,&user.dmgen);CHKERRQ(ierr);
  ierr = DMSetOptionsPrefix(user.dmgen,"dmgen_");CHKERRQ(ierr);
  ierr = DMSetFromOptions(user.dmgen);CHKERRQ(ierr);
  ierr = DMSetUp(user.dmgen);CHKERRQ(ierr);
  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,user.neqs_net,1,1,NULL,&user.dmnet);CHKERRQ(ierr);
  ierr = DMSetOptionsPrefix(user.dmnet,"dmnet_");CHKERRQ(ierr);
  ierr = DMSetFromOptions(user.dmnet);CHKERRQ(ierr);
  ierr = DMSetUp(user.dmnet);CHKERRQ(ierr);

  /* Create a composite DM packer and add the two DMs */
  ierr = DMCompositeCreate(PETSC_COMM_WORLD,&user.dmpgrid);CHKERRQ(ierr);
  ierr = DMSetOptionsPrefix(user.dmpgrid,"pgrid_");CHKERRQ(ierr);
  ierr = DMCompositeAddDM(user.dmpgrid,user.dmgen);CHKERRQ(ierr);
  ierr = DMCompositeAddDM(user.dmpgrid,user.dmnet);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(user.dmpgrid,&X);CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD,&J);CHKERRQ(ierr);
  ierr = MatSetSizes(J,PETSC_DECIDE,PETSC_DECIDE,user.neqs_pgrid,user.neqs_pgrid);CHKERRQ(ierr);
  ierr = MatSetFromOptions(J);CHKERRQ(ierr);
  ierr = PreallocateJacobian(J,&user);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set initial conditions
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = SetInitialGuess(X,&user);CHKERRQ(ierr);

  /* Form Jacobian */
  ierr = ResidualJacobian(X,J,(void*)&user);CHKERRQ(ierr);
  ierr = MatScale(J,-1);CHKERRQ(ierr);
  is0 = user.is_diff;
  is1 = user.is_alg;

  ierr = MatGetSchurComplement(J,is1,is1,is0,is0,MAT_IGNORE_MATRIX,NULL,MAT_SCHUR_COMPLEMENT_AINV_DIAG,MAT_INITIAL_MATRIX,&Jred);CHKERRQ(ierr);

  if (!terse) {
    ierr = MatView(Jred,NULL);CHKERRQ(ierr);
  }

  ierr = MatCreateVecs(Jred,NULL,&Xr);CHKERRQ(ierr);
  ierr = MatCreateVecs(Jred,NULL,&Xi);CHKERRQ(ierr);

  /* Create the eigensolver and set the various options */
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
  ierr = EPSSetOperators(eps,Jred,NULL);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_NHEP);CHKERRQ(ierr);
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

  /* Solve the eigenvalue problem */
  ierr = EPSSolve(eps);CHKERRQ(ierr);

  ierr = EPSGetIterationNumber(eps,&its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the eigensolver: %D\n",its);CHKERRQ(ierr);
  ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n", type);CHKERRQ(ierr);
  ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if (terse) {
    ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL);CHKERRQ(ierr);
  } else {
    /* Get number of converged approximate eigenpairs */
    ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %D\n\n",nconv);CHKERRQ(ierr);

    if (nconv>0) {
      /* Display eigenvalues and relative errors */
      ierr = PetscPrintf(PETSC_COMM_WORLD,
           "           k          ||Ax-kx||/||kx||\n"
           "   ----------------- ------------------\n");CHKERRQ(ierr);

      for (i=0;i<nconv;i++) {
        /* Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and
          ki (imaginary part) */
        ierr = EPSGetEigenpair(eps,i,&kr,&ki,Xr,Xi);CHKERRQ(ierr);
        /* Compute the relative error associated to each eigenpair */
        ierr = EPSComputeError(eps,i,EPS_ERROR_RELATIVE,&error);CHKERRQ(ierr);

#if defined(PETSC_USE_COMPLEX)
        re = PetscRealPart(kr);
        im = PetscImaginaryPart(kr);
#else
        re = kr;
        im = ki;
#endif
        if (im!=0.0) {
          ierr = PetscPrintf(PETSC_COMM_WORLD," %9f%+9fi %12g\n",(double)re,(double)im,(double)error);CHKERRQ(ierr);
        } else {
          ierr = PetscPrintf(PETSC_COMM_WORLD,"   %12f       %12g\n",(double)re,(double)error);CHKERRQ(ierr);
        }
      }
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
    }
  }

  /* Free work space */
  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = MatDestroy(&J);CHKERRQ(ierr);
  ierr = MatDestroy(&Jred);CHKERRQ(ierr);
  ierr = MatDestroy(&user.Ybus);CHKERRQ(ierr);
  ierr = VecDestroy(&X);CHKERRQ(ierr);
  ierr = VecDestroy(&Xr);CHKERRQ(ierr);
  ierr = VecDestroy(&Xi);CHKERRQ(ierr);
  ierr = VecDestroy(&user.V0);CHKERRQ(ierr);
  ierr = DMDestroy(&user.dmgen);CHKERRQ(ierr);
  ierr = DMDestroy(&user.dmnet);CHKERRQ(ierr);
  ierr = DMDestroy(&user.dmpgrid);CHKERRQ(ierr);
  ierr = ISDestroy(&user.is_diff);CHKERRQ(ierr);
  ierr = ISDestroy(&user.is_alg);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   build:
      requires: !complex

   test:
      suffix: 1
      args: -terse
      requires: double !complex !define(PETSC_USE_64BIT_INDICES)
      localrunfiles: X.bin Ybus.bin

TEST*/
