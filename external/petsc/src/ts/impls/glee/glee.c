/*
  Code for time stepping with the General Linear with Error Estimation method


  Notes:
  The general system is written as

  Udot = F(t,U)

*/
#include <petsc/private/tsimpl.h>                /*I   "petscts.h"   I*/
#include <petscdm.h>

static PetscBool  cited = PETSC_FALSE;
static const char citation[] =
  "@ARTICLE{Constantinescu_TR2016b,\n"
  " author = {Constantinescu, E.M.},\n"
  " title = {Estimating Global Errors in Time Stepping},\n"
  " journal = {ArXiv e-prints},\n"
  " year = 2016,\n"
  " adsurl = {http://adsabs.harvard.edu/abs/2015arXiv150305166C}\n}\n";

static TSGLEEType   TSGLEEDefaultType = TSGLEE35;
static PetscBool    TSGLEERegisterAllCalled;
static PetscBool    TSGLEEPackageInitialized;
static PetscInt     explicit_stage_time_id;

typedef struct _GLEETableau *GLEETableau;
struct _GLEETableau {
  char      *name;
  PetscInt   order;               /* Classical approximation order of the method i*/
  PetscInt   s;                   /* Number of stages */
  PetscInt   r;                   /* Number of steps */
  PetscReal  gamma;               /* LTE ratio */
  PetscReal *A,*B,*U,*V,*S,*F,*c; /* Tableau */
  PetscReal *Fembed;              /* Embedded final method coefficients */
  PetscReal *Ferror;              /* Coefficients for computing error   */
  PetscReal *Serror;              /* Coefficients for initializing the error   */
  PetscInt   pinterp;             /* Interpolation order */
  PetscReal *binterp;             /* Interpolation coefficients */
  PetscReal  ccfl;                /* Placeholder for CFL coefficient relative to forward Euler  */
};
typedef struct _GLEETableauLink *GLEETableauLink;
struct _GLEETableauLink {
  struct _GLEETableau tab;
  GLEETableauLink     next;
};
static GLEETableauLink GLEETableauList;

typedef struct {
  GLEETableau  tableau;
  Vec          *Y;         /* Solution vector (along with auxiliary solution y~ or eps) */
  Vec          *X;         /* Temporary solution vector */
  Vec          *YStage;    /* Stage values */
  Vec          *YdotStage; /* Stage right hand side */
  Vec          W;          /* Right-hand-side for implicit stage solve */
  Vec          Ydot;       /* Work vector holding Ydot during residual evaluation */
  Vec          yGErr;      /* Vector holding the global error after a step is completed */
  PetscScalar  *swork;     /* Scalar work (size of the number of stages)*/
  PetscScalar  *rwork;     /* Scalar work (size of the number of steps)*/
  PetscReal    scoeff;     /* shift = scoeff/dt */
  PetscReal    stage_time;
  TSStepStatus status;
} TS_GLEE;

/*MC
     TSGLEE23 - Second order three stage GLEE method

     This method has three stages.
     s = 3, r = 2

     Level: advanced

.seealso: TSGLEE
M*/
/*MC
     TSGLEE24 - Second order four stage GLEE method

     This method has four stages.
     s = 4, r = 2

     Level: advanced

.seealso: TSGLEE
M*/
/*MC
     TSGLEE25i - Second order five stage GLEE method

     This method has five stages.
     s = 5, r = 2

     Level: advanced

.seealso: TSGLEE
M*/
/*MC
     TSGLEE35  - Third order five stage GLEE method

     This method has five stages.
     s = 5, r = 2

     Level: advanced

.seealso: TSGLEE
M*/
/*MC
     TSGLEEEXRK2A  - Second order six stage GLEE method

     This method has six stages.
     s = 6, r = 2

     Level: advanced

.seealso: TSGLEE
M*/
/*MC
     TSGLEERK32G1  - Third order eight stage GLEE method

     This method has eight stages.
     s = 8, r = 2

     Level: advanced

.seealso: TSGLEE
M*/
/*MC
     TSGLEERK285EX  - Second order nine stage GLEE method

     This method has nine stages.
     s = 9, r = 2

     Level: advanced

.seealso: TSGLEE
M*/

/*@C
  TSGLEERegisterAll - Registers all of the General Linear with Error Estimation methods in TSGLEE

  Not Collective, but should be called by all processes which will need the schemes to be registered

  Level: advanced

.seealso:  TSGLEERegisterDestroy()
@*/
PetscErrorCode TSGLEERegisterAll(void)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (TSGLEERegisterAllCalled) PetscFunctionReturn(0);
  TSGLEERegisterAllCalled = PETSC_TRUE;

  {
#define GAMMA 0.5
    /* y-eps form */
    const PetscInt
      p = 1,
      s = 3,
      r = 2;
    const PetscReal
      A[3][3]   = {{1.0,0,0},{0,0.5,0},{0,0.5,0.5}},
      B[2][3]   = {{1.0,0,0},{-2.0,1.0,1.0}},
      U[3][2]   = {{1.0,0},{1.0,0.5},{1.0,0.5}},
      V[2][2]   = {{1,0},{0,1}},
      S[2]      = {1,0},
      F[2]      = {1,0},
      Fembed[2] = {1,1-GAMMA},
      Ferror[2] = {0,1},
      Serror[2] = {1,0};
    ierr = TSGLEERegister(TSGLEEi1,p,s,r,GAMMA,&A[0][0],&B[0][0],&U[0][0],&V[0][0],S,F,NULL,Fembed,Ferror,Serror,0,NULL);CHKERRQ(ierr);
  }
  {
#undef GAMMA
#define GAMMA 0.0
    /* y-eps form */
    const PetscInt
      p = 2,
      s = 3,
      r = 2;
    const PetscReal
      A[3][3]   = {{0,0,0},{1,0,0},{0.25,0.25,0}},
      B[2][3]   = {{1.0/12.0,1.0/12.0,5.0/6.0},{1.0/12.0,1.0/12.0,-1.0/6.0}},
      U[3][2]   = {{1,0},{1,10},{1,-1}},
      V[2][2]   = {{1,0},{0,1}},
      S[2]      = {1,0},
      F[2]      = {1,0},
      Fembed[2] = {1,1-GAMMA},
      Ferror[2] = {0,1},
      Serror[2] = {1,0};
    ierr = TSGLEERegister(TSGLEE23,p,s,r,GAMMA,&A[0][0],&B[0][0],&U[0][0],&V[0][0],S,F,NULL,Fembed,Ferror,Serror,0,NULL);CHKERRQ(ierr);
  }
  {
#undef GAMMA
#define GAMMA 0.0
    /* y-y~ form */
    const PetscInt
      p = 2,
      s = 4,
      r = 2;
    const PetscReal
      A[4][4]   = {{0,0,0,0},{0.75,0,0,0},{0.25,29.0/60.0,0,0},{-21.0/44.0,145.0/44.0,-20.0/11.0,0}},
      B[2][4]   = {{109.0/275.0,58.0/75.0,-37.0/110.0,1.0/6.0},{3.0/11.0,0,75.0/88.0,-1.0/8.0}},
      U[4][2]   = {{0,1},{75.0/58.0,-17.0/58.0},{0,1},{0,1}},
      V[2][2]   = {{1,0},{0,1}},
      S[2]      = {1,1},
      F[2]      = {1,0},
      Fembed[2] = {0,1},
      Ferror[2] = {-1.0/(1.0-GAMMA),1.0/(1.0-GAMMA)},
      Serror[2] = {1.0-GAMMA,1.0};
      ierr = TSGLEERegister(TSGLEE24,p,s,r,GAMMA,&A[0][0],&B[0][0],&U[0][0],&V[0][0],S,F,NULL,Fembed,Ferror,Serror,0,NULL);CHKERRQ(ierr);
  }
  {
#undef GAMMA
#define GAMMA 0.0
    /* y-y~ form */
    const PetscInt
      p = 2,
      s = 5,
      r = 2;
    const PetscReal
      A[5][5]   = {{0,0,0,0,0},
                   {-0.94079244066783383269,0,0,0,0},
                   {0.64228187778301907108,0.10915356933958500042,0,0,0},
                   {-0.51764297742287450812,0.74414270351096040738,-0.71404164927824538121,0,0},
                   {-0.44696561556825969206,-0.76768425657590196518,0.20111608138142987881,0.93828186737840469796,0}},
      B[2][5]   = {{-0.029309178948150356153,-0.49671981884013874923,0.34275801517650053274,0.32941112623949194988,0.85385985637229662276},
                   {0.78133219686062535272,0.074238691892675897635,0.57957363498384957966,-0.24638502829674959968,-0.18875949544040123033}},
      U[5][2]   = {{0.16911424754448327735,0.83088575245551672265},
                   {0.53638465733199574340,0.46361534266800425660},
                   {0.39901579167169582526,0.60098420832830417474},
                   {0.87689005530618575480,0.12310994469381424520},
                   {0.99056100455550913009,0.0094389954444908699092}},
      V[2][2]   = {{1,0},{0,1}},
      S[2]      = {1,1},
      F[2]      = {1,0},
      Fembed[2] = {0,1},
      Ferror[2] = {-1.0/(1.0-GAMMA),1.0/(1.0-GAMMA)},
      Serror[2] = {1.0-GAMMA,1.0};
    ierr = TSGLEERegister(TSGLEE25I,p,s,r,GAMMA,&A[0][0],&B[0][0],&U[0][0],&V[0][0],S,F,NULL,Fembed,Ferror,Serror,0,NULL);CHKERRQ(ierr);
  }
  {
#undef GAMMA
#define GAMMA 0.0
    /* y-y~ form */
    const PetscInt
      p = 3,
      s = 5,
      r = 2;
    const PetscReal
      A[5][5]   = {{0,0,0,0,0},
                   {- 2169604947363702313.0 /  24313474998937147335.0,0,0,0,0},
                   {46526746497697123895.0 /  94116917485856474137.0,-10297879244026594958.0 /  49199457603717988219.0,0,0,0},
                   {23364788935845982499.0 /  87425311444725389446.0,-79205144337496116638.0 / 148994349441340815519.0,40051189859317443782.0 /  36487615018004984309.0,0,0},
                   {42089522664062539205.0 / 124911313006412840286.0,-15074384760342762939.0 / 137927286865289746282.0,-62274678522253371016.0 / 125918573676298591413.0,13755475729852471739.0 /  79257927066651693390.0,0}},
      B[2][5]   = {{61546696837458703723.0 /  56982519523786160813.0,-55810892792806293355.0 / 206957624151308356511.0,24061048952676379087.0 / 158739347956038723465.0,3577972206874351339.0 /   7599733370677197135.0,-59449832954780563947.0 / 137360038685338563670.0},
                   {- 9738262186984159168.0 /  99299082461487742983.0,-32797097931948613195.0 /  61521565616362163366.0,42895514606418420631.0 /  71714201188501437336.0,22608567633166065068.0 /  55371917805607957003.0,94655809487476459565.0 / 151517167160302729021.0}},
      U[5][2]   = {{70820309139834661559.0 /  80863923579509469826.0,10043614439674808267.0 /  80863923579509469826.0},
                   {161694774978034105510.0 / 106187653640211060371.0,-55507121337823045139.0 / 106187653640211060371.0},
                   {78486094644566264568.0 /  88171030896733822981.0,9684936252167558413.0 /  88171030896733822981.0},
                   {65394922146334854435.0 /  84570853840405479554.0,19175931694070625119.0 /  84570853840405479554.0},
                   {8607282770183754108.0 / 108658046436496925911.0,100050763666313171803.0 / 108658046436496925911.0}},
      V[2][2]   = {{1,0},{0,1}},
      S[2]      = {1,1},
      F[2]      = {1,0},
      Fembed[2] = {0,1},
      Ferror[2] = {-1.0/(1.0-GAMMA),1.0/(1.0-GAMMA)},
      Serror[2] = {1.0-GAMMA,1.0};
    ierr = TSGLEERegister(TSGLEE35,p,s,r,GAMMA,&A[0][0],&B[0][0],&U[0][0],&V[0][0],S,F,NULL,Fembed,Ferror,Serror,0,NULL);CHKERRQ(ierr);
  }
  {
#undef GAMMA
#define GAMMA 0.25
    /* y-eps form */
    const PetscInt
      p = 2,
      s = 6,
      r = 2;
    const PetscReal
      A[6][6]   = {{0,0,0,0,0,0},
                   {1,0,0,0,0,0},
                   {0,0,0,0,0,0},
                   {0,0,0.5,0,0,0},
                   {0,0,0.25,0.25,0,0},
                   {0,0,0.25,0.25,0.5,0}},
      B[2][6]   = {{0.5,0.5,0,0,0,0},
                   {-2.0/3.0,-2.0/3.0,1.0/3.0,1.0/3.0,1.0/3.0,1.0/3.0}},
      U[6][2]   = {{1,0},{1,0},{1,0.75},{1,0.75},{1,0.75},{1,0.75}},
      V[2][2]   = {{1,0},{0,1}},
      S[2]      = {1,0},
      F[2]      = {1,0},
      Fembed[2] = {1,1-GAMMA},
      Ferror[2] = {0,1},
      Serror[2] = {1,0};
    ierr = TSGLEERegister(TSGLEEEXRK2A,p,s,r,GAMMA,&A[0][0],&B[0][0],&U[0][0],&V[0][0],S,F,NULL,Fembed,Ferror,Serror,0,NULL);CHKERRQ(ierr);
  }
  {
#undef GAMMA
#define GAMMA 0.0
    /* y-eps form */
    const PetscInt
      p = 3,
      s = 8,
      r = 2;
    const PetscReal
      A[8][8]   = {{0,0,0,0,0,0,0,0},
                   {0.5,0,0,0,0,0,0,0},
                   {-1,2,0,0,0,0,0,0},
                   {1.0/6.0,2.0/3.0,1.0/6.0,0,0,0,0,0},
                   {0,0,0,0,0,0,0,0},
                   {-7.0/24.0,1.0/3.0,1.0/12.0,-1.0/8.0,0.5,0,0,0},
                   {7.0/6.0,-4.0/3.0,-1.0/3.0,0.5,-1.0,2.0,0,0},
                   {0,0,0,0,1.0/6.0,2.0/3.0,1.0/6.0,0}},
      B[2][8]   = {{1.0/6.0,2.0/3.0,1.0/6.0,0,0,0,0,0},
                   {-1.0/6.0,-2.0/3.0,-1.0/6.0,0,1.0/6.0,2.0/3.0,1.0/6.0,0}},
      U[8][2]   = {{1,0},{1,0},{1,0},{1,0},{1,1},{1,1},{1,1},{1,1}},
      V[2][2]   = {{1,0},{0,1}},
      S[2]      = {1,0},
      F[2]      = {1,0},
      Fembed[2] = {1,1-GAMMA},
      Ferror[2] = {0,1},
      Serror[2] = {1,0};
    ierr = TSGLEERegister(TSGLEERK32G1,p,s,r,GAMMA,&A[0][0],&B[0][0],&U[0][0],&V[0][0],S,F,NULL,Fembed,Ferror,Serror,0,NULL);CHKERRQ(ierr);
  }
  {
#undef GAMMA
#define GAMMA 0.25
    /* y-eps form */
    const PetscInt
      p = 2,
      s = 9,
      r = 2;
    const PetscReal
      A[9][9]   = {{0,0,0,0,0,0,0,0,0},
                   {0.585786437626904966,0,0,0,0,0,0,0,0},
                   {0.149999999999999994,0.849999999999999978,0,0,0,0,0,0,0},
                   {0,0,0,0,0,0,0,0,0},
                   {0,0,0,0.292893218813452483,0,0,0,0,0},
                   {0,0,0,0.0749999999999999972,0.424999999999999989,0,0,0,0},
                   {0,0,0,0.176776695296636893,0.176776695296636893,0.146446609406726241,0,0,0},
                   {0,0,0,0.176776695296636893,0.176776695296636893,0.146446609406726241,0.292893218813452483,0,0},
                   {0,0,0,0.176776695296636893,0.176776695296636893,0.146446609406726241,0.0749999999999999972,0.424999999999999989,0}},
      B[2][9]   = {{0.353553390593273786,0.353553390593273786,0.292893218813452483,0,0,0,0,0,0},
                   {-0.471404520791031678,-0.471404520791031678,-0.390524291751269959,0.235702260395515839,0.235702260395515839,0.195262145875634979,0.235702260395515839,0.235702260395515839,0.195262145875634979}},
      U[9][2]   = {{1,0},{1,0},{1,0},{1,0.75},{1,0.75},{1,0.75},{1,0.75},{1,0.75},{1,0.75}},
      V[2][2]   = {{1,0},{0,1}},
      S[2]      = {1,0},
      F[2]      = {1,0},
      Fembed[2] = {1,1-GAMMA},
      Ferror[2] = {0,1},
      Serror[2] = {1,0};
    ierr = TSGLEERegister(TSGLEERK285EX,p,s,r,GAMMA,&A[0][0],&B[0][0],&U[0][0],&V[0][0],S,F,NULL,Fembed,Ferror,Serror,0,NULL);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

/*@C
   TSGLEERegisterDestroy - Frees the list of schemes that were registered by TSGLEERegister().

   Not Collective

   Level: advanced

.seealso: TSGLEERegister(), TSGLEERegisterAll()
@*/
PetscErrorCode TSGLEERegisterDestroy(void)
{
  PetscErrorCode ierr;
  GLEETableauLink link;

  PetscFunctionBegin;
  while ((link = GLEETableauList)) {
    GLEETableau t = &link->tab;
    GLEETableauList = link->next;
    ierr = PetscFree5(t->A,t->B,t->U,t->V,t->c);CHKERRQ(ierr);
    ierr = PetscFree2(t->S,t->F);               CHKERRQ(ierr);
    ierr = PetscFree (t->Fembed);               CHKERRQ(ierr);
    ierr = PetscFree (t->Ferror);               CHKERRQ(ierr);
    ierr = PetscFree (t->Serror);               CHKERRQ(ierr);
    ierr = PetscFree (t->binterp);              CHKERRQ(ierr);
    ierr = PetscFree (t->name);                 CHKERRQ(ierr);
    ierr = PetscFree (link);                    CHKERRQ(ierr);
  }
  TSGLEERegisterAllCalled = PETSC_FALSE;
  PetscFunctionReturn(0);
}

/*@C
  TSGLEEInitializePackage - This function initializes everything in the TSGLEE package. It is called
  from TSInitializePackage().

  Level: developer

.seealso: PetscInitialize()
@*/
PetscErrorCode TSGLEEInitializePackage(void)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (TSGLEEPackageInitialized) PetscFunctionReturn(0);
  TSGLEEPackageInitialized = PETSC_TRUE;
  ierr = TSGLEERegisterAll();CHKERRQ(ierr);
  ierr = PetscObjectComposedDataRegister(&explicit_stage_time_id);CHKERRQ(ierr);
  ierr = PetscRegisterFinalize(TSGLEEFinalizePackage);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
  TSGLEEFinalizePackage - This function destroys everything in the TSGLEE package. It is
  called from PetscFinalize().

  Level: developer

.seealso: PetscFinalize()
@*/
PetscErrorCode TSGLEEFinalizePackage(void)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  TSGLEEPackageInitialized = PETSC_FALSE;
  ierr = TSGLEERegisterDestroy();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   TSGLEERegister - register an GLEE scheme by providing the entries in the Butcher tableau

   Not Collective, but the same schemes should be registered on all processes on which they will be used

   Input Parameters:
+  name   - identifier for method
.  order  - order of method
.  s - number of stages
.  r - number of steps
.  gamma - LTE ratio
.  A - stage coefficients (dimension s*s, row-major)
.  B - step completion coefficients (dimension r*s, row-major)
.  U - method coefficients (dimension s*r, row-major)
.  V - method coefficients (dimension r*r, row-major)
.  S - starting coefficients
.  F - finishing coefficients
.  c - abscissa (dimension s; NULL to use row sums of A)
.  Fembed - step completion coefficients for embedded method
.  Ferror - error computation coefficients
.  Serror - error initialization coefficients
.  pinterp - order of interpolation (0 if unavailable)
-  binterp - array of interpolation coefficients (NULL if unavailable)

   Notes:
   Several GLEE methods are provided, this function is only needed to create new methods.

   Level: advanced

.seealso: TSGLEE
@*/
PetscErrorCode TSGLEERegister(TSGLEEType name,PetscInt order,PetscInt s, PetscInt r,
                              PetscReal gamma,
                              const PetscReal A[],const PetscReal B[],
                              const PetscReal U[],const PetscReal V[],
                              const PetscReal S[],const PetscReal F[],
                              const PetscReal c[],
                              const PetscReal Fembed[],const PetscReal Ferror[],
                              const PetscReal Serror[],
                              PetscInt pinterp, const PetscReal binterp[])
{
  PetscErrorCode    ierr;
  GLEETableauLink   link;
  GLEETableau       t;
  PetscInt          i,j;

  PetscFunctionBegin;
  ierr     = TSGLEEInitializePackage();CHKERRQ(ierr);
  ierr     = PetscNew(&link);CHKERRQ(ierr);
  t        = &link->tab;
  ierr     = PetscStrallocpy(name,&t->name);CHKERRQ(ierr);
  t->order = order;
  t->s     = s;
  t->r     = r;
  t->gamma = gamma;
  ierr     = PetscMalloc5(s*s,&t->A,r*r,&t->V,s,&t->c,r*s,&t->B,s*r,&t->U);CHKERRQ(ierr);
  ierr     = PetscMalloc2(r,&t->S,r,&t->F);CHKERRQ(ierr);
  ierr     = PetscArraycpy(t->A,A,s*s);CHKERRQ(ierr);
  ierr     = PetscArraycpy(t->B,B,r*s);CHKERRQ(ierr);
  ierr     = PetscArraycpy(t->U,U,s*r);CHKERRQ(ierr);
  ierr     = PetscArraycpy(t->V,V,r*r);CHKERRQ(ierr);
  ierr     = PetscArraycpy(t->S,S,r  );CHKERRQ(ierr);
  ierr     = PetscArraycpy(t->F,F,r  );CHKERRQ(ierr);
  if (c) {
    ierr   = PetscArraycpy(t->c,c,s);CHKERRQ(ierr);
  } else {
    for (i=0; i<s; i++) for (j=0,t->c[i]=0; j<s; j++) t->c[i] += A[i*s+j];
  }
  ierr = PetscMalloc1(r,&t->Fembed);CHKERRQ(ierr);
  ierr = PetscMalloc1(r,&t->Ferror);CHKERRQ(ierr);
  ierr = PetscMalloc1(r,&t->Serror);CHKERRQ(ierr);
  ierr = PetscArraycpy(t->Fembed,Fembed,r);CHKERRQ(ierr);
  ierr = PetscArraycpy(t->Ferror,Ferror,r);CHKERRQ(ierr);
  ierr = PetscArraycpy(t->Serror,Serror,r);CHKERRQ(ierr);
  t->pinterp = pinterp;
  ierr       = PetscMalloc1(s*pinterp,&t->binterp);CHKERRQ(ierr);
  ierr       = PetscArraycpy(t->binterp,binterp,s*pinterp);CHKERRQ(ierr);

  link->next      = GLEETableauList;
  GLEETableauList = link;
  PetscFunctionReturn(0);
}

static PetscErrorCode TSEvaluateStep_GLEE(TS ts,PetscInt order,Vec X,PetscBool *done)
{
  TS_GLEE         *glee = (TS_GLEE*) ts->data;
  GLEETableau     tab = glee->tableau;
  PetscReal       h, *B = tab->B, *V = tab->V,
                  *F = tab->F,
                  *Fembed = tab->Fembed;
  PetscInt        s = tab->s, r = tab->r, i, j;
  Vec             *Y = glee->Y, *YdotStage = glee->YdotStage;
  PetscScalar     *ws = glee->swork, *wr = glee->rwork;
  PetscErrorCode  ierr;

  PetscFunctionBegin;

  switch (glee->status) {
    case TS_STEP_INCOMPLETE:
    case TS_STEP_PENDING:
      h = ts->time_step; break;
    case TS_STEP_COMPLETE:
      h = ts->ptime - ts->ptime_prev; break;
    default: SETERRQ(PetscObjectComm((PetscObject)ts),PETSC_ERR_PLIB,"Invalid TSStepStatus");
  }

  if (order == tab->order) {

    /* Note: Irrespective of whether status is TS_STEP_INCOMPLETE
             or TS_STEP_COMPLETE, glee->X has the solution at the
             beginning of the time step. So no need to roll-back.
    */
    if (glee->status == TS_STEP_INCOMPLETE) {
      for (i=0; i<r; i++) {
        ierr = VecZeroEntries(Y[i]);CHKERRQ(ierr);
        for (j=0; j<r; j++) wr[j] = V[i*r+j];
        ierr = VecMAXPY(Y[i],r,wr,glee->X);CHKERRQ(ierr);
        for (j=0; j<s; j++) ws[j] = h*B[i*s+j];
        ierr = VecMAXPY(Y[i],s,ws,YdotStage);CHKERRQ(ierr);
      }
      ierr = VecZeroEntries(X);CHKERRQ(ierr);
      for (j=0; j<r; j++) wr[j] = F[j];
      ierr = VecMAXPY(X,r,wr,Y);CHKERRQ(ierr);
    } else {ierr = VecCopy(ts->vec_sol,X);CHKERRQ(ierr);}
    PetscFunctionReturn(0);

  } else if (order == tab->order-1) {

    /* Complete with the embedded method (Fembed) */
    for (i=0; i<r; i++) {
      ierr = VecZeroEntries(Y[i]);CHKERRQ(ierr);
      for (j=0; j<r; j++) wr[j] = V[i*r+j];
      ierr = VecMAXPY(Y[i],r,wr,glee->X);CHKERRQ(ierr);
      for (j=0; j<s; j++) ws[j] = h*B[i*s+j];
      ierr = VecMAXPY(Y[i],s,ws,YdotStage);CHKERRQ(ierr);
    }
    ierr = VecZeroEntries(X);CHKERRQ(ierr);
    for (j=0; j<r; j++) wr[j] = Fembed[j];
    ierr = VecMAXPY(X,r,wr,Y);CHKERRQ(ierr);

    if (done) *done = PETSC_TRUE;
    PetscFunctionReturn(0);
  }
  if (done) *done = PETSC_FALSE;
  else SETERRQ3(PetscObjectComm((PetscObject)ts),PETSC_ERR_SUP,"GLEE '%s' of order %D cannot evaluate step at order %D",tab->name,tab->order,order);
  PetscFunctionReturn(0);
}

static PetscErrorCode TSStep_GLEE(TS ts)
{
  TS_GLEE         *glee = (TS_GLEE*)ts->data;
  GLEETableau     tab = glee->tableau;
  const PetscInt  s = tab->s, r = tab->r;
  PetscReal       *A = tab->A, *U = tab->U,
                  *F = tab->F,
                  *c = tab->c;
  Vec             *Y = glee->Y, *X = glee->X,
                  *YStage = glee->YStage, 
                  *YdotStage = glee->YdotStage,
                  W = glee->W;
  SNES            snes;
  PetscScalar     *ws = glee->swork, *wr = glee->rwork;
  TSAdapt         adapt;
  PetscInt        i,j,reject,next_scheme,its,lits;
  PetscReal       next_time_step;
  PetscReal       t;
  PetscBool       accept;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  ierr = PetscCitationsRegister(citation,&cited);CHKERRQ(ierr);

  for (i=0; i<r; i++) { ierr = VecCopy(Y[i],X[i]);CHKERRQ(ierr); }

  ierr           = TSGetSNES(ts,&snes);CHKERRQ(ierr);
  next_time_step = ts->time_step;
  t              = ts->ptime;
  accept         = PETSC_TRUE;
  glee->status   = TS_STEP_INCOMPLETE;

  for (reject=0; reject<ts->max_reject && !ts->reason; reject++,ts->reject++) {

    PetscReal h = ts->time_step;
    ierr = TSPreStep(ts);CHKERRQ(ierr);

    for (i=0; i<s; i++) {

      glee->stage_time = t + h*c[i];
      ierr = TSPreStage(ts,glee->stage_time);CHKERRQ(ierr);

      if (A[i*s+i] == 0) {  /* Explicit stage */
        ierr = VecZeroEntries(YStage[i]);CHKERRQ(ierr);
        for (j=0; j<r; j++) wr[j] = U[i*r+j];
        ierr = VecMAXPY(YStage[i],r,wr,X);CHKERRQ(ierr);
        for (j=0; j<i; j++) ws[j] = h*A[i*s+j];
        ierr = VecMAXPY(YStage[i],i,ws,YdotStage);CHKERRQ(ierr);
      } else {              /* Implicit stage */
        glee->scoeff = 1.0/A[i*s+i];
        /* compute right-hand-side */
        ierr = VecZeroEntries(W);CHKERRQ(ierr);
        for (j=0; j<r; j++) wr[j] = U[i*r+j];
        ierr = VecMAXPY(W,r,wr,X);CHKERRQ(ierr);
        for (j=0; j<i; j++) ws[j] = h*A[i*s+j];
        ierr = VecMAXPY(W,i,ws,YdotStage);CHKERRQ(ierr);
        ierr = VecScale(W,glee->scoeff/h);CHKERRQ(ierr);
        /* set initial guess */
        ierr = VecCopy(i>0 ? YStage[i-1] : ts->vec_sol,YStage[i]);CHKERRQ(ierr);
        /* solve for this stage */
        ierr = SNESSolve(snes,W,YStage[i]);CHKERRQ(ierr);
        ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr);
        ierr = SNESGetLinearSolveIterations(snes,&lits);CHKERRQ(ierr);
        ts->snes_its += its; ts->ksp_its += lits;
      }
      ierr = TSGetAdapt(ts,&adapt);CHKERRQ(ierr);
      ierr = TSAdaptCheckStage(adapt,ts,glee->stage_time,YStage[i],&accept);CHKERRQ(ierr);
      if (!accept) goto reject_step;
      ierr = TSPostStage(ts,glee->stage_time,i,YStage);CHKERRQ(ierr);
      ierr = TSComputeRHSFunction(ts,t+h*c[i],YStage[i],YdotStage[i]);CHKERRQ(ierr);
    }
    ierr = TSEvaluateStep(ts,tab->order,ts->vec_sol,NULL);CHKERRQ(ierr);
    glee->status = TS_STEP_PENDING;

    /* Register only the current method as a candidate because we're not supporting multiple candidates yet. */
    ierr = TSGetAdapt(ts,&adapt);CHKERRQ(ierr);
    ierr = TSAdaptCandidatesClear(adapt);CHKERRQ(ierr);
    ierr = TSAdaptCandidateAdd(adapt,tab->name,tab->order,1,tab->ccfl,(PetscReal)tab->s,PETSC_TRUE);CHKERRQ(ierr);
    ierr = TSAdaptChoose(adapt,ts,ts->time_step,&next_scheme,&next_time_step,&accept);CHKERRQ(ierr);
    if (accept) {
      /* ignore next_scheme for now */
      ts->ptime     += ts->time_step;
      ts->time_step  = next_time_step;
      glee->status = TS_STEP_COMPLETE;
      /* compute and store the global error */
      /* Note: this is not needed if TSAdaptGLEE is not used */
      ierr = TSGetTimeError(ts,0,&(glee->yGErr));CHKERRQ(ierr);
      ierr = PetscObjectComposedDataSetReal((PetscObject)ts->vec_sol,explicit_stage_time_id,ts->ptime);CHKERRQ(ierr);
      break;
    } else {                    /* Roll back the current step */
      for (j=0; j<r; j++) wr[j] = F[j];
      ierr = VecMAXPY(ts->vec_sol,r,wr,X);CHKERRQ(ierr);
      ts->time_step = next_time_step;
      glee->status  = TS_STEP_INCOMPLETE;
    }
reject_step: continue;
  }
  if (glee->status != TS_STEP_COMPLETE && !ts->reason) ts->reason = TS_DIVERGED_STEP_REJECTED;
  PetscFunctionReturn(0);
}

static PetscErrorCode TSInterpolate_GLEE(TS ts,PetscReal itime,Vec X)
{
  TS_GLEE         *glee = (TS_GLEE*)ts->data;
  PetscInt        s=glee->tableau->s, pinterp=glee->tableau->pinterp,i,j;
  PetscReal       h,tt,t;
  PetscScalar     *b;
  const PetscReal *B = glee->tableau->binterp;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  if (!B) SETERRQ1(PetscObjectComm((PetscObject)ts),PETSC_ERR_SUP,"TSGLEE %s does not have an interpolation formula",glee->tableau->name);
  switch (glee->status) {
    case TS_STEP_INCOMPLETE:
    case TS_STEP_PENDING:
      h = ts->time_step;
      t = (itime - ts->ptime)/h;
      break;
    case TS_STEP_COMPLETE:
      h = ts->ptime - ts->ptime_prev;
      t = (itime - ts->ptime)/h + 1; /* In the interval [0,1] */
      break;
    default: SETERRQ(PetscObjectComm((PetscObject)ts),PETSC_ERR_PLIB,"Invalid TSStepStatus");
  }
  ierr = PetscMalloc1(s,&b);CHKERRQ(ierr);
  for (i=0; i<s; i++) b[i] = 0;
  for (j=0,tt=t; j<pinterp; j++,tt*=t) {
    for (i=0; i<s; i++) {
      b[i]  += h * B[i*pinterp+j] * tt;
    }
  }
  ierr = VecCopy(glee->YStage[0],X);CHKERRQ(ierr);
  ierr = VecMAXPY(X,s,b,glee->YdotStage);CHKERRQ(ierr);
  ierr = PetscFree(b);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*------------------------------------------------------------*/
static PetscErrorCode TSReset_GLEE(TS ts)
{
  TS_GLEE        *glee = (TS_GLEE*)ts->data;
  PetscInt       s, r;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!glee->tableau) PetscFunctionReturn(0);
  s    = glee->tableau->s;
  r    = glee->tableau->r;
  ierr = VecDestroyVecs(r,&glee->Y);CHKERRQ(ierr);
  ierr = VecDestroyVecs(r,&glee->X);CHKERRQ(ierr);
  ierr = VecDestroyVecs(s,&glee->YStage);CHKERRQ(ierr);
  ierr = VecDestroyVecs(s,&glee->YdotStage);CHKERRQ(ierr);
  ierr = VecDestroy(&glee->Ydot);CHKERRQ(ierr);
  ierr = VecDestroy(&glee->yGErr);CHKERRQ(ierr);
  ierr = VecDestroy(&glee->W);CHKERRQ(ierr);
  ierr = PetscFree2(glee->swork,glee->rwork);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode TSGLEEGetVecs(TS ts,DM dm,Vec *Ydot)
{
  TS_GLEE     *glee = (TS_GLEE*)ts->data;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (Ydot) {
    if (dm && dm != ts->dm) {
      ierr = DMGetNamedGlobalVector(dm,"TSGLEE_Ydot",Ydot);CHKERRQ(ierr);
    } else *Ydot = glee->Ydot;
  }
  PetscFunctionReturn(0);
}


static PetscErrorCode TSGLEERestoreVecs(TS ts,DM dm,Vec *Ydot)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (Ydot) {
    if (dm && dm != ts->dm) {
      ierr = DMRestoreNamedGlobalVector(dm,"TSGLEE_Ydot",Ydot);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

/*
  This defines the nonlinear equation that is to be solved with SNES
*/
static PetscErrorCode SNESTSFormFunction_GLEE(SNES snes,Vec X,Vec F,TS ts)
{
  TS_GLEE       *glee = (TS_GLEE*)ts->data;
  DM             dm,dmsave;
  Vec            Ydot;
  PetscReal      shift = glee->scoeff / ts->time_step;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr   = SNESGetDM(snes,&dm);CHKERRQ(ierr);
  ierr   = TSGLEEGetVecs(ts,dm,&Ydot);CHKERRQ(ierr);
  /* Set Ydot = shift*X */
  ierr   = VecCopy(X,Ydot);CHKERRQ(ierr);
  ierr   = VecScale(Ydot,shift);CHKERRQ(ierr);
  dmsave = ts->dm;
  ts->dm = dm;

  ierr = TSComputeIFunction(ts,glee->stage_time,X,Ydot,F,PETSC_FALSE);CHKERRQ(ierr);

  ts->dm = dmsave;
  ierr   = TSGLEERestoreVecs(ts,dm,&Ydot);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode SNESTSFormJacobian_GLEE(SNES snes,Vec X,Mat A,Mat B,TS ts)
{
  TS_GLEE        *glee = (TS_GLEE*)ts->data;
  DM             dm,dmsave;
  Vec            Ydot;
  PetscReal      shift = glee->scoeff / ts->time_step;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = SNESGetDM(snes,&dm);CHKERRQ(ierr);
  ierr = TSGLEEGetVecs(ts,dm,&Ydot);CHKERRQ(ierr);
  /* glee->Ydot has already been computed in SNESTSFormFunction_GLEE (SNES guarantees this) */
  dmsave = ts->dm;
  ts->dm = dm;

  ierr = TSComputeIJacobian(ts,glee->stage_time,X,Ydot,shift,A,B,PETSC_FALSE);CHKERRQ(ierr);

  ts->dm = dmsave;
  ierr   = TSGLEERestoreVecs(ts,dm,&Ydot);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode DMCoarsenHook_TSGLEE(DM fine,DM coarse,void *ctx)
{
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

static PetscErrorCode DMRestrictHook_TSGLEE(DM fine,Mat restrct,Vec rscale,Mat inject,DM coarse,void *ctx)
{
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}


static PetscErrorCode DMSubDomainHook_TSGLEE(DM dm,DM subdm,void *ctx)
{
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

static PetscErrorCode DMSubDomainRestrictHook_TSGLEE(DM dm,VecScatter gscat,VecScatter lscat,DM subdm,void *ctx)
{
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

static PetscErrorCode TSSetUp_GLEE(TS ts)
{
  TS_GLEE        *glee = (TS_GLEE*)ts->data;
  GLEETableau    tab;
  PetscInt       s,r;
  PetscErrorCode ierr;
  DM             dm;

  PetscFunctionBegin;
  if (!glee->tableau) {
    ierr = TSGLEESetType(ts,TSGLEEDefaultType);CHKERRQ(ierr);
  }
  tab  = glee->tableau;
  s    = tab->s;
  r    = tab->r;
  ierr = VecDuplicateVecs(ts->vec_sol,r,&glee->Y);CHKERRQ(ierr);
  ierr = VecDuplicateVecs(ts->vec_sol,r,&glee->X);CHKERRQ(ierr);
  ierr = VecDuplicateVecs(ts->vec_sol,s,&glee->YStage);CHKERRQ(ierr);
  ierr = VecDuplicateVecs(ts->vec_sol,s,&glee->YdotStage);CHKERRQ(ierr);
  ierr = VecDuplicate(ts->vec_sol,&glee->Ydot);CHKERRQ(ierr);
  ierr = VecDuplicate(ts->vec_sol,&glee->yGErr);CHKERRQ(ierr);
  ierr = VecZeroEntries(glee->yGErr);CHKERRQ(ierr);
  ierr = VecDuplicate(ts->vec_sol,&glee->W);CHKERRQ(ierr);
  ierr = PetscMalloc2(s,&glee->swork,r,&glee->rwork);CHKERRQ(ierr);
  ierr = TSGetDM(ts,&dm);CHKERRQ(ierr);
  ierr = DMCoarsenHookAdd(dm,DMCoarsenHook_TSGLEE,DMRestrictHook_TSGLEE,ts);CHKERRQ(ierr);
  ierr = DMSubDomainHookAdd(dm,DMSubDomainHook_TSGLEE,DMSubDomainRestrictHook_TSGLEE,ts);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode TSStartingMethod_GLEE(TS ts)
{
  TS_GLEE        *glee = (TS_GLEE*)ts->data;
  GLEETableau    tab  = glee->tableau;
  PetscInt       r=tab->r,i;
  PetscReal      *S=tab->S;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  for (i=0; i<r; i++) { 
    ierr = VecZeroEntries(glee->Y[i]);CHKERRQ(ierr); 
    ierr = VecAXPY(glee->Y[i],S[i],ts->vec_sol);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}


/*------------------------------------------------------------*/

static PetscErrorCode TSSetFromOptions_GLEE(PetscOptionItems *PetscOptionsObject,TS ts)
{
  PetscErrorCode ierr;
  char           gleetype[256];

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"GLEE ODE solver options");CHKERRQ(ierr);
  {
    GLEETableauLink link;
    PetscInt        count,choice;
    PetscBool       flg;
    const char      **namelist;

    ierr = PetscStrncpy(gleetype,TSGLEEDefaultType,sizeof(gleetype));CHKERRQ(ierr);
    for (link=GLEETableauList,count=0; link; link=link->next,count++) ;
    ierr = PetscMalloc1(count,(char***)&namelist);CHKERRQ(ierr);
    for (link=GLEETableauList,count=0; link; link=link->next,count++) namelist[count] = link->tab.name;
    ierr = PetscOptionsEList("-ts_glee_type","Family of GLEE method","TSGLEESetType",(const char*const*)namelist,count,gleetype,&choice,&flg);CHKERRQ(ierr);
    ierr = TSGLEESetType(ts,flg ? namelist[choice] : gleetype);CHKERRQ(ierr);
    ierr = PetscFree(namelist);CHKERRQ(ierr);
  }
  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode TSView_GLEE(TS ts,PetscViewer viewer)
{
  TS_GLEE        *glee   = (TS_GLEE*)ts->data;
  GLEETableau    tab  = glee->tableau;
  PetscBool      iascii;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&iascii);CHKERRQ(ierr);
  if (iascii) {
    TSGLEEType gleetype;
    char       buf[512];
    ierr = TSGLEEGetType(ts,&gleetype);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  GLEE type %s\n",gleetype);CHKERRQ(ierr);
    ierr = PetscFormatRealArray(buf,sizeof(buf),"% 8.6f",tab->s,tab->c);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  Abscissa     c = %s\n",buf);CHKERRQ(ierr);
    /* Note: print out r as well */
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode TSLoad_GLEE(TS ts,PetscViewer viewer)
{
  PetscErrorCode ierr;
  SNES           snes;
  TSAdapt        tsadapt;

  PetscFunctionBegin;
  ierr = TSGetAdapt(ts,&tsadapt);CHKERRQ(ierr);
  ierr = TSAdaptLoad(tsadapt,viewer);CHKERRQ(ierr);
  ierr = TSGetSNES(ts,&snes);CHKERRQ(ierr);
  ierr = SNESLoad(snes,viewer);CHKERRQ(ierr);
  /* function and Jacobian context for SNES when used with TS is always ts object */
  ierr = SNESSetFunction(snes,NULL,NULL,ts);CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,NULL,NULL,NULL,ts);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
  TSGLEESetType - Set the type of GLEE scheme

  Logically collective

  Input Parameter:
+  ts - timestepping context
-  gleetype - type of GLEE-scheme

  Level: intermediate

.seealso: TSGLEEGetType(), TSGLEE
@*/
PetscErrorCode TSGLEESetType(TS ts,TSGLEEType gleetype)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ts,TS_CLASSID,1);
  PetscValidCharPointer(gleetype,2);
  ierr = PetscTryMethod(ts,"TSGLEESetType_C",(TS,TSGLEEType),(ts,gleetype));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
  TSGLEEGetType - Get the type of GLEE scheme

  Logically collective

  Input Parameter:
.  ts - timestepping context

  Output Parameter:
.  gleetype - type of GLEE-scheme

  Level: intermediate

.seealso: TSGLEESetType()
@*/
PetscErrorCode TSGLEEGetType(TS ts,TSGLEEType *gleetype)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ts,TS_CLASSID,1);
  ierr = PetscUseMethod(ts,"TSGLEEGetType_C",(TS,TSGLEEType*),(ts,gleetype));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode  TSGLEEGetType_GLEE(TS ts,TSGLEEType *gleetype)
{
  TS_GLEE     *glee = (TS_GLEE*)ts->data;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!glee->tableau) {
    ierr = TSGLEESetType(ts,TSGLEEDefaultType);CHKERRQ(ierr);
  }
  *gleetype = glee->tableau->name;
  PetscFunctionReturn(0);
}
PetscErrorCode  TSGLEESetType_GLEE(TS ts,TSGLEEType gleetype)
{
  TS_GLEE         *glee = (TS_GLEE*)ts->data;
  PetscErrorCode  ierr;
  PetscBool       match;
  GLEETableauLink link;

  PetscFunctionBegin;
  if (glee->tableau) {
    ierr = PetscStrcmp(glee->tableau->name,gleetype,&match);CHKERRQ(ierr);
    if (match) PetscFunctionReturn(0);
  }
  for (link = GLEETableauList; link; link=link->next) {
    ierr = PetscStrcmp(link->tab.name,gleetype,&match);CHKERRQ(ierr);
    if (match) {
      ierr = TSReset_GLEE(ts);CHKERRQ(ierr);
      glee->tableau = &link->tab;
      PetscFunctionReturn(0);
    }
  }
  SETERRQ1(PetscObjectComm((PetscObject)ts),PETSC_ERR_ARG_UNKNOWN_TYPE,"Could not find '%s'",gleetype);
  PetscFunctionReturn(0);
}

static PetscErrorCode  TSGetStages_GLEE(TS ts,PetscInt *ns,Vec **Y)
{
  TS_GLEE *glee = (TS_GLEE*)ts->data;

  PetscFunctionBegin;
  if (ns) *ns = glee->tableau->s;
  if (Y)  *Y  = glee->YStage;
  PetscFunctionReturn(0);
}

PetscErrorCode TSGetSolutionComponents_GLEE(TS ts,PetscInt *n,Vec *Y)
{
  TS_GLEE         *glee = (TS_GLEE*)ts->data;
  GLEETableau     tab   = glee->tableau;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  if (!Y) *n = tab->r;
  else {
    if ((*n >= 0) && (*n < tab->r)) { 
      ierr = VecCopy(glee->Y[*n],*Y);CHKERRQ(ierr);
    } else SETERRQ3(PetscObjectComm((PetscObject)ts),PETSC_ERR_ARG_OUTOFRANGE,"Second argument (%d) out of range[%d,%d].",*n,0,tab->r-1);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode TSGetAuxSolution_GLEE(TS ts,Vec *X)
{
  TS_GLEE         *glee = (TS_GLEE*)ts->data;
  GLEETableau     tab   = glee->tableau;
  PetscReal       *F    = tab->Fembed;
  PetscInt        r     = tab->r;
  Vec             *Y    = glee->Y;
  PetscScalar     *wr   = glee->rwork;
  PetscInt        i;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  ierr = VecZeroEntries(*X);CHKERRQ(ierr);
  for (i=0; i<r; i++) wr[i] = F[i];
  ierr = VecMAXPY((*X),r,wr,Y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode TSGetTimeError_GLEE(TS ts,PetscInt n,Vec *X)
{
  TS_GLEE         *glee = (TS_GLEE*)ts->data;
  GLEETableau     tab   = glee->tableau;
  PetscReal       *F    = tab->Ferror;
  PetscInt        r     = tab->r;
  Vec             *Y    = glee->Y;
  PetscScalar     *wr   = glee->rwork;
  PetscInt        i;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  ierr = VecZeroEntries(*X);CHKERRQ(ierr);
  if(n==0){
    for (i=0; i<r; i++) wr[i] = F[i];
    ierr = VecMAXPY((*X),r,wr,Y);CHKERRQ(ierr);
  } else if(n==-1) {
    *X=glee->yGErr;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode TSSetTimeError_GLEE(TS ts,Vec X)
{
  TS_GLEE         *glee = (TS_GLEE*)ts->data;
  GLEETableau     tab   = glee->tableau;
  PetscReal       *S    = tab->Serror;
  PetscInt        r     = tab->r,i;
  Vec             *Y    = glee->Y;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  if (r != 2) SETERRQ2(PetscObjectComm((PetscObject)ts),PETSC_ERR_SUP,"TSSetTimeError_GLEE not supported for '%s' with r=%D.",tab->name,tab->r);
  for (i=1; i<r; i++) {
    ierr = VecCopy(ts->vec_sol,Y[i]);CHKERRQ(ierr);
    ierr = VecAXPBY(Y[i],S[0],S[1],X);CHKERRQ(ierr);
    ierr = VecCopy(X,glee->yGErr);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode TSDestroy_GLEE(TS ts)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = TSReset_GLEE(ts);CHKERRQ(ierr);
  if (ts->dm) {
    ierr = DMCoarsenHookRemove(ts->dm,DMCoarsenHook_TSGLEE,DMRestrictHook_TSGLEE,ts);CHKERRQ(ierr);
    ierr = DMSubDomainHookRemove(ts->dm,DMSubDomainHook_TSGLEE,DMSubDomainRestrictHook_TSGLEE,ts);CHKERRQ(ierr);
  }
  ierr = PetscFree(ts->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)ts,"TSGLEEGetType_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)ts,"TSGLEESetType_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------ */
/*MC
      TSGLEE - ODE and DAE solver using General Linear with Error Estimation schemes

  The user should provide the right hand side of the equation 
  using TSSetRHSFunction().

  Notes:
  The default is TSGLEE35, it can be changed with TSGLEESetType() or -ts_glee_type

  Level: beginner

.seealso:  TSCreate(), TS, TSSetType(), TSGLEESetType(), TSGLEEGetType(),
           TSGLEE23, TTSGLEE24, TSGLEE35, TSGLEE25I, TSGLEEEXRK2A, 
           TSGLEERK32G1, TSGLEERK285EX, TSGLEEType, TSGLEERegister()

M*/
PETSC_EXTERN PetscErrorCode TSCreate_GLEE(TS ts)
{
  TS_GLEE         *th;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
#if !defined(PETSC_USE_DYNAMIC_LIBRARIES)
  ierr = TSGLEEInitializePackage();CHKERRQ(ierr);
#endif

  ts->ops->reset                  = TSReset_GLEE;
  ts->ops->destroy                = TSDestroy_GLEE;
  ts->ops->view                   = TSView_GLEE;
  ts->ops->load                   = TSLoad_GLEE;
  ts->ops->setup                  = TSSetUp_GLEE;
  ts->ops->step                   = TSStep_GLEE;
  ts->ops->interpolate            = TSInterpolate_GLEE;
  ts->ops->evaluatestep           = TSEvaluateStep_GLEE;
  ts->ops->setfromoptions         = TSSetFromOptions_GLEE;
  ts->ops->getstages              = TSGetStages_GLEE;
  ts->ops->snesfunction           = SNESTSFormFunction_GLEE;
  ts->ops->snesjacobian           = SNESTSFormJacobian_GLEE;
  ts->ops->getsolutioncomponents  = TSGetSolutionComponents_GLEE;
  ts->ops->getauxsolution         = TSGetAuxSolution_GLEE;
  ts->ops->gettimeerror           = TSGetTimeError_GLEE;
  ts->ops->settimeerror           = TSSetTimeError_GLEE;
  ts->ops->startingmethod         = TSStartingMethod_GLEE;
  ts->default_adapt_type          = TSADAPTGLEE;

  ts->usessnes = PETSC_TRUE;

  ierr = PetscNewLog(ts,&th);CHKERRQ(ierr);
  ts->data = (void*)th;

  ierr = PetscObjectComposeFunction((PetscObject)ts,"TSGLEEGetType_C",TSGLEEGetType_GLEE);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)ts,"TSGLEESetType_C",TSGLEESetType_GLEE);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
