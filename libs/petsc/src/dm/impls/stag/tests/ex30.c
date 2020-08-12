static char help[] = "Test DMStagMatGetValuesStencil() in 3D\n\n";

#include <petscdm.h>
#include <petscksp.h>
#include <petscdmstag.h>

/* Shorter, more convenient names for DMStagStencilLocation entries */
#define BACK_DOWN_LEFT   DMSTAG_BACK_DOWN_LEFT
#define BACK_DOWN        DMSTAG_BACK_DOWN
#define BACK_DOWN_RIGHT  DMSTAG_BACK_DOWN_RIGHT
#define BACK_LEFT        DMSTAG_BACK_LEFT
#define BACK             DMSTAG_BACK
#define BACK_RIGHT       DMSTAG_BACK_RIGHT
#define BACK_UP_LEFT     DMSTAG_BACK_UP_LEFT
#define BACK_UP          DMSTAG_BACK_UP
#define BACK_UP_RIGHT    DMSTAG_BACK_UP_RIGHT
#define DOWN_LEFT        DMSTAG_DOWN_LEFT
#define DOWN             DMSTAG_DOWN
#define DOWN_RIGHT       DMSTAG_DOWN_RIGHT
#define LEFT             DMSTAG_LEFT
#define ELEMENT          DMSTAG_ELEMENT
#define RIGHT            DMSTAG_RIGHT
#define UP_LEFT          DMSTAG_UP_LEFT
#define UP               DMSTAG_UP
#define UP_RIGHT         DMSTAG_UP_RIGHT
#define FRONT_DOWN_LEFT  DMSTAG_FRONT_DOWN_LEFT
#define FRONT_DOWN       DMSTAG_FRONT_DOWN
#define FRONT_DOWN_RIGHT DMSTAG_FRONT_DOWN_RIGHT
#define FRONT_LEFT       DMSTAG_FRONT_LEFT
#define FRONT            DMSTAG_FRONT
#define FRONT_RIGHT      DMSTAG_FRONT_RIGHT
#define FRONT_UP_LEFT    DMSTAG_FRONT_UP_LEFT
#define FRONT_UP         DMSTAG_FRONT_UP
#define FRONT_UP_RIGHT   DMSTAG_FRONT_UP_RIGHT

static PetscErrorCode CreateMat(DM,Mat*);
static PetscErrorCode CheckMat(DM,Mat);

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  DM             dmSol;
  Mat            A;

  ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  {
    const PetscInt dof0 = 0, dof1 = 0,dof2 = 1, dof3 = 1; /* 1 dof on each face and element center */
    const PetscInt stencilWidth = 1;
    ierr = DMStagCreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,4,5,6,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,dof0,dof1,dof2,dof3,DMSTAG_STENCIL_BOX,stencilWidth,NULL,NULL,NULL,&dmSol);CHKERRQ(ierr);
    ierr = DMSetFromOptions(dmSol);CHKERRQ(ierr);
    ierr = DMSetUp(dmSol);CHKERRQ(ierr);
    ierr = DMStagSetUniformCoordinatesExplicit(dmSol,0.0,1.0,0.0,1.0,0.0,1.0);CHKERRQ(ierr);
  }
  ierr = CreateMat(dmSol,&A);CHKERRQ(ierr);
  ierr = CheckMat(dmSol,A);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = DMDestroy(&dmSol);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}

static PetscErrorCode CreateMat(DM dmSol,Mat *pA)
{
  PetscErrorCode    ierr;
  Vec               coordLocal;
  Mat               A;
  PetscInt          startx,starty,startz,N[3],nx,ny,nz,ex,ey,ez,d;
  PetscInt          icp[3],icux[3],icuy[3],icuz[3],icux_right[3],icuy_up[3],icuz_front[3];
  PetscBool         isLastRankx,isLastRanky,isLastRankz,isFirstRankx,isFirstRanky,isFirstRankz;
  PetscReal         hx,hy,hz;
  DM                dmCoord;
  PetscScalar       ****arrCoord;

  PetscFunctionBeginUser;
  ierr = DMCreateMatrix(dmSol,pA);CHKERRQ(ierr);
  A = *pA;
  ierr = DMStagGetCorners(dmSol,&startx,&starty,&startz,&nx,&ny,&nz,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DMStagGetGlobalSizes(dmSol,&N[0],&N[1],&N[2]);CHKERRQ(ierr);
  if (N[0] < 2 || N[1] < 2 || N[2] < 2) SETERRQ(PetscObjectComm((PetscObject)dmSol),PETSC_ERR_ARG_SIZ,"This example requires at least two elements in each dimensions");
  ierr = DMStagGetIsLastRank(dmSol,&isLastRankx,&isLastRanky,&isLastRankz);CHKERRQ(ierr);
  ierr = DMStagGetIsFirstRank(dmSol,&isFirstRankx,&isFirstRanky,&isFirstRankz);CHKERRQ(ierr);
  hx = 1.0/N[0]; hy = 1.0/N[1]; hz = 1.0/N[2];
  ierr = DMGetCoordinateDM(dmSol,&dmCoord);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(dmSol,&coordLocal);CHKERRQ(ierr);
  ierr = DMStagVecGetArrayRead(dmCoord,coordLocal,&arrCoord);CHKERRQ(ierr);
  for (d=0; d<3; ++d) {
    ierr = DMStagGetLocationSlot(dmCoord,ELEMENT,d,&icp[d]       );CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dmCoord,LEFT,   d,&icux[d]      );CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dmCoord,DOWN,   d,&icuy[d]      );CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dmCoord,BACK,   d,&icuz[d]      );CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dmCoord,RIGHT,  d,&icux_right[d]);CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dmCoord,UP,     d,&icuy_up[d]   );CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dmCoord,FRONT,  d,&icuz_front[d]);CHKERRQ(ierr);
  }

  for (ez = startz; ez<startz+nz; ++ez) {
    for (ey = starty; ey<starty+ny; ++ey) {
      for (ex = startx; ex<startx+nx; ++ex) {
        if (ex == N[0]-1) {
          /* Right Boundary velocity Dirichlet */
          DMStagStencil row;
          const PetscScalar valA = 1.0;
          row.i = ex; row.j = ey; row.k = ez; row.loc = RIGHT; row.c = 0;
          ierr = DMStagMatSetValuesStencil(dmSol,A,1,&row,1,&row,&valA,INSERT_VALUES);CHKERRQ(ierr);
        }
        if (ey == N[1]-1) {
          /* Top boundary velocity Dirichlet */
          DMStagStencil row;
          const PetscScalar valA = 1.0;
          row.i = ex; row.j = ey; row.k = ez; row.loc = UP; row.c = 0;
          ierr = DMStagMatSetValuesStencil(dmSol,A,1,&row,1,&row,&valA,INSERT_VALUES);CHKERRQ(ierr);
        }
        if (ez == N[2]-1) {
          /* Top boundary velocity Dirichlet */
          DMStagStencil row;
          const PetscScalar valA = 1.0;
          row.i = ex; row.j = ey; row.k = ez; row.loc = FRONT; row.c = 0;
          ierr = DMStagMatSetValuesStencil(dmSol,A,1,&row,1,&row,&valA,INSERT_VALUES);CHKERRQ(ierr);
        }

        /* Equation on left face of this element */
        if (ex == 0) {
          /* Left velocity Dirichlet */
          DMStagStencil row;
          const PetscScalar valA = 1.0;
          row.i = ex; row.j = ey; row.k = ez; row.loc = LEFT; row.c = 0;
          ierr = DMStagMatSetValuesStencil(dmSol,A,1,&row,1,&row,&valA,INSERT_VALUES);CHKERRQ(ierr);
        } else {
          /* X-momentum interior equation : (u_xx + u_yy + u_zz) - p_x = f^x */
          DMStagStencil row,col[9];
          PetscScalar   valA[9];
          PetscInt      nEntries;

          row.i = ex; row.j = ey; row.k = ez; row.loc = LEFT; row.c = 0;
          if (ey == 0) {
            if (ez == 0) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = LEFT;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -1.0 / (hy*hy) -1.0 / (hz*hz);
              /* Missing down term */
              col[1].i = ex  ; col[1].j = ey+1;  col[1].k = ez  ; col[1].loc = LEFT;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex-1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = LEFT;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
              col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = LEFT;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              /* Missing back term */
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez+1; col[4].loc = LEFT;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex-1; col[5].j = ey  ;  col[5].k = ez  ; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hx;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hx;
            } else if (ez == N[2]-1) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = LEFT;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -1.0 / (hy*hy) -1.0 / (hz*hz);
              /* Missing down term */
              col[1].i = ex  ; col[1].j = ey+1;  col[1].k = ez  ; col[1].loc = LEFT;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex-1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = LEFT;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
              col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = LEFT;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = LEFT;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              /* Missing front term */
              col[5].i = ex-1; col[5].j = ey  ;  col[5].k = ez  ; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hx;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hx;
            } else {
              nEntries = 8;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = LEFT;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -1.0 / (hy*hy) -2.0 / (hz*hz);
              /* Missing down term */
              col[1].i = ex  ; col[1].j = ey+1;  col[1].k = ez  ; col[1].loc = LEFT;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex-1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = LEFT;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
              col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = LEFT;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = LEFT;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez+1; col[5].loc = LEFT;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
              col[6].i = ex-1; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hx;
              col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hx;
            }
          } else if (ey == N[1]-1) {
            if (ez == 0) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = LEFT;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -2.0 / (hy*hy) -1.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = LEFT;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              /* Missing up term */
              col[2].i = ex-1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = LEFT;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
              col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = LEFT;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              /* Missing back entry */
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez+1; col[4].loc = LEFT;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex-1; col[5].j = ey  ;  col[5].k = ez  ; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hx;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hx;
            } else if (ez == N[2]-1) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = LEFT;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -2.0 / (hy*hy) -1.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = LEFT;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              /* Missing up term */
              col[2].i = ex-1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = LEFT;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
              col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = LEFT;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = LEFT;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              /* Missing front term */
              col[5].i = ex-1; col[5].j = ey  ;  col[5].k = ez  ; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hx;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hx;
            } else {
              nEntries = 8;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = LEFT;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -2.0 / (hy*hy) -2.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = LEFT;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              /* Missing up term */
              col[2].i = ex-1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = LEFT;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
              col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = LEFT;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = LEFT;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez+1; col[5].loc = LEFT;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
              col[6].i = ex-1; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hx;
              col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hx;
            }
          } else if (ez == 0) {
            nEntries = 8;
            col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = LEFT;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -2.0 / (hy*hy) -1.0 / (hz*hz);
            col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = LEFT;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
            col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = LEFT;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
            col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = LEFT;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
            col[4].i = ex+1; col[4].j = ey  ;  col[4].k = ez  ; col[4].loc = LEFT;    col[4].c  = 0; valA[4] =  1.0 / (hx*hx);
            /* Missing back term */
            col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez+1; col[5].loc = LEFT;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
            col[6].i = ex-1; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hx;
            col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hx;
          } else if (ez == N[2]-1) {
            nEntries = 8;
            col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = LEFT;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -2.0 / (hy*hy) -1.0 / (hz*hz);
            col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = LEFT;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
            col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = LEFT;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
            col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = LEFT;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
            col[4].i = ex+1; col[4].j = ey  ;  col[4].k = ez  ; col[4].loc = LEFT;    col[4].c  = 0; valA[4] =  1.0 / (hx*hx);
            col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez-1; col[5].loc = LEFT;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
            /* Missing front term */
            col[6].i = ex-1; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hx;
            col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hx;
          } else {
            nEntries = 9;
            col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = LEFT;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -2.0 / (hy*hy) -2.0 / (hz*hz);
            col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = LEFT;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
            col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = LEFT;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
            col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = LEFT;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
            col[4].i = ex+1; col[4].j = ey  ;  col[4].k = ez  ; col[4].loc = LEFT;    col[4].c  = 0; valA[4] =  1.0 / (hx*hx);
            col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez-1; col[5].loc = LEFT;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
            col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez+1; col[6].loc = LEFT;    col[6].c  = 0; valA[6] =  1.0 / (hz*hz);
            col[7].i = ex-1; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] =  1.0 / hx;
            col[8].i = ex  ; col[8].j = ey  ;  col[8].k = ez  ; col[8].loc = ELEMENT; col[8].c  = 0; valA[8] = -1.0 / hx;
          }
          ierr = DMStagMatSetValuesStencil(dmSol,A,1,&row,nEntries,col,valA,INSERT_VALUES);CHKERRQ(ierr);
        }

        /* Equation on bottom face of this element */
        if (ey == 0) {
          /* Bottom boundary velocity Dirichlet */
          DMStagStencil row;
          const PetscScalar valA = 1.0;
          row.i = ex; row.j = ey; row.k = ez; row.loc = DOWN; row.c = 0;
          ierr = DMStagMatSetValuesStencil(dmSol,A,1,&row,1,&row,&valA,INSERT_VALUES);CHKERRQ(ierr);
        } else {
          /* Y-momentum equation, (v_xx + v_yy + v_zz) - p_y = f^y */
          DMStagStencil row,col[9];
          PetscScalar   valA[9];
          PetscInt      nEntries;

          row.i = ex; row.j = ey; row.k = ez; row.loc = DOWN; row.c = 0;
          if (ex ==0) {
            if (ez == 0) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = DOWN;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -2.0 / (hy*hy) -1.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = DOWN;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = DOWN;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
              /* Left term missing */
              col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = DOWN;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              /* Back term missing */
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez+1; col[4].loc = DOWN;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey-1;  col[5].k = ez  ; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hy;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hy;
            } else if (ez == N[2]-1) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = DOWN;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -2.0 / (hy*hy) -1.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = DOWN;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = DOWN;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
              /* Left term missing */
              col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = DOWN;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = DOWN;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              /* Front term missing */
              col[5].i = ex  ; col[5].j = ey-1;  col[5].k = ez  ; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hy;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hy;
            } else {
              nEntries = 8;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = DOWN;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -2.0 / (hy*hy) -2.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = DOWN;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = DOWN;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
              /* Left term missing */
              col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = DOWN;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = DOWN;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez+1; col[5].loc = DOWN;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
              col[6].i = ex  ; col[6].j = ey-1;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hy;
              col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hy;
            }
          } else if (ex == N[0]-1) {
            if (ez == 0) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = DOWN;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -2.0 / (hy*hy) -1.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = DOWN;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = DOWN;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
              col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = DOWN;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              /* Right term missing */
              /* Back term missing */
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez+1; col[4].loc = DOWN;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey-1;  col[5].k = ez  ; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hy;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hy;
            } else if (ez == N[2]-1) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = DOWN;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -2.0 / (hy*hy) -1.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = DOWN;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = DOWN;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
              col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = DOWN;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              /* Right term missing */
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = DOWN;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              /* Front term missing */
              col[5].i = ex  ; col[5].j = ey-1;  col[5].k = ez  ; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hy;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hy;
            } else {
              nEntries = 8;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = DOWN;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -2.0 / (hy*hy) -2.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = DOWN;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = DOWN;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
              col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = DOWN;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              /* Right term missing */
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = DOWN;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez+1; col[5].loc = DOWN;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
              col[6].i = ex  ; col[6].j = ey-1;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hy;
              col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hy;
            }
          } else if (ez == 0) {
            nEntries = 8;
            col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = DOWN;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -2.0 / (hy*hy) -1.0 / (hz*hz);
            col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = DOWN;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
            col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = DOWN;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
            col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = DOWN;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
            col[4].i = ex+1; col[4].j = ey  ;  col[4].k = ez  ; col[4].loc = DOWN;    col[4].c  = 0; valA[4] =  1.0 / (hx*hx);
            /* Back term missing */
            col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez+1; col[5].loc = DOWN;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
            col[6].i = ex  ; col[6].j = ey-1;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hy;
            col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hy;
          } else if (ez == N[2]-1) {
            nEntries = 8;
            col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = DOWN;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -2.0 / (hy*hy) -1.0 / (hz*hz);
            col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = DOWN;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
            col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = DOWN;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
            col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = DOWN;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
            col[4].i = ex+1; col[4].j = ey  ;  col[4].k = ez  ; col[4].loc = DOWN;    col[4].c  = 0; valA[4] =  1.0 / (hx*hx);
            col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez-1; col[5].loc = DOWN;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
            /* Front term missing */
            col[6].i = ex  ; col[6].j = ey-1;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hy;
            col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hy;
          } else {
            nEntries = 9;
            col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = DOWN;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -2.0 / (hy*hy) -2.0 / (hz*hz);
            col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = DOWN;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
            col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = DOWN;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
            col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = DOWN;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
            col[4].i = ex+1; col[4].j = ey  ;  col[4].k = ez  ; col[4].loc = DOWN;    col[4].c  = 0; valA[4] =  1.0 / (hx*hx);
            col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez-1; col[5].loc = DOWN;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
            col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez+1; col[6].loc = DOWN;    col[6].c  = 0; valA[6] =  1.0 / (hz*hz);
            col[7].i = ex  ; col[7].j = ey-1;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] =  1.0 / hy;
            col[8].i = ex  ; col[8].j = ey  ;  col[8].k = ez  ; col[8].loc = ELEMENT; col[8].c  = 0; valA[8] = -1.0 / hy;
          }
          ierr = DMStagMatSetValuesStencil(dmSol,A,1,&row,nEntries,col,valA,INSERT_VALUES);CHKERRQ(ierr);
        }

        /* Equation on back face of this element */
        if (ez == 0) {
          /* Back boundary velocity Dirichlet */
          DMStagStencil row;
          const PetscScalar valA = 1.0;
          row.i = ex; row.j = ey; row.k = ez; row.loc = BACK; row.c = 0;
          ierr = DMStagMatSetValuesStencil(dmSol,A,1,&row,1,&row,&valA,INSERT_VALUES);CHKERRQ(ierr);
        } else {
          /* Z-momentum equation, (w_xx + w_yy + w_zz) - p_z = f^z */
          DMStagStencil row,col[9];
          PetscScalar   valA[9];
          PetscInt      nEntries;

          row.i = ex; row.j = ey; row.k = ez; row.loc = BACK; row.c = 0;
          if (ex == 0) {
            if (ey == 0) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = BACK;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -1.0 / (hy*hy) -2.0 / (hz*hz);
              /* Down term missing */
              col[1].i = ex  ; col[1].j = ey+1;  col[1].k = ez  ; col[1].loc = BACK;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              /* Left term missing */
              col[2].i = ex+1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = BACK;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
              col[3].i = ex  ; col[3].j = ey  ;  col[3].k = ez-1; col[3].loc = BACK;    col[3].c  = 0; valA[3] =  1.0 / (hz*hz);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez+1; col[4].loc = BACK;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez-1; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hz;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hz;
            } else if (ey == N[1]-1) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = BACK;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -1.0 / (hy*hy) -2.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = BACK;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              /* Up term missing */
              /* Left term missing */
              col[2].i = ex+1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = BACK;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
              col[3].i = ex  ; col[3].j = ey  ;  col[3].k = ez-1; col[3].loc = BACK;    col[3].c  = 0; valA[3] =  1.0 / (hz*hz);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez+1; col[4].loc = BACK;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez-1; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hz;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hz;
            } else {
              nEntries = 8;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = BACK;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -2.0 / (hy*hy) -2.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = BACK;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = BACK;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
              /* Left term missing */
              col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = BACK;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = BACK;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez+1; col[5].loc = BACK;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez-1; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hz;
              col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hz;
            }
          } else if (ex == N[0]-1) {
            if (ey == 0) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = BACK;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -1.0 / (hy*hy) -2.0 / (hz*hz);
              /* Down term missing */
              col[1].i = ex  ; col[1].j = ey+1;  col[1].k = ez  ; col[1].loc = BACK;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex-1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = BACK;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
              /* Right term missing */
              col[3].i = ex  ; col[3].j = ey  ;  col[3].k = ez-1; col[3].loc = BACK;    col[3].c  = 0; valA[3] =  1.0 / (hz*hz);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez+1; col[4].loc = BACK;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez-1; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hz;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hz;
            } else if (ey == N[1]-1) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = BACK;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -1.0 / (hy*hy) -2.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = BACK;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              /* Up term missing */
              col[2].i = ex-1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = BACK;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
              /* Right term missing */
              col[3].i = ex  ; col[3].j = ey  ;  col[3].k = ez-1; col[3].loc = BACK;    col[3].c  = 0; valA[3] =  1.0 / (hz*hz);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez+1; col[4].loc = BACK;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez-1; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hz;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hz;
            } else {
              nEntries = 8;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = BACK;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -2.0 / (hy*hy) -2.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = BACK;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = BACK;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
              col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = BACK;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              /* Right term missing */
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = BACK;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez+1; col[5].loc = BACK;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez-1; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hz;
              col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hz;
            }
          } else if (ey == 0) {
            nEntries = 8;
            col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = BACK;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -1.0 / (hy*hy) -2.0 / (hz*hz);
            /* Down term missing */
            col[1].i = ex  ; col[1].j = ey+1;  col[1].k = ez  ; col[1].loc = BACK;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
            col[2].i = ex-1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = BACK;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
            col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = BACK;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
            col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = BACK;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
            col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez+1; col[5].loc = BACK;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
            col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez-1; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hz;
            col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hz;
          } else if (ey == N[1]-1) {
            nEntries = 8;
            col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = BACK;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -1.0 / (hy*hy) -2.0 / (hz*hz);
            col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = BACK;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
            /* Up term missing */
            col[2].i = ex-1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = BACK;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
            col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = BACK;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
            col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = BACK;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
            col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez+1; col[5].loc = BACK;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
            col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez-1; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hz;
            col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hz;
          } else {
            nEntries = 9;
            col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = BACK;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -2.0 / (hy*hy) -2.0 / (hz*hz);
            col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = BACK;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
            col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = BACK;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
            col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = BACK;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
            col[4].i = ex+1; col[4].j = ey  ;  col[4].k = ez  ; col[4].loc = BACK;    col[4].c  = 0; valA[4] =  1.0 / (hx*hx);
            col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez-1; col[5].loc = BACK;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
            col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez+1; col[6].loc = BACK;    col[6].c  = 0; valA[6] =  1.0 / (hz*hz);
            col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez-1; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] =  1.0 / hz;
            col[8].i = ex  ; col[8].j = ey  ;  col[8].k = ez  ; col[8].loc = ELEMENT; col[8].c  = 0; valA[8] = -1.0 / hz;
          }
          ierr = DMStagMatSetValuesStencil(dmSol,A,1,&row,nEntries,col,valA,INSERT_VALUES);CHKERRQ(ierr);
        }

        /* P equation : u_x + v_y + w_z = g
           Note that this includes an explicit zero on the diagonal. This is only needed for
           direct solvers (not required if using an iterative solver and setting the constant-pressure nullspace) */
        {
          DMStagStencil row,col[7];
          PetscScalar   valA[7];

          row.i    = ex; row.j    = ey; row.k    = ez; row.loc    = ELEMENT; row.c    = 0;
          col[0].i = ex; col[0].j = ey; col[0].k = ez; col[0].loc = LEFT;    col[0].c = 0; valA[0] = -1.0 / hx;
          col[1].i = ex; col[1].j = ey; col[1].k = ez; col[1].loc = RIGHT;   col[1].c = 0; valA[1] =  1.0 / hx;
          col[2].i = ex; col[2].j = ey; col[2].k = ez; col[2].loc = DOWN;    col[2].c = 0; valA[2] = -1.0 / hy;
          col[3].i = ex; col[3].j = ey; col[3].k = ez; col[3].loc = UP;      col[3].c = 0; valA[3] =  1.0 / hy;
          col[4].i = ex; col[4].j = ey; col[4].k = ez; col[4].loc = BACK;    col[4].c = 0; valA[4] = -1.0 / hz;
          col[5].i = ex; col[5].j = ey; col[5].k = ez; col[5].loc = FRONT;   col[5].c = 0; valA[5] =  1.0 / hz;
          col[6]   = row;                                                                  valA[6] =  0.0;
          ierr = DMStagMatSetValuesStencil(dmSol,A,1,&row,7,col,valA,INSERT_VALUES);CHKERRQ(ierr);
        }
      }
    }
  }
  ierr = DMStagVecRestoreArrayRead(dmCoord,coordLocal,&arrCoord);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/* A helper function to check values */
static PetscErrorCode check_vals(PetscInt ex, PetscInt ey, PetscInt ez, PetscInt n,const PetscScalar *ref,const PetscScalar *computed)
{
  PetscInt i;

  PetscFunctionBeginUser;
  for (i=0; i<n; ++i) {
    if (ref[i] != computed[i]) SETERRQ7(PETSC_COMM_SELF,PETSC_ERR_PLIB,"(%D,%D,%D) Assertion Failure. (ref[%D]) %g != %g (computed)[%D]",ex,ey,ez,i,(double)PetscRealPart(ref[i]),(double)PetscRealPart(computed[i]),i);
  }
  PetscFunctionReturn(0);
}

/* The same function as above, but getting and checking values, instead of setting them */
static PetscErrorCode CheckMat(DM dmSol,Mat A)
{
  PetscErrorCode    ierr;
  Vec               coordLocal;
  PetscInt          startx,starty,startz,N[3],nx,ny,nz,ex,ey,ez,d;
  PetscInt          icp[3],icux[3],icuy[3],icuz[3],icux_right[3],icuy_up[3],icuz_front[3];
  PetscBool         isLastRankx,isLastRanky,isLastRankz,isFirstRankx,isFirstRanky,isFirstRankz;
  PetscReal         hx,hy,hz;
  DM                dmCoord;
  PetscScalar       ****arrCoord;
  PetscScalar       computed[1024];

  PetscFunctionBeginUser;
  ierr = DMStagGetCorners(dmSol,&startx,&starty,&startz,&nx,&ny,&nz,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DMStagGetGlobalSizes(dmSol,&N[0],&N[1],&N[2]);CHKERRQ(ierr);
  if (N[0] < 2 || N[1] < 2 || N[2] < 2) SETERRQ(PetscObjectComm((PetscObject)dmSol),PETSC_ERR_ARG_SIZ,"This example requires at least two elements in each dimensions");
  ierr = DMStagGetIsLastRank(dmSol,&isLastRankx,&isLastRanky,&isLastRankz);CHKERRQ(ierr);
  ierr = DMStagGetIsFirstRank(dmSol,&isFirstRankx,&isFirstRanky,&isFirstRankz);CHKERRQ(ierr);
  hx = 1.0/N[0]; hy = 1.0/N[1]; hz = 1.0/N[2];
  ierr = DMGetCoordinateDM(dmSol,&dmCoord);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(dmSol,&coordLocal);CHKERRQ(ierr);
  ierr = DMStagVecGetArrayRead(dmCoord,coordLocal,&arrCoord);CHKERRQ(ierr);
  for (d=0; d<3; ++d) {
    ierr = DMStagGetLocationSlot(dmCoord,ELEMENT,d,&icp[d]       );CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dmCoord,LEFT,   d,&icux[d]      );CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dmCoord,DOWN,   d,&icuy[d]      );CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dmCoord,BACK,   d,&icuz[d]      );CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dmCoord,RIGHT,  d,&icux_right[d]);CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dmCoord,UP,     d,&icuy_up[d]   );CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dmCoord,FRONT,  d,&icuz_front[d]);CHKERRQ(ierr);
  }

  for (ez = startz; ez<startz+nz; ++ez) {
    for (ey = starty; ey<starty+ny; ++ey) {
      for (ex = startx; ex<startx+nx; ++ex) {
        if (ex == N[0]-1) {
          /* Right Boundary velocity Dirichlet */
          DMStagStencil row;
          const PetscScalar valA = 1.0;
          row.i = ex; row.j = ey; row.k = ez; row.loc = RIGHT; row.c = 0;
          ierr = DMStagMatGetValuesStencil(dmSol,A,1,&row,1,&row,computed);CHKERRQ(ierr);
          ierr = check_vals(ex,ey,ez,1,&valA,computed);CHKERRQ(ierr);
        }
        if (ey == N[1]-1) {
          /* Top boundary velocity Dirichlet */
          DMStagStencil row;
          const PetscScalar valA = 1.0;
          row.i = ex; row.j = ey; row.k = ez; row.loc = UP; row.c = 0;
          ierr = DMStagMatGetValuesStencil(dmSol,A,1,&row,1,&row,computed);CHKERRQ(ierr);
          ierr = check_vals(ex,ey,ez,1,&valA,computed);CHKERRQ(ierr);
        }
        if (ez == N[2]-1) {
          /* Top boundary velocity Dirichlet */
          DMStagStencil row;
          const PetscScalar valA = 1.0;
          row.i = ex; row.j = ey; row.k = ez; row.loc = FRONT; row.c = 0;
          ierr = DMStagMatGetValuesStencil(dmSol,A,1,&row,1,&row,computed);CHKERRQ(ierr);
          ierr = check_vals(ex,ey,ez,1,&valA,computed);CHKERRQ(ierr);
        }

        /* Equation on left face of this element */
        if (ex == 0) {
          /* Left velocity Dirichlet */
          DMStagStencil row;
          const PetscScalar valA = 1.0;
          row.i = ex; row.j = ey; row.k = ez; row.loc = LEFT; row.c = 0;
          ierr = DMStagMatGetValuesStencil(dmSol,A,1,&row,1,&row,computed);CHKERRQ(ierr);
          ierr = check_vals(ex,ey,ez,1,&valA,computed);CHKERRQ(ierr);
        } else {
          /* X-momentum interior equation : (u_xx + u_yy + u_zz) - p_x = f^x */
          DMStagStencil row,col[9];
          PetscScalar   valA[9];
          PetscInt      nEntries;

          row.i = ex; row.j = ey; row.k = ez; row.loc = LEFT; row.c = 0;
          if (ey == 0) {
            if (ez == 0) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = LEFT;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -1.0 / (hy*hy) -1.0 / (hz*hz);
              /* Missing down term */
              col[1].i = ex  ; col[1].j = ey+1;  col[1].k = ez  ; col[1].loc = LEFT;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex-1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = LEFT;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
              col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = LEFT;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              /* Missing back term */
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez+1; col[4].loc = LEFT;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex-1; col[5].j = ey  ;  col[5].k = ez  ; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hx;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hx;
            } else if (ez == N[2]-1) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = LEFT;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -1.0 / (hy*hy) -1.0 / (hz*hz);
              /* Missing down term */
              col[1].i = ex  ; col[1].j = ey+1;  col[1].k = ez  ; col[1].loc = LEFT;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex-1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = LEFT;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
              col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = LEFT;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = LEFT;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              /* Missing front term */
              col[5].i = ex-1; col[5].j = ey  ;  col[5].k = ez  ; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hx;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hx;
            } else {
              nEntries = 8;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = LEFT;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -1.0 / (hy*hy) -2.0 / (hz*hz);
              /* Missing down term */
              col[1].i = ex  ; col[1].j = ey+1;  col[1].k = ez  ; col[1].loc = LEFT;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex-1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = LEFT;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
              col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = LEFT;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = LEFT;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez+1; col[5].loc = LEFT;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
              col[6].i = ex-1; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hx;
              col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hx;
            }
          } else if (ey == N[1]-1) {
            if (ez == 0) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = LEFT;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -2.0 / (hy*hy) -1.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = LEFT;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              /* Missing up term */
              col[2].i = ex-1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = LEFT;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
              col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = LEFT;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              /* Missing back entry */
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez+1; col[4].loc = LEFT;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex-1; col[5].j = ey  ;  col[5].k = ez  ; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hx;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hx;
            } else if (ez == N[2]-1) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = LEFT;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -2.0 / (hy*hy) -1.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = LEFT;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              /* Missing up term */
              col[2].i = ex-1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = LEFT;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
              col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = LEFT;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = LEFT;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              /* Missing front term */
              col[5].i = ex-1; col[5].j = ey  ;  col[5].k = ez  ; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hx;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hx;
            } else {
              nEntries = 8;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = LEFT;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -2.0 / (hy*hy) -2.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = LEFT;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              /* Missing up term */
              col[2].i = ex-1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = LEFT;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
              col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = LEFT;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = LEFT;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez+1; col[5].loc = LEFT;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
              col[6].i = ex-1; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hx;
              col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hx;
            }
          } else if (ez == 0) {
            nEntries = 8;
            col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = LEFT;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -2.0 / (hy*hy) -1.0 / (hz*hz);
            col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = LEFT;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
            col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = LEFT;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
            col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = LEFT;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
            col[4].i = ex+1; col[4].j = ey  ;  col[4].k = ez  ; col[4].loc = LEFT;    col[4].c  = 0; valA[4] =  1.0 / (hx*hx);
            /* Missing back term */
            col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez+1; col[5].loc = LEFT;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
            col[6].i = ex-1; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hx;
            col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hx;
          } else if (ez == N[2]-1) {
            nEntries = 8;
            col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = LEFT;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -2.0 / (hy*hy) -1.0 / (hz*hz);
            col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = LEFT;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
            col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = LEFT;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
            col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = LEFT;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
            col[4].i = ex+1; col[4].j = ey  ;  col[4].k = ez  ; col[4].loc = LEFT;    col[4].c  = 0; valA[4] =  1.0 / (hx*hx);
            col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez-1; col[5].loc = LEFT;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
            /* Missing front term */
            col[6].i = ex-1; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hx;
            col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hx;
          } else {
            nEntries = 9;
            col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = LEFT;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -2.0 / (hy*hy) -2.0 / (hz*hz);
            col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = LEFT;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
            col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = LEFT;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
            col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = LEFT;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
            col[4].i = ex+1; col[4].j = ey  ;  col[4].k = ez  ; col[4].loc = LEFT;    col[4].c  = 0; valA[4] =  1.0 / (hx*hx);
            col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez-1; col[5].loc = LEFT;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
            col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez+1; col[6].loc = LEFT;    col[6].c  = 0; valA[6] =  1.0 / (hz*hz);
            col[7].i = ex-1; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] =  1.0 / hx;
            col[8].i = ex  ; col[8].j = ey  ;  col[8].k = ez  ; col[8].loc = ELEMENT; col[8].c  = 0; valA[8] = -1.0 / hx;
          }
          ierr = DMStagMatGetValuesStencil(dmSol,A,1,&row,nEntries,col,computed);CHKERRQ(ierr);
          ierr = check_vals(ex,ey,ez,nEntries,valA,computed);CHKERRQ(ierr);
        }

        /* Equation on bottom face of this element */
        if (ey == 0) {
          /* Bottom boundary velocity Dirichlet */
          DMStagStencil row;
          const PetscScalar valA = 1.0;
          row.i = ex; row.j = ey; row.k = ez; row.loc = DOWN; row.c = 0;
          ierr = DMStagMatGetValuesStencil(dmSol,A,1,&row,1,&row,computed);CHKERRQ(ierr);
          ierr = check_vals(ex,ey,ez,1,&valA,computed);CHKERRQ(ierr);
        } else {
          /* Y-momentum equation, (v_xx + v_yy + v_zz) - p_y = f^y */
          DMStagStencil row,col[9];
          PetscScalar   valA[9];
          PetscInt      nEntries;

          row.i = ex; row.j = ey; row.k = ez; row.loc = DOWN; row.c = 0;
          if (ex ==0) {
            if (ez == 0) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = DOWN;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -2.0 / (hy*hy) -1.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = DOWN;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = DOWN;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
              /* Left term missing */
              col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = DOWN;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              /* Back term missing */
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez+1; col[4].loc = DOWN;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey-1;  col[5].k = ez  ; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hy;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hy;
            } else if (ez == N[2]-1) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = DOWN;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -2.0 / (hy*hy) -1.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = DOWN;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = DOWN;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
              /* Left term missing */
              col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = DOWN;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = DOWN;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              /* Front term missing */
              col[5].i = ex  ; col[5].j = ey-1;  col[5].k = ez  ; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hy;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hy;
            } else {
              nEntries = 8;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = DOWN;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -2.0 / (hy*hy) -2.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = DOWN;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = DOWN;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
              /* Left term missing */
              col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = DOWN;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = DOWN;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez+1; col[5].loc = DOWN;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
              col[6].i = ex  ; col[6].j = ey-1;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hy;
              col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hy;
            }
          } else if (ex == N[0]-1) {
            if (ez == 0) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = DOWN;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -2.0 / (hy*hy) -1.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = DOWN;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = DOWN;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
              col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = DOWN;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              /* Right term missing */
              /* Back term missing */
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez+1; col[4].loc = DOWN;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey-1;  col[5].k = ez  ; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hy;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hy;
            } else if (ez == N[2]-1) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = DOWN;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -2.0 / (hy*hy) -1.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = DOWN;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = DOWN;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
              col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = DOWN;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              /* Right term missing */
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = DOWN;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              /* Front term missing */
              col[5].i = ex  ; col[5].j = ey-1;  col[5].k = ez  ; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hy;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hy;
            } else {
              nEntries = 8;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = DOWN;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -2.0 / (hy*hy) -2.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = DOWN;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = DOWN;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
              col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = DOWN;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              /* Right term missing */
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = DOWN;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez+1; col[5].loc = DOWN;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
              col[6].i = ex  ; col[6].j = ey-1;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hy;
              col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hy;
            }
          } else if (ez == 0) {
            nEntries = 8;
            col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = DOWN;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -2.0 / (hy*hy) -1.0 / (hz*hz);
            col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = DOWN;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
            col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = DOWN;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
            col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = DOWN;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
            col[4].i = ex+1; col[4].j = ey  ;  col[4].k = ez  ; col[4].loc = DOWN;    col[4].c  = 0; valA[4] =  1.0 / (hx*hx);
            /* Back term missing */
            col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez+1; col[5].loc = DOWN;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
            col[6].i = ex  ; col[6].j = ey-1;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hy;
            col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hy;
          } else if (ez == N[2]-1) {
            nEntries = 8;
            col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = DOWN;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -2.0 / (hy*hy) -1.0 / (hz*hz);
            col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = DOWN;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
            col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = DOWN;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
            col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = DOWN;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
            col[4].i = ex+1; col[4].j = ey  ;  col[4].k = ez  ; col[4].loc = DOWN;    col[4].c  = 0; valA[4] =  1.0 / (hx*hx);
            col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez-1; col[5].loc = DOWN;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
            /* Front term missing */
            col[6].i = ex  ; col[6].j = ey-1;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hy;
            col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hy;
          } else {
            nEntries = 9;
            col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = DOWN;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -2.0 / (hy*hy) -2.0 / (hz*hz);
            col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = DOWN;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
            col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = DOWN;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
            col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = DOWN;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
            col[4].i = ex+1; col[4].j = ey  ;  col[4].k = ez  ; col[4].loc = DOWN;    col[4].c  = 0; valA[4] =  1.0 / (hx*hx);
            col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez-1; col[5].loc = DOWN;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
            col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez+1; col[6].loc = DOWN;    col[6].c  = 0; valA[6] =  1.0 / (hz*hz);
            col[7].i = ex  ; col[7].j = ey-1;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] =  1.0 / hy;
            col[8].i = ex  ; col[8].j = ey  ;  col[8].k = ez  ; col[8].loc = ELEMENT; col[8].c  = 0; valA[8] = -1.0 / hy;
          }
          ierr = DMStagMatGetValuesStencil(dmSol,A,1,&row,nEntries,col,computed);CHKERRQ(ierr);
          ierr = check_vals(ex,ey,ez,nEntries,valA,computed);CHKERRQ(ierr);
        }

        /* Equation on back face of this element */
        if (ez == 0) {
          /* Back boundary velocity Dirichlet */
          DMStagStencil row;
          const PetscScalar valA = 1.0;
          row.i = ex; row.j = ey; row.k = ez; row.loc = BACK; row.c = 0;
          ierr = DMStagMatGetValuesStencil(dmSol,A,1,&row,1,&row,computed);CHKERRQ(ierr);
          ierr = check_vals(ex,ey,ez,1,&valA,computed);CHKERRQ(ierr);
        } else {
          /* Z-momentum equation, (w_xx + w_yy + w_zz) - p_z = f^z */
          DMStagStencil row,col[9];
          PetscScalar   valA[9];
          PetscInt      nEntries;

          row.i = ex; row.j = ey; row.k = ez; row.loc = BACK; row.c = 0;
          if (ex == 0) {
            if (ey == 0) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = BACK;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -1.0 / (hy*hy) -2.0 / (hz*hz);
              /* Down term missing */
              col[1].i = ex  ; col[1].j = ey+1;  col[1].k = ez  ; col[1].loc = BACK;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              /* Left term missing */
              col[2].i = ex+1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = BACK;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
              col[3].i = ex  ; col[3].j = ey  ;  col[3].k = ez-1; col[3].loc = BACK;    col[3].c  = 0; valA[3] =  1.0 / (hz*hz);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez+1; col[4].loc = BACK;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez-1; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hz;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hz;
            } else if (ey == N[1]-1) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = BACK;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -1.0 / (hy*hy) -2.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = BACK;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              /* Up term missing */
              /* Left term missing */
              col[2].i = ex+1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = BACK;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
              col[3].i = ex  ; col[3].j = ey  ;  col[3].k = ez-1; col[3].loc = BACK;    col[3].c  = 0; valA[3] =  1.0 / (hz*hz);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez+1; col[4].loc = BACK;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez-1; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hz;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hz;
            } else {
              nEntries = 8;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = BACK;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -2.0 / (hy*hy) -2.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = BACK;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = BACK;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
              /* Left term missing */
              col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = BACK;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = BACK;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez+1; col[5].loc = BACK;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez-1; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hz;
              col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hz;
            }
          } else if (ex == N[0]-1) {
            if (ey == 0) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = BACK;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -1.0 / (hy*hy) -2.0 / (hz*hz);
              /* Down term missing */
              col[1].i = ex  ; col[1].j = ey+1;  col[1].k = ez  ; col[1].loc = BACK;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex-1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = BACK;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
              /* Right term missing */
              col[3].i = ex  ; col[3].j = ey  ;  col[3].k = ez-1; col[3].loc = BACK;    col[3].c  = 0; valA[3] =  1.0 / (hz*hz);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez+1; col[4].loc = BACK;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez-1; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hz;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hz;
            } else if (ey == N[1]-1) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = BACK;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -1.0 / (hy*hy) -2.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = BACK;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              /* Up term missing */
              col[2].i = ex-1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = BACK;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
              /* Right term missing */
              col[3].i = ex  ; col[3].j = ey  ;  col[3].k = ez-1; col[3].loc = BACK;    col[3].c  = 0; valA[3] =  1.0 / (hz*hz);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez+1; col[4].loc = BACK;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez-1; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hz;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hz;
            } else {
              nEntries = 8;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = BACK;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -2.0 / (hy*hy) -2.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = BACK;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = BACK;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
              col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = BACK;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              /* Right term missing */
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = BACK;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez+1; col[5].loc = BACK;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez-1; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hz;
              col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hz;
            }
          } else if (ey == 0) {
            nEntries = 8;
            col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = BACK;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -1.0 / (hy*hy) -2.0 / (hz*hz);
            /* Down term missing */
            col[1].i = ex  ; col[1].j = ey+1;  col[1].k = ez  ; col[1].loc = BACK;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
            col[2].i = ex-1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = BACK;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
            col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = BACK;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
            col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = BACK;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
            col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez+1; col[5].loc = BACK;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
            col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez-1; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hz;
            col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hz;
          } else if (ey == N[1]-1) {
            nEntries = 8;
            col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = BACK;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -1.0 / (hy*hy) -2.0 / (hz*hz);
            col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = BACK;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
            /* Up term missing */
            col[2].i = ex-1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = BACK;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
            col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = BACK;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
            col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = BACK;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
            col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez+1; col[5].loc = BACK;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
            col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez-1; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hz;
            col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hz;
          } else {
            nEntries = 9;
            col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = BACK;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -2.0 / (hy*hy) -2.0 / (hz*hz);
            col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = BACK;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
            col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = BACK;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
            col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = BACK;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
            col[4].i = ex+1; col[4].j = ey  ;  col[4].k = ez  ; col[4].loc = BACK;    col[4].c  = 0; valA[4] =  1.0 / (hx*hx);
            col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez-1; col[5].loc = BACK;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
            col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez+1; col[6].loc = BACK;    col[6].c  = 0; valA[6] =  1.0 / (hz*hz);
            col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez-1; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] =  1.0 / hz;
            col[8].i = ex  ; col[8].j = ey  ;  col[8].k = ez  ; col[8].loc = ELEMENT; col[8].c  = 0; valA[8] = -1.0 / hz;
          }
          ierr = DMStagMatGetValuesStencil(dmSol,A,1,&row,nEntries,col,computed);CHKERRQ(ierr);
          ierr = check_vals(ex,ey,ez,nEntries,valA,computed);CHKERRQ(ierr);
        }

        /* P equation : u_x + v_y + w_z = g
           Note that this includes an explicit zero on the diagonal. This is only needed for
           direct solvers (not required if using an iterative solver and setting the constant-pressure nullspace) */
        {
          DMStagStencil row,col[7];
          PetscScalar   valA[7];

          row.i    = ex; row.j    = ey; row.k    = ez; row.loc    = ELEMENT; row.c    = 0;
          col[0].i = ex; col[0].j = ey; col[0].k = ez; col[0].loc = LEFT;    col[0].c = 0; valA[0] = -1.0 / hx;
          col[1].i = ex; col[1].j = ey; col[1].k = ez; col[1].loc = RIGHT;   col[1].c = 0; valA[1] =  1.0 / hx;
          col[2].i = ex; col[2].j = ey; col[2].k = ez; col[2].loc = DOWN;    col[2].c = 0; valA[2] = -1.0 / hy;
          col[3].i = ex; col[3].j = ey; col[3].k = ez; col[3].loc = UP;      col[3].c = 0; valA[3] =  1.0 / hy;
          col[4].i = ex; col[4].j = ey; col[4].k = ez; col[4].loc = BACK;    col[4].c = 0; valA[4] = -1.0 / hz;
          col[5].i = ex; col[5].j = ey; col[5].k = ez; col[5].loc = FRONT;   col[5].c = 0; valA[5] =  1.0 / hz;
          col[6]   = row;                                                                  valA[6] =  0.0;
          ierr = DMStagMatGetValuesStencil(dmSol,A,1,&row,7,col,computed);CHKERRQ(ierr);
          ierr = check_vals(ex,ey,ez,7,valA,computed);CHKERRQ(ierr);
        }
      }
    }
  }
  ierr = DMStagVecRestoreArrayRead(dmCoord,coordLocal,&arrCoord);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*TEST

   test:
      suffix: 1
      nsize: 1

   test:
      suffix: 2
      nsize: 4

TEST*/
