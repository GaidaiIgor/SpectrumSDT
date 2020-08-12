static char help[] = "Test DMStag ghosted boundaries in 3d\n\n";
/* This solves a very contrived problem - the "pressure" terms are set to a constant function
   and the "velocity" terms are just the sum of neighboring values of these, hence twice the
   constant */
#include <petscdm.h>
#include <petscksp.h>
#include <petscdmstag.h>

#define PRESSURE_CONST 2.0

PetscErrorCode ApplyOperator(Mat,Vec,Vec);

int main(int argc,char **argv)
{
  PetscErrorCode  ierr;
  DM              dmSol;
  Vec             sol,solRef,solRefLocal,rhs,rhsLocal;
  Mat             A;
  KSP             ksp;
  PC              pc;
  PetscInt        startx,starty,startz,nx,ny,nz,ex,ey,ez,nExtrax,nExtray,nExtraz;
  PetscInt        iux,iuy,iuz,ip;
  PetscInt        dof0,dof1,dof2,dof3;
  PetscScalar     ****arrSol,****arrRHS;

  ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  /* Note: these defaults are chosen to suit the problem. We allow adjusting
     them to check that things still work when you add unused extra dof */
  dof0 = 0;
  dof1 = 0;
  dof2 = 1;
  dof3 = 1;
  ierr = DMStagCreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_GHOSTED,DM_BOUNDARY_GHOSTED,DM_BOUNDARY_GHOSTED,3,3,3,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,dof0,dof1,dof2,dof3,DMSTAG_STENCIL_BOX,1,NULL,NULL,NULL,&dmSol);CHKERRQ(ierr);
  ierr = DMSetFromOptions(dmSol);CHKERRQ(ierr);
  ierr = DMSetUp(dmSol);CHKERRQ(ierr);

  /* Compute reference solution on the grid, using direct array access */
  ierr = DMCreateGlobalVector(dmSol,&rhs);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(dmSol,&solRef);CHKERRQ(ierr);
  ierr = DMGetLocalVector(dmSol,&solRefLocal);CHKERRQ(ierr);
  ierr = DMGetLocalVector(dmSol,&rhsLocal);CHKERRQ(ierr);
  ierr = DMStagVecGetArray(dmSol,solRefLocal,&arrSol);CHKERRQ(ierr);

  ierr = DMStagGetCorners(dmSol,&startx,&starty,&startz,&nx,&ny,&nz,&nExtrax,&nExtray,&nExtraz);CHKERRQ(ierr);
  ierr = DMStagVecGetArray(dmSol,rhsLocal,&arrRHS);CHKERRQ(ierr);

  /* Get the correct entries for each of our variables in local element-wise storage */
  ierr = DMStagGetLocationSlot(dmSol,DMSTAG_LEFT,0,&iux);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(dmSol,DMSTAG_DOWN,0,&iuy);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(dmSol,DMSTAG_BACK,0,&iuz);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(dmSol,DMSTAG_ELEMENT,0,&ip);CHKERRQ(ierr);
  for (ez=startz; ez<startz+nz+nExtraz; ++ez) {
    for (ey=starty; ey<starty+ny+nExtray; ++ey) {
      for (ex=startx; ex<startx+nx+nExtrax; ++ex) {
        arrSol[ez][ey][ex][iux] = 2*PRESSURE_CONST;
        arrRHS[ez][ey][ex][iux] = 0.0;
        arrSol[ez][ey][ex][iuy] = 2*PRESSURE_CONST;
        arrRHS[ez][ey][ex][iuy] = 0.0;
        arrSol[ez][ey][ex][iuz] = 2*PRESSURE_CONST;
        arrRHS[ez][ey][ex][iuz] = 0.0;
        if (ex < startx+nx && ey < starty+ny && ez < startz+nz) {
          arrSol[ez][ey][ex][ip] = PRESSURE_CONST;
          arrRHS[ez][ey][ex][ip] = PRESSURE_CONST;
        }
      }
    }
  }
  ierr = DMStagVecRestoreArray(dmSol,rhsLocal,&arrRHS);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(dmSol,rhsLocal,INSERT_VALUES,rhs);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(dmSol,rhsLocal,INSERT_VALUES,rhs);CHKERRQ(ierr);
  ierr = DMStagVecRestoreArray(dmSol,solRefLocal,&arrSol);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(dmSol,solRefLocal,INSERT_VALUES,solRef);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(dmSol,solRefLocal,INSERT_VALUES,solRef);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dmSol,&solRefLocal);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dmSol,&rhsLocal);CHKERRQ(ierr);

  /* Matrix-free Operator */
  ierr = DMSetMatType(dmSol,MATSHELL);CHKERRQ(ierr);
  ierr = DMCreateMatrix(dmSol,&A);CHKERRQ(ierr);
  ierr = MatShellSetOperation(A,MATOP_MULT,(void(*) (void)) ApplyOperator);CHKERRQ(ierr);

  /* Solve */
  ierr = DMCreateGlobalVector(dmSol,&sol);CHKERRQ(ierr);
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCNONE);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = KSPSolve(ksp,rhs,sol);CHKERRQ(ierr);

  /* Check Solution */
  {
    Vec       diff;
    PetscReal normsolRef,errAbs,errRel;

    ierr = VecDuplicate(sol,&diff);CHKERRQ(ierr);
    ierr = VecCopy(sol,diff);CHKERRQ(ierr);
    ierr = VecAXPY(diff,-1.0,solRef);CHKERRQ(ierr);
    ierr = VecNorm(diff,NORM_2,&errAbs);CHKERRQ(ierr);
    ierr = VecNorm(solRef,NORM_2,&normsolRef);CHKERRQ(ierr);
    errRel = errAbs/normsolRef;
    if (errAbs > 1e14 || errRel > 1e14) {
      ierr = PetscPrintf(PetscObjectComm((PetscObject)dmSol),"Error (abs): %g\nError (rel): %g\n",(double)errAbs,(double)errRel);CHKERRQ(ierr);
      ierr = PetscPrintf(PetscObjectComm((PetscObject)dmSol),"Non-zero error. Probable failure.\n");CHKERRQ(ierr);
    }
    ierr = VecDestroy(&diff);CHKERRQ(ierr);
  }

  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = VecDestroy(&sol);CHKERRQ(ierr);
  ierr = VecDestroy(&solRef);CHKERRQ(ierr);
  ierr = VecDestroy(&rhs);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = DMDestroy(&dmSol);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}

PetscErrorCode ApplyOperator(Mat A,Vec in,Vec out)
{
  PetscErrorCode    ierr;
  DM                dm;
  Vec               inLocal,outLocal;
  PetscScalar       ****arrIn;
  PetscScalar       ****arrOut;
  PetscInt          startx,starty,startz,nx,ny,nz,nExtrax,nExtray,nExtraz,ex,ey,ez,idxP,idxUx,idxUy,idxUz,startGhostx,startGhosty,startGhostz,nGhostx,nGhosty,nGhostz;
  PetscBool         isFirstx,isFirsty,isFirstz,isLastx,isLasty,isLastz;

  PetscFunctionBeginUser;
  ierr = MatGetDM(A,&dm);CHKERRQ(ierr);
  ierr = DMGetLocalVector(dm,&inLocal);CHKERRQ(ierr);
  ierr = DMGetLocalVector(dm,&outLocal);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(dm,in,INSERT_VALUES,inLocal);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(dm,in,INSERT_VALUES,inLocal);CHKERRQ(ierr);
  ierr = DMStagGetCorners(dm,&startx,&starty,&startz,&nx,&ny,&nz,&nExtrax,&nExtray,&nExtraz);CHKERRQ(ierr);
  ierr = DMStagGetGhostCorners(dm,&startGhostx,&startGhosty,&startGhostz,&nGhostx,&nGhosty,&nGhostz);CHKERRQ(ierr);
  ierr = DMStagVecGetArrayRead(dm,inLocal,&arrIn);CHKERRQ(ierr);
  ierr = DMStagVecGetArray(dm,outLocal,&arrOut);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(dm,DMSTAG_LEFT,0,&idxUx);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(dm,DMSTAG_DOWN,0,&idxUy);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(dm,DMSTAG_BACK,0,&idxUz);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(dm,DMSTAG_ELEMENT,0,&idxP);CHKERRQ(ierr);
  ierr = DMStagGetIsFirstRank(dm,&isFirstx,&isFirsty,&isFirstz);CHKERRQ(ierr);
  ierr = DMStagGetIsLastRank(dm,&isLastx,&isLasty,&isLastz);CHKERRQ(ierr);

  /* Set "pressures" on ghost boundaries by copying neighboring values*/
  if (isFirstx) {
    for (ez=startz; ez<startz+nz+nExtraz; ++ez) {
      for (ey=starty; ey<starty+ny+nExtray; ++ey) {
        arrIn[ez][ey][-1][idxP] = arrIn[ez][ey][0][idxP];
      }
    }
  }
  if (isLastx){
    for (ez=startz; ez<startz+nz+nExtraz; ++ez) {
      for (ey=starty; ey<starty+ny+nExtray; ++ey) {
        arrIn[ez][ey][startx + nx][idxP] = arrIn[ez][ey][startx + nx - 1][idxP];
      }
    }
  }
  if (isFirsty) {
    for (ez=startz; ez<startz+nz+nExtraz; ++ez) {
      for (ex=startx; ex<startx+nx+nExtrax; ++ex) {
        arrIn[ez][-1][ex][idxP] = arrIn[ez][0][ex][idxP];
      }
    }
  }
  if  (isLasty) {
    for (ez=startz; ez<startz+nz+nExtraz; ++ez) {
      for (ex=startx; ex<startx+nx+nExtrax; ++ex) {
        arrIn[ez][starty + ny][ex][idxP] = arrIn[ez][starty + ny - 1][ex][idxP];
      }
    }
  }

  if (isFirstz) {
    for (ey=starty; ey<starty+ny+nExtray; ++ey) {
      for (ex=startx; ex<startx+nx+nExtrax; ++ex) {
        arrIn[-1][ey][ex][idxP] = arrIn[0][ey][ex][idxP];
      }
    }
  }
  if (isLastz) {
    for (ey=starty; ey<starty+ny+nExtray; ++ey) {
      for (ex=startx; ex<startx+nx+nExtrax; ++ex) {
        arrIn[startz + nz][ey][ex][idxP] = arrIn[startz + nz -1][ey][ex][idxP];
      }
    }
  }

  /* Apply operator on physical points */
  for (ez=startz; ez<startz+nz+nExtraz; ++ez) {
    for (ey=starty; ey<starty+ny+nExtray; ++ey) {
      for (ex=startx; ex<startx+nx+nExtrax; ++ex) {
        if (ex < startx+nx && ey < starty+ny && ez < startz+nz) {/* Don't compute pressure outside domain */
          arrOut[ez][ey][ex][idxP] = arrIn[ez][ey][ex][idxP];
        }
        arrOut[ez][ey][ex][idxUx] = arrIn[ez][ey][ex][idxP] + arrIn[ez][ey][ex-1][idxP] - arrIn[ez][ey][ex][idxUx];
        arrOut[ez][ey][ex][idxUy] = arrIn[ez][ey][ex][idxP] + arrIn[ez][ey-1][ex][idxP] - arrIn[ez][ey][ex][idxUy];
        arrOut[ez][ey][ex][idxUz] = arrIn[ez][ey][ex][idxP] + arrIn[ez-1][ey][ex][idxP] - arrIn[ez][ey][ex][idxUz];
      }
    }
  }
  ierr = DMStagVecRestoreArrayRead(dm,inLocal,&arrIn);CHKERRQ(ierr);
  ierr = DMStagVecRestoreArray(dm,outLocal,&arrOut);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(dm,outLocal,INSERT_VALUES,out);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(dm,outLocal,INSERT_VALUES,out);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dm,&inLocal);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dm,&outLocal);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*TEST

   test:
      suffix: 1
      nsize: 1

   test:
      suffix: 2
      nsize: 8

   test:
      suffix: 3
      nsize: 1
      args: -stag_stencil_width 2

   test:
      suffix: 4
      nsize: 8
      args: -stag_grid_x 4 -stag_grid_y 5 -stag_grid_z 4 -stag_stencil_width 2

   test:
      suffix: 5
      nsize: 8
      args: -stag_dof_0 3 -stag_dof_1 2 -stag_dof_2 4 -stag_dof_3 2 -stag_stencil_width 3 -stag_grid_x 6 -stag_grid_y 6 -stag_grid_z 6

TEST*/
