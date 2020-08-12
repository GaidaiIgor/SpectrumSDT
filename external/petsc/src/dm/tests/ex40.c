
static char help[] = "Tests mirror boundary conditions in 2-d.\n\n";

#include <petscdm.h>
#include <petscdmda.h>

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  PetscInt       M = 8, N = 8,stencil_width = 1, dof = 1,m,n,xstart,ystart,i,j,c;
  DM             da;
  Vec            global,local;
  PetscScalar    ***vglobal;
  PetscViewer    sview;

  ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,0,"-stencil_width",&stencil_width,0);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,0,"-dof",&dof,0);CHKERRQ(ierr);

  ierr = DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_MIRROR,DM_BOUNDARY_MIRROR,DMDA_STENCIL_STAR,M,N,PETSC_DECIDE,PETSC_DECIDE,dof,stencil_width,NULL,NULL,&da);CHKERRQ(ierr);
  ierr = DMSetFromOptions(da);CHKERRQ(ierr);
  ierr = DMSetUp(da);CHKERRQ(ierr);  
  ierr = DMDAGetCorners(da,&xstart,&ystart,0,&m,&n,0);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(da,&global);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(da,global,&vglobal);CHKERRQ(ierr);
  for (j=ystart; j<ystart+n; j++) {
    for (i=xstart; i<xstart+m; i++) {
      for (c=0; c<dof; c++) {
        vglobal[j][i][c] = 100*j + 10*(i+1) + c;
      }
    }
  }
  ierr = DMDAVecRestoreArrayDOF(da,global,&vglobal);CHKERRQ(ierr);

  ierr = DMCreateLocalVector(da,&local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(da,global,INSERT_VALUES,local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da,global,INSERT_VALUES,local);CHKERRQ(ierr);

  ierr = PetscViewerGetSubViewer(PETSC_VIEWER_STDOUT_WORLD,PETSC_COMM_SELF,&sview);CHKERRQ(ierr);
  ierr = VecView(local,sview);CHKERRQ(ierr);
  ierr = PetscViewerRestoreSubViewer(PETSC_VIEWER_STDOUT_WORLD,PETSC_COMM_SELF,&sview);CHKERRQ(ierr);
  ierr = PetscViewerFlush(PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = VecView(global,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  ierr = DMDestroy(&da);CHKERRQ(ierr);
  ierr = VecDestroy(&local);CHKERRQ(ierr);
  ierr = VecDestroy(&global);CHKERRQ(ierr);

  ierr = PetscFinalize();
  return ierr;
}







/*TEST

   test:

   test:
      suffix: 2
      nsize: 4
      filter: grep -v "Vec Object"

TEST*/
