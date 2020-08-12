static char help[] = "Test resource recycling and MPI_Comm and keyval creation in mpi or mpiuni\n";

#include <petscsys.h>

#define  CHKMPIERR(ierr)  do {if (ierr) MPI_Abort(MPI_COMM_WORLD, ierr);} while(0)

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscMPIInt    key1,key2,attr1=100,attr2=200,*attr,flag;
  MPI_Comm       newcomm;

  ierr = MPI_Init(&argc,&argv);CHKMPIERR(ierr);

  /* Repeated keyval or comm create/free should not blow up MPI */
  for (i=0; i<500; i++) {
    ierr = MPI_Comm_create_keyval(MPI_COMM_NULL_COPY_FN,MPI_COMM_NULL_DELETE_FN,&key1,NULL);CHKMPIERR(ierr);
    ierr = MPI_Comm_free_keyval(&key1);CHKMPIERR(ierr);
    ierr = MPI_Comm_dup(MPI_COMM_WORLD,&newcomm);CHKMPIERR(ierr);
    ierr = MPI_Comm_free(&newcomm);CHKMPIERR(ierr);
  }

  /* The following keyval/attr code exposes a bug in old mpiuni code, where it had wrong newcomm returned in MPI_Comm_dup. */
  ierr = MPI_Comm_create_keyval(MPI_COMM_NULL_COPY_FN,MPI_COMM_NULL_DELETE_FN,&key1,NULL);CHKMPIERR(ierr);
  ierr = MPI_Comm_create_keyval(MPI_COMM_NULL_COPY_FN,MPI_COMM_NULL_DELETE_FN,&key2,NULL);CHKMPIERR(ierr);
  ierr = MPI_Comm_dup(MPI_COMM_WORLD,&newcomm);CHKMPIERR(ierr);
  if (MPI_COMM_WORLD == newcomm) printf("Error: wrong newcomm returned by MPI_Comm_dup()\n");

  ierr = MPI_Comm_set_attr(MPI_COMM_WORLD,key1,&attr1);CHKMPIERR(ierr);
  ierr = MPI_Comm_set_attr(newcomm,key2,&attr2);CHKMPIERR(ierr);
  ierr = MPI_Comm_get_attr(newcomm,key1,&attr,&flag);CHKMPIERR(ierr);
  if (flag) printf("Error: newcomm should not have attribute for keyval %d\n", (int)key1);
  ierr = MPI_Comm_get_attr(MPI_COMM_WORLD,key1,&attr,&flag);CHKMPIERR(ierr);
  if (*attr != attr1) printf("Error: expected attribute %d, but got %d\n", (int)attr1, (int)*attr);
  ierr = MPI_Comm_get_attr(newcomm,key2,&attr,&flag);CHKMPIERR(ierr);
  if (*attr != attr2) printf("Error: expected attribute %d, but got %d\n", (int)attr2, (int)*attr);

  ierr = MPI_Comm_delete_attr(MPI_COMM_WORLD,key1);CHKMPIERR(ierr);
  ierr = MPI_Comm_delete_attr(newcomm,key2);CHKMPIERR(ierr);
  ierr = MPI_Comm_free_keyval(&key1);CHKMPIERR(ierr);
  ierr = MPI_Comm_free_keyval(&key2);CHKMPIERR(ierr);
  ierr = MPI_Comm_free(&newcomm);CHKMPIERR(ierr);

  /* Init/Finalize PETSc multiple times when MPI is initialized */
  for (i=0; i<500; i++) {
    ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
    ierr = PetscFinalize();if (ierr) return ierr;
  }

  ierr = MPI_Finalize();
  return ierr;
}

/*TEST
   # Elemental in debug mode has bugs that it can not be repeatedly init/finalize'd for more than 300 times
   testset:
    output_file: output/ex53_1.out
    test:
      suffix: 1
      requires: !elemental

    test:
      suffix: 2
      requires: elemental !define(PETSC_USE_DEBUG)
TEST*/
