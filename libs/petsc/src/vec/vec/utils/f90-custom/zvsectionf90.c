#include <petscvec.h>
#include <petsc/private/f90impl.h>

#if defined(PETSC_HAVE_FORTRAN_CAPS)
#define vecsetvaluessectionf90_ VECSETVALUESSECTIONF90
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define vecsetvaluessectionf90_ vecsetvaluessectionf90
#endif

/* Definitions of Fortran Wrapper routines */

PETSC_EXTERN void vecsetvaluessectionf90_(Vec *v, PetscSection *section, PetscInt *point, F90Array1d *ptr, InsertMode *mode, int *__ierr PETSC_F90_2PTR_PROTO(ptrd))
{
  PetscScalar *array;

  *__ierr = F90Array1dAccess(ptr, MPIU_SCALAR, (void**) &array PETSC_F90_2PTR_PARAM(ptrd));if (*__ierr) return;
  *__ierr = VecSetValuesSection(*v, *section, *point, array, *mode);
}
