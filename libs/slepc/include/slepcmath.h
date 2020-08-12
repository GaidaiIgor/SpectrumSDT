/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc mathematics include file. Defines basic operations and functions.
   This file is included by slepcsys.h and should not be used directly.
*/

#if !defined(SLEPCMATH_H)
#define SLEPCMATH_H

/*
    Default tolerance for the different solvers, depending on the precision
*/
#if defined(PETSC_USE_REAL_SINGLE)
#  define SLEPC_DEFAULT_TOL   1e-6
#elif defined(PETSC_USE_REAL_DOUBLE)
#  define SLEPC_DEFAULT_TOL   1e-8
#elif defined(PETSC_USE_REAL___FLOAT128)
#  define SLEPC_DEFAULT_TOL   1e-16
#else
#  define SLEPC_DEFAULT_TOL   1e-7
#endif

/*@C
   SlepcAbs - Returns sqrt(x**2+y**2), taking care not to cause unnecessary
   overflow. It is based on LAPACK's DLAPY2.

   Not Collective

   Input parameters:
.  x,y - the real numbers

   Output parameter:
.  return - the result

   Note:
   This function is not available from Fortran.

   Level: developer
@*/
PETSC_STATIC_INLINE PetscReal SlepcAbs(PetscReal x,PetscReal y)
{
  PetscReal w,z,t,xabs=PetscAbs(x),yabs=PetscAbs(y);

  w = PetscMax(xabs,yabs);
  z = PetscMin(xabs,yabs);
  if (z == 0.0) return w;
  t = z/w;
  return w*PetscSqrtReal(1.0+t*t);
}

/*MC
   SlepcAbsEigenvalue - Returns the absolute value of a complex number given
   its real and imaginary parts.

   Synopsis:
   PetscReal SlepcAbsEigenvalue(PetscScalar x,PetscScalar y)

   Not Collective

   Input parameters:
+  x  - the real part of the complex number
-  y  - the imaginary part of the complex number

   Notes:
   This function computes sqrt(x**2+y**2), taking care not to cause unnecessary
   overflow. It is based on LAPACK's DLAPY2.

   This function is not available from Fortran.

   Level: developer
M*/
#if !defined(PETSC_USE_COMPLEX)
#define SlepcAbsEigenvalue(x,y) SlepcAbs(x,y)
#else
#define SlepcAbsEigenvalue(x,y) PetscAbsScalar(x)
#endif

#endif


/*
   SlepcSetFlushToZero - Set the FTZ flag in floating-point arithmetic.
*/
PETSC_STATIC_INLINE PetscErrorCode SlepcSetFlushToZero(unsigned int *state)
{
  PetscFunctionBegin;
#if defined(PETSC_HAVE_XMMINTRIN_H)
#if defined(_MM_FLUSH_ZERO_ON) && defined(__SSE__)
  *state = _MM_GET_FLUSH_ZERO_MODE();
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#else
  *state = 0;
#endif
#endif
  PetscFunctionReturn(0);
}

/*
   SlepcResetFlushToZero - Reset the FTZ flag in floating-point arithmetic.
*/
PETSC_STATIC_INLINE PetscErrorCode SlepcResetFlushToZero(unsigned int *state)
{
  PetscFunctionBegin;
#if defined(PETSC_HAVE_XMMINTRIN_H)
#if defined(_MM_FLUSH_ZERO_MASK) && defined(__SSE__)
  _MM_SET_FLUSH_ZERO_MODE(*state & _MM_FLUSH_ZERO_MASK);
#endif
#endif
  PetscFunctionReturn(0);
}

