/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Private data structure used by the ARPACK interface
*/

#if !defined(SLEPC_ARPACK_H)
#define SLEPC_ARPACK_H

typedef struct {
  PetscInt     *select;
  PetscScalar  *workev;
  PetscScalar  *workd;
  PetscScalar  *workl;
  PetscInt     lworkl;
  PetscReal    *rwork;
} EPS_ARPACK;

/*
   Definition of routines from the ARPACK package
*/

#if defined(PETSC_HAVE_MPIUNI)

#if defined(PETSC_USE_COMPLEX)

#if defined(PETSC_USE_REAL_SINGLE)

#if defined(SLEPC_ARPACK_HAVE_UNDERSCORE)
#define SLEPC_ARPACK(lcase,ucase) c##lcase##_
#elif defined(SLEPC_ARPACK_HAVE_CAPS)
#define SLEPC_ARPACK(lcase,ucase) C##ucase
#else
#define SLEPC_ARPACK(lcase,ucase) c##lcase
#endif

#else

#if defined(SLEPC_ARPACK_HAVE_UNDERSCORE)
#define SLEPC_ARPACK(lcase,ucase) z##lcase##_
#elif defined(SLEPC_ARPACK_HAVE_CAPS)
#define SLEPC_ARPACK(lcase,ucase) Z##ucase
#else
#define SLEPC_ARPACK(lcase,ucase) z##lcase
#endif

#endif

#else

#if defined(PETSC_USE_REAL_SINGLE)

#if defined(SLEPC_ARPACK_HAVE_UNDERSCORE)
#define SLEPC_ARPACK(lcase,ucase) s##lcase##_
#elif defined(SLEPC_ARPACK_HAVE_CAPS)
#define SLEPC_ARPACK(lcase,ucase) S##ucase
#else
#define SLEPC_ARPACK(lcase,ucase) s##lcase
#endif

#else

#if defined(SLEPC_ARPACK_HAVE_UNDERSCORE)
#define SLEPC_ARPACK(lcase,ucase) d##lcase##_
#elif defined(SLEPC_ARPACK_HAVE_CAPS)
#define SLEPC_ARPACK(lcase,ucase) D##ucase
#else
#define SLEPC_ARPACK(lcase,ucase) d##lcase
#endif

#endif

#endif

#else  /* not MPIUNI */

#if defined(PETSC_USE_COMPLEX)

#if defined(PETSC_USE_REAL_SINGLE)

#if defined(SLEPC_ARPACK_HAVE_UNDERSCORE)
#define SLEPC_ARPACK(lcase,ucase) pc##lcase##_
#elif defined(SLEPC_ARPACK_HAVE_CAPS)
#define SLEPC_ARPACK(lcase,ucase) PC##ucase
#else
#define SLEPC_ARPACK(lcase,ucase) pc##lcase
#endif

#else

#if defined(SLEPC_ARPACK_HAVE_UNDERSCORE)
#define SLEPC_ARPACK(lcase,ucase) pz##lcase##_
#elif defined(SLEPC_ARPACK_HAVE_CAPS)
#define SLEPC_ARPACK(lcase,ucase) PZ##ucase
#else
#define SLEPC_ARPACK(lcase,ucase) pz##lcase
#endif

#endif

#else

#if defined(PETSC_USE_REAL_SINGLE)

#if defined(SLEPC_ARPACK_HAVE_UNDERSCORE)
#define SLEPC_ARPACK(lcase,ucase) ps##lcase##_
#elif defined(SLEPC_ARPACK_HAVE_CAPS)
#define SLEPC_ARPACK(lcase,ucase) PS##ucase
#else
#define SLEPC_ARPACK(lcase,ucase) ps##lcase
#endif

#else

#if defined(SLEPC_ARPACK_HAVE_UNDERSCORE)
#define SLEPC_ARPACK(lcase,ucase) pd##lcase##_
#elif defined(SLEPC_ARPACK_HAVE_CAPS)
#define SLEPC_ARPACK(lcase,ucase) PD##ucase
#else
#define SLEPC_ARPACK(lcase,ucase) pd##lcase
#endif

#endif

#endif

#endif

#if defined(PETSC_HAVE_MPIUNI)

#define COMM_ARG

#if !defined(PETSC_USE_COMPLEX)

#define ARPACKnaupd_(comm,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p) SLEPC_ARPACK(naupd,NAUPD) ((a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o),(p),1,2)
#define ARPACKneupd_(comm,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y) SLEPC_ARPACK(neupd,NEUPD) ((a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o),(p),(q),(r),(s),(t),(u),(v),(w),(x),(y),1,1,2)
#define ARPACKsaupd_(comm,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p) SLEPC_ARPACK(saupd,SAUPD) ((a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o),(p),1,2)
#define ARPACKseupd_(comm,a,b,c,d,e,f,g,h,i,j,k,l,m,o,p,q,r,s,t,u,v,w) SLEPC_ARPACK(seupd,SEUPD) ((a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(o),(p),(q),(r),(s),(t),(u),(v),(w),1,1,2)

#else

#define ARPACKnaupd_(comm,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q) SLEPC_ARPACK(naupd,NAUPD) ((a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o),(p),(q),1,2)
#define ARPACKneupd_(comm,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x) SLEPC_ARPACK(neupd,NEUPD) ((a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o),(p),(q),(r),(s),(t),(u),(v),(w),(x),1,1,2)

#endif

#else /* not MPIUNI */

#define COMM_ARG MPI_Fint*,

#if !defined(PETSC_USE_COMPLEX)

#define ARPACKnaupd_(comm,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p) SLEPC_ARPACK(naupd,NAUPD) ((comm),(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o),(p),1,2)
#define ARPACKneupd_(comm,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y) SLEPC_ARPACK(neupd,NEUPD) ((comm),(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o),(p),(q),(r),(s),(t),(u),(v),(w),(x),(y),1,1,2)
#define ARPACKsaupd_(comm,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p) SLEPC_ARPACK(saupd,SAUPD) ((comm),(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o),(p),1,2)
#define ARPACKseupd_(comm,a,b,c,d,e,f,g,h,i,j,k,l,m,o,p,q,r,s,t,u,v,w) SLEPC_ARPACK(seupd,SEUPD) ((comm),(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(o),(p),(q),(r),(s),(t),(u),(v),(w),1,1,2)

#else

#define ARPACKnaupd_(comm,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q) SLEPC_ARPACK(naupd,NAUPD) ((comm),(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o),(p),(q),1,2)
#define ARPACKneupd_(comm,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x) SLEPC_ARPACK(neupd,NEUPD) ((comm),(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o),(p),(q),(r),(s),(t),(u),(v),(w),(x),1,1,2)

#endif

#endif

SLEPC_EXTERN void   SLEPC_ARPACK(saupd,SAUPD)(COMM_ARG PetscInt*,char*,PetscInt*,const char*,PetscInt*,PetscReal*,PetscScalar*,PetscInt*,PetscScalar*,PetscInt*,PetscInt*,PetscInt*,PetscScalar*,PetscScalar*,PetscInt*,PetscInt*,int,int);
SLEPC_EXTERN void   SLEPC_ARPACK(seupd,SEUPD)(COMM_ARG PetscInt*,char*,PetscInt*,PetscReal*,PetscReal*,PetscInt*,PetscReal*,char*,PetscInt*,const char*,PetscInt*,PetscReal*,PetscScalar*,PetscInt*,PetscScalar*,PetscInt*,PetscInt*,PetscInt*,PetscScalar*,PetscScalar*,PetscInt*,PetscInt*,int,int,int);

#if !defined(PETSC_USE_COMPLEX)
SLEPC_EXTERN void   SLEPC_ARPACK(naupd,NAUPD)(COMM_ARG PetscInt*,char*,PetscInt*,const char*,PetscInt*,PetscReal*,PetscScalar*,PetscInt*,PetscScalar*,PetscInt*,PetscInt*,PetscInt*,PetscScalar*,PetscScalar*,PetscInt*,PetscInt*,int,int);
SLEPC_EXTERN void   SLEPC_ARPACK(neupd,NEUPD)(COMM_ARG PetscInt*,char*,PetscInt*,PetscReal*,PetscReal*,PetscReal*,PetscInt*,PetscReal*,PetscReal*,PetscReal*,char*,PetscInt*,const char*,PetscInt*,PetscReal*,PetscScalar*,PetscInt*,PetscScalar*,PetscInt*,PetscInt*,PetscInt*,PetscScalar*,PetscScalar*,PetscInt*,PetscInt*,int,int,int);
#else
SLEPC_EXTERN void   SLEPC_ARPACK(naupd,NAUPD)(COMM_ARG PetscInt*,char*,PetscInt*,const char*,PetscInt*,PetscReal*,PetscScalar*,PetscInt*,PetscScalar*,PetscInt*,PetscInt*,PetscInt*,PetscScalar*,PetscScalar*,PetscInt*,PetscReal*,PetscInt*,int,int);
SLEPC_EXTERN void   SLEPC_ARPACK(neupd,NEUPD)(COMM_ARG PetscInt*,char*,PetscInt*,PetscScalar*,PetscScalar*,PetscInt*,PetscScalar*,PetscScalar*,char*,PetscInt*,const char*,PetscInt*,PetscReal*,PetscScalar*,PetscInt*,PetscScalar*,PetscInt*,PetscInt*,PetscInt*,PetscScalar*,PetscScalar*,PetscInt*,PetscReal*,PetscInt*,int,int,int);
#endif

#endif

