#
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  SLEPc - Scalable Library for Eigenvalue Problem Computations
#  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
#
#  This file is part of SLEPc.
#  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#

CFLAGS     =
FFLAGS     =
CPPFLAGS   =
FPPFLAGS   =
LOCDIR     = src/eps/tests/
EXAMPLESC  = test1.c test2.c test3.c test4.c test5.c test6.c \
             test8.c test9.c test10.c test11.c test12.c test13.c \
             test14.c test16.c test17.c test18.c test19.c test20.c \
             test21.c test22.c test23.c test24.c test25.c test26.c test27.c \
             test28.c test29.c test30.c test31.c test32.c test33.c test34.c \
             test35.c test36.c test37.c
EXAMPLESF  = test7f.F test14f.F test15f.F test17f.F90
MANSEC     = EPS
TESTS      = test1 test2 test3 test4 test5 test6 test7f test8 test9 test10 \
             test11 test12 test13 test14 test14f test15f test16 test17 test17f \
             test18 test19 test20 test21 test22 test23 test24 test25 test26 test27 \
             test28 test29 test30 test31 test32 test33 test34 test35 test36 test37

include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

#------------------------------------------------------------------------------------

testtest5_blopex: test5.PETSc
	@fail=0; if [ "${PETSC_WITH_BATCH}" != "" -o "${MPIEXEC}" = "/bin/false" ]; then \
	  echo "Skipping BLOPEX test"; \
	elif [ -f test5 ]; then \
	  test=test5_2_blopex; check=test5_2; \
	  if [ "${PETSC_SCALAR}" = "complex" ]; then check=$${check}_complex; fi; \
	  ${MPIEXEC} -n 1 ./test5 -symm -eps_type blopex -eps_nev 4 > $${test}.tmp 2>&1; \
	  if (${DIFF} output/$${check}.out $${test}.tmp > /dev/null 2>&1) then \
	    echo "BLOPEX example src/eps/tests/test5 run successfully with 1 MPI process"; \
	  else \
	    echo "Possible error running BLOPEX src/eps/tests/test5 with 1 MPI process"; \
	    cat $${test}.tmp; fail=1; \
	  fi; \
	  ${RM} $${test}.tmp; \
	  ${MAKE} SLEPC_DIR=${SLEPC_DIR} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} test5.rm ; \
	else fail=1; fi; \
	exit $$fail

testtest7f: test7f.PETSc
	@fail=0; if [ "${PETSC_WITH_BATCH}" != "" ]; then \
	  echo "Running with batch filesystem; to test run src/eps/tests/test7f " ; \
	  echo "with your systems batch system"; \
	elif [ "${MPIEXEC}" = "/bin/false" ]; then \
	  echo "*mpiexec not found*. Please run src/eps/tests/test7f manually"; \
	elif [ -f test7f ]; then \
	  test=test7f_1; check=test7f_1; \
	  GFORTRAN_UNBUFFERED_ALL=1 ${MPIEXEC} -n 1 ./test7f -eps_nev 4 -eps_ncv 22 > $${test}.tmp 2>&1; \
	  if (${DIFF} output/$${check}.out $${test}.tmp > /dev/null 2>&1) then \
	    echo "Fortran example src/eps/tests/test7f run successfully with 1 MPI process"; \
	  else \
	    echo "Possible error running Fortran src/eps/tests/test7f with 1 MPI process"; \
	    cat $${test}.tmp; fail=1; \
	  fi; \
	  ${RM} $${test}.tmp; \
	  ${MAKE} SLEPC_DIR=${SLEPC_DIR} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} test7f.rm ; \
	else fail=1; fi; \
	exit $$fail

testtest10: test10.PETSc
	@fail=0; if [ "${PETSC_WITH_BATCH}" != "" ]; then \
	  echo "Running with batch filesystem; to test run src/eps/tests/test10" ; \
	  echo "with your systems batch system"; \
	elif [ "${MPIEXEC}" = "/bin/false" ]; then \
	  echo "*mpiexec not found*. Please run src/eps/tests/test10 manually"; \
	elif [ -f test10 ]; then \
	  test=test10_1; check=test10_1_ks; \
	  options="-eps_nev 4 -eps_ncv 14 -m 11 -eps_largest_magnitude"; \
	  ${MPIEXEC} -n 1 ./test10 $$options > $${test}.tmp 2>&1; \
	  if (${DIFF} output/$${check}.out $${test}.tmp > /dev/null 2>&1) then \
	    echo "C/C++ example src/eps/tests/test10 run successfully with 1 MPI process"; \
	  else \
	    echo "Possible error running C/C++ src/eps/tests/test10 with 1 MPI process"; \
	    cat $${test}.tmp; fail=1; \
	  fi; \
	  if [ "${MPIEXEC}" != "${PETSC_DIR}/lib/petsc/bin/petsc-mpiexec.uni" ]; then \
	    ${MPIEXEC} -n 2 ./test10 $$options > $${test}.tmp 2>&1; \
	    if (${DIFF} output/$${check}.out $${test}.tmp > /dev/null 2>&1) then \
	      echo "C/C++ example src/eps/tests/test10 run successfully with 2 MPI process"; \
	    else \
	      echo "Possible error running C/C++ src/eps/tests/test10 with 2 MPI process"; \
	      cat $${test}.tmp; fail=1; \
	    fi; \
	  fi; \
	  ${RM} $${test}.tmp; \
	  ${MAKE} SLEPC_DIR=${SLEPC_DIR} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} test10.rm; \
	else fail=1; fi; \
	exit $$fail

