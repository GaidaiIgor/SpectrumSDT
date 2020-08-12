#
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  SLEPc - Scalable Library for Eigenvalue Problem Computations
#  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
#
#  This file is part of SLEPc.
#  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#  This is the makefile for installing SLEPc. See the Users Manual
#  for directions on installing SLEPc.
#

ALL: all
LOCDIR = ./
DIRS   = src include docs

# Include the rest of makefiles
include ./${PETSC_ARCH}/lib/slepc/conf/slepcvariables
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

# This makefile doesn't really do any work. Sub-makes still benefit from parallelism.
.NOTPARALLEL:

OMAKE_SELF = $(OMAKE) -f makefile

#
# Basic targets to build SLEPc library
all:
	+@${OMAKE_SELF} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR} chk_petscdir chk_slepcdir | tee ${PETSC_ARCH}/lib/slepc/conf/make.log
	@ln -sf ${PETSC_ARCH}/lib/slepc/conf/make.log make.log
	+@${OMAKE_SELF} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR} all-local 2>&1 | tee -a ${PETSC_ARCH}/lib/slepc/conf/make.log;
	@egrep -i "( error | error: |no such file or directory)" ${PETSC_ARCH}/lib/slepc/conf/make.log | tee ./${PETSC_ARCH}/lib/slepc/conf/error.log > /dev/null
	+@if test -s ./${PETSC_ARCH}/lib/slepc/conf/error.log; then \
           printf ${PETSC_TEXT_HILIGHT}"*******************************ERROR************************************\n" 2>&1 | tee -a ${PETSC_ARCH}/lib/slepc/conf/make.log; \
           echo "  Error during compile, check ${PETSC_ARCH}/lib/slepc/conf/make.log" 2>&1 | tee -a ${PETSC_ARCH}/lib/slepc/conf/make.log; \
           echo "  Send all contents of ./${PETSC_ARCH}/lib/slepc/conf to slepc-maint@upv.es" 2>&1 | tee -a ${PETSC_ARCH}/lib/slepc/conf/make.log;\
           printf "************************************************************************"${PETSC_TEXT_NORMAL}"\n" 2>&1 | tee -a ${PETSC_ARCH}/lib/slepc/conf/make.log; \
	 elif [ "${SLEPC_INSTALLDIR}" = "${SLEPC_DIR}/${PETSC_ARCH}" ]; then \
           echo "Now to check if the library is working do:";\
           echo "make SLEPC_DIR=${SLEPC_DIR} PETSC_DIR=${PETSC_DIR} check";\
           echo "=========================================";\
	 else \
	   echo "Now to install the library do:";\
	   echo "make SLEPC_DIR=${SLEPC_DIR} PETSC_DIR=${PETSC_DIR} install";\
	   echo "=========================================";\
	 fi
	@echo "Finishing make run at `date +'%a, %d %b %Y %H:%M:%S %z'`" >> ${PETSC_ARCH}/lib/slepc/conf/make.log
	@if test -s ./${PETSC_ARCH}/lib/slepc/conf/error.log; then exit 1; fi

all-local: info slepc_libs

#
# Prints information about the system and version of SLEPc being compiled
#
info:
	-@echo "=========================================="
	-@echo Starting make run on `hostname` at `date +'%a, %d %b %Y %H:%M:%S %z'`
	-@echo Machine characteristics: `uname -a`
	-@echo "-----------------------------------------"
	-@echo "Using SLEPc directory: ${SLEPC_DIR}"
	-@echo "Using PETSc directory: ${PETSC_DIR}"
	-@echo "Using PETSc arch: ${PETSC_ARCH}"
	-@echo "-----------------------------------------"
	-@grep "define SLEPC_VERSION" ${SLEPC_DIR}/include/slepcversion.h | ${SED} "s/........//"
	-@echo "-----------------------------------------"
	-@echo "Using SLEPc configure options: ${SLEPC_CONFIGURE_OPTIONS}"
	-@echo "Using SLEPc configuration flags:"
	-@grep "\#define " ${SLEPCCONF_H}
	-@echo "-----------------------------------------"
	-@grep "define PETSC_VERSION" ${PETSC_DIR}/include/petscversion.h | ${SED} "s/........//"
	-@echo "-----------------------------------------"
	-@echo "Using PETSc configure options: ${CONFIGURE_OPTIONS}"
	-@echo "Using PETSc configuration flags:"
	-@grep "\#define " ${PETSCCONF_H}
	-@echo "-----------------------------------------"
	-@echo "Using C/C++ include paths: ${SLEPC_CC_INCLUDES}"
	-@echo "Using C/C++ compiler: ${PCC} ${PCC_FLAGS} ${COPTFLAGS} ${CFLAGS}"
	-@if [ "${FC}" != "" ]; then \
	   echo "Using Fortran include/module paths: ${SLEPC_FC_INCLUDES}";\
	   echo "Using Fortran compiler: ${FC} ${FC_FLAGS} ${FFLAGS} ${FPP_FLAGS}";\
         fi
	-@echo "-----------------------------------------"
	-@echo "Using C/C++ linker: ${PCC_LINKER}"
	-@echo "Using C/C++ flags: ${PCC_LINKER_FLAGS}"
	-@if [ "${FC}" != "" ]; then \
	   echo "Using Fortran linker: ${FC_LINKER}";\
	   echo "Using Fortran flags: ${FC_LINKER_FLAGS}";\
         fi
	-@echo "-----------------------------------------"
	-@echo "Using libraries: ${SLEPC_LIB}"
	-@echo "------------------------------------------"
	-@echo "Using mpiexec: ${MPIEXEC}"
	-@echo "------------------------------------------"
	-@echo "Using MAKEFLAGS: -j$(MAKE_NP) -l$(MAKE_LOAD) $(MAKEFLAGS)"
	-@echo "=========================================="

# Simple test examples for checking a correct installation
check_install: check
check:
	-+@${OMAKE_SELF} PATH="${PETSC_DIR}/${PETSC_ARCH}/lib:${SLEPC_DIR}/${PETSC_ARCH}/lib:${PATH}" PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR} check_build 2>&1 | tee ./${PETSC_ARCH}/lib/slepc/conf/check.log
check_build:
	-@echo "Running test examples to verify correct installation"
	-@echo "Using SLEPC_DIR=${SLEPC_DIR}, PETSC_DIR=${PETSC_DIR} and PETSC_ARCH=${PETSC_ARCH}"
	+@cd src/eps/tests && \
         ${OMAKE_SELF} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR} testtest10 && \
	 egrep "^#define PETSC_HAVE_FORTRAN 1" ${PETSCCONF_H} | tee .ftn.log > /dev/null; \
         if test -s .ftn.log; then \
           ${OMAKE_SELF} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR} testtest7f; \
         fi ; ${RM} .ftn.log && \
	 if [ "${BLOPEX_LIB}" != "" ]; then \
           ${OMAKE_SELF} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR} testtest5_blopex; \
         fi
	-@echo "Completed test examples"

# Compare ABI/API of two versions of PETSc library with the old one defined by PETSC_{DIR,ARCH}_ABI_OLD
abitest:
	@if [ "x${SLEPC_DIR_ABI_OLD}" = "x" ] || [ "x${PETSC_ARCH_ABI_OLD}" = "x" ] || [ "x${PETSC_DIR_ABI_OLD}" = "x" ]; \
		then printf "You must set environment variables SLEPC_DIR_ABI_OLD, PETSC_ARCH_ABI_OLD and PETSC_DIR_ABI_OLD to run abitest\n"; \
		exit 1; \
	fi;
	-@echo "Comparing ABI/API of the following two SLEPc versions (you must have already configured and built them using GCC and with -g):"
	-@echo "========================================================================================="
	-@echo "    Old: SLEPC_DIR_ABI_OLD  = ${SLEPC_DIR_ABI_OLD}"
	-@echo "         PETSC_ARCH_ABI_OLD = ${PETSC_ARCH_ABI_OLD}"
	-@echo "         PETSC_DIR_ABI_OLD  = ${PETSC_DIR_ABI_OLD}"
	-@cd ${SLEPC_DIR_ABI_OLD}; echo "         Branch             = "`git rev-parse --abbrev-ref HEAD`
	-@echo "    New: SLEPC_DIR          = ${SLEPC_DIR}"
	-@echo "         PETSC_ARCH         = ${PETSC_ARCH}"
	-@echo "         PETSC_DIR          = ${PETSC_DIR}"
	-@echo "         Branch             = "`git rev-parse --abbrev-ref HEAD`
	-@echo "========================================================================================="
	-@$(PYTHON)	${SLEPC_DIR}/lib/slepc/bin/maint/abicheck.py -old_dir ${SLEPC_DIR_ABI_OLD} -old_arch ${PETSC_ARCH_ABI_OLD} -old_petsc_dir ${PETSC_DIR_ABI_OLD} -new_dir ${SLEPC_DIR} -new_arch ${PETSC_ARCH} -new_petsc_dir ${PETSC_DIR} -report_format html

# Deletes SLEPc library
deletelibs:
	-${RM} -r ${SLEPC_LIB_DIR}/libslepc*.*
deletemods:
	-${RM} -f ${SLEPC_DIR}/${PETSC_ARCH}/include/slepc*.mod

allclean:
	-@${OMAKE} -f gmakefile clean

clean:: allclean

reconfigure:
	@unset MAKEFLAGS && ${PYTHON} ${PETSC_ARCH}/lib/slepc/conf/reconfigure-${PETSC_ARCH}.py

#
# Install relevant files in the prefix directory
#
install:
	@${PYTHON} ./config/install.py ${SLEPC_DIR} ${PETSC_DIR} ${SLEPC_INSTALLDIR} -destDir=${DESTDIR} ${PETSC_ARCH} ${AR_LIB_SUFFIX} ${RANLIB};

# ------------------------------------------------------------------
#
# All remaining actions are intended for SLEPc developers only.
# SLEPc users should not generally need to use these commands.
#

# Builds all the documentation
alldoc: allcite allpdf alldoc1 alldoc2 docsetdate

# Build just citations
allcite: chk_loc deletemanualpages
	-${OMAKE_SELF} ACTION=manualpages_buildcite tree_basic LOC=${LOC}
	-@sed -e s%man+../%man+manualpages/% ${LOC}/docs/manualpages/manualpages.cit > ${LOC}/docs/manualpages/htmlmap
	-@cat ${PETSC_DIR}/src/docs/mpi.www.index >> ${LOC}/docs/manualpages/htmlmap

# Build just PDF manual + prerequisites
allpdf:
	-cd docs/manual; ${OMAKE_SELF} slepc.pdf clean; mv slepc.pdf ../../docs

# Build just manual pages + prerequisites
allmanpages: chk_loc allcite
	-${OMAKE_SELF} ACTION=slepc_manualpages tree_basic LOC=${LOC}

# Build just manual examples + prerequisites
allmanexamples: chk_loc allmanpages
	-${OMAKE_SELF} ACTION=slepc_manexamples tree_basic LOC=${LOC}

# Build everything that goes into 'doc' dir except html sources
alldoc1: chk_loc allcite allmanpages allmanexamples
	-${PYTHON} ${PETSC_DIR}/lib/petsc/bin/maint/wwwindex.py ${SLEPC_DIR} ${LOC}
	-@echo "<html>" > singleindex.html
	-@echo "<head>" >> singleindex.html
	-@echo "  <title>Subroutine Index</title>" >> singleindex.html
	-@echo "  <meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\">" >> singleindex.html
	-@echo "  <link rel=\"stylesheet\" href=\"/slepc.css\" type=\"text/css\">" >> singleindex.html
	-@echo "</head>" >> singleindex.html
	-@echo "<body>" >> singleindex.html
	-@cat ${LOC}/docs/manualpages/singleindex.html >> singleindex.html
	-@sed -e 's/CC3333/883300/' singleindex.html > ${LOC}/docs/manualpages/singleindex.html
	-@${RM} singleindex.html

# Builds .html versions of the source
alldoc2: chk_loc allcite
	-${OMAKE_SELF} ACTION=slepc_html PETSC_DIR=${PETSC_DIR} alltree LOC=${LOC}
	cp ${LOC}/docs/manual.html ${LOC}/docs/index.html

# modify all generated html files and add in version number, date, canonical URL info.
docsetdate: chk_petscdir
	@echo "Updating generated html files with slepc version, date, canonical URL info";\
        version_release=`grep '^#define SLEPC_VERSION_RELEASE ' include/slepcversion.h |tr -s ' ' | cut -d ' ' -f 3`; \
        version_major=`grep '^#define SLEPC_VERSION_MAJOR ' include/slepcversion.h |tr -s ' ' | cut -d ' ' -f 3`; \
        version_minor=`grep '^#define SLEPC_VERSION_MINOR ' include/slepcversion.h |tr -s ' ' | cut -d ' ' -f 3`; \
        version_subminor=`grep '^#define SLEPC_VERSION_SUBMINOR ' include/slepcversion.h |tr -s ' ' | cut -d ' ' -f 3`; \
        if  [ $${version_release} = 0 ]; then \
          slepcversion=slepc-master; \
          export slepcversion; \
        elif [ $${version_release} = 1 ]; then \
          slepcversion=slepc-$${version_major}.$${version_minor}.$${version_subminor}; \
          export slepcversion; \
        else \
          echo "Unknown SLEPC_VERSION_RELEASE: $${version_release}"; \
          exit; \
        fi; \
        datestr=`git log -1 --pretty=format:%ci | cut -d ' ' -f 1`; \
        export datestr; \
        gitver=`git describe --match "v*"`; \
        export gitver; \
        find include src docs/manualpages -type f -name \*.html \
          -exec perl -pi -e 's^(<body.*>)^$$1\n   <div id=\"version\" align=right><b>$$ENV{slepcversion} $$ENV{datestr}</b></div>\n   <div id="bugreport" align=right><a href="mailto:slepc-maint\@upv.es?subject=Typo or Error in Documentation &body=Please describe the typo or error in the documentation: $$ENV{slepcversion} $$ENV{gitver} {} "><small>Report Typos and Errors</small></a></div>^i' {} \; \
          -exec perl -pi -e 's^(<head>)^$$1 <link rel="canonical" href="https://slepc.upv.es/documentation/current/{}" />^i' {} \; ; \
        echo "Done fixing version number, date, canonical URL info"

# Deletes documentation
alldocclean: deletemanualpages allcleanhtml
deletemanualpages: chk_loc
	-@if [ -d ${LOC} -a -d ${LOC}/docs/manualpages ]; then \
          ${RM} -rf ${LOC}/docs/manualpages ;\
          ${RM} -f ${LOC}/docs/slepc.pdf ;\
        fi
allcleanhtml:
	-${OMAKE_SELF} ACTION=cleanhtml PETSC_DIR=${PETSC_DIR} alltree

# Builds Fortran stub files
allfortranstubs:
	-@${RM} -rf ${PETSC_ARCH}/include/slepc/finclude/ftn-auto/*-tmpdir
	@${PYTHON} ${SLEPC_DIR}/lib/slepc/bin/maint/generatefortranstubs.py ${BFORT} ${VERBOSE}
	-@${PYTHON} ${SLEPC_DIR}/lib/slepc/bin/maint/generatefortranstubs.py -merge ${VERBOSE}
	-@${RM} -rf include/slepc/finclude/ftn-auto/*-tmpdir
deletefortranstubs:
	-@find . -type d -name ftn-auto | xargs rm -rf

check_output:
	-@${OMAKE} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR} -f gmakefile.test check_output

check_ascii:
	@ ! git --no-pager grep -l -I -P "[^\x00-\x7F]"

checkbadSource_slepc:
	-@${OMAKE_SELF} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR} checkbadSource 2>&1
	@ exit `grep -c 'files with errors' checkbadSource.out`

# -------------------------------------------------------------------------------
#
# Some macros to check if the Fortran interface is up-to-date.
#
countfortranfunctions:
	-@for D in `find ${SLEPC_DIR}/src -name ftn-auto` \
	`find ${SLEPC_DIR}/src -name ftn-custom`; do cd $$D; \
	egrep '^void' *.c | \
	cut -d'(' -f1 | tr -s  ' ' | cut -d' ' -f3 | uniq | egrep -v "(^$$|Petsc)" | \
	sed "s/_$$//"; done | sort > /tmp/countfortranfunctions

countcfunctions:
	-@ ls ${SLEPC_DIR}/include/*.h | grep -v slepcblaslapack.h | \
	xargs grep extern | grep "(" | tr -s ' ' | \
	cut -d'(' -f1 | cut -d' ' -f3 | grep -v "\*" | tr -s '\012' |  \
	tr 'A-Z' 'a-z' |  sort > /tmp/countcfunctions

difffortranfunctions: countfortranfunctions countcfunctions
	-@echo -------------- Functions missing in the Fortran interface ---------------------
	-@${DIFF} /tmp/countcfunctions /tmp/countfortranfunctions | grep "^<" | cut -d' ' -f2
	-@echo ----------------- Functions missing in the C interface ------------------------
	-@${DIFF} /tmp/countcfunctions /tmp/countfortranfunctions | grep "^>" | cut -d' ' -f2
	-@${RM}  /tmp/countcfunctions /tmp/countfortranfunctions

checkbadfortranstubs:
	-@echo "========================================="
	-@echo "Functions with MPI_Comm as an Argument"
	-@echo "========================================="
	-@for D in `find ${SLEPC_DIR}/src -name ftn-auto`; do cd $$D; \
	grep '^void' *.c | grep 'MPI_Comm' | \
	tr -s ' ' | tr -s ':' ' ' |cut -d'(' -f1 | cut -d' ' -f1,3; done
	-@echo "========================================="
	-@echo "Functions with a String as an Argument"
	-@echo "========================================="
	-@for D in `find ${SLEPC_DIR}/src -name ftn-auto`; do cd $$D; \
	grep '^void' *.c | grep 'char \*' | \
	tr -s ' ' | tr -s ':' ' ' |cut -d'(' -f1 | cut -d' ' -f1,3; done
	-@echo "========================================="
	-@echo "Functions with Pointers to PETSc Objects as Argument"
	-@echo "========================================="
	-@_p_OBJ=`grep _p_ ${PETSC_DIR}/include/*.h | tr -s ' ' | \
	cut -d' ' -f 3 | tr -s '\012' | grep -v '{' | cut -d'*' -f1 | \
	sed "s/_p_//g" | tr -s '\012 ' ' *|' ` ; \
	_p_OBJS=`grep _p_ ${SLEPC_DIR}/include/*.h | tr -s ' ' | \
	cut -d' ' -f 3 | tr -s '\012' | grep -v '{' | cut -d'*' -f1 | \
	sed "s/_p_//g" | tr -s '\012 ' ' *|' ` ; \
	for D in `find ${SLEPC_DIR}/src -name ftn-auto`; do cd $$D; \
	for OBJ in $$_p_OBJ $$_p_OBJS; do \
	grep "$$OBJ \*" *.c | tr -s ' ' | tr -s ':' ' ' | \
	cut -d'(' -f1 | cut -d' ' -f1,4; \
	done; done

# Generate tags
alletags:
	-@${PYTHON} ${SLEPC_DIR}/lib/slepc/bin/maint/generateetags.py
	-@find config -type f -name "*.py" |grep -v SCCS | xargs etags -o TAGS_PYTHON

.PHONY: info all deletelibs allclean alletags alldoc allcleanhtml countfortranfunctions install

