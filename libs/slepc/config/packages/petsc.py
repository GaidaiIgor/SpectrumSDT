#
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  SLEPc - Scalable Library for Eigenvalue Problem Computations
#  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
#
#  This file is part of SLEPc.
#  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#

import package, os, sys

class PETSc(package.Package):

  def __init__(self,argdb,log):
    package.Package.__init__(self,argdb,log)
    self.packagename = 'petsc'

  def Check(self):
    (result, output) = self.Link([],[],[])
    self.havepackage = result

  def InitDir(self,prefixdir):
    if 'PETSC_DIR' in os.environ:
      self.dir = os.environ['PETSC_DIR']
      if not os.path.exists(self.dir):
        self.log.Exit('PETSC_DIR enviroment variable is not valid')
    else:
      if prefixdir:
        self.dir = prefixdir
        os.environ['PETSC_DIR'] = self.dir
      else:
        self.log.Exit('PETSC_DIR enviroment variable is not set')

  def LoadVersion(self):
    try:
      f = open(os.path.join(self.dir,'include','petscversion.h'))
      for l in f.readlines():
        l = l.split()
        if len(l) == 3:
          if l[1] == 'PETSC_VERSION_RELEASE':
            self.release = l[2]
          if l[1] == 'PETSC_VERSION_MAJOR':
            major = l[2]
          elif l[1] == 'PETSC_VERSION_MINOR':
            minor = l[2]
          elif l[1] == 'PETSC_VERSION_SUBMINOR':
            subminor = l[2]
      f.close()
      self.version = major + '.' + minor
      self.lversion = major + '.' + minor + '.' + subminor
      self.nversion = int(major)*100 + int(minor)
    except:
      self.log.Exit('File error while reading PETSc version')

    # Check whether this is a working copy of the repository
    self.isrepo = False
    if os.path.exists(os.path.join(self.dir,'.git')):
      (status, output) = self.RunCommand('cd '+self.dir+';git rev-parse')
      if not status:
        self.isrepo = True
        (status, self.gitrev)  = self.RunCommand('cd '+self.dir+';git log -1 --pretty=format:%H')
        (status, self.gitdate) = self.RunCommand('cd '+self.dir+';git log -1 --pretty=format:%ci')
        (status, self.branch)  = self.RunCommand('cd '+self.dir+';git describe --contains --all HEAD')

  def LoadConf(self):
    if 'PETSC_ARCH' in os.environ and os.environ['PETSC_ARCH']:
      self.isinstall = False
      self.arch = os.environ['PETSC_ARCH']
      if os.path.basename(self.arch) != self.arch:
        suggest = os.path.basename(self.arch)
        if not suggest: suggest = os.path.basename(self.arch[0:-1])
        self.log.Exit('Variable PETSC_ARCH must not be a full path\nYou set PETSC_ARCH=%s, maybe you meant PETSC_ARCH=%s'% (self.arch,suggest))
      petscvariables = os.path.join(self.dir,self.arch,'lib','petsc','conf','petscvariables')
      petscconf_h = os.path.join(self.dir,self.arch,'include','petscconf.h')
    else:
      self.isinstall = True
      petscvariables = os.path.join(self.dir,'lib','petsc','conf','petscvariables')
      petscconf_h = os.path.join(self.dir,'include','petscconf.h')

    self.buildsharedlib = False
    self.bfort = 'nobfortinpetsc'
    try:
      f = open(petscvariables)
      for l in f.readlines():
        r = l.split('=',1)
        if len(r)!=2: continue
        k = r[0].strip()
        v = r[1].strip()
        v = v.replace('${PETSC_DIR}',self.dir)  # needed in some Cray installations
        if k == 'PETSC_SCALAR':
          self.scalar = v
        elif k == 'PETSC_PRECISION':
          self.precision = v
        elif k == 'MAKE':
          self.make = v
        elif k == 'PREFIXDIR':
          self.prefixdir = v
        elif k == 'BFORT':
          self.bfort = v
        elif k == 'CC':
          self.cc = v
        elif k == 'CC_FLAGS':
          self.cc_flags = v
        elif k == 'CXX':
          self.cxx = v
        elif k == 'CXX_FLAGS':
          self.cxx_flags = v
        elif k == 'FC' and not v=='':
          self.fc = v
        elif k == 'FC_FLAGS':
          self.fc_flags = v
        elif k == 'AR':
          self.ar = v
        elif k == 'AR_FLAGS':
          self.ar_flags = v
        elif k == 'AR_LIB_SUFFIX':
          self.ar_lib_suffix = v
        elif k == 'BUILDSHAREDLIB' and v=='yes':
          self.buildsharedlib = True
        elif k == 'CC_LINKER_SLFLAG':
          self.slflag = v
        elif k == 'SL_LINKER_SUFFIX':
          self.sl_suffix = v
        elif k == 'RANLIB':
          self.ranlib = v
      f.close()
    except:
      self.log.Exit('Cannot process file ' + petscvariables)

    self.ind64 = False
    self.mpiuni = False
    self.debug = False
    self.singlelib = False
    self.blaslapackmangling = ''
    self.blaslapackint64 = False
    self.fortran = False
    self.language = 'c'
    self.cxxdialectcxx11 = False
    self.hpddm = False
    try:
      f = open(petscconf_h)
      for l in f.readlines():
        l = l.split()
        if len(l)==3 and l[0]=='#define' and l[1]=='PETSC_USE_64BIT_INDICES' and l[2]=='1':
          self.ind64 = True
        elif len(l)==3 and l[0]=='#define' and l[1]=='PETSC_HAVE_MPIUNI' and l[2]=='1':
          self.mpiuni = True
        elif len(l)==3 and l[0]=='#define' and l[1]=='PETSC_USE_DEBUG' and l[2]=='1':
          self.debug = True
        elif len(l)==3 and l[0]=='#define' and l[1]=='PETSC_USE_SINGLE_LIBRARY' and l[2]=='1':
          self.singlelib = True
        elif len(l)==3 and l[0]=='#define' and l[1]=='PETSC_BLASLAPACK_UNDERSCORE' and l[2]=='1':
          self.blaslapackmangling = 'underscore'
        elif len(l)==3 and l[0]=='#define' and l[1]=='PETSC_BLASLAPACK_CAPS' and l[2]=='1':
          self.blaslapackmangling = 'caps'
        elif len(l)==3 and l[0]=='#define' and l[1]=='PETSC_HAVE_64BIT_BLAS_INDICES' and l[2]=='1':
          self.blaslapackint64 = True
        elif len(l)==3 and l[0]=='#define' and l[1]=='PETSC_HAVE_FORTRAN' and l[2]=='1':
          self.fortran = True
        elif len(l)==3 and l[0]=='#define' and l[1]=='PETSC_CLANGUAGE_CXX' and l[2]=='1':
          self.language = 'c++'
        elif len(l)==3 and l[0]=='#define' and l[1]=='PETSC_HAVE_CXX_DIALECT_CXX11' and l[2]=='1':
          self.cxxdialectcxx11 = True
        elif len(l)==3 and l[0]=='#define' and l[1]=='PETSC_HAVE_HPDDM' and l[2]=='1':
          self.hpddm = True
        elif self.isinstall and len(l)==3 and l[0]=='#define' and l[1]=='PETSC_ARCH':
          self.arch = l[2].strip('"')
      f.close()
    except:
      if self.isinstall:
        self.log.Exit('Cannot process file ' + petscconf_h + ', maybe you forgot to set PETSC_ARCH')
      else:
        self.log.Exit('Cannot process file ' + petscconf_h)

