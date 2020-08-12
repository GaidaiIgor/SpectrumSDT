#
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  SLEPc - Scalable Library for Eigenvalue Problem Computations
#  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
#
#  This file is part of SLEPc.
#  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#

import os, shutil, log, package

class Primme(package.Package):

  def __init__(self,argdb,log):
    package.Package.__init__(self,argdb,log)
    self.packagename    = 'primme'
    self.installable    = True
    self.downloadable   = True
    self.version        = '3.1.1'
    self.url            = 'https://github.com/primme/primme/archive/v'+self.version+'.tar.gz'
    self.archive        = 'primme-'+self.version+'.tar.gz'
    self.dirname        = 'primme-'+self.version
    self.supportssingle = True
    self.supports64bint = True
    self.hasheaders     = True
    self.hasdloadflags  = True
    self.ProcessArgs(argdb)

  def SampleCode(self,petsc):
    if petsc.scalar == 'real':
      if petsc.precision == 'single':
        function = 'sprimme'
        rdtype = 'float'
      else:
        function = 'dprimme'
        rdtype = 'double'
      cdtype = rdtype
    else:
      if petsc.precision == 'single':
        function = 'cprimme'
        rdtype = 'float'
        cdtype = 'PRIMME_COMPLEX_FLOAT'
      else:
        function = 'zprimme'
        rdtype = 'double'
        cdtype = 'PRIMME_COMPLEX_DOUBLE'

    code = '#include "primme.h"\n'
    code += 'int main() {\n'
    code += '  ' + rdtype + ' *a=NULL,*c=NULL;\n'
    code += '  ' + cdtype + ' *b=NULL;\n'
    code += '  primme_params primme;\n'
    code += '  primme_initialize(&primme);\n'
    code += '  primme_set_method(PRIMME_DYNAMIC,&primme);\n'
    code += '  ' + function + '(a,b,c,&primme);\n'
    code += '  primme_free(&primme);\n'
    code += '  return 0;\n}\n'
    return code


  def Check(self,conf,vars,petsc,archdir):
    code = self.SampleCode(petsc)
    if self.packagedir:
      dirs = [os.path.join(self.packagedir,'lib'),self.packagedir]
      incdirs = [os.path.join(self.packagedir,'include'),self.packagedir]
    else:
      dirs = self.GenerateGuesses('Primme',archdir)
      incdirs = self.GenerateGuesses('Primme',archdir,'include')

    libs = self.packagelibs
    if not libs:
      libs = ['-lprimme']
    includes = self.packageincludes
    if not includes:
      includes = ['.']

    for (d,i) in zip(dirs,incdirs):
      if d:
        if petsc.buildsharedlib:
          l = [petsc.slflag + d] + ['-L' + d] + libs
        else:
          l = ['-L' + d] + libs
        f = ['-I' + i]
      else:
        l = libs
        f = ['-I' + includes[0]]
      result = self.Link([],[],l+f,code,' '.join(f),petsc.language)
      if result:
        conf.write('#define SLEPC_HAVE_PRIMME 1\n')
        vars.write('PRIMME_LIB = ' + ' '.join(l) + '\n')
        vars.write('PRIMME_INCLUDE = ' + ' '.join(f) + '\n')
        self.havepackage = True
        self.packageflags = l+f
        self.location = includes[0] if self.packageincludes else i
        return

    self.log.Exit('Unable to link with PRIMME library in directories'+' '.join(dirs)+' with libraries and link flags '+' '.join(libs)+' [NOTE: make sure PRIMME version is 2.0 at least]')


  def DownloadAndInstall(self,conf,vars,slepc,petsc,archdir,prefixdir):
    externdir = slepc.CreateDir(archdir,'externalpackages')
    builddir  = self.Download(externdir,slepc.downloaddir)

    # Configure
    g = open(os.path.join(builddir,'mymake_flags'),'w')
    g.write('export LIBRARY     = libprimme.'+petsc.ar_lib_suffix+'\n')
    g.write('export SOLIBRARY   = libprimme.'+petsc.sl_suffix+'\n')
    g.write('export SONAMELIBRARY = libprimme.'+petsc.sl_suffix+self.version+'\n')
    g.write('export CC          = '+petsc.cc+'\n')
    if hasattr(petsc,'fc'):
      g.write('export F77         = '+petsc.fc+'\n')
    g.write('export DEFINES     = ')
    if petsc.blaslapackmangling == 'underscore':
      g.write('-DF77UNDERSCORE ')
    if petsc.blaslapackint64:
      g.write('-DPRIMME_BLASINT_SIZE=64')
    g.write('\n')
    g.write('export INCLUDE     = \n')
    g.write('export CFLAGS      = '+petsc.cc_flags.replace('-Wall','').replace('-Wshadow','').replace('-fvisibility=hidden','')+' '+self.buildflags+'\n')
    g.write('export RANLIB      = '+petsc.ranlib+'\n')
    g.write('export PREFIX      = '+prefixdir+'\n')
    g.write('include makefile\n')
    g.close()

    # Build package
    target = ' install' if petsc.buildsharedlib else ' lib'
    mymake = petsc.make + ' -f mymake_flags '
    (result, output) = self.RunCommand('cd '+builddir+'&&'+mymake+' clean && '+mymake+target)
    if result:
      self.log.Exit('Installation of PRIMME failed')

    # Move files
    incdir,libdir = slepc.CreatePrefixDirs(prefixdir)
    if not petsc.buildsharedlib:
      os.rename(os.path.join(builddir,'lib','libprimme.'+petsc.ar_lib_suffix),os.path.join(libdir,'libprimme.'+petsc.ar_lib_suffix))
      for root, dirs, files in os.walk(os.path.join(builddir,'include')):
        for name in files:
          shutil.copyfile(os.path.join(builddir,'include',name),os.path.join(incdir,name))

    if petsc.buildsharedlib:
      l = petsc.slflag + libdir + ' -L' + libdir + ' -lprimme'
    else:
      l = '-L' + libdir + ' -lprimme'
    f = '-I' + incdir

    # Check build
    code = self.SampleCode(petsc)
    (result, output) = self.Link([],[],[l]+[f],code,f,petsc.language)
    if not result:
      self.log.Exit('Unable to link with downloaded PRIMME')

    # Write configuration files
    conf.write('#define SLEPC_HAVE_PRIMME 1\n')
    vars.write('PRIMME_LIB = ' + l + '\n')
    vars.write('PRIMME_INCLUDE = ' + f + '\n')

    self.location = incdir
    self.havepackage = True
    self.packageflags = [l] + [f]


  def LoadVersion(self,conf):
    try:
      f = open(os.path.join(self.location,'primme.h'))
      for l in f.readlines():
        l = l.split()
        if len(l) == 3:
          if l[1] == 'PRIMME_VERSION_MAJOR':
            major = l[2]
          elif l[1] == 'PRIMME_VERSION_MINOR':
            minor = l[2]
      f.close()
      self.iversion = major + '.' + minor
      if major=='3':
        conf.write('#define SLEPC_HAVE_PRIMME3 1\n')
    except: pass

