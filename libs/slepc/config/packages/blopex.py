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

class Blopex(package.Package):

  def __init__(self,argdb,log):
    package.Package.__init__(self,argdb,log)
    self.packagename  = 'blopex'
    self.installable  = True
    self.downloadable = True
    self.gitcommit    = '6eba31f0e071f134a6e4be8eccfb8d9d7bdd5ac7'  #master dec-2019
    self.url          = 'https://github.com/lobpcg/blopex/archive/'+self.gitcommit+'.tar.gz'
    self.archive      = 'blopex-'+self.gitcommit+'.tar.gz'
    self.dirname      = 'blopex-'+self.gitcommit
    self.hasheaders   = True
    self.ProcessArgs(argdb)

  def SampleCode(self,petsc):
    if petsc.scalar == 'real':
      function = 'lobpcg_solve_double'
    else:
      function = 'lobpcg_solve_complex'
    code =  '#include <stdlib.h>\n'
    code += '#include <lobpcg.h>\n'
    code += 'int main() {\n'
    code += '  lobpcg_BLASLAPACKFunctions fn;\n'
    code += '  lobpcg_Tolerance           tol;\n'
    code += '  tol.absolute=1e-20; tol.relative=1e-7;\n'
    code += '  fn.dpotrf=NULL; fn.dsygv=NULL;\n'
    code += '  ' + function + '(NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,fn,tol,0,0,NULL,NULL,NULL,0,NULL,NULL,0);\n'
    code += '  return 0;\n}\n'
    return code

  def Check(self,conf,vars,petsc,archdir):
    code = self.SampleCode(petsc)
    if self.packagedir:
      dirs = [os.path.join(self.packagedir,'blopex_abstract','lib')]
      incdirs = [os.path.join(self.packagedir,'blopex_abstract','include')]
    else:
      dirs = self.GenerateGuesses('blopex',archdir)
      incdirs = self.GenerateGuesses('blopex',archdir,'include')

    libs = self.packagelibs
    if not libs:
      libs = ['-lBLOPEX']
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
      (result, output) = self.Link([],[],l+f,code,' '.join(f),petsc.language)
      if result:
        conf.write('#define SLEPC_HAVE_BLOPEX 1\n')
        vars.write('BLOPEX_LIB = ' + ' '.join(l) + '\n')
        vars.write('BLOPEX_INCLUDE = ' + ' '.join(f) + '\n')
        self.havepackage = True
        self.packageflags = l+f
        self.location = includes[0] if self.packageincludes else i
        return

    self.log.Exit('Unable to link with BLOPEX library in directories'+' '.join(dirs)+' with libraries and link flags '+' '.join(libs))

  def DownloadAndInstall(self,conf,vars,slepc,petsc,archdir,prefixdir):
    externdir = slepc.CreateDir(archdir,'externalpackages')
    downloadd = self.Download(externdir,slepc.downloaddir)
    builddir  = os.path.join(downloadd,'blopex_abstract')

    # Configure
    g = open(os.path.join(builddir,'Makefile.inc'),'w')
    g.write('CC          = '+petsc.cc+'\n')
    g.write('CFLAGS      = '+petsc.cc_flags.replace('-Wall','').replace('-Wshadow','')+'\n')
    g.write('AR          = '+petsc.ar+' '+petsc.ar_flags+'\n')
    g.write('AR_LIB_SUFFIX = '+petsc.ar_lib_suffix+'\n')
    g.write('RANLIB      = '+petsc.ranlib+'\n')
    g.write('TARGET_ARCH = \n')
    g.close()

    # Build package
    (result,output) = self.RunCommand('cd '+builddir+'&&'+petsc.make+' clean &&'+petsc.make)
    if result:
      self.log.Exit('Installation of BLOPEX failed')

    # Move files
    incdir,libDir = slepc.CreatePrefixDirs(prefixdir)
    incblopexdir  = slepc.CreateDir(incdir,'blopex')
    os.rename(os.path.join(builddir,'lib','libBLOPEX.'+petsc.ar_lib_suffix),os.path.join(libDir,'libBLOPEX.'+petsc.ar_lib_suffix))
    for root, dirs, files in os.walk(os.path.join(builddir,'include')):
      for name in files:
        shutil.copyfile(os.path.join(builddir,'include',name),os.path.join(incblopexdir,name))

    if petsc.buildsharedlib:
      l = petsc.slflag + libDir + ' -L' + libDir + ' -lBLOPEX'
    else:
      l = '-L' + libDir + ' -lBLOPEX'
    f = '-I' + incdir + ' -I' + incblopexdir

    # Check build
    code = self.SampleCode(petsc)
    (result, output) = self.Link([],[],[l]+[f],code,f,petsc.language)
    if not result:
      self.log.Exit('Unable to link with downloaded BLOPEX')

    # Write configuration files
    conf.write('#define SLEPC_HAVE_BLOPEX 1\n')
    vars.write('BLOPEX_LIB = ' + l + '\n')
    vars.write('BLOPEX_INCLUDE = ' + f + '\n')

    self.havepackage = True
    self.packageflags = [l] + [f]

