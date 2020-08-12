#
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  SLEPc - Scalable Library for Eigenvalue Problem Computations
#  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
#
#  This file is part of SLEPc.
#  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#

import os, log, package

class Slicot(package.Package):

  def __init__(self,argdb,log):
    package.Package.__init__(self,argdb,log)
    self.packagename    = 'slicot'
    self.installable    = True
    self.version        = '4.5'
    self.url            = 'http://slicot.org/objects/software/shared/slicot45.tar.gz'
    self.archive        = 'slicot45.tar.gz'
    self.dirname        = 'slicot'
    self.supportsscalar = ['real']
    self.fortran        = True
    self.ProcessArgs(argdb)


  def Check(self,conf,vars,petsc,archdir):
    functions = ['sb03od','sb03md']
    if self.packagelibs:
      libs = [self.packagelibs]
    else:
      libs = [['-lslicot']]

    if self.packagedir:
      dirs = [os.path.join(self.packagedir,'lib'),self.packagedir]
    else:
      dirs = self.GenerateGuesses('slicot',archdir)

    self.FortranLib(conf,vars,dirs,libs,functions)


  def DownloadAndInstall(self,conf,vars,slepc,petsc,archdir,prefixdir):
    externdir = slepc.CreateDir(archdir,'externalpackages')
    builddir  = self.Download(externdir,slepc.downloaddir)
    libname = 'libslicot.a'

    # Configure
    g = open(os.path.join(builddir,'make.inc'),'w')
    g.write('FORTRAN   = '+petsc.fc+'\n')
    g.write('OPTS      = '+petsc.fc_flags.replace('-Wall','').replace('-Wshadow','')+'\n')
    g.write('ARCH      = '+petsc.ar+'\n')
    g.write('ARCHFLAGS = '+petsc.ar_flags+'\n')
    g.write('SLICOTLIB = ../'+libname+'\n')
    g.close()

    # Build package
    target = 'lib'
    (result,output) = self.RunCommand('cd '+builddir+'&&'+petsc.make+' clean &&'+petsc.make+' '+target)
    if result:
      self.log.Exit('Installation of SLICOT failed')

    # Move files
    incdir,libdir = slepc.CreatePrefixDirs(prefixdir)
    os.rename(os.path.join(builddir,libname),os.path.join(libdir,libname))

    # Check build
    functions = ['sb03od']
    libs = [['-lslicot']]
    dirs = [libdir]
    self.FortranLib(conf,vars,dirs,libs,functions)

