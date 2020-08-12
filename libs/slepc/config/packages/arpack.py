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

class Arpack(package.Package):

  def __init__(self,argdb,log):
    package.Package.__init__(self,argdb,log)
    self.packagename    = 'arpack'
    self.installable    = True
    self.downloadable   = True
    self.version        = '3.7.0'
    self.url            = 'https://github.com/opencollab/arpack-ng/archive/'+self.version+'.tar.gz'
    self.archive        = 'arpack-ng-'+self.version+'.tar.gz'
    self.dirname        = 'arpack-ng-'+self.version
    self.supportssingle = True
    self.supports64bint = True
    self.fortran        = True
    self.ProcessArgs(argdb)

  def Functions(self,petsc):
    if petsc.mpiuni:
      if petsc.scalar == 'real':
        if petsc.precision == 'single':
          functions = ['snaupd','sneupd','ssaupd','sseupd']
        else:
          functions = ['dnaupd','dneupd','dsaupd','dseupd']
      else:
        if petsc.precision == 'single':
          functions = ['cnaupd','cneupd']
        else:
          functions = ['znaupd','zneupd']
    else:
      if petsc.scalar == 'real':
        if petsc.precision == 'single':
          functions = ['psnaupd','psneupd','pssaupd','psseupd']
        else:
          functions = ['pdnaupd','pdneupd','pdsaupd','pdseupd']
      else:
        if petsc.precision == 'single':
          functions = ['pcnaupd','pcneupd']
        else:
          functions = ['pznaupd','pzneupd']
    return functions


  def Check(self,conf,vars,petsc,archdir):
    functions = self.Functions(petsc)
    if self.packagelibs:
      libs = [self.packagelibs]
    else:
      if petsc.mpiuni:
        libs = [['-larpack'],['-larpack_LINUX'],['-larpack_SUN4']]
      else:
        libs = [['-lparpack','-larpack'],['-lparpack_MPI','-larpack'],['-lparpack_MPI-LINUX','-larpack_LINUX'],['-lparpack_MPI-SUN4','-larpack_SUN4']]

    if self.packagedir:
      dirs = [os.path.join(self.packagedir,'lib'),self.packagedir]
    else:
      dirs = self.GenerateGuesses('Arpack',archdir)
    self.FortranLib(conf,vars,dirs,libs,functions)


  def DownloadAndInstall(self,conf,vars,slepc,petsc,archdir,prefixdir):
    externdir = slepc.CreateDir(archdir,'externalpackages')
    builddir  = self.Download(externdir,slepc.downloaddir)

    # Check for autoreconf
    (result,output) = self.RunCommand('autoreconf --help')
    if result:
      self.log.Exit('--download-arpack requires that the command autoreconf is available on your PATH')

    # Build package
    confopt = '--prefix='+prefixdir+' CC="'+petsc.cc+'" CFLAGS="'+petsc.cc_flags+'" F77="'+petsc.fc+'" FFLAGS="'+petsc.fc_flags.replace('-Wall','').replace('-Wshadow','')+'" FC="'+petsc.fc+'" FCFLAGS="'+petsc.fc_flags.replace('-Wall','').replace('-Wshadow','')+'"'
    if not petsc.mpiuni:
      confopt = confopt+' --enable-mpi MPICC="'+petsc.cc+'" MPIF77="'+petsc.fc+'" MPIFC="'+petsc.fc+'"'
    if not petsc.buildsharedlib:
      confopt = confopt+' --disable-shared'
    if petsc.ind64:
      if not petsc.blaslapackint64:
        self.log.Exit('To install ARPACK with 64-bit integers you also need a BLAS with 64-bit integers')
      confopt = confopt+' INTERFACE64=1'
    (result,output) = self.RunCommand('cd '+builddir+'&& sh bootstrap && ./configure '+confopt+' && '+petsc.make+' && '+petsc.make+' install')
    if result:
      self.log.Exit('Installation of ARPACK failed')

    # Check build
    functions = self.Functions(petsc)
    if petsc.mpiuni:
      libs = [['-larpack']]
    else:
      libs = [['-lparpack','-larpack']]
    libdir = os.path.join(prefixdir,'lib')
    dirs = [libdir]
    self.FortranLib(conf,vars,dirs,libs,functions)

