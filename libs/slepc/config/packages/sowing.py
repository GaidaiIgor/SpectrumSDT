#
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  SLEPc - Scalable Library for Eigenvalue Problem Computations
#  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
#
#  This file is part of SLEPc.
#  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#

import os, sys, log, package

class Sowing(package.Package):

  def __init__(self,argdb,log):
    package.Package.__init__(self,argdb,log)
    self.packagename  = 'sowing'
    self.downloadable = True
    self.url          = 'https://bitbucket.org/petsc/pkg-sowing.git'
    self.ProcessArgs(argdb)

  def DownloadAndInstall(self,slepc,petsc,archdir):
    name = self.packagename.upper()
    self.log.NewSection('Installing '+name+'...')

    # Create externalpackages directory
    externdir = slepc.CreateDir(archdir,'externalpackages')

    # Check if source is already available
    builddir = os.path.join(externdir,'pkg-sowing')
    if os.path.exists(builddir):
      self.log.write('Using '+builddir)
    else: # clone Sowing repo
      url = self.packageurl
      if url=='':
        url = self.url
      try:
        (result,output) = self.RunCommand('cd '+externdir+'&& git clone '+url)
      except RuntimeError as e:
        self.log.Exit('Cannot clone '+url+': '+str(e))

    # Configure, build and install package
    (result,output) = self.RunCommand('cd '+builddir+'&& ./configure --prefix='+archdir+'&&'+petsc.make+'&&'+petsc.make+' install')

    self.havepackage = True
    return os.path.join(archdir,'bin','bfort')
