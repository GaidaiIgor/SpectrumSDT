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

class Blzpack(package.Package):

  def __init__(self,argdb,log):
    package.Package.__init__(self,argdb,log)
    self.packagename    = 'blzpack'
    self.installable    = True
    self.supportsscalar = ['real']
    self.supportssingle = True
    self.ProcessArgs(argdb)

  def Check(self,conf,vars,petsc,archdir):
    if petsc.precision == 'single':
      functions = ['blzdrs']
    else:
      functions = ['blzdrd']

    if self.packagelibs:
      libs = [self.packagelibs]
    else:
      libs = [['-lblzpack']]

    if self.packagedir:
      dirs = [os.path.join(self.packagedir,'lib'),self.packagedir]
    else:
      dirs = self.GenerateGuesses('Blzpack',archdir)

    self.FortranLib(conf,vars,dirs,libs,functions)

