#
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  SLEPc - Scalable Library for Eigenvalue Problem Computations
#  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
#
#  This file is part of SLEPc.
#  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#

import log, package

class Lapack(package.Package):

  def __init__(self,argdb,log):
    package.Package.__init__(self,argdb,log)
    self.packagename = 'lapack'

  def ShowInfo(self):
    if hasattr(self,'missing'):
      self.log.Warn('Missing LAPACK functions: %s.\nSome SLEPc functionality will not be available.\nPlease reconfigure and recompile PETSc with a full LAPACK implementation'%(' '.join(self.missing)))

  def Process(self,conf,vars,slepc,petsc,archdir=''):
    self.make = petsc.make
    self.mangling = petsc.blaslapackmangling
    if petsc.buildsharedlib:
      self.slflag = petsc.slflag
    self.log.NewSection('Checking LAPACK library...')
    self.Check(conf,vars,petsc)

  def LinkBlasLapack(self,functions,callbacks,flags,petsc):
    code = ''
    for f in functions:
      if petsc.language == 'c++':
        code += 'extern "C" void ' + f + '();\n'
      else:
        code += 'extern void ' + f + '();\n'
    code += 'int main() {\n'
    for f in functions:
      code += f + '();\n'
    code += 'return 0;\n}\n'
    (result, output) = self.Link(functions,callbacks,flags,code)
    return result

  def Check(self,conf,vars,petsc):

    # LAPACK standard functions
    l = ['laev2','gehrd','lanhs','trexc','trevc','tgexc','tgevc','stedc','hsein','larfg','larf','lascl','trsyl']

    # LAPACK functions with different real and complex names
    if petsc.scalar == 'real':
      l += ['orghr','sytrd','orgtr']
      if petsc.precision == 'single':
        prefix = 's'
      elif petsc.precision == '__float128':
        prefix = 'q'
      else:
        prefix = 'd'
    else:
      l += ['unghr','hetrd','ungtr']
      if petsc.precision == 'single':
        prefix = 'c'
      elif petsc.precision == '__float128':
        prefix = 'w'
      else:
        prefix = 'z'

    # add prefix to LAPACK names
    functions = []
    for i in l:
      functions.append(prefix + i)

    # in this case, the real name represents both versions
    namesubst = {'unghr':'orghr','hetrd':'sytrd','ungtr':'orgtr'}

    # LAPACK functions which are always used in real version
    l = ['stevr','bdsdc','lag2','lasv2','lartg','laed4','lamrg']
    if petsc.precision == 'single':
      prefix = 's'
    elif petsc.precision == '__float128':
      prefix = 'q'
    else:
      prefix = 'd'
    for i in l:
      functions.append(prefix + i)

    # check for all functions at once
    allf = []
    for i in functions:
      if self.mangling == 'underscore':
        f = i + '_'
      elif self.mangling == 'caps':
        f = i.upper()
      else:
        f = i
      allf.append(f)

    self.log.write('=== Checking all LAPACK functions...')
    if self.LinkBlasLapack(allf,[],[],petsc):
      return

    # check functions one by one
    self.missing = []
    for i in functions:
      if self.mangling == 'underscore':
        f = i + '_\n'
      elif self.mangling == 'caps':
        f = i.upper() + '\n'
      else:
        f = i + '\n'

      self.log.write('=== Checking LAPACK '+i+' function...')
      if not self.LinkBlasLapack([f],[],[],petsc):
        self.missing.append(i)
        # some complex functions are represented by their real names
        if i[1:] in namesubst:
          nf = namesubst[i[1:]]
        else:
          nf = i[1:]
        conf.write('#define SLEPC_MISSING_LAPACK_' + nf.upper() + ' 1\n')

