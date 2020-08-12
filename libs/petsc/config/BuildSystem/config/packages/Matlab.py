from __future__ import generators
import config.package

import os

class Configure(config.package.Package):
  def __init__(self, framework):
    config.package.Package.__init__(self, framework)
    self.hastests       = 1
    self.executablename = 'matlab'
    return

  def setupHelp(self, help):
    import nargs
    help.addArgument('MATLAB', '-with-matlab=<bool>',         nargs.ArgBool(None, 0, 'Activate Matlab'))
    help.addArgument('MATLAB', '-with-matlab-socket=<bool>',  nargs.ArgBool(None, 1, 'Build socket code for Matlab'))
    help.addArgument('MATLAB', '-with-matlab-dir=<root dir>', nargs.ArgDir(None, None, 'Specify the root directory of the Matlab installation'))
    help.addArgument('MATLAB', '-with-matlab-arch=<string>',  nargs.ArgString(None, None, 'Use Matlab Architecture (default use first-found)'))
    return

  def generateGuesses(self):
    '''Generate list of possible locations of Matlab'''
    if 'with-matlab-dir' in self.argDB:
      yield self.argDB['with-matlab-dir']
      raise RuntimeError('You set a value for --with-matlab-dir, but '+self.argDB['with-matlab-dir']+' cannot be used\n')
    if self.getExecutable('matlab', getFullPath = 1):
      # follow any symbolic link of this path
      self.matlab = os.path.realpath(self.matlab)
      yield os.path.dirname(os.path.dirname(self.matlab))
    if os.path.isdir('/Applications'):
      for dir in os.listdir('/Applications'):
        if dir.startswith('MATLAB'):
          if os.path.isfile(os.path.join('/Applications',dir,'bin','matlab')):
            yield os.path.join('/Applications',dir)
    return

  def configureLibrary(self):
    '''Find a Matlab installation and check if it can work with PETSc'''
    import re

    versionPattern = re.compile('Version ([0-9]*.[0-9]*)')
    for matlab in self.generateGuesses():
      self.log.write('Testing Matlab at '+matlab+'\n')
      interpreter = os.path.join(matlab,'bin','matlab')
      if 'with-matlab-arch' in self.argDB:
        interpreter = interpreter+' -'+self.argDB['with-matlab-arch']

      output      = ''
      try:
        output,err,ret = config.package.Package.executeShellCommand(interpreter+' -nodisplay -r "display([\'Version \' version]); exit"', log = self.log)
      except  RuntimeError as e:
        self.log.write('WARNING: Found Matlab at '+matlab+' but unable to run'+str(e)+'\n')
        continue

      match  = versionPattern.search(output)
      r = float(match.group(1))
      if r < 6.0:
        self.log.write('WARNING: Matlab version must be at least 6; yours is '+str(r))
        continue
      # make sure this is true root of Matlab
      if not os.path.isdir(os.path.join(matlab,'extern','lib')):
        self.log.write('WARNING:'+matlab+' is not the root directory for Matlab\n')
        self.log.write('        Run with --with-matlab-dir=Matlabrootdir if you know where it is\n')
      else:
        self.matlab      = matlab
        ls = os.listdir(os.path.join(matlab,'extern','lib'))
        if ls:
          if 'with-matlab-arch' in self.argDB:
            self.matlab_arch = self.argDB['with-matlab-arch']
            if not self.matlab_arch in ls:
              self.log.write('WARNING: You indicated --with-matlab-arch='+self.matlab_arch+' but that arch does not exist;\n possibilities are '+str(ls))
              continue
          else:
            self.matlab_arch = ls[0]
          self.log.write('Configuring PETSc to use the Matlab at '+matlab+' Matlab arch '+self.matlab_arch+'\n')
          self.mex = os.path.join(matlab,'bin','mex')
          if 'with-matlab-arch' in self.argDB:
            self.mex = self.mex+' -'+self.argDB['with-matlab-arch']

          self.command = os.path.join(matlab,'bin','matlab -'+self.matlab_arch)
          self.include = [os.path.join(matlab,'extern','include')]
          self.framework.packages.append(self)
          self.addMakeMacro('MATLAB_MEX',self.mex)
          self.addMakeMacro('MATLAB_COMMAND',self.command)
          self.addDefine('MATLAB_COMMAND','"'+self.command+'"')
          self.found = 1
          if not 'with-matlab-socket' in self.argDB or self.argDB['with-matlab-socket']:
            self.addDefine('USE_MATLAB_SOCKET','1')
            self.addMakeMacro('MATLAB_SOCKET','yes')
          return
        else:
          self.log.write('WARNING:Unable to use Matlab because cannot locate Matlab external libraries at '+os.path.join(matlab,'extern','lib')+'\n')
    raise RuntimeError('Could not find a functional Matlab\nRun with --with-matlab-dir=Matlabrootdir if you know where it is\n')
    return
