import config.package

class Configure(config.package.Package):
  def __init__(self, framework):
    config.package.Package.__init__(self, framework)
    self.cython = 0
    self.numpy = 0
    return

  def __str__(self):
    return ''

  def setupHelp(self,help):
    config.package.Package.setupHelp(self,help)
    import nargs
    help.addArgument('PYTHON', '-with-python-exec=<executable>',                nargs.Arg(None, None, 'Alternate Python executable to use for mpi4py/petsc4py'))
    return

  def configure(self):
    '''determine python binary to use'''
    if 'with-python-exec' in self.argDB:
      self.getExecutable(self.argDB['with-python-exec'], getFullPath=1, resultName='pyexe', setMakeMacro = 0)
    else:
      import sys
      self.pyexe = sys.executable
    self.addDefine('PYTHON_EXE','"'+self.pyexe+'"')
    try:
      output1,err1,ret1  = config.package.Package.executeShellCommand(self.pyexe + ' -c "import Cython"',timeout=60, log = self.log)
      self.cython = 1
    except: pass
    try:
      output1,err1,ret1  = config.package.Package.executeShellCommand(self.pyexe + ' -c "import numpy"',timeout=60, log = self.log)
      self.numpy = 1
    except: pass
    return
