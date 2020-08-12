import config.package
import os

class Configure(config.package.Package):
  def __init__(self, framework):
    config.package.Package.__init__(self, framework)
    self.gitcommit         = 'b09c65a'
    self.download          = ['git://https://github.com/viennacl/viennacl-dev']
    self.downloaddirname   = [str('viennacl-dev')]
    self.includes          = ['viennacl/forwards.h']
    self.cxx               = 1
    self.downloadonWindows = 1
    self.complex           = 0
    return

  def setupDependencies(self, framework):
    config.package.Package.setupDependencies(self, framework)
    self.cuda = framework.require('config.packages.cuda',self)
    self.opencl = framework.require('config.packages.opencl',self)
    self.openmp = framework.require('config.packages.openmp',self)
    self.deps = []
    self.odeps = [self.cuda,self.opencl,self.openmp]
    return

  def Install(self):
    import shutil
    import os
    self.log.write('ViennaCLDir = '+self.packageDir+' installDir '+self.installDir+'\n')
    #includeDir = self.packageDir
    srcdir     = os.path.join(self.packageDir, 'viennacl')
    destdir    = os.path.join(self.installDir, 'include', 'viennacl')
    if self.installSudo:
      self.installDirProvider.printSudoPasswordMessage()
      try:
        output,err,ret  = config.package.Package.executeShellCommand(self.installSudo+'mkdir -p '+destdir+' && '+self.installSudo+'rm -rf '+destdir+'  && '+self.installSudo+'cp -rf '+srcdir+' '+destdir, timeout=6000, log = self.log)
      except RuntimeError as e:
        raise RuntimeError('Error copying ViennaCL include files from '+os.path.join(self.packageDir, 'ViennaCL')+' to '+packageDir)
    else:
      try:
        if os.path.isdir(destdir): shutil.rmtree(destdir)
        shutil.copytree(srcdir,destdir)
      except RuntimeError as e:
        raise RuntimeError('Error installing ViennaCL include files: '+str(e))
    return self.installDir

  def configureLibrary(self):
    config.package.Package.configureLibrary(self)
    #check for CUDA:
    if not self.cuda.found:
      self.addDefine('HAVE_VIENNACL_NO_CUDA', 1)
