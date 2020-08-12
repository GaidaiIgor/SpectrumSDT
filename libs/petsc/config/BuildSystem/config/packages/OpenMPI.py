import config.package
import os

class Configure(config.package.GNUPackage):
  def __init__(self, framework):
    config.package.GNUPackage.__init__(self, framework)
    self.download               = ['https://download.open-mpi.org/release/open-mpi/v3.1/openmpi-3.1.4.tar.gz',
                                   'http://ftp.mcs.anl.gov/pub/petsc/externalpackages/openmpi-3.1.4.tar.gz']
    self.downloaddirnames       = ['openmpi']
    self.skippackagewithoptions = 1
    self.isMPI                  = 1
    return

  def formGNUConfigureArgs(self):
    args = config.package.GNUPackage.formGNUConfigureArgs(self)
    args.append('--with-rsh=ssh')
    args.append('MAKE='+self.make.make)
    if not hasattr(self.compilers, 'CXX'):
      raise RuntimeError('Error: OpenMPI requires C++ compiler. None specified')
    if hasattr(self.compilers, 'FC'):
      self.pushLanguage('FC')
      if not self.fortran.fortranIsF90:
        args.append('--disable-mpi-f90')
        args.append('FC=""')
      self.popLanguage()
    else:
      args.append('--disable-mpi-f77')
      args.append('--disable-mpi-f90')
      args.append('F77=""')
      args.append('FC=""')
      args.append('--enable-mpi-fortran=no')
    if not self.argDB['with-shared-libraries']:
      args.append('--enable-shared=no')
      args.append('--enable-static=yes')
    args.append('--disable-vt')
    # have OpenMPI build its own private copy of hwloc to prevent possible conflict with one used by PETSc
    args.append('--with-hwloc=internal')
    return args

  def checkDownload(self):
    if self.argDB['download-'+self.downloadname.lower()] and  'package-prefix-hash' in self.argDB and self.argDB['package-prefix-hash'] == 'reuse':
      self.logWrite('Reusing package prefix install of '+self.defaultInstallDir+' for OpenMPI')
      self.installDir = self.defaultInstallDir
      self.updateCompilers(self.installDir,'mpicc','mpicxx','mpif77','mpif90')
      return self.installDir
    if self.argDB['download-'+self.downloadname.lower()]:
      return self.getInstallDir()
    return ''

  def Install(self):
    '''After downloading and installing OpenMPI we need to reset the compilers to use those defined by the OpenMPI install'''
    if 'package-prefix-hash' in self.argDB and self.argDB['package-prefix-hash'] == 'reuse':
      return self.defaultInstallDir
    installDir = config.package.GNUPackage.Install(self)
    self.updateCompilers(installDir,'mpicc','mpicxx','mpif77','mpif90')
    return installDir

