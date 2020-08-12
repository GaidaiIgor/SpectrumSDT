import config.package
import os

class Configure(config.package.CMakePackage):
  def __init__(self, framework):
    config.package.CMakePackage.__init__(self, framework)
    self.minversion       = '6.1.1'
    self.version          = '6.3.0'
    self.versionname      = 'SUPERLU_DIST_MAJOR_VERSION.SUPERLU_DIST_MINOR_VERSION.SUPERLU_DIST_PATCH_VERSION'
    self.gitcommit        = 'v'+self.version
    self.download         = ['git://https://github.com/xiaoyeli/superlu_dist','https://github.com/xiaoyeli/superlu_dist/archive/'+self.gitcommit+'.tar.gz']
    self.downloaddirnames = ['SuperLU_DIST','superlu_dist']
    self.functions        = ['set_default_options_dist']
    self.includes         = ['superlu_ddefs.h']
    self.liblist          = [['libsuperlu_dist.a']]
    # SuperLU_Dist does not work with --download-fblaslapack with Compaqf90 compiler on windows.
    # However it should work with intel ifort.
    self.downloadonWindows= 1
    self.hastests         = 1
    self.hastestsdatafiles= 1
    self.requirec99flag   = 1 # SuperLU_Dist uses C99 features
    self.precisions       = ['double']
    self.cxx              = 1
    self.requirescxx11    = 1
    return

  def setupDependencies(self, framework):
    config.package.CMakePackage.setupDependencies(self, framework)
    self.blasLapack     = framework.require('config.packages.BlasLapack',self)
    self.parmetis       = framework.require('config.packages.parmetis',self)
    self.mpi            = framework.require('config.packages.MPI',self)
    self.cuda           = framework.require('config.packages.cuda',self)
    self.openmp         = framework.require('config.packages.openmp',self)
    self.odeps          = [self.parmetis,self.cuda,self.openmp]
    self.deps           = [self.mpi,self.blasLapack]
    return

  def formCMakeConfigureArgs(self):
    args = config.package.CMakePackage.formCMakeConfigureArgs(self)
    if self.openmp.found:
      self.usesopenmp = 'yes'
    else:
      args.append('-DCMAKE_DISABLE_FIND_PACKAGE_OpenMP=TRUE')
    if self.cuda.found:
      if not self.openmp.found:
        raise RuntimeError('SuperLU_DIST GPU code currently requires OpenMP. Use --with-openmp=1')
      # SuperLU_DIST CMake doesn't know about GPU builds
      for place,item in enumerate(args):
        if item.find('CMAKE_C_FLAGS') >= 0:
          args[place]=item[:-1]+' -DGPU_ACC '+self.headers.toString(self.cuda.include)+'"'
    args.append('-DUSE_XSDK_DEFAULTS=YES')
    args.append('-DTPL_BLAS_LIBRARIES="'+self.libraries.toString(self.blasLapack.dlib)+'"')
    args.append('-DTPL_LAPACK_LIBRARIES="'+self.libraries.toString(self.blasLapack.dlib)+'"')
    if self.parmetis.found:
      args.append('-DTPL_PARMETIS_INCLUDE_DIRS="'+';'.join(self.parmetis.dinclude)+'"')
      args.append('-DTPL_PARMETIS_LIBRARIES="'+self.libraries.toString(self.parmetis.dlib)+'"')
    else:
      args.append('-Denable_parmetislib=FALSE')
      args.append('-DTPL_ENABLE_PARMETISLIB=FALSE')

    if self.getDefaultIndexSize() == 64:
      args.append('-DXSDK_INDEX_SIZE=64')

    if not hasattr(self.compilers, 'FC'):
      args.append('-DXSDK_ENABLE_Fortran=OFF')

    args.append('-Denable_tests=0')
    args.append('-Denable_examples=0')
    #  CMake in SuperLU should set this; but like many other packages it does not
    args.append('-DCMAKE_INSTALL_NAME_DIR:STRING="'+os.path.join(self.installDir,self.libdir)+'"')
    args.append('-DMPI_C_COMPILER:STRING="'+self.framework.getCompiler()+'"')
    args.append('-DMPI_C_COMPILE_FLAGS:STRING=""')
    args.append('-DMPI_C_INCLUDE_PATH:STRING=""')
    args.append('-DMPI_C_HEADER_DIR:STRING=""')
    args.append('-DMPI_C_LIBRARIES:STRING=""')
    args.append('-DCMAKE_INSTALL_LIBDIR:STRING="'+os.path.join(self.installDir,self.libdir)+'"')

    # Add in fortran mangling flag
    if self.blasLapack.mangling == 'underscore':
      mangledef = '-DAdd_'
    elif self.blasLapack.mangling == 'caps':
      mangledef = '-DUpCase'
    else:
      mangledef = '-DNoChange'
    for place,item in enumerate(args):
      if item.find('CMAKE_C_FLAGS') >= 0 or item.find('CMAKE_CXX_FLAGS') >= 0:
        args[place]=item[:-1]+' '+mangledef+'"'

    return args



