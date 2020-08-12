#!/usr/bin/env python

"""
SLEPc: Scalable Library for Eigenvalue Problem Computations
===========================================================

SLEPc is a software library for the solution of large scale sparse
eigenvalue problems on parallel computers. It is an extension of PETSc
and can be used for either standard or generalized eigenproblems, with
real or complex arithmetic. It can also be used for computing a
partial SVD of a large, sparse, rectangular matrix, and to solve
nonlinear eigenvalue problems

.. note::

   To install ``PETSc``, ``SLEPc``, ``petsc4py``, and ``slepc4py``
   (``mpi4py`` is optional but highly recommended) use::

     $ pip install numpy mpi4py
     $ pip install petsc petsc4py
     $ pip install slepc slepc4py

.. tip::

  You can also install the in-development versions with::

    $ pip install Cython numpy mpi4py
    $ pip install --no-deps https://gitlab.com/petsc/petsc/-/archive/master/petsc-master.tar.gz
    $ pip install --no-deps https://bitbucket.org/petsc/petsc4py/get/master.tar.gz
    $ pip install --no-deps https://bitbucket.org/slepc/slepc/get/master.tar.gz
    $ pip install --no-deps https://bitbucket.org/slepc/slepc4py/get/master.tar.gz

"""

import sys, os
from setuptools import setup
from setuptools.command.install import install as _install
from distutils.util import get_platform, split_quoted
from distutils.spawn import find_executable
from distutils import log

init_py = """\
# Author:  SLEPc Team
# Contact: slepc-maint@upv.es

def get_slepc_dir():
    import os
    return os.path.dirname(__file__)

def get_config():
    conf = {}
    conf['SLEPC_DIR'] = get_slepc_dir()
    return conf
"""

metadata = {
    'provides' : ['slepc'],
    'zip_safe' : False,
}

CONFIGURE_OPTIONS = []

def bootstrap():
    from os.path import join, isdir, abspath
    # Set SLEPC_DIR
    SLEPC_DIR  = abspath(os.getcwd())
    os.environ['SLEPC_DIR']  = SLEPC_DIR
    # Check PETSC_DIR/PETSC_ARCH
    PETSC_DIR  = os.environ.get('PETSC_DIR',  "")
    PETSC_ARCH = os.environ.get('PETSC_ARCH', "")
    if not (PETSC_DIR and isdir(PETSC_DIR)):
        PETSC_DIR = None
        try: del os.environ['PETSC_DIR']
        except KeyError: pass
        PETSC_ARCH = None
        try: del os.environ['PETSC_ARCH']
        except KeyError: pass
    elif not (PETSC_ARCH and isdir(join(PETSC_DIR, PETSC_ARCH))):
        PETSC_ARCH = None
        try: del os.environ['PETSC_ARCH']
        except KeyError: pass
    # Generate package __init__.py file
    from distutils.dir_util import mkpath
    pkgdir = os.path.join(SLEPC_DIR, 'pypi')
    pkgfile = os.path.join(pkgdir, '__init__.py')
    if not os.path.exists(pkgdir): mkpath(pkgdir)
    fh = open(pkgfile, 'wt')
    fh.write(init_py)
    fh.close()
    # Configure options
    options = os.environ.get('PETSC_CONFIGURE_OPTIONS', '')
    CONFIGURE_OPTIONS.extend(split_quoted(options))
    #
    if not PETSC_DIR:
        vstr = version().split('.')[:2]
        x, y = int(vstr[0]), int(vstr[1])
        reqs = ">=%s.%s,<%s.%s" % (x, y, x, y+1)
        metadata['install_requires'] = ['petsc'+reqs]

def get_petsc_dir():
    PETSC_DIR = os.environ.get('PETSC_DIR')
    if PETSC_DIR: return PETSC_DIR
    try:
        import petsc
        PETSC_DIR = petsc.get_petsc_dir()
    except ImportError:
        log.warn("PETSC_DIR not specified")
        PETSC_DIR = os.path.join(os.path.sep, 'usr', 'local', 'petsc')
    return PETSC_DIR

def get_petsc_arch():
    PETSC_ARCH = os.environ.get('PETSC_ARCH', "")
    return PETSC_ARCH

def config(prefix, dry_run=False):
    log.info('SLEPc: configure')
    options = [
        '--prefix=' + prefix,
        ]
    options.extend(CONFIGURE_OPTIONS)
    #
    log.info('configure options:')
    for opt in options:
        log.info(' '*4 + opt)
    # Run SLEPc configure
    if dry_run: return
    os.environ['PETSC_DIR'] = get_petsc_dir()
    python = find_executable('python2') or find_executable('python')
    command = [python, './configure', '--prefix='+prefix]
    status = os.system(" ".join(command))
    if status != 0: raise RuntimeError(status)

def build(dry_run=False):
    log.info('SLEPc: build')
    # Run SLEPc build
    if dry_run: return
    PETSC_ARCH = get_petsc_arch()
    if PETSC_ARCH: PETSC_ARCH = 'PETSC_ARCH=' + PETSC_ARCH
    make = find_executable('make')
    command = [make, 'all',
               'PETSC_DIR='+get_petsc_dir(), PETSC_ARCH]
    status = os.system(" ".join(command))
    if status != 0: raise RuntimeError(status)

def install(dest_dir, dry_run=False):
    log.info('SLEPc: install')
    # Run SLEPc install
    if dry_run: return
    PETSC_ARCH = get_petsc_arch()
    if PETSC_ARCH: PETSC_ARCH = 'PETSC_ARCH=' + PETSC_ARCH
    make = find_executable('make')
    command = [make, 'install',
               'PETSC_DIR='+get_petsc_dir(), PETSC_ARCH]
    status = os.system(" ".join(command))
    if status != 0: raise RuntimeError(status)

class context(object):
    def __init__(self):
        self.sys_argv = sys.argv[:]
        self.wdir = os.getcwd()
    def enter(self):
        del sys.argv[1:]
        pdir = os.environ['SLEPC_DIR']
        os.chdir(pdir)
        return self
    def exit(self):
        sys.argv[:] = self.sys_argv
        os.chdir(self.wdir)

class cmd_install(_install):

    def initialize_options(self):
        _install.initialize_options(self)
        self.optimize = 1

    def finalize_options(self):
        _install.finalize_options(self)
        self.install_lib = self.install_platlib
        self.install_libbase = self.install_lib

    def run(self):
        root_dir = os.path.abspath(self.install_lib)
        dest_dir = prefix = os.path.join(root_dir, 'slepc')
        #
        #
        ctx = context().enter()
        try:
            config(prefix, self.dry_run)
            build(self.dry_run)
            install(dest_dir, self.dry_run)
        finally:
            ctx.exit()
        #
        self.outputs = []
        for dirpath, _, filenames in os.walk(dest_dir):
            for fn in filenames:
                self.outputs.append(os.path.join(dirpath, fn))
        #
        _install.run(self)

    def get_outputs(self):
        outputs = getattr(self, 'outputs', [])
        outputs += _install.get_outputs(self)
        return outputs

def version():
    import re
    version_re = {
        'major'  : re.compile(r"#define\s+SLEPC_VERSION_MAJOR\s+(\d+)"),
        'minor'  : re.compile(r"#define\s+SLEPC_VERSION_MINOR\s+(\d+)"),
        'micro'  : re.compile(r"#define\s+SLEPC_VERSION_SUBMINOR\s+(\d+)"),
        'release': re.compile(r"#define\s+SLEPC_VERSION_RELEASE\s+(\d+)"),
        }
    slepcversion_h = os.path.join('include','slepcversion.h')
    data = open(slepcversion_h, 'r').read()
    major = int(version_re['major'].search(data).groups()[0])
    minor = int(version_re['minor'].search(data).groups()[0])
    micro = int(version_re['micro'].search(data).groups()[0])
    release = int(version_re['release'].search(data).groups()[0])
    if release:
        v = "%d.%d.%d" % (major, minor, micro)
    else:
        v = "%d.%d.0.dev%d" % (major, minor+1, 0)
    return v

def tarball():
    VERSION = version()
    if '.dev' in VERSION: return None
    return ('https://slepc.upv.es/download/distrib/'
            'slepc-%s.tar.gz#egg=slepc-%s' % (VERSION, VERSION))

description = __doc__.split('\n')[1:-1]; del description[1:3]
classifiers = """
Development Status :: 5 - Production/Stable
Intended Audience :: Developers
Intended Audience :: Science/Research
License :: OSI Approved :: BSD License
Operating System :: POSIX
Programming Language :: C
Programming Language :: C++
Programming Language :: Fortran
Programming Language :: Python
Topic :: Scientific/Engineering
Topic :: Software Development :: Libraries
"""

if 'bdist_wheel' in sys.argv:
    sys.stderr.write("slepc: this package cannot be built as a wheel\n")
    sys.exit(1)

bootstrap()
setup(name='slepc',
      version=version(),
      description=description.pop(0),
      long_description='\n'.join(description),
      classifiers= classifiers.split('\n')[1:-1],
      keywords = ['SLEPc','PETSc', 'MPI'],
      platforms=['POSIX'],
      license='BSD',

      url='https://slepc.upv.es/',
      download_url=tarball(),

      author='SLEPc Team',
      author_email='slepc-maint@upv.es',
      maintainer='Lisandro Dalcin',
      maintainer_email='dalcinl@gmail.com',

      packages = ['slepc'],
      package_dir = {'slepc': 'pypi'},
      cmdclass={'install': cmd_install},
      **metadata)
