#!/usr/bin/env python
#
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  SLEPc - Scalable Library for Eigenvalue Problem Computations
#  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
#
#  This file is part of SLEPc.
#  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#

from __future__ import print_function
import os, sys, time, shutil

def AddDefine(conffile,name,value,prefix='SLEPC_'):
  conffile.write('#define '+prefix+name+' "'+value+'"\n')

def WriteModulesFile(modules,version,sdir):
  ''' Write the contents of the Modules file '''
  modules.write('#%Module\n\n')
  modules.write('proc ModulesHelp { } {\n')
  modules.write('    puts stderr "This module sets the path and environment variables for slepc-%s"\n' % version)
  modules.write('    puts stderr "     see https://slepc.upv.es/ for more information"\n')
  modules.write('    puts stderr ""\n}\n')
  modules.write('module-whatis "SLEPc - Scalable Library for Eigenvalue Problem Computations"\n\n')
  modules.write('module load petsc\n')
  modules.write('set slepc_dir %s\n' % sdir)
  modules.write('setenv SLEPC_DIR $slepc_dir\n')

def WritePkgconfigFile(pkgconfig,version,pversion,sdir,isinstall,prefixdir,singlelib):
  ''' Write the contents of the pkg-config file '''
  pkgconfig.write('prefix=%s\n' % prefixdir)
  pkgconfig.write('exec_prefix=${prefix}\n')
  pkgconfig.write('includedir=${prefix}/include\n')
  pkgconfig.write('libdir=${prefix}/lib\n\n')
  pkgconfig.write('Name: SLEPc\n')
  pkgconfig.write('Description: the Scalable Library for Eigenvalue Problem Computations\n')
  pkgconfig.write('Version: %s\n' % version)
  pkgconfig.write('Requires: PETSc >= %s\n' % pversion)
  pkgconfig.write('Cflags: -I${includedir}')
  if not isinstall:
    pkgconfig.write(' -I'+os.path.join(sdir,'include'))
  pkgconfig.write('\nLibs:')
  if singlelib:
    pkgconfig.write(' -L${libdir} -lslepc\n')
  else:
    pkgconfig.write(' -L${libdir} -lslepcnep -lslepcpep -lslepcsvd -lslepceps -lslepcmfn -lslepclme -lslepcsys\n')

def WriteReconfigScript(reconfig,slepcdir,usedargs):
  ''' Write the contents of the reconfigure script '''
  reconfig.write('#!/usr/bin/env python\n\n')
  reconfig.write('import os, sys\n')
  if usedargs:
    reconfig.write('sys.argv.extend(\''+usedargs+'\'.split())\n')
  reconfig.write('execfile(os.path.join(\''+slepcdir+'\',\'config\',\'configure.py\'))\n')

# Use en_US as language so that compiler messages are in English
def fixLang(lang):
  if lang in os.environ and os.environ[lang] != '':
    lv = os.environ[lang]
    enc = ''
    try: lv,enc = lv.split('.')
    except: pass
    if lv not in ['en_US','C']: lv = 'en_US'
    if enc: lv = lv+'.'+enc
    os.environ[lang] = lv

fixLang('LC_LOCAL')
fixLang('LANG')

# Set python path
configdir = os.path.abspath('config')
if not os.path.isdir(configdir):
  sys.exit('ERROR: Run configure from $SLEPC_DIR, not '+os.path.abspath('.'))
sys.path.insert(0,configdir)
sys.path.insert(0,os.path.join(configdir,'packages'))

# Load auxiliary classes
import argdb, log
argdb = argdb.ArgDB(sys.argv)
log   = log.Log()

# Load classes for packages and process command-line options
import slepc, petsc, arpack, blzpack, trlan, primme, blopex, sowing, lapack, slicot, hpddm
slepc   = slepc.SLEPc(argdb,log)
petsc   = petsc.PETSc(argdb,log)
arpack  = arpack.Arpack(argdb,log)
blopex  = blopex.Blopex(argdb,log)
blzpack = blzpack.Blzpack(argdb,log)
primme  = primme.Primme(argdb,log)
trlan   = trlan.Trlan(argdb,log)
sowing  = sowing.Sowing(argdb,log)
lapack  = lapack.Lapack(argdb,log)
slicot  = slicot.Slicot(argdb,log)
hpddm   = hpddm.HPDDM(argdb,log)

externalpackages = [arpack, blopex, blzpack, hpddm, primme, slicot, trlan]
optionspackages  = [slepc, sowing] + externalpackages
checkpackages    = externalpackages + [lapack]

# Print help if requested and check for wrong command-line options
if argdb.PopHelp():
  print('SLEPc Configure Help')
  print('-'*80)
  for pkg in optionspackages:
    pkg.ShowHelp()
  sys.exit(0)
argdb.ErrorPetscOptions()
argdb.ErrorIfNotEmpty()

# Check enviroment and PETSc version
log.Print('Checking environment...')
petsc.InitDir(slepc.prefixdir)
slepc.InitDir()
petsc.LoadVersion()
slepc.LoadVersion()

# Load PETSc configuration
petsc.LoadConf()

# Check for empty PETSC_ARCH
emptyarch = not ('PETSC_ARCH' in os.environ and os.environ['PETSC_ARCH'])
if emptyarch:
  pseudoarch = 'arch-' + sys.platform.replace('cygwin','mswin')+ '-' + petsc.language.lower().replace('+','x')
  if petsc.debug:
    pseudoarch += '-debug'
  else:
    pseudoarch += '-opt'
  if not 'real' in petsc.scalar:
    pseudoarch += '-' + petsc.scalar
  archname = 'installed-'+pseudoarch.replace('linux-','linux2-')
else:
  archname = petsc.arch

# Create directories for configuration files
archdir, archdirexisted = slepc.CreateDirTest(slepc.dir,archname)
libdir  = slepc.CreateDir(archdir,'lib')
confdir = slepc.CreateDirTwo(libdir,'slepc','conf')

# Open log file
log.Open(slepc.dir,confdir,'configure.log')
log.write('='*80)
log.write('Starting Configure Run at '+time.ctime(time.time()))
log.write('Configure Options: '+' '.join(sys.argv[1:]))
log.write('Working directory: '+os.getcwd())
log.write('Python version:\n'+sys.version)
log.write('make: '+petsc.make)

# Some checks related to PETSc configuration
if petsc.nversion < slepc.nversion:
  log.Exit('This SLEPc version is not compatible with PETSc version '+petsc.version)
if not petsc.precision in ['double','single','__float128']:
  log.Exit('This SLEPc version does not work with '+petsc.precision+' precision')

# Display versions and paths
log.write('PETSc source directory: '+petsc.dir)
log.write('PETSc install directory: '+petsc.prefixdir)
log.write('PETSc version: '+petsc.lversion)
if not emptyarch:
  log.write('PETSc architecture: '+petsc.arch)
log.write('SLEPc source directory: '+slepc.dir)
if slepc.isinstall:
  log.write('SLEPc install directory: '+slepc.prefixdir)
log.write('SLEPc version: '+slepc.lversion)

# Clean previous configuration if needed
if archdirexisted:
  if slepc.isinstall and not slepc.clean:
    log.Exit('You are requesting a prefix install but the arch directory '+archdir+' already exists and may contain files from previous builds; consider adding option --with-clean')
  if slepc.clean:
    log.Println('\nCleaning arch dir '+archdir+'...')
    try:
      for root, dirs, files in os.walk(archdir,topdown=False):
        for name in files:
          if name!='configure.log':
            os.remove(os.path.join(root,name))
    except:
      log.Exit('Cannot remove existing files in '+archdir)
    for rdir in ['obj','externalpackages']:
      try:
        shutil.rmtree(os.path.join(archdir,rdir))
      except: pass

# Create other directories and configuration files
if not slepc.prefixdir:
  slepc.prefixdir = archdir
includedir = slepc.CreateDir(archdir,'include')
modulesdir = slepc.CreateDirTwo(confdir,'modules','slepc')
pkgconfdir = slepc.CreateDir(libdir,'pkgconfig')
slepcvars  = slepc.CreateFile(confdir,'slepcvariables')
slepcconf  = slepc.CreateFile(includedir,'slepcconf.h')
pkgconfig  = slepc.CreateFile(pkgconfdir,'SLEPc.pc')
if slepc.isinstall:
  modules  = slepc.CreateFile(modulesdir,slepc.lversion)
else:
  modules  = slepc.CreateFile(modulesdir,slepc.lversion+'-'+archname)
  reconfig = slepc.CreateFile(confdir,'reconfigure-'+archname+'.py')
  reconfigpath = os.path.join(confdir,'reconfigure-'+archname+'.py')

# Write initial part of file slepcvariables
slepcvars.write('SLEPC_CONFIGURE_OPTIONS = '+argdb.UsedArgs()+'\n')
slepcvars.write('SLEPC_INSTALLDIR = '+slepc.prefixdir+'\n')
if emptyarch:
  slepcvars.write('INSTALLED_PETSC = 1\n')
if slepc.datadir:
  slepcvars.write('DATAFILESPATH = '+slepc.datadir+'\n')

# Write initial part of file slepcconf.h
slepcconf.write('#if !defined(SLEPCCONF_H)\n')
slepcconf.write('#define SLEPCCONF_H\n\n')
AddDefine(slepcconf,'PETSC_DIR',petsc.dir)
AddDefine(slepcconf,'PETSC_ARCH',petsc.arch)
AddDefine(slepcconf,'DIR',slepc.dir)
AddDefine(slepcconf,'LIB_DIR',os.path.join(slepc.prefixdir,'lib'))
if slepc.isrepo:
  AddDefine(slepcconf,'VERSION_GIT',slepc.gitrev)
  AddDefine(slepcconf,'VERSION_DATE_GIT',slepc.gitdate)
  AddDefine(slepcconf,'VERSION_BRANCH_GIT',slepc.branch)

# Create global configuration file for the case of empty PETSC_ARCH
if emptyarch:
  globconf = slepc.CreateFile(os.path.join(slepc.dir,'lib','slepc','conf'),'slepcvariables')
  globconf.write('SLEPC_DIR = '+slepc.dir+'\n')
  globconf.write('PETSC_ARCH = '+archname+'\n')
  globconf.close()

# Check if PETSc is working
log.NewSection('Checking PETSc installation...')
if petsc.nversion > slepc.nversion:
  log.Warn('PETSc version '+petsc.version+' is newer than SLEPc version '+slepc.version)
if slepc.release=='1' and not petsc.release=='1':
  log.Exit('A release version of SLEPc requires a release version of PETSc, not a development version')
if slepc.release=='0' and petsc.release=='1':
  log.Exit('A development version of SLEPc cannot be built with a release version of PETSc')
if petsc.isinstall:
  if os.path.realpath(petsc.prefixdir) != os.path.realpath(petsc.dir):
    log.Warn('PETSC_DIR does not point to PETSc installation path')
petsc.Check()
if not petsc.havepackage:
  log.Exit('Unable to link with PETSc')

# Single library installation
if petsc.singlelib:
  slepcvars.write('SHLIBS = libslepc\n')
  slepcvars.write('LIBNAME = '+os.path.join('${INSTALL_LIB_DIR}','libslepc.${AR_LIB_SUFFIX}')+'\n')
  for module in ['SYS','EPS','SVD','PEP','NEP','MFN','LME']:
    slepcvars.write('SLEPC_'+module+'_LIB = ${CC_LINKER_SLFLAG}${SLEPC_LIB_DIR} -L${SLEPC_LIB_DIR} -lslepc ${SLEPC_EXTERNAL_LIB} ${PETSC_SNES_LIB}\n')
  slepcvars.write('SLEPC_LIB = ${CC_LINKER_SLFLAG}${SLEPC_LIB_DIR} -L${SLEPC_LIB_DIR} -lslepc ${SLEPC_EXTERNAL_LIB} ${PETSC_SNES_LIB}\n')

# Check for external packages and for missing LAPACK functions
for pkg in checkpackages:
  pkg.Process(slepcconf,slepcvars,slepc,petsc,archdir)

# Write Modules and pkg-config configuration files
log.NewSection('Writing various configuration files...')
log.write('Modules file in '+modulesdir)
if slepc.isinstall:
  WriteModulesFile(modules,slepc.lversion,slepc.prefixdir)
else:
  WriteModulesFile(modules,slepc.lversion,slepc.dir)
log.write('pkg-config file in '+pkgconfdir)
WritePkgconfigFile(pkgconfig,slepc.lversion,petsc.version,slepc.dir,slepc.isinstall,slepc.prefixdir,petsc.singlelib)

# Write reconfigure file
if not slepc.isinstall:
  WriteReconfigScript(reconfig,slepc.dir,argdb.UsedArgs())
  try:
    os.chmod(reconfigpath,0o775)
  except OSError as e:
    log.Exit('Unable to make reconfigure script executable:\n'+str(e))

# Finish with configuration files (except slepcvars)
slepcconf.write('\n#endif\n')
slepcconf.close()
pkgconfig.close()
modules.close()
if not slepc.isinstall: reconfig.close()

# Download sowing if requested and make Fortran stubs if necessary
bfort = petsc.bfort
if sowing.downloadpackage:
  bfort = sowing.DownloadAndInstall(slepc,petsc,archdir)

if slepc.isrepo and petsc.fortran:
  try:
    if not os.path.exists(bfort):
      bfort = os.path.join(archdir,'bin','bfort')
    if not os.path.exists(bfort):
      bfort = sowing.DownloadAndInstall(slepc,petsc,archdir)
    log.NewSection('Generating Fortran stubs...')
    log.write('Using BFORT='+bfort)
    sys.path.insert(0, os.path.abspath(os.path.join('lib','slepc','bin','maint')))
    import generatefortranstubs
    generatefortranstubs.main(slepc.dir,bfort,os.getcwd(),0)
    generatefortranstubs.processf90interfaces(slepc.dir,0)
  except:
    log.Exit('Try configuring with --download-sowing or use a git version of PETSc')

if bfort != petsc.bfort:
  slepcvars.write('BFORT = '+bfort+'\n')

# Finally we can close the slepcvariables file
slepcvars.close()

# Print summary
log.NewSection('')
log.Println('')
log.Println('='*80)
log.Println('SLEPc Configuration')
log.Println('='*80)
log.Println('\nSLEPc directory:\n '+slepc.dir)
if slepc.isrepo:
  log.Println('  It is a git repository on branch: '+slepc.branch)
if slepc.isinstall:
  log.Println('SLEPc prefix directory:\n '+slepc.prefixdir)
log.Println('PETSc directory:\n '+petsc.dir)
if petsc.isrepo:
  log.Println('  It is a git repository on branch: '+petsc.branch)
  if slepc.isrepo and petsc.branch!='maint' and slepc.branch!='maint':
    try:
      import dateutil.parser, datetime
      petscdate = dateutil.parser.parse(petsc.gitdate)
      slepcdate = dateutil.parser.parse(slepc.gitdate)
      if abs(petscdate-slepcdate)>datetime.timedelta(days=30):
        log.Warn('Your PETSc and SLEPc repos may not be in sync (more than 30 days apart)')
    except ImportError: pass
if emptyarch and slepc.isinstall:
  log.Println('Prefix install with '+petsc.precision+' precision '+petsc.scalar+' numbers')
else:
  log.Println('Architecture "'+archname+'" with '+petsc.precision+' precision '+petsc.scalar+' numbers')
for pkg in checkpackages:
  pkg.ShowInfo()
log.write('\nFinishing Configure Run at '+time.ctime(time.time()))
log.write('='*80)
print()
print('xxx'+'='*74+'xxx')
print(' Configure stage complete. Now build the SLEPc library with:')
if emptyarch:
  print('   make SLEPC_DIR='+slepc.dir+' PETSC_DIR='+petsc.dir)
else:
  print('   make SLEPC_DIR='+slepc.dir+' PETSC_DIR='+petsc.dir+' PETSC_ARCH='+archname)
print('xxx'+'='*74+'xxx')
print()
