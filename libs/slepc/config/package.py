from __future__ import print_function
#
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  SLEPc - Scalable Library for Eigenvalue Problem Computations
#  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
#
#  This file is part of SLEPc.
#  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#

import os, sys, tempfile, shutil, tarfile
import log, argdb
try:
  from urllib import urlretrieve
except ImportError:
  from urllib.request import urlretrieve
try:
  import urlparse as urlparse_local
except ImportError:
  from urllib import parse as urlparse_local
if sys.version_info < (3,):
  import commands
else:
  import subprocess
import socket

# Fix parsing for nonstandard schemes
urlparse_local.uses_netloc.extend(['bk', 'ssh', 'svn'])

class Package:

  def __init__(self,argdb,log):
    self.installable     = False  # an already installed package can be picked --with-xxx-dir
    self.downloadable    = False  # package can be downloaded and installed with --download-xxx
    self.downloadpackage = 0
    self.packagedir      = ''
    self.packagelibs     = []
    self.packageincludes = []
    self.packageurl      = ''
    self.buildflags      = ''
    self.log             = log
    self.supportsscalar  = ['real', 'complex']
    self.supportssingle  = False
    self.supports64bint  = False
    self.fortran         = False
    self.hasheaders      = False
    self.hasdloadflags   = False

  def RunCommand(self,instr):
    try:
      self.log.write('- '*35+'\nRunning command:\n'+instr+'\n'+'- '*35)
    except AttributeError: pass
    if sys.version_info < (3,):
      (result,output) = commands.getstatusoutput(instr)
    else:
      try:
        output = subprocess.check_output(instr,shell=True,stderr=subprocess.STDOUT)
        result = 0
      except subprocess.CalledProcessError as ex:
        output = ex.output
        result = ex.returncode
      output = output.decode(encoding='UTF-8',errors='replace').rstrip()
    try:
      self.log.write(output)
    except AttributeError: pass
    return (result,output)

  def ProcessArgs(self,argdb):
    self.requested = False
    self.havepackage = False
    if self.installable:
      string,found = argdb.PopPath('with-'+self.packagename+'-dir',exist=True)
      if found:
        self.requested = True
        self.packagedir = string
      string,found = argdb.PopString('with-'+self.packagename+'-lib')
      if found:
        self.requested = True
        self.packagelibs = string.split(',')
      if self.hasheaders:
        string,found = argdb.PopString('with-'+self.packagename+'-include')
        if found:
          self.requested = True
          self.packageincludes = string.split(',')
      value,found = argdb.PopBool('with-'+self.packagename)
      if found:
        self.requested = value
    if self.downloadable:
      string,cflagsfound = argdb.PopString('download-'+self.packagename+'-cflags')
      url,flag,found = argdb.PopUrl('download-'+self.packagename)
      if found:
        if self.requested:
          self.log.Exit('Cannot request both download and install simultaneously')
        self.requested = True
        self.download = True
        self.packageurl = url
        self.downloadpackage = flag
      if cflagsfound:
        if not hasattr(self,'download') or not self.download:
          self.log.Exit('--download-'+self.packagename+'-cflags must be used together with --download-'+self.packagename)
        self.buildflags = string

  def Process(self,conf,vars,slepc,petsc,archdir=''):
    self.make = petsc.make
    if petsc.buildsharedlib:
      self.slflag = petsc.slflag
    if self.requested:
      name = self.packagename.upper()
      if self.downloadpackage:
        if hasattr(self,'version'):
          self.log.NewSection('Installing '+name+' version '+self.version+'...')
        else:
          self.log.NewSection('Installing '+name+'...')
        self.Precondition(petsc)
        self.DownloadAndInstall(conf,vars,slepc,petsc,archdir,slepc.prefixdir)
      elif self.installable:
        self.log.NewSection('Checking '+name+'...')
        self.Precondition(petsc)
        self.Check(conf,vars,petsc,archdir)
      try:
        self.LoadVersion(conf)
        self.log.write('Version number for '+name+' is '+self.iversion)
      except AttributeError:
        pass

  def Precondition(self,petsc):
    package = self.packagename.upper()
    if petsc.scalar == 'complex':
      if 'complex' not in self.supportsscalar:
        self.log.Exit(package+' does not support complex scalars')
    elif petsc.scalar == 'real':
      if 'real' not in self.supportsscalar:
        self.log.Exit(package+' is supported only with complex scalars')
    if petsc.precision == 'single':
      if not self.supportssingle:
        self.log.Exit(package+' is supported only in double precision')
    elif petsc.precision != 'double':
      self.log.Exit('Precision '+petsc.precision+' is not supported for external packages')
    if petsc.ind64 and not self.supports64bint:
      self.log.Exit(package+' cannot be used with 64-bit integers')
    if self.downloadpackage and self.fortran and not hasattr(petsc,'fc'):
      self.log.Exit('Option --download-'+self.packagename+' requires a Fortran compiler')

  def Download(self,externdir,downloaddir):
    # Check if source is already available
    if os.path.exists(os.path.join(externdir,self.dirname)):
      self.log.write('Using '+os.path.join(externdir,self.dirname))
      dirname = self.dirname
    else:

      if downloaddir:
        # Get tarball from download dir
        if self.packageurl=='':
          localFile = os.path.join(downloaddir,self.archive)
        else:
          localFile = os.path.join(downloaddir,self.packageurl)
        if not os.path.exists(localFile):
          self.log.Exit('Could not find file '+self.archive+' under '+downloaddir)
        url = localFile
        filename = os.path.basename(url)
      else:
        # Download tarball
        url = self.packageurl
        if url=='':
          url = self.url
        if os.path.exists(url):
          url = 'file:'+url
        filename = os.path.basename(urlparse_local.urlparse(url)[2])
        localFile = os.path.join(externdir,self.archive)
        self.log.write('Downloading '+url+' to '+localFile)

        if os.path.exists(localFile):
          os.remove(localFile)
        try:
          sav_timeout = socket.getdefaulttimeout()
          socket.setdefaulttimeout(30)
          urlretrieve(url, localFile)
          socket.setdefaulttimeout(sav_timeout)
        except Exception as e:
          socket.setdefaulttimeout(sav_timeout)
          failureMessage = '''\
Unable to download package %s from: %s
* If URL specified manually - perhaps there is a typo?
* If your network is disconnected - please reconnect and rerun ./configure
* Or perhaps you have a firewall blocking the download
* You can run with --with-packages-download-dir=/adirectory and ./configure will instruct you what packages to download manually
* or you can download the above URL manually, to /yourselectedlocation/%s
  and use the configure option:
  --download-%s=/yourselectedlocation/%s
''' % (self.packagename.upper(), url, filename, self.packagename, filename)
          self.log.Exit(failureMessage)

      # Uncompress tarball
      extractdir = os.path.join(externdir,self.dirname)
      self.log.write('Uncompressing '+localFile+' to directory '+extractdir)
      if os.path.exists(extractdir):
        for root, dirs, files in os.walk(extractdir, topdown=False):
          for name in files:
            os.remove(os.path.join(root,name))
          for name in dirs:
            os.rmdir(os.path.join(root,name))
      failureMessage = '''\
Downloaded package %s from: %s is not a tarball.
[or installed python cannot process compressed files]
* If you are behind a firewall - please fix your proxy and rerun ./configure
  For example at LANL you may need to set the environmental variable http_proxy (or HTTP_PROXY?) to  http://proxyout.lanl.gov
* You can run with --with-packages-download-dir=/adirectory and ./configure will instruct you what packages to download manually
* or you can download the above URL manually, to /yourselectedlocation/%s
  and use the configure option:
  --download-%s=/yourselectedlocation/%s
''' % (self.packagename.upper(), url, filename, self.packagename, filename)
      try:
        tf = tarfile.open(localFile)
      except tarfile.ReadError as e:
        self.log.Exit(str(e)+'\n'+failureMessage)
      if not tf: self.log.Exit(failureMessage)
      #git puts 'pax_global_header' as the first entry and some tar utils process this as a file
      firstname = tf.getnames()[0]
      if firstname == 'pax_global_header':
        firstmember = tf.getmembers()[1]
      else:
        firstmember = tf.getmembers()[0]
      # some tarfiles list packagename/ but some list packagename/filename in the first entry
      if firstmember.isdir():
        dirname = firstmember.name
      else:
        dirname = os.path.dirname(firstmember.name)
      tf.extractall(path=externdir)
      tf.close()

      # fix file permissions for the untared tarballs
      try:
        # check if 'dirname' is set'
        if dirname:
          (result,output) = self.RunCommand('cd '+externdir+'; chmod -R a+r '+dirname+'; find '+dirname+' -type d -name "*" -exec chmod a+rx {} \;')
        else:
          self.log.Warn('Could not determine dirname extracted by '+localFile+' to fix file permissions')
      except RuntimeError as e:
        self.log.Exit('Error changing permissions for '+dirname+' obtained from '+localFile+ ' : '+str(e))
      if not downloaddir:
        os.remove(localFile)
    return os.path.join(externdir,dirname)

  wd = 36

  def ShowHelp(self):
    wd = Package.wd
    if self.downloadable or self.installable:
      print(self.packagename.upper()+':')
    if self.downloadable:
      print(('  --download-'+self.packagename+'[=<fname>]').ljust(wd)+': Download and install '+self.packagename.upper())
      if self.hasdloadflags:
        print(('  --download-'+self.packagename+'-cflags=<flags>').ljust(wd)+': Indicate extra flags to compile '+self.packagename.upper())
    if self.installable:
      print(('  --with-'+self.packagename+'=<bool>').ljust(wd)+': Indicate if you wish to test for '+self.packagename.upper())
      print(('  --with-'+self.packagename+'-dir=<dir>').ljust(wd)+': Indicate the root directory of the '+self.packagename.upper()+' installation')
      print(('  --with-'+self.packagename+'-lib=<libraries>').ljust(wd)+': Indicate comma-separated libraries and link flags for '+self.packagename.upper())
      if self.hasheaders:
        print(('  --with-'+self.packagename+'-include=<dirs>').ljust(wd)+': Indicate the directory of the '+self.packagename.upper()+' include files')

  def ShowInfo(self):
    if self.havepackage:
      self.log.Println(self.packagename.upper()+' library flags:')
      self.log.Println(' '+' '.join(self.packageflags))

  def TestRuns(self,petsc):
    if self.havepackage:
      return [self.packagename.upper()]
    else:
      return []

  def Link(self,functions,callbacks,flags,givencode='',cflags='',clanguage='c',logdump=True):

    # Create temporary directory and makefile
    try:
      tmpdir = tempfile.mkdtemp(prefix='slepc-')
      if not os.path.isdir(tmpdir): os.mkdir(tmpdir)
    except:
      self.log.Exit('Cannot create temporary directory')
    try:
      makefile = open(os.path.join(tmpdir,'makefile'),'w')
      if cflags!='':
        if clanguage=='c++': makefile.write('CXXFLAGS='+cflags+'\n')
        else: makefile.write('CFLAGS='+cflags+'\n')
      makefile.write('checklink: checklink.o\n')
      makefile.write('\t${CLINKER} -o checklink checklink.o ${LINKFLAGS} ${PETSC_SNES_LIB}\n')
      makefile.write('\t@${RM} -f checklink checklink.o\n')
      makefile.write('LOCDIR = ./\n')
      makefile.write('include '+os.path.join('${PETSC_DIR}','lib','petsc','conf','variables')+'\n')
      makefile.write('include '+os.path.join('${PETSC_DIR}','lib','petsc','conf','rules')+'\n')
      makefile.close()
    except:
      self.log.Exit('Cannot create makefile in temporary directory')

    # Create source file
    if givencode == '':
      code = '#include "petscsnes.h"\n'
      for f in functions:
        code += 'PETSC_EXTERN int\n' + f + '();\n'

      for c in callbacks:
        code += 'int '+ c + '() { return 0; } \n'

      code += 'int main() {\n'
      code += 'Vec v; Mat m; KSP k;\n'
      code += 'PetscInitializeNoArguments();\n'
      code += 'VecCreate(PETSC_COMM_WORLD,&v);\n'
      code += 'MatCreate(PETSC_COMM_WORLD,&m);\n'
      code += 'KSPCreate(PETSC_COMM_WORLD,&k);\n'
      for f in functions:
        code += f + '();\n'
      code += 'return 0;\n}\n'
    else:
      code = givencode

    cfile = open(os.path.join(tmpdir,'checklink.c'),'w')
    cfile.write(code)
    cfile.close()
    if logdump:
      try:
        self.log.write('- '*35+'\nChecking link with code:\n')
        self.log.write(code)
      except AttributeError: pass

    # Try to compile test program
    (result, output) = self.RunCommand('cd ' + tmpdir + ';' + self.make + ' checklink LINKFLAGS="'+' '.join(flags)+'"')
    shutil.rmtree(tmpdir)

    if result:
      return (0,code + output)
    else:
      return (1,code + output)

  def FortranLink(self,functions,callbacks,flags):
    output = '\n=== With linker flags: '+' '.join(flags)

    f = []
    for i in functions:
      f.append(i+'_')
    c = []
    for i in callbacks:
      c.append(i+'_')
    (result, output1) = self.Link(f,c,flags,logdump=False)
    output1 = '\n====== With underscore Fortran names\n' + output1
    if result: return ('UNDERSCORE',output1)

    f = []
    for i in functions:
      f.append(i.upper())
    c = []
    for i in callbacks:
      c.append(i.upper())
    (result, output2) = self.Link(f,c,flags,logdump=False)
    output2 = '\n====== With capital Fortran names\n' + output2
    if result: return ('CAPS',output2)

    return ('',output + output1 + output2)

  def GenerateGuesses(self,name,archdir,word='lib'):
    installdirs = [os.path.join(os.path.sep,'usr','local'),os.path.join(os.path.sep,'opt')]
    if 'HOME' in os.environ:
      installdirs.insert(0,os.environ['HOME'])

    dirs = []
    for i in installdirs:
      dirs = dirs + [os.path.join(i,word)]
      for d in [name,name.upper(),name.lower()]:
        dirs = dirs + [os.path.join(i,d)]
        dirs = dirs + [os.path.join(i,d,word)]
        dirs = dirs + [os.path.join(i,word,d)]

    for d in dirs[:]:
      if not os.path.exists(d):
        dirs.remove(d)
    dirs = [''] + dirs + [os.path.join(archdir,word)]
    return dirs

  def FortranLib(self,conf,vars,dirs,libs,functions,callbacks = []):
    name = self.packagename.upper()
    error = ''
    mangling = ''
    for d in dirs:
      for l in libs:
        if d:
          if hasattr(self,'slflag'):
            flags = [self.slflag + d] + ['-L' + d] + l
          else:
            flags = ['-L' + d] + l
        else:
          flags = l
        (mangling, output) = self.FortranLink(functions,callbacks,flags)
        error += output
        if mangling: break
      if mangling: break

    if mangling:
      self.log.write(output)
    else:
      self.log.write(error)
      self.log.Exit('Unable to link with '+name+' library in directories '+' '.join(dirs)+' with libraries and link flags '+' '.join(flags))

    conf.write('#define SLEPC_HAVE_' + name + ' 1\n#define SLEPC_' + name + '_HAVE_'+mangling+' 1\n')
    vars.write(name + '_LIB = '+' '.join(flags)+'\n')
    self.havepackage = True
    self.packageflags = flags

