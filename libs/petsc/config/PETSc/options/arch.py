import config.base
import os
import re

class Configure(config.base.Configure):
  def __init__(self, framework):
    config.base.Configure.__init__(self, framework)
    self.headerPrefix = 'PETSC'
    self.substPrefix  = 'PETSC'
    return

  def __str1__(self):
    if not hasattr(self, 'arch'):
      return ''
    desc = ['PETSc:']
    desc.append('  PETSC_ARCH: '+str(self.arch))
    return '\n'.join(desc)+'\n'

  def setupHelp(self, help):
    import nargs
    help.addArgument('PETSc', '-PETSC_ARCH=<string>',     nargs.Arg(None, None, 'The configuration name'))
    help.addArgument('PETSc', '-with-petsc-arch=<string>',nargs.Arg(None, None, 'The configuration name'))
    help.addArgument('PETSc', '-force=<bool>',            nargs.ArgBool(None, 0, 'Bypass configure hash caching, and run to completion'))
    return

  def setupDependencies(self, framework):
    self.sourceControl = framework.require('config.sourceControl',self)
    self.petscdir = framework.require('PETSc.options.petscdir', self)
    return

  def setNativeArchitecture(self):
    import sys
    arch = 'arch-' + sys.platform.replace('cygwin','mswin')
    # use opt/debug, c/c++ tags.s
    arch+= '-'+self.framework.argDB['with-clanguage'].lower().replace('+','x')
    if self.framework.argDB['with-debugging']:
      arch += '-debug'
    else:
      arch += '-opt'
    self.nativeArch = arch
    return

  def configureArchitecture(self):
    '''Checks PETSC_ARCH and sets if not set'''
    # Warn if PETSC_ARCH doesn't match env variable
    if 'PETSC_ARCH' in self.framework.argDB and 'PETSC_ARCH' in os.environ and self.framework.argDB['PETSC_ARCH'] != os.environ['PETSC_ARCH']:
      self.logPrintBox('''\
Warning: PETSC_ARCH from environment does not match command-line or name of script.
Warning: Using from command-line or name of script: %s, ignoring environment: %s''' % (str(self.framework.argDB['PETSC_ARCH']), str(os.environ['PETSC_ARCH'])))
      os.environ['PETSC_ARCH'] = self.framework.argDB['PETSC_ARCH']
    if 'with-petsc-arch' in self.framework.argDB:
      self.arch = self.framework.argDB['with-petsc-arch']
      msg = 'option -with-petsc-arch='+str(self.arch)
    elif 'PETSC_ARCH' in self.framework.argDB:
      self.arch = self.framework.argDB['PETSC_ARCH']
      msg = 'option PETSC_ARCH='+str(self.arch)
    elif 'PETSC_ARCH' in os.environ:
      self.arch = os.environ['PETSC_ARCH']
      msg = 'environment variable PETSC_ARCH='+str(self.arch)
    else:
      self.arch = self.nativeArch
    if self.arch.find('/') >= 0 or self.arch.find('\\') >= 0:
      raise RuntimeError('PETSC_ARCH should not contain path characters, but you have specified with '+msg)
    if self.arch.startswith('-'):
      raise RuntimeError('PETSC_ARCH should not start with "-", but you have specified with '+msg)
    if self.arch.startswith('.'):
      raise RuntimeError('PETSC_ARCH should not start with ".", but you have specified with '+msg)
    if not len(self.arch):
      raise RuntimeError('PETSC_ARCH cannot be empty string. Use a valid string or do not set one. Currently set with '+msg)
    self.archBase = re.sub(r'^(\w+)[-_]?.*$', r'\1', self.arch)
    return

  def makeDependency(self,hash,hashfile,hashfilepackages):
    '''Deletes the current hashfile and saves the hashfile names and its value in framework so that'''
    '''framework.Configure can create the file upon success of configure'''
    import os
    if hash:
      self.framework.hash = hash
      self.framework.hashfile = hashfile
      self.logPrint('Setting hashfile: '+hashfile)
      if hashfilepackages: self.framework.hashfilepackages = hashfilepackages
    try:
      self.logPrint('Deleting configure hash file: '+hashfile)
      os.remove(hashfile)
      self.logPrint('Deleted configure hash file: '+hashfile)
    except:
      self.logPrint('Unable to delete configure hash file: '+hashfile)


  def checkDependency(self):
    '''Checks if files in config have changed, the command line options have changed or the PATH has changed'''
    '''  By default - checks if configure needs to be run'''
    '''  If --arch-hash it manages the same information but it:'''
    '''     * computes a short hash for the configuration <hashvalue>'''
    '''     * sets self.arch and PETSC_ARCH to arch-<hashvalue>'''
    '''       This results in the downloaded packages being installed once to the arch-<hasvalue> directory'''
    '''       and a new directory with a different hash is created if the configuration changes.'''
    '''     This mode is intended mostly for testing to reduce reconfigure and recompile times (not currently used)'''
    '''  If --package-prefix-hash=directory is provided'''
    '''     * computes a short hash for the configuration <hashvalue>'''
    '''     * puts the downloaded external packages into location directory/hash'''
    '''       This results in the downloaded packages being installed once'''
    '''       and a new directory with a different hash is created if the configuration changes.'''
    '''     This mode is intended mostly for testing to reduce time of reinstalling external packages'''
    import os
    import sys
    import hashlib
    import platform
    if sys.version_info < (3,):
      hash = 'Uname: '+platform.uname()[0]+' '+platform.uname()[4]+'\n'
    else:
      hash = 'Uname: '+platform.uname().system+' '+platform.uname().processor+'\n'
    hash += 'PATH=' + os.environ.get('PATH', '') + '\n'
    args = sorted(set(filter(lambda x: not (x.startswith('PETSC_ARCH') or x == '--force'),sys.argv[1:])))
    hash += 'args:\n' + '\n'.join('    '+a for a in args) + '\n'
    chash=''
    try:
      for root, dirs, files in os.walk('config'):
        if root == 'config':
          dirs.remove('examples')
        for f in files:
          if not f.endswith('.py') or f.startswith('.') or f.startswith('#'):
            continue
          fname = os.path.join(root, f)
          with open(fname,'rb') as f:
            chash += hashlib.sha256(f.read()).hexdigest() + '  ' + fname + '\n'
    except:
      self.logPrint('Error generating file list/hash from config directory for configure hash, forcing new configuration')
      return
    hash += '\n'.join(sorted(chash.splitlines()))
    hashfilepackages = None
    # Generate short hash to use for the arch so the same arch can be reused if the configuration files don't change
    if 'arch-hash' in self.argDB:
      if self.argDB['prefix']:
        raise RuntimeError('Cannot provide --prefix and --arch-hash')
      if hasattr(self.argDB,'PETSC_ARCH'):
        raise RuntimeError('Cannot provide PETSC_ARCH and --arch-hash')
      if 'package-prefix-hash' in self.argDB:
        raise RuntimeError('Cannot provide --arch-hash and --package-prefix-hash')
      if os.getenv('PETSC_ARCH'):
        raise RuntimeError('Do not set the environmental variable PETSC_ARCH and use --arch-hash')
    if 'arch-hash' in self.argDB or 'package-prefix-hash' in self.argDB:
      import hashlib
      m = hashlib.md5()
      m.update(hash.encode('utf-8'))
      hprefix = m.hexdigest()
      if 'arch-hash' in self.argDB:
        self.argDB['PETSC_ARCH'] = 'arch-'+hprefix[0:6]
        self.arch = 'arch-'+hprefix[0:6]
      else:
        if not os.path.isdir(self.argDB['package-prefix-hash']):
          self.logPrintBox('Specified package-prefix-hash location %s not found! Attemping to create this dir!' % self.argDB['package-prefix-hash'])
          try:
            os.makedirs(self.argDB['package-prefix-hash'])
          except Exception as e:
            self.logPrint('Error creating package-prefix-hash directory '+self.argDB['package-prefix-hash']+': '+str(e))
            raise RuntimeError('You must have write permission to create this directory!')
        status = False
        for idx in range(6,len(hprefix)):
          hashdirpackages = os.path.join(self.argDB['package-prefix-hash'],hprefix[0:idx])
          hashfilepackages = os.path.join(hashdirpackages,'configure-hash')
          if os.path.isdir(hashdirpackages):
            if os.path.exists(hashfilepackages):
              self.argDB['package-prefix-hash'] = 'reuse' # indicates prefix libraries already built, no need to rebuild
              status = True
              break
            else: continue # perhaps an incomplete build? use a longer hash
          else:
            try:
              os.mkdir(hashdirpackages)
            except Exception as e:
              self.logPrint('Error creating package-prefix-hash directory '+hashdirpackages+': '+str(e))
              raise RuntimeError('You must have write permission on --package-prefix-hash='+self.argDB['package-prefix-hash']+' directory')
            status = True
            break
        if not status:
          raise RuntimeError('Unable to create package-prefix-hash dir! Suggest cleaning up %s* !' % os.path.join(self.argDB['package-prefix-hash'],hprefix[0:6]) )
        self.argDB['prefix'] = hashdirpackages

    hashfile = os.path.join(self.arch,'lib','petsc','conf','configure-hash')

    if self.argDB['force']:
      self.logPrint('Forcing a new configuration requested by use')
      self.makeDependency(hash,hashfile,hashfilepackages)
      return
    a = ''
    try:
      with open(hashfile, 'r') as f:
        a = f.read()
    except:
      self.logPrint('No previous hashfile found')
      self.makeDependency(hash,hashfile,hashfilepackages)
      return
    if a == hash:
      try:
        self.logPrint('Attempting to save lib/petsc/conf/petscvariables file')
        with open(os.path.join('lib','petsc','conf','petscvariables'), 'w') as g:
          g.write('PETSC_ARCH='+self.arch+'\n')
          g.write('PETSC_DIR='+self.petscdir.dir+'\n')
          g.write('include $(PETSC_DIR)/$(PETSC_ARCH)/lib/petsc/conf/petscvariables\n')
          self.logPrint('Saved lib/petsc/conf/petscvariables file')
      except:
        self.logPrint('Unable to save lib/petsc/conf/petscvariables file')
      self.logPrint('configure hash file: '+hashfile+' matches; no need to run configure.')
      print('Your configure options and state has not changed; no need to run configure')
      print('However you can force a configure run using the option: --force')
      sys.exit()
    self.logPrint('configure hash file: '+hashfile+' does not match\n'+a+'\n---\n'+hash+'\n need to run configure')
    self.makeDependency(hash,hashfile,hashfilepackages)

  def configure(self):
    self.executeTest(self.setNativeArchitecture)
    self.executeTest(self.configureArchitecture)
    # required by top-level configure.py
    self.framework.arch = self.arch
    self.checkDependency()
    return
