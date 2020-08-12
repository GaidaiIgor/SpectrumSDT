#!/usr/bin/env python

# Used by the makefile lib/petsc/conf/rules html to replace #includes in the source files with links when possible

import sys, os, re

if __name__ == "__main__":
  if (len(sys.argv) < 3): sys.exit(1)
  filename = sys.argv[1]
  petscdir = sys.argv[2]
  root = os.path.relpath(os.path.realpath(petscdir),os.path.realpath(os.getcwd()))
  froot = os.path.relpath(os.path.realpath(petscdir),os.path.dirname(os.path.realpath(filename)))

  for line in sys.stdin:
    if 'include' in line:
      ins = re.search('#include [ ]*&lt;',line)
      if ins:
        includename = line[ins.end():re.search('&gt;[a-zA-Z0-9/<>#*"=. ]*',line).start()]
        ln = re.search('<a name="line[0-9]*">[ 0-9]*: </a>',line)
        linenumber = line[ln.start():ln.end()]
        if os.path.isfile(includename):
          sys.stdout.write(linenumber+' #include <A href="'+includename+'.html">&lt;'+includename+'&gt;</A>\n')
        elif os.path.isfile(os.path.join(root,'include',includename)):
          sys.stdout.write(linenumber+' #include <A href="'+os.path.relpath(os.path.join(froot,'include',includename))+'.html">&lt;'+includename+'&gt;</A>\n')
        elif os.path.isfile(os.path.join(root,includename)):
          sys.stdout.write(linenumber+' #include <A href="'+os.path.relpath(os.path.join(froot,includename))+'.html">&lt;'+includename+'&gt;</A>\n')
        else:
          sys.stdout.write(line)
      else:
        sys.stdout.write(line)
    else:
      sys.stdout.write(line)

