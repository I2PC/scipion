#!/usr/bin/env python

import sys, os

from pyworkflow.manager import Manager
import pyworkflow.utils as pwutils
import pyworkflow.em as em


def usage(error):
    print """
    ERROR: %s
    
    Usage: fixlinks.py PROJECT SEARCH_DIR
        PROJECT: provide the project name to fix broken links in the imports.
        SEARCH_DIR: provide a directory where to look for the files.
        and fix the links.    
    """ % error
    sys.exit(1)    


if len(sys.argv) != 3:
    usage("Incorrect number of input parameters")

projName = sys.argv[1]
searchDir = os.path.abspath(sys.argv[2])

# Create a new project
manager = Manager()

if not manager.hasProject(projName):
    usage("Unexistent project: %s" % pwutils.red(projName))
    
if not os.path.exists(searchDir):
    usage("Unexistent SEARCH_DIR: %s" % pwutils.red(searchDir))
    
project = manager.loadProject(projName)

runs = project.getRuns()

for prot in runs:
    broken = False
    if isinstance(prot, em.ProtImport):
        for _, attr in prot.iterOutputEM():
            fn = attr.getFiles()
            for f in attr.getFiles():
                if ':' in f:
                    f = f.split(':')[0]

                if not os.path.exists(f):                    
                    if not broken:
                        broken = True
                        print "Found broken links in run: ", pwutils.magenta(prot.getRunName())
                    print "  Missing: ", pwutils.magenta(f)
                    if os.path.islink(f):
                        print "    -> ", pwutils.red(os.path.realpath(f))
                    newFile = pwutils.findFile(os.path.basename(f), searchDir, recursive=True)
                    if newFile:
                        print "  Found file %s, creating link..." % newFile
                        print pwutils.green("   %s -> %s" % (f, newFile))
                        pwutils.createAbsLink(newFile, f)
