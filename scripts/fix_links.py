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

project.fixLinks(searchDir)
