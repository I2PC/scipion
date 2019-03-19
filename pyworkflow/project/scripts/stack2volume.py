#!/usr/bin/env python

import sys, os
from glob import glob

from pyworkflow.em import Ccp4Header
from pyworkflow.project import Manager
import pyworkflow.utils as pwutils
import pyworkflow.em as em


def usage(error):
    print """
    ERROR: %s
    
    stack2Volume will swap the dimension in the header of stack to make them 
    volumes. Something like 10 x 1 x 10 x 10 will be change to 1 x 10 x 10 x 10
    Usage: stack2volume.py PATH
        PATH: path to look for stack files    
    """ % error
    sys.exit(1)    

if len(sys.argv) != 2:
    usage("Incorrect number of input parameters")

path = sys.argv[1]

print ("Looking for files like: %s" % path)


for file in glob(path):

    print ("Changing header of %s" % file)
    try:
        header = Ccp4Header(file, readHeader=True)
        # Flag it as volume.
        header.setISPG(401)
        header.writeHeader()

    except Exception as e:
        print(pwutils.red("Failed to change header: % s" % e))
