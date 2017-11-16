#!/usr/bin/env python

"""
This script will simply create symbolic links from files in an input
folder. The links will be created in an output folder and every some
time the last link will be 'touched' simulating an image that is being
acquired by the microscope.
"""
# target = /home/scipionuser/OffloadData/2017_11_09_administrator_test_sample/GRID_01/DATA/Images-Disc1/GridSquare_27783809/Data
#origin = ~/data/md1/01 
import sys, os
import time
import shutil
from glob import glob

import contextlib


def usage(prog, error):
    print """
    ERROR: %s
    
    Usage: simulate_acquisition.py INPUT_PATTERN OUTPUT_FOLDER
        INPUT_PATTERN: input pattern matching input files.
        OUTPUT_FOLDER: where to create the output links.
    example: %s "/home/scipionuser/OffloadData/Data/01/*.mrc" /home/scipionuser/OffloadData/2017_11_09_administrator_test_sample/GRID_01/DATA/Images-Disc1/GridSquare_27783809/Data
    """ % (error, prog)
    sys.exit(1)    


if len(sys.argv) != 3:
    usage(sys.argv[0], "Incorrect number of input parameters")

inputPattern = sys.argv[1]
outputDir = sys.argv[2]

inputFiles = glob(inputPattern)
inputFiles.sort()

print "Cleaning output directory: ", outputDir
if os.path.exists(outputDir):
    shutil.rmtree(outputDir)
os.makedirs(outputDir)

aTime = 30  # sec
n = 5  # number of touchs 
t = aTime / n  # touch image each aTime/n sec

for f in inputFiles:
    outputPath = os.path.join(outputDir, os.path.basename(f))
    print "Linking %s -> %s" % (f -> outputPath)

    for i in range(n):
        open(outputPath, 'w').close()
        time.sleep(t)
    try:
        os.remove(outputPath)  # remove file if exists        
    except:
        pass
    os.symlink(f, outputPath)

