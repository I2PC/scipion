#!/usr/bin/env python

"""
This script will simply create symbolic links from files in an input
folder. The links will be created in an output folder and every some
time the last link will be 'touched' simulating an image that is being
acquired by the microscope.
"""

import sys, os
import time
import pyworkflow.utils as pwutils


def usage(error):
    print """
    ERROR: %s
    
    Usage: simulate_acquisition.py INPUT_FOLDER OUTPUT_FOLDER
        INPUT_FOLDER: input folder from where to read input files
        OUTPUT_FOLDER: where to create the output links.
        time: the time that the last link will be touched.        
    """ % error
    sys.exit(1)    


if len(sys.argv) < 3:
    usage("Incorrect number of input parameters")

inputDir = sys.argv[1]
outputDir = sys.argv[2]

inputFiles = os.listdir(inputDir)
inputFiles.sort() 

print "Cleaning output directory: ", outputDir
pwutils.cleanPath(outputDir)
pwutils.makePath(outputDir)

aTime = 30
n = 5
t = aTime / n
print "t=%s" % t

for f in inputFiles:
	inputPath = os.path.join(inputDir, f)
	outputPath = os.path.join(outputDir, f)
	print "Linking %s -> %s" % (outputPath, inputPath)
	
	for i in range(n):		
		time.sleep(t)		
		pwutils.copyFile(inputPath, outputPath)
		
	
	
