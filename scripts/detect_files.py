#!/usr/bin/env python

"""
This script will detect new files that are ben output
in a folder. We will wait for some time that a file
have not changed to said is have been output.
"""

import sys, os
import time
from datetime import timedelta, datetime
from glob import glob

import pyworkflow.utils as pwutils


def usage(error):
    print """
    ERROR: %s
    
    Usage: detect_files.py INPUT_FOLDER 
        INPUT_FOLDER: input folder from where to read input files
    """ % error
    sys.exit(1)    
    
    
class FilesWatcher(object):
    """ Class to check a folder for new files. """
    
    def __init__(self, filesPattern, fileWaitTime=30, waitTime=300, 
                 foundFiles={}):
        """ Initialize the FilesWatcher 
        Params:
            filesPattern: the pattern of files to match (using glob).
            fileWaitTime: how long to wait to mark a file of no-changed
            waitTime: how long to wait for files before stop.
            foundFiles: pass some previous found files.
        """
        self.filesPattern = filesPattern
        self.fileWait = timedelta(seconds=fileWaitTime)
        self.wait = timedelta(seconds=waitTime)        
        self.foundFiles = dict(foundFiles)
        self.finished = False
        self.startTime = datetime.now()
        
    def getNewFiles(self):
        """ Check if there are new files matching the pattern.
        Return an empty list if there is no files but some were found.
        Return None if no more files found anymore.
        """
        if self.finished:
            return None
            
        inputFiles = glob(self.filesPattern)
        newFiles = []
        someNew = False
        
        for f in inputFiles:            
            if f not in self.foundFiles:
                someNew = True
                mTime = datetime.fromtimestamp(os.path.getmtime(f))
                delta = datetime.now() - mTime
                if pwutils.envVarOn('SCIPION_DEBUG'):
                    print "Checking file: '%s' (%s) wait (%s)" % (f, delta.seconds, self.fileWait.seconds)
                if delta > self.fileWait:
                    newFiles.append(f)
                    self.foundFiles[f] = True         

        if not someNew and datetime.now() - self.startTime > self.wait:
            self.finished = True
            return None
        
        return newFiles
        

if len(sys.argv) < 2:
    usage("Incorrect number of input parameters")

inputDir = sys.argv[1]


print "Monitoring directory: ", inputDir

# 60 seconds to wait for files changes
wTime = 20 
# 5 mins to wait for new files
nTime = 60

fwatcher = FilesWatcher(os.path.join(inputDir, '*'), wTime, nTime)

while True:
    print "=" * 70
    newFiles = fwatcher.getNewFiles()

    if newFiles is None:
        print "Waited enough time, QUITTING..."
        break
        
    if newFiles: # non empty result
        print "New files: "
        for nf in newFiles:
            print " - '%s'" % nf
    
    time.sleep(5) 
