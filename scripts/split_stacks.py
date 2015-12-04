#!/usr/bin/env python

"""
Split one or several stack files into other stacks with less elements 
or individual images.

Examples:

- Extract individual .mrc micrographs from a movie .mrc file:
scipion python scripts/split_stacks.py --files movie.mrcs --ext mrc 

- Regroup a movie into stacks of 4 elements:
scipion python scripts/split_stacks.py --files movie.mrcs --ext mrcs -n 4

- Same as before but for all movies in the current dir and output in spider format:
scipion python scripts/split_stacks.py --files "movie*.mrcs" --ext stk -n 4
  
"""
import glob
import os
import argparse
from pyworkflow.em import ImageHandler


class Main():

    def getNumOfElements(self, fileN):
        _, _, z, n = self.imgh.getDimensions(fileN)
        if z > 1:
            return z
        else:
            return n
        
    def getOutputLoc(self, filePrefix, i):
        if self.n == 1:
            index = 1
            fn = "%s_%03d.%s" % (filePrefix, i, self.args.ext)
        else:
            index = self.outIndex
            fn = "%s_%03d.%s" % (filePrefix, self.groupCount, self.args.ext)
            self.outIndex += 1
            if self.outIndex > self.n:
                self.outIndex = 1
                self.groupCount += 1
                
        return (index, fn)
        
    def __init__(self):
        parser = argparse.ArgumentParser(description="Extract images from a stack and write them into "
                                                     "smaller stacks or even individual images.")
        add = parser.add_argument  # shortcut
        add('--files', default='',
            help='Pattern to match input files.')
        add('--ext', default='mrc',
            help='Define extension of output files.')
        add('-n', default=1,
            help='Group output images in stacks of this size.')
            
        self.args = parser.parse_args()
        
        self.n = int(self.args.n)            
        filePaths = glob.glob(self.args.files)
        filePaths.sort()
        self.imgh = ImageHandler()
        
        for fileN in filePaths:
            self.outIndex = 1
            self.groupCount = 1
            n = self.getNumOfElements(fileN)
            
            if fileN.endswith('.mrc'):
                fileN += ':mrcs'
            
            filePrefix = os.path.splitext(fileN)[0]
            
            for i in range(1, n+1):
                outputLoc = self.getOutputLoc(filePrefix, i)
                self.imgh.convert((i, fileN), outputLoc)
    

if __name__ == "__main__":
    Main()
