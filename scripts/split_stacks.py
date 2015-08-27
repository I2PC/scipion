#!/usr/bin/env python

import glob
import os
import argparse
from pyworkflow.em import ImageHandler


def main():
    parser = argparse.ArgumentParser(description="make individual micrographs from a stack")
    add = parser.add_argument  # shortcut
    add('--files', default='',
        help='list of files to split on individual micrographs')
    add('--ext', default='mrc',
        help='define extension of individual micrographs')
    args = parser.parse_args()

    def getNumOfElements(fileN):
        _, _, z, n = imgh.getDimensions(fileN)
        if z > 1:
            return z
        else:
            return n
    
    filePaths = glob.glob(args.files)
    filePaths.sort()
    imgh = ImageHandler()
    
    for fileN in filePaths:
        n = getNumOfElements(fileN)
        if fileN.endswith('.mrc'):
            fileN += ':mrcs'
        for i in range(1, n+1):
            outputFn = "%s_%03d.%s" % (os.path.splitext(fileN)[0], i, args.ext)
            imgh.convert((i, fileN), outputFn)


if __name__ == "__main__":
    main()
