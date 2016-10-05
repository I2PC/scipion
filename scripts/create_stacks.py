#!/usr/bin/env python

"""
Create movie stacks from the individual frame files.
"""

import glob
import os
import re
import argparse

from pyworkflow.em import ImageHandler
import pyworkflow.utils as pwutils


class Main():
    def __init__(self):
        parser = argparse.ArgumentParser(
            description="Create movie stacks from the individual "
                        "frame files.")
        add = parser.add_argument  # shortcut
        add('--files', default='',
            help='Pattern to match input frame files.')        
        add('-n', 
            help='Number of frames per movie.')
        add('--suffix',
            help='Provide suffix added to create movie file. '
                 'e.g. _frames.mrcs')

        add('--delete_frames', action='store_true',
            help='Provide this option if you want to delete individual frame '
                 'files after the movie stack is created. ')
            
        args = parser.parse_args()
        n = int(args.n)

        frameRegex = re.compile("(?P<prefix>.+[^\d]+)(?P<frameid>\d+)")
        #  Group all frames for each movie
        # Key of the dictionary will be the common prefix and the value
        # will be a list with all frames in that movie
        frameDict = {}
        filePaths = glob.glob(args.files)
        filePaths.sort()

        for fileName in filePaths:
            fnNoExt = pwutils.removeExt(fileName)

            match = frameRegex.match(fnNoExt)

            if match is None:
                raise Exception("Incorrect match of frame files pattern!")

            d = match.groupdict()
            prefix = d['prefix']
            frameid = int(d['frameid'])

            if prefix not in frameDict:
                frameDict[prefix] = []

            frameDict[prefix].append((frameid, fileName))

        suffix = args.suffix
        ih = ImageHandler()

        for k, v in frameDict.iteritems():
            if len(v) != n:
                raise Exception("Incorrect number of frames!")

            movieFn = k + suffix
            movieOut = movieFn

            if movieOut.endswith("mrc"):
                movieOut += ":mrcs"

            print "Writing movie stack: ", movieFn

            for i, frame in enumerate(sorted(v, key=lambda x: x[0])):
                frameFn = frame[1] # Frame name stored previously
                ih.convert(frameFn, (i+1, movieFn))

                if args.delete_frames:
                    pwutils.cleanPath(frameFn)



if __name__ == "__main__":
    Main()
