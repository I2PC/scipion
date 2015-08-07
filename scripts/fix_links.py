#!/usr/bin/env python

import sys, os
from os.path import exists, join


def usage(error):
    print """
    ERROR: %s
    
    Usage: fixlinks.py LINKS_DIR SEARCH_DIR
        LINKS_DIR: provide the folder where the broken links are.
        SEARCH_DIR: provide a directory where to look for the files
        and fix the links.    
    """ % error
    sys.exit(1)    


if len(sys.argv) != 3:
    usage("Incorrect number of input parameters")

linksDir = sys.argv[1]
searchDir = sys.argv[2]
    
if not exists(linksDir):
    usage("Unexistent LINKS_DIR: %s" % linksDir)
    
if not exists(searchDir):
    usage("Unexistent SEARCH_DIR: %s" % searchDir)


# Find the broken links under the LINKS_DIR

files = os.listdir(linksDir)
brokenLinks = []

for f in files:
    fn = join(linksDir, f)
    realPath = os.path.realpath(fn)
    if (os.path.islink(fn) and 
        not exists(realPath)):
        brokenLinks.append((fn, realPath.split(os.sep))) # store as list 
        
print "Found %d broken links" % len(brokenLinks)

if len(brokenLinks) == 0:
    print "  No work to do. EXITING..."
    sys.exit(0)
    
first = brokenLinks[0][1]

for i in range(1, len(first)):
    print "Trying from path: %s" % join(searchDir, os.sep.join(first[i:]))
    if all(exists(join(searchDir, os.sep.join(parts[i:]))) for _, parts in brokenLinks):
        # We found all links from a common root, let's FIX the links
        print "FOUND ALL LINKS!!! Let's fix them. "
        for fn, parts in brokenLinks:
            os.remove(fn)
            os.symlink(os.path.abspath(join(searchDir, os.sep.join(parts[i:]))), fn)
        sys.exit(0)
        
print "ERROR: Links where not fixed."
