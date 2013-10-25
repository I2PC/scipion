#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

# use: python eigendoc.py cor_EIG.dat > eigenvalues.dat
#
#    2012-01-26 -- Prints only the lines for the factors
#    2012-01-26 -- Adapted from JSL

import string,sys

length = len(sys.argv[1:])

if length == 0:
    print "Usage: python eigendoc.py PREFIX_EIG_file"
    sys.exit(0)

file = sys.argv[1]

fp = open(file,'r')
B = fp.readlines()
fp.close()

firstline = string.split(B[0])
numfactors = int(firstline[0])
#print numfactors, type(numfactors)

B = B[1:numfactors+1]  # drop the 1st line

k = 1
for line in B:
    s = string.split(line)
    percent = string.atof(s[1])
    print "%5d 1%12f" % (k, percent)
    k += 1
