#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
Generate a Set with file-md5sum information for all files
under a given directory.
"""

import sys
import os
from os.path import join, exists, relpath
import hashlib


def main():
    datasets = sys.argv[1:-1]
    datasetPath = sys.argv[-1]
    for path in datasets:
        print "Formatting %s" % path
        absPath = join(datasetPath, path)
        if not exists(absPath):
            sys.exit("ERROR: %s folder doesn't exist in datasets folder %s." %
                     (path, datasetPath))
        with open(join(absPath, 'MANIFEST'), 'w+') as manifest:
            for root, dirs, files in os.walk(absPath):
                for filename in files:
                    if filename == 'MANIFEST':
                        continue  # do not include ourselves
                    fn = join(root, filename)  # file to check
                    manifest.write('%s %s\n' % (relpath(fn, absPath), md5(fn)))


def md5(fname):
    """ Return the md5 hash of file fname """
    mhash = hashlib.md5()
    with open(fname) as f:
        for chunk in iter(lambda: f.read(128 * mhash.block_size), ''):
            mhash.update(chunk)
    return mhash.hexdigest()



if __name__ == '__main__':
    main()
