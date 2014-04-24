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

import os, sys
import hashlib
from os.path import join
from pyworkflow.utils import getExt

if __name__ == '__main__':
    datasets = sys.argv[1:-1]
    datasetPath = sys.argv[-1]
    for path in datasets:
        print "Formatting %(dataset)s" % ({'dataset': path})
        absPath = os.path.join(datasetPath, path)
        if not os.path.exists(absPath):
            print "ERROR: %(path)s folder doesn't exist in datasets folder %(datasetsFolder)s." % ({'path': path, 'datasetsFolder': datasetPath})
            sys.exit(1)
        manifestPath = join(absPath, 'MANIFEST')
        manifestFile = open(manifestPath, 'w+')
        for root, dirs, files in os.walk(absPath):
            for filename in files:
                fn = join(root, filename)
                if fn != manifestPath:
                    md5sum = 0
                    md5 = hashlib.md5()
                    with open(fn, 'r+') as fileToCheck:
                        for chunk in iter(lambda: fileToCheck.read(128*md5.block_size), b''):
                            md5.update(chunk)
                    md5sum = md5.hexdigest()
                    manifestFile.write('%(fn)s %(md5sum)s\n' % ({'md5sum': md5sum, 'fn': os.path.relpath(fn, absPath)}))
        manifestFile.close()
