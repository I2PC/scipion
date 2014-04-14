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

if __name__ == '__main__':
    path = sys.argv[1]
    
    md5Fn = join(path, 'data.md5')
    f = open(md5Fn, 'w+')
    
    for root, dirs, files in os.walk(path):
        for filename in files:
            fn = join(root, filename)
            if fn != md5Fn:
                with open(fn) as fileToCheck:
                    data = fileToCheck.read()
                    md5sum = hashlib.md5(data).hexdigest()
                f.write('%s %s\n' % (md5sum, fn))

    f.close()                

        
    