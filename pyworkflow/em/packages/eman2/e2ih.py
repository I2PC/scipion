# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
This script should be launched using the EMAN2 python interpreter 
and its own environment.
This scripts will be used from ImageHandler to manipulate image
format that Xmipp does not support right now.
This way should be deprecated when we read more formats.
"""

import os
import sys
import EMAN2 as eman


if __name__ == '__main__':
    n = len(sys.argv)

    if 1 < n < 4:
        inputFile = sys.argv[1]

        if n > 2:  # convert from input to output
            outputFile = sys.argv[2]
            emd = eman.EMData()
            if '@' in inputFile or '@' in outputFile:
                def _getLoc(fn):
                    if '@' in fn:
                        i, fn = fn.split('@')
                        return [fn, int(i) - 1]  # 0-index in eman
                    return [fn, 0]
                emd.read_image(*_getLoc(inputFile))
                emd.write_image(*_getLoc(outputFile))
            else:
                nsize = eman.EMUtil.get_image_count(inputFile)
                for i in range(nsize):
                    emd.read_image(inputFile, i)
                    emd.write_image(outputFile, i)
        else:  # just print dimensions
            nsize = eman.EMUtil.get_image_count(inputFile)
            emd = eman.EMData(inputFile, 0, True)
            print emd.get_xsize(), emd.get_ysize(), emd.get_zsize(), nsize
    else:
        print "usage: %s inputFile [outputFile]" % os.path.basename(sys.argv[0])
