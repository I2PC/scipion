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
This script should be launched using the EMAN2 python interpreter 
and its own environment.
This scripts will convert any SetOfImages to an EMAN2 .hdf stack
It will read from the stdin the number of images and then
for each image will read: index, filename and a possible trasformation.
As parameters will receive the output filename for the hdf stack
"""

import os, sys

if __name__ == '__main__':
    if len(sys.argv) > 1:
        outputFile = sys.argv[1]
        
        #print 'PATH', os.environ['PATH'], '\n'
        #print 'PYTHONPATH', os.environ['PYTHONPATH'], '\n'
        #print 'LD_LIBRARY_PATH', os.environ['LD_LIBRARY_PATH']
            
        from EMAN2 import EMData, EMUtil
        
        i = 0
        line = sys.stdin.readline()
        while line:
        #for line in sys.stdin:
            imgId, index, filename = line.split()
            if index: # NO_INDEX is zero, otherwise remove one due EMAN2 is zero based index
                index = int(index) - 1
            imageData = EMData(filename, index, False)
            imageData['item_id'] = imgId
            imageData.write_image(outputFile, i, EMUtil.ImageType.IMAGE_HDF, False)
            i += 1
            print "EMAN2: ", line.strip()
            sys.stdout.flush()
            line = sys.stdin.readline()
        print "DONE"
    else:
        print "usage: %s outputFile" % os.path.basename(sys.argv[0])


    
