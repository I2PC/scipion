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
import json
import numpy as np
import EMAN2 as eman


if __name__ == '__main__':
    if len(sys.argv) > 1:
        outputFile = sys.argv[1]
        
        #print 'PATH', os.environ['PATH'], '\n'
        #print 'PYTHONPATH', os.environ['PYTHONPATH'], '\n'
        #print 'LD_LIBRARY_PATH', os.environ['LD_LIBRARY_PATH']
            
        i = 0
        line = sys.stdin.readline()
        while line:
        #for line in sys.stdin:
            objDict=json.loads(line)
            ###imgId, index, filename = line.split()
            print >> sys.stderr,"objDict", objDict
            if '_index' in objDict.keys():
                index = int(objDict['_index']) - 1
            if '_filename' in objDict.keys():
                filename = str(objDict['_filename'])
            else:
                raise Exception('ERROR (e2converter): Cannot process a particle without filename')
            imageData = eman.EMData()
            transformation = None
            if '_angles' in objDict.keys():
                #TODO: convert to vector not matrix
                angles = objDict['_angles']
                shifts = objDict['_shifts']
                transformation = eman.Transform({"type":"spider",
                               "phi":angles[0],
                               "theta":angles[1],
                               "psi":angles[2],
                               "tx":shifts[0],
                               "ty":shifts[1],
                               "tz":shifts[2],
                               "mirror":0,####TODO: test flip
                               "scale":1.0})

            imageData.read_image(filename, index)
            if transformation is not None:
                imageData.set_attr('xform.projection', transformation)

#            if '_itemId' in objDict.keys():
#                itemId = objDict['_itemId']
#            else:
#                raise Exception('ERROR (e2converter): Cannot process a particle without _itemId')
#            imageData['item_id'] = itemId

            imageData.write_image(outputFile, i, eman.EMUtil.ImageType.IMAGE_HDF, False)
            i += 1
            print "OK"
            sys.stdout.flush()
            line = sys.stdin.readline()
        print "DONE"
    else:
        print "usage: %s outputFile" % os.path.basename(sys.argv[0])


