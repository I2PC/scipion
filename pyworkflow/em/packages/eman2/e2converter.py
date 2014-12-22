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
from EMAN2 import Transform

def geometryFromMatrix(matrix, inverseTransform):
    from pyworkflow.em.transformations import translation_from_matrix, euler_from_matrix
    if inverseTransform:
        from numpy.linalg import inv
        matrix = inv(matrix)
        shifts = -translation_from_matrix(matrix)
    else:
        shifts = translation_from_matrix(matrix)
    rad_to_ang = 180./np.pi
    angles = -1.*rad_to_ang* euler_from_matrix(matrix, axes='szyz')
    return shifts, angles

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
            objDict=json.loads(line)
            ###imgId, index, filename = line.split()
            print >> sys.stderr,"objDict", objDict
            if '_index' in objDict.keys():
                index = int(objDict['_index']) - 1
            if '_filename' in objDict.keys():
                filename = str(objDict['_filename'])
            else:
                raise Exception('ERROR (e2converter): Cannot process a particle without filename')
            imageData = EMData()
            t= None
            if '_transform._matrix' in objDict.keys():
                transform_matrix = np.matrix(str(objDict['_transform._matrix']))
                print >> sys.stderr, "transform_matrix", transform_matrix[0,0]
                shifts, angles = geometryFromMatrix(transform_matrix, False)
                t = Transform({"type":"spider",
                               "phi":angles[0],
                               "theta":angles[1],
                               "psi":angles[2],
                               "tx":shifts[0],
                               "ty":shifts[1],
                               "tz":shifts[2],
                               "mirror":0,####TODO: test flip
                               "scale":1.0})

            #EMData.read_image(imageData, filename,index)
            imageData.read_image(filename,index)
            if t is not None:
                imageData.set_attr('xform.projection', t)

#            if '_itemId' in objDict.keys():
#                itemId = objDict['_itemId']
#            else:
#                raise Exception('ERROR (e2converter): Cannot process a particle without _itemId')
#            imageData['item_id'] = itemId

            imageData.write_image(outputFile, i, EMUtil.ImageType.IMAGE_HDF, False)
            i += 1
            print "EMAN2: "
            sys.stdout.flush()
            line = sys.stdin.readline()
        print "DONE"
    else:
        print "usage: %s outputFile" % os.path.basename(sys.argv[0])


