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

MODE_WRITE = 'write'
MODE_READ = 'read'


def writeParticles(outputFile):
    i = 0
    line = sys.stdin.readline()
    
    while line:
    #for line in sys.stdin:
        objDict=json.loads(line)
        ###imgId, index, filename = line.split()
        if '_index' in objDict.keys():
            index = int(objDict['_index'] - 1)
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
                                             "mirror":0,  ####TODO: test flip
                                             "scale":1.0})
        imageData.read_image(filename, index)
        if transformation is not None:
            imageData.set_attr('xform.projection', transformation)
    
        imageData.write_image(outputFile, i, eman.EMUtil.ImageType.IMAGE_HDF, False)
        i += 1
        print "OK"
        sys.stdout.flush()
        line = sys.stdin.readline()
    print "DONE"
    
    
def readParticles(inputParts, inputCls, inputClasses, outputTxt):
    clsImgsEven = eman.EMData.read_images(inputCls+"_even.hdf")
    clsImgsOdd = eman.EMData.read_images(inputCls+"_odd.hdf")
    classesEven = eman.EMData.read_images(inputClasses+"_even.hdf")
    classesOdd = eman.EMData.read_images(inputClasses+"_odd.hdf")
    imgs = eman.EMData.read_images(inputParts)

    clsClassListEven = clsImgsEven[0]
    clsClassListOdd = clsImgsOdd[0]
    
    clsClassList = {}
    shiftXList = {}
    shiftYList = {}
    dAlphaList = {}
    flipList = {}
    
    for i in range(clsClassListEven.get_attr_dict()['ny']):
        clsClassList[2*i] = clsClassListEven[0,i]
        shiftXList[2*i] = clsImgsEven[2][0,i]
        shiftYList[2*i] = clsImgsEven[3][0,i]
        dAlphaList[2*i] = clsImgsEven[4][0,i]
        flipList[2*i] = clsImgsEven[5][0,i]
    
    for i in range(clsClassListOdd.get_attr_dict()['ny']):
        clsClassList[2*i+1] = clsClassListOdd[0,i]
        shiftXList[2*i+1] = clsImgsOdd[2][0,i]
        shiftYList[2*i+1] = clsImgsOdd[3][0,i]
        dAlphaList[2*i+1] = clsImgsOdd[4][0,i]
        flipList[2*i+1] = clsImgsOdd[5][0,i]
    
    f = open(outputTxt, 'w')
#     index = 1
    
    for index in range(len(imgs)):
        classNum = clsClassList[index]
        if index%2 == 0:
            classes = classesEven
        else:
            classes = classesOdd
        az = classes[int(classNum)].get_attr_dict()['xform.projection'].get_rotation("eman")['az']
        alt = classes[int(classNum)].get_attr_dict()['xform.projection'].get_rotation("eman")['alt']
        
        if flipList[index] == 1:
            transform = eman.Transform({"type": "eman",
                                        "az": az + 180,
                                        "alt": 180 - alt,
                                        "phi": dAlphaList[index]*-1,
                                        "tx": shiftXList[index],
                                        "ty": shiftYList[index]
                                        })
        else:
            transform = eman.Transform({"type": "eman",
                                        "az": az,
                                        "alt": alt,
                                        "phi": dAlphaList[index],
                                        "tx": shiftXList[index], 
                                        "ty":shiftYList[index]
                                        })
        transform = transform.inverse()
        
        rot = transform.get_rotation("spider")['psi']
        tilt = transform.get_rotation("spider")['theta']
        psi = transform.get_rotation("spider")['phi']
        shiftX = transform.get_trans()[0]
        shiftY = transform.get_trans()[1]
        
        print >> f, index, rot, tilt, psi, shiftX, shiftY
    f.close()
    
    
if __name__ == '__main__':
    if len(sys.argv) > 1:
        mode = sys.argv[1]
        
        if mode == MODE_WRITE:
            outputFile = sys.argv[2]
            writeParticles(outputFile)
        elif mode == MODE_READ:
            inputParts = sys.argv[2]
            inputCls = sys.argv[3]
            inputClasses = sys.argv[4]
            outputTxt = sys.argv[5]
            readParticles(inputParts, inputCls, inputClasses, outputTxt)
        else:
            raise Exception("e2converter: Unknown mode '%s'" % mode)
    else:
        print "usage: %s outputFile" % os.path.basename(sys.argv[0])

