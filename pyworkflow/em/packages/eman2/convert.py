# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Laura del Cano (ldelcano@cnb.csic.es)
# *              Josue Gomez Blanco (ldelcano@cnb.csic.es)
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
This module contains converter functions that will serve to:
1. Write from base classes to Eman specific files
2. Read from Eman files to base classes
"""

import os, glob
import json
import numpy
import subprocess

from os.path import join, exists

import pyworkflow as pw
import pyworkflow.em as em
from pyworkflow.em.data import Coordinate
from pyworkflow.em.packages.eman2 import getEmanCommand, getEnviron
from pyworkflow.utils.path import createLink, removeBaseExt, replaceBaseExt, cleanPath


# LABEL_TYPES = { 
#                xmipp.LABEL_SIZET: long,
#                xmipp.LABEL_DOUBLE: float,
#                xmipp.LABEL_INT: int,
#                xmipp.LABEL_BOOL: bool              
#                }

#     @static
def loadJson(jsonFn):
    """ This function loads the Json dictionary into memory """
    jsonFile = open(jsonFn)
    jsonDict = json.load(jsonFile)
    jsonFile.close()
    return jsonDict

 
def writeJson(jsonDict, jsonFn):
    """ This function write a Json dictionary """
    with open(jsonFn, 'w') as outfile:
        json.dump(jsonDict, outfile)


def objectToRow(obj, row, attrDict):
    pass


def _rowToObject(row, obj, attrDict):
    pass

    
def rowToObject(md, objId, obj, attrDict):
    pass

    
def readCTFModel(filename):
    pass


def writeCTFModel(ctfObj, filename):
    """ Write a CTFModel object as Xmipp .ctfparam"""
    pass


def readMicrograph(md, objId):
    """ Create a Micrograph object from a row of Xmipp metadata. """
    pass


def locationToEman(index, filename):
    pass


def micrographToRow(mic, micRow):
    pass


def rowToCoordinate(md, objId):
    """ Create a Coordinate from a json. """
    pass


def readSetOfMicrographs(filename):
    pass


def writeSetOfMicrographs(micSet, filename, rowFunc=None):
    pass

    
def readPosCoordinates(posFile):
    pass
            

def readSetOfCoordinates(workDir, micSet, coordSet):
    """ Read from Eman .json files.
    It is expected a file named: base.json under the workDir.
    Params:
        workDir: where the Eman boxer output files are located.
        micSet: the SetOfMicrographs to associate the .json, which 
            name should be the same of the micrographs.
        coordSet: the SetOfCoordinates that will be populated.
    """
    # Read the boxSize from the e2boxercache/base.json
    jsonFnbase = join(workDir, 'e2boxercache', 'base.json')
    jsonBoxDict = loadJson(jsonFnbase)
    size = int(jsonBoxDict["box_size"])
    jsonFninfo = join(workDir, 'info/')
    
    for mic in micSet:
        micBase = removeBaseExt(mic.getFileName())
        micPosFn = ''.join(glob.glob(jsonFninfo + '*' + micBase + '_info.json'))
        readCoordinates(mic, micPosFn, coordSet)
    coordSet.setBoxSize(size)


def readCoordinates(mic, fileName, coordsSet):
    if exists(fileName):
        jsonPosDict = loadJson(fileName)

        if jsonPosDict.has_key("boxes"):
            boxes = jsonPosDict["boxes"]

            for box in boxes:
                x, y = box[:2]
                coord = Coordinate()
                coord.setPosition(x, y)
                coord.setMicrograph(mic)
                coordsSet.append(coord)
        


def writeSetOfCoordinates():
    pass


def createEmanProcess(script='e2converter.py', args=None, direc="."):
    """ Open a new Process with all EMAN environment (python...etc)
    that will server as an adaptor to use EMAN library
    """
    program = pw.join('em', 'packages', 'eman2', script)        
    cmd = getEmanCommand(program, args)
    
#    gcmd = greenStr(cmd)
    print "** Running: '%s'" % cmd
    proc = subprocess.Popen(cmd, shell=True, env=getEnviron(), 
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE, cwd=direc)
    
    return proc
    
    
def writeSetOfParticles(partSet, path, **kwargs):
    """ Convert the imgSet particles to a single .hdf file as expected by Eman. 
    This function should be called from a current dir where
    the images in the set are available.
    """
#     from math import fabs
    pathMrcs = path.split("/particles")[0]
    isMrcs, fnDict = convertBinaryFiles(partSet, pathMrcs)
    
    
    partSet.getFirstItem().getFileName()
    fileName = ""
    a = 0
#     listHdf = []
    proc = createEmanProcess(args='write')
    for i, part in iterParticlesByMic(partSet):
        micId = part.getMicId()
        objDict = part.getObjDict()
        
        if not micId:
            micId = 0
        
        objDict['hdfFn'] = os.path.join(path, "mic_%0.6d.hdf" % micId)
#             listHdf.append(basename(hdfFn))
        
        alignType = kwargs.get('alignType')
#         print "objDict, ", objDict
    
        if alignType != em.ALIGN_NONE:
            shift, angles = alignmentToRow(part.getTransform(), alignType)
            #json cannot encode arrays so I convert them to lists
            #json fail if has -0 as value
            objDict['_shifts'] = shift.tolist()
            objDict['_angles'] = angles.tolist()
#         fn = objDict['_filename']
        # check if the is a file mapping
#         objDict['_filename'] = filesDict.get(fn, fn)
        objDict['_itemId'] = part.getObjId()
        if isMrcs:
            objDict['_filename'] = fnDict[objDict['_filename']]
        
        # the index in EMAN begin in 0.
        if fileName != objDict['_filename']:
            fileName = objDict['_filename']
            if objDict['_index'] == 0:
                a = 0
            else:
                a = 1
        objDict['_index'] = int(objDict['_index'] - a)
        # Write the e2converter.py process from where to read the image
        print >> proc.stdin, json.dumps(objDict)
        proc.stdin.flush()
        proc.stdout.readline()
    proc.kill()


def readSetOfParticles(filename, partSet, **kwargs):
    pass


def geometryFromMatrix(matrix, inverseTransform):
    from pyworkflow.em.transformations import translation_from_matrix, euler_from_matrix
    if inverseTransform:
        from numpy.linalg import inv
        matrix = inv(matrix)
        shifts = -translation_from_matrix(matrix)
    else:
        shifts = translation_from_matrix(matrix)
    rad_to_ang = 180./numpy.pi
    angles = -numpy.rad2deg(euler_from_matrix(matrix, axes='szyz'))
    return shifts, angles


def matrixFromGeometry(shifts, angles, inverseTransform):
    """ Create the transformation matrix from a given
    2D shifts in X and Y...and the 3 euler angles.
    """
    from pyworkflow.em.transformations import euler_matrix
    from numpy import deg2rad
    radAngles = -deg2rad(angles)

    M = euler_matrix(radAngles[0], radAngles[1], radAngles[2], 'szyz')
    if inverseTransform:
        from numpy.linalg import inv
        M[:3, 3] = -shifts[:3]
        M = inv(M)
    else:
        M[:3, 3] = shifts[:3]

    return M


def alignmentToRow(alignment, alignType):
    """
    is2D == True-> matrix is 2D (2D images alignment)
            otherwise matrix is 3D (3D volume alignment or projection)
    invTransform == True  -> for xmipp implies projection
                          -> for xmipp implies alignment
    """
#     is2D = alignType == em.ALIGN_2D
#     inverseTransform = alignType == em.ALIGN_PROJ
    
    #tranformation matrix is procesed here because
    #it uses routines available thrrough scipion python
    matrix = alignment.getMatrix()
    return geometryFromMatrix(matrix, True)


def rowToAlignment(alignmentList, alignType):
    """
    is2D == True-> matrix is 2D (2D images alignment)
            otherwise matrix is 3D (3D volume alignment or projection)
    invTransform == True  -> for xmipp implies projection
        """
    is2D = alignType == em.ALIGN_2D
    inverseTransform = alignType == em.ALIGN_PROJ
     
    alignment = em.Transform()
    angles = numpy.zeros(3)
    shifts = numpy.zeros(3)
    shifts[0] = alignmentList[3]
    shifts[1] = alignmentList[4]
    if not is2D:
        angles[0] = alignmentList[0]
        angles[1] = alignmentList[1]
        shifts[2] = 0
        angles[2] = alignmentList[2]
    else:
        psi = alignmentList[0]
        rot = alignmentList[1]
        if rot !=0. and psi !=0:
            print "HORROR rot and psi are different from zero"
        angles[0] = alignmentList[0] + alignmentList[1]
    #if alignment
    matrix = matrixFromGeometry(shifts, angles, inverseTransform)
    alignment.setMatrix(matrix)
     
    #FIXME: now are also storing the alignment parameters since
    # the conversions to the Transform matrix have not been extensively tested.
    # After this, we should only keep the matrix 
    #for paramName, label in ALIGNMENT_DICT.iteritems():
    #    if alignmentRow.hasLabel(label):
    #        setattr(alignment, paramName, alignmentRow.getValueAsObject(label))    
     
    return alignment


def convertBinaryFiles(imgSet, outputDir):
    """ Convert binary images files to a format read by EMAN.
    Params:
        imgSet: input image set to be converted.
        outputDir: where to put the converted file(s)
    Return:
        A dictionary with old-file as key and new-file as value
        If empty, not conversion was done.
    """
    filesDict = {}
    ih = em.ImageHandler()
    # This approach can be extended when
    # converting from a binary file format that
    # is not read from Relion
    def linkMrcToMrcs(fn):
        """ Just create a link named .mrc to Eman understand 
        that it is a mrc binary stack.
        """
        newFn = join(outputDir, replaceBaseExt(fn, 'mrcs'))
        createLink(fn, newFn)
        return newFn
        
    def convertStack(fn):
        """ Convert from a format that is not read by Relion
        to an spider stack.
        """
        newFn = join(outputDir, replaceBaseExt(fn, 'stk'))
        ih.convertStack(fn, newFn)
        return newFn
        
    ext = imgSet.getFirstItem().getFileName()
    if ext.endswith('.mrc'):
        mapFunc = linkMrcToMrcs
        ismrcs = True
#     elif ext.endswith('.hdf'): # assume eman .hdf format
#         mapFunc = convertStack
    else:
        mapFunc = None
        ismrcs = False
        
    if mapFunc is not None:
        for fn in imgSet.getFiles():
            filesDict[fn] = mapFunc(fn) # convert and map new filename

    return ismrcs, filesDict


def iterParticlesByMic(partSet):
    """ Iterate the particles ordered by micrograph """
    for i, part in enumerate(partSet.iterItems(orderBy=['_micId', 'id'],
                                                                 direction='ASC')):
        yield i, part
