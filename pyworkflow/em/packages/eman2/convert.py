# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              Laura del Cano (ldelcano@cnb.csic.es) [1]
# *              Josue Gomez Blanco (josue.gomez-blanco@mcgill.ca) [1]
# *              Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [2]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# * [2] MRC Laboratory of Molecular Biology (MRC-LMB)
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

import glob
import json
import numpy
import subprocess

import pyworkflow as pw
import pyworkflow.em as em
import pyworkflow.utils as pwutils
from pyworkflow.em.data import Coordinate
from eman2 import getEmanCommand, getEnviron


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


def rowToObject(md, objId, obj, attrDict):
    pass


def readCTFModel(ctfModel, filename):
    jsonDict = loadJson(filename)
    keyPos = None
    ctfPhaseShift = 0.0

    if jsonDict.has_key('ctf_frame'):
        keyPos = jsonDict['ctf_frame'][1]
    elif jsonDict.has_key('ctf'):
        keyPos = jsonDict['ctf'][0]
    else:
        setWrongDefocus(ctfModel)

    if keyPos:
        defocus = float(keyPos['defocus'])
        defocusAngle = float(keyPos['dfang'])
        dfdiff = float(keyPos['dfdiff'])
        ampcont = float(keyPos['ampcont'])
        defocusU = 10000.0 * defocus + 5000.0 * dfdiff
        defocusV = 20000.0 * defocus - defocusU

        # calculate phase shift as in EMAN2 ctf.cpp
        if ampcont > -100.0 and ampcont <= 100.0:
            PhaseShift = numpy.arcsin(ampcont / 100.0)
        elif (ampcont > 100.0):
            PhaseShift = numpy.pi - numpy.arcsin(2.0 - ampcont / 100.0)
        else:
            PhaseShift = -numpy.pi - numpy.arcsin(-2.0 - ampcont / 100.0)
        ctfPhaseShift = numpy.rad2deg(PhaseShift)

        ctfModel.setStandardDefocus(defocusU, defocusV, defocusAngle)
        if jsonDict.has_key('ctf_im2d'):
            # psdFile = jsonDict['ctf_im2d']['__image__'][0]
            fnBase = pwutils.removeExt(filename) + '_jsonimg'
            psdFile = "1@%s.hdf" % fnBase
            if pwutils.exists(psdFile):
                ctfModel.setPsdFile(psdFile)
    ctfModel.setPhaseShift(float(ctfPhaseShift))


def setWrongDefocus(ctfModel):
    ctfModel.setDefocusU(-999)
    ctfModel.setDefocusV(-1)
    ctfModel.setDefocusAngle(-999)


def writeCTFModel(ctfObj, filename):
    """ Write a CTFModel object as Xmipp .ctfparam"""
    pass


def jsonToCtfModel(ctfJsonFn, ctfModel):
    """ Create a CTFModel from a json file """
    mdFn = str(ctfJsonFn).replace('particles', 'info')
    mdFn = mdFn.split('__ctf_flip')[0] + '_info.json'
    if pwutils.exists(mdFn):
        readCTFModel(ctfModel, mdFn)


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


def readSetOfCoordinates(workDir, micSet, coordSet, invertY=False, newBoxer=False):
    """ Read from Eman .json files.
    Params:
        workDir: where the Eman boxer output files are located.
        micSet: the SetOfMicrographs to associate the .json, which
            name should be the same of the micrographs.
        coordSet: the SetOfCoordinates that will be populated.
    """
    if newBoxer:
        # read boxSize from info/project.json
        jsonFnbase = pwutils.join(workDir, 'info', 'project.json')
        jsonBoxDict = loadJson(jsonFnbase)
        size = int(jsonBoxDict["global.boxsize"])
    else:
        # read boxSize from e2boxercache/base.json
        jsonFnbase = pwutils.join(workDir, 'e2boxercache', 'base.json')
        jsonBoxDict = loadJson(jsonFnbase)
        size = int(jsonBoxDict["box_size"])

    jsonFninfo = pwutils.join(workDir, 'info/')

    for mic in micSet:
        micBase = pwutils.removeBaseExt(mic.getFileName())
        micPosFn = ''.join(glob.glob(jsonFninfo + '*' + micBase + '_info.json'))
        readCoordinates(mic, micPosFn, coordSet, invertY)
    coordSet.setBoxSize(size)


def readCoordinates(mic, fileName, coordsSet, invertY=False):
    if pwutils.exists(fileName):
        jsonPosDict = loadJson(fileName)

        if jsonPosDict.has_key("boxes"):
            boxes = jsonPosDict["boxes"]

            for box in boxes:
                x, y = box[:2]

                if invertY:
                    y = mic.getYDim() - y

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
    firstCoord = partSet.getFirstItem().getCoordinate() or None
    hasMicName = False
    if firstCoord:
        hasMicName = firstCoord.getMicName() or False

    fileName = ""
    a = 0
    proc = createEmanProcess(args='write')

    for i, part in iterParticlesByMic(partSet):
        micName = micId = part.getMicId()
        if hasMicName:
            micName = pwutils.removeBaseExt(part.getCoordinate().getMicName())
        objDict = part.getObjDict()

        if not micId:
            micId = 0

        suffix = kwargs.get('suffix', '')
        if hasMicName and (micName != str(micId)):
            objDict['hdfFn'] = pwutils.join(path, "%s%s.hdf" % (micName, suffix))
        else:
            objDict['hdfFn'] = pwutils.join(path, "mic_%06d%s.hdf" % (micId, suffix))

        alignType = kwargs.get('alignType')

        if alignType != em.ALIGN_NONE:
            shift, angles = alignmentToRow(part.getTransform(), alignType)
            # json cannot encode arrays so I convert them to lists
            # json fail if has -0 as value
            objDict['_shifts'] = shift.tolist()
            objDict['_angles'] = angles.tolist()
        objDict['_itemId'] = part.getObjId()

        # the index in EMAN begins with 0
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


def getImageDimensions(imageFile):
    """ This function will allow us to use EMAN2 to read some formats
     not currently supported by the native image library (Xmipp).
     Underneath, it will call a script to do the job.
    """
    proc = createEmanProcess('e2ih.py', args=imageFile)
    return tuple(map(int, proc.stdout.readline().split()))


def convertImage(inputLoc, outputLoc):
    """ This function will allow us to use EMAN2 to write some formats
     not currently supported by the native image library (Xmipp).
     Underneath, it will call an script to do the job.
    """

    def _getFn(loc):
        """ Use similar naming convention as in Xmipp.
        This does not works for EMAN out of here.
        """
        if isinstance(loc, tuple):
            if loc[0] != em.NO_INDEX:
                return "%06d@%s" % loc
            return loc[1]
        else:
            return loc

    proc = createEmanProcess('e2ih.py', args='%s %s' % (_getFn(inputLoc),
                                                        _getFn(outputLoc)))
    proc.wait()


def iterLstFile(filename):
    f = open(filename)
    for line in f:
        if '#' not in line:
            # Decompose Eman filename
            index, filename = int(line.split()[0]) + 1, line.split()[1]
            yield (index, filename)
    f.close()


def geometryFromMatrix(matrix, inverseTransform):
    from pyworkflow.em.transformations import translation_from_matrix, euler_from_matrix
    if inverseTransform:
        from numpy.linalg import inv
        matrix = inv(matrix)
        shifts = -translation_from_matrix(matrix)
    else:
        shifts = translation_from_matrix(matrix)
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

    # tranformation matrix is procesed here because
    # it uses routines available thrrough scipion python
    matrix = alignment.getMatrix()
    return geometryFromMatrix(matrix, True)


def rowToAlignment(alignmentList, alignType):
    """
    is2D == True-> matrix is 2D (2D images alignment)
            otherwise matrix is 3D (3D volume alignment or projection)
    invTransform == True  -> for xmipp implies projection
        """
    # use all angles in 2D since we might have mirrors
    # is2D = alignType == em.ALIGN_2D
    inverseTransform = alignType == em.ALIGN_PROJ

    alignment = em.Transform()
    angles = numpy.zeros(3)
    shifts = numpy.zeros(3)
    shifts[0] = alignmentList[3]
    shifts[1] = alignmentList[4]
    shifts[2] = 0
    angles[0] = alignmentList[0]
    angles[1] = alignmentList[1]
    angles[2] = alignmentList[2]

    matrix = matrixFromGeometry(shifts, angles, inverseTransform)
    alignment.setMatrix(matrix)

    return alignment


def iterParticlesByMic(partSet):
    """ Iterate the particles ordered by micrograph """
    for i, part in enumerate(partSet.iterItems(orderBy=['_micId', 'id'],
                                               direction='ASC')):
        yield i, part


def convertReferences(refSet, outputFn):
    """ Simplified version of writeSetOfParticles function.
    Writes out an hdf stack.
    """
    fileName = ""
    a = 0
    proc = createEmanProcess(args='write')

    for part in refSet:
        objDict = part.getObjDict()
        objDict['hdfFn'] = outputFn
        objDict['_itemId'] = part.getObjId()

        # the index in EMAN begins with 0
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
