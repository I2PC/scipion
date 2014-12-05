# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Laura del Cano (ldelcano@cnb.csic.es)
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
from pyworkflow.em.constants import ALIGN_NONE, ALIGN_2D, ALIGN_PROJ
"""
This module contains converter functions that will serve to:
1. Write from base classes to Relion specific files
2. Read from Relion files to base classes
"""

import os
from os.path import join
import numpy
from collections import OrderedDict

from pyworkflow.utils import Environ
from pyworkflow.utils.path import (createLink, cleanPath, copyFile,
                                   findRootFrom)
from pyworkflow.em import ImageHandler, ObjectWrap, NO_INDEX, String
from pyworkflow.em.data import (Particle, SetOfClasses2D, SetOfClasses3D,
                                SetOfParticles, CTFModel, Transform,
                                Acquisition, Coordinate)
from pyworkflow.em.packages.xmipp3 import XmippMdRow, getLabelPythonType
from pyworkflow.em.relion.constants import *
# Since Relion share several conventions with Xmipp, we will reuse 
# the xmipp3 tools implemented for Scipion here
import xmipp


# This dictionary will be used to map
# between CTFModel properties and Xmipp labels
ACQUISITION_DICT = OrderedDict([ 
       ("_amplitudeContrast", xmipp.RLN_CTF_Q0),
       ("_sphericalAberration", xmipp.RLN_CTF_CS),
       ("_voltage", xmipp.RLN_CTF_VOLTAGE)
       ])

COOR_DICT = OrderedDict([
             ("_x", xmipp.RLN_IMAGE_COORD_X),
             ("_y", xmipp.RLN_IMAGE_COORD_Y)
             ])

CTF_DICT = OrderedDict([
       ("_defocusU", xmipp.RLN_CTF_DEFOCUSU),
       ("_defocusV", xmipp.RLN_CTF_DEFOCUSV),
       ("_defocusAngle", xmipp.RLN_CTF_DEFOCUS_ANGLE)
       ])

CTF_PSD_DICT = OrderedDict([
       ("_psdFile", xmipp.RLN_CTF_IMAGE)
       ])

CTF_EXTRA_LABELS = [   
    xmipp.RLN_CTF_FOM,
    # In Relion the ctf also contains acquisition information
    xmipp.RLN_CTF_Q0,
    xmipp.RLN_CTF_CS,
    xmipp.RLN_CTF_VOLTAGE,
    xmipp.RLN_CTF_MAGNIFICATION,
    xmipp.RLN_CTF_DETECTOR_PIXEL_SIZE
    ]

# Some extra labels to take into account the zscore
IMAGE_EXTRA_LABELS = [
    xmipp.RLN_SELECT_PARTICLES_ZSCORE,
    xmipp.RLN_PARTICLE_ID,
    xmipp.RLN_FRAME_ID,
    ]
 
# ANGLES_DICT = OrderedDict([
#        ("_angleY", xmipp.RLN_ANGLE_Y),
#        ("_angleY2", xmipp.RLN_ANGLE_Y2),
#        ("_angleTilt", xmipp.RLN_ANGLE_TILT)
#        ])
 
ALIGNMENT_DICT = OrderedDict([ 
       ("_rlnOriginX", xmipp.RLN_ORIENT_ORIGIN_X),
       ("_rlnOriginY", xmipp.RLN_ORIENT_ORIGIN_Y),
       ("_rlnOriginZ", xmipp.RLN_ORIENT_ORIGIN_Z),
       ("_rlnAngleRot", xmipp.RLN_ORIENT_ROT),
       ("_rlnAngleTilt", xmipp.RLN_ORIENT_TILT),
       ("_rlnAnglePsi", xmipp.RLN_ORIENT_PSI),
       ])


def getEnviron():
    """ Setup the environment variables needed to launch Relion. """
    environ = Environ(os.environ)
    environ.update({
            'PATH': join(os.environ['RELION_HOME'], 'bin'),
            'LD_LIBRARY_PATH': join(os.environ['RELION_HOME'], 'lib') + ":" + join(os.environ['RELION_HOME'], 'lib64'),
#             'LD_LIBRARY_PATH': join(os.environ['RELION_HOME'], 'lib64'),
            }, position=Environ.BEGIN)
    return environ


def locationToRelion(index, filename):
    """ Convert an index and filename location
    to a string with @ as expected in Relion.
    """
    #TODO: Maybe we need to add more logic dependent of the format
    if index != NO_INDEX:
        return "%06d@%s" % (index, filename)
    
    return filename


def relionToLocation(filename):
    """ Return a location (index, filename) given
    a Relion filename with the index@filename structure. """
    if '@' in filename:
        return xmipp.FileName(filename).decompose()
    else:
        return NO_INDEX, str(filename)


class ParticleAdaptor():
    """ Class used to convert a set of particles for Relion.
    It will write an stack in Spider format and also
    modify the output star file to point to the new stack.
    """
    def __init__(self, imgSet, filesMapping={}, originalSet=None,
                 preprocessImageRow=None,
                 postprocessImageRow=None):
        self._rowCount = 1
        self._ih = ImageHandler()
        self._imgSet = imgSet
        self._filesMapping = filesMapping
        self._originalSet = originalSet
        self._preprocessImageRow = preprocessImageRow
        self._postprocessImageRow = postprocessImageRow
        
    def preprocessImageRow(self, img, imgRow):
        """ If binary images where converted or just 
        created links, replace the orignal names.
        """
        # preprocessImageRow is a way relion protocol have to
        # do some preprocessing on image object or row 
        # before the normal setup
        if self._preprocessImageRow:
            self._preprocessImageRow(img, imgRow)
            
        # Re-write the row with the new location
        if self._filesMapping:
            # If the filesMapping is not empty
            # we assume that the stack files have been converted
            # to be properly read by relion and also that
            # each particle preserve its index and only change the filename
            newFile = self._filesMapping[img.getFileName()]
            newLoc = (img.getIndex(), newFile)
            img.setLocation(newLoc)

    def postprocessImageRow(self, img, imgRow):
        """ If binary images where converted or just 
        created links, replace the orignal names.
        """            
        if img.hasMicId():
            imgRow.setValue('rlnMicrographName', 'fake_micrograph_%06d.mrc' % img.getMicId())
            imgRow.setValue(xmipp.MDL_MICROGRAPH_ID, long(img.getMicId()))
            
        if imgRow.hasLabel(xmipp.MDL_PARTICLE_ID):
            # Write stuff for particle polishing
            movieId = imgRow.getValue(xmipp.MDL_MICROGRAPH_ID)
            movieName = 'fake_micrograph_%06d_movie.mrcs' % movieId 
            frameId = imgRow.getValue(xmipp.MDL_FRAME_ID) + 1
            micName = '%06d@%s' % (frameId, movieName)
            particleId = imgRow.getValue(xmipp.MDL_PARTICLE_ID)
            particle = self._originalSet[particleId]
            if particle is None:
                particleName = 'None'
                raise Exception("Particle with id %d not found!!!" % particleId)
            else:
                particleName = locationToRelion(*particle.getLocation())
                
            imgRow.setValue('rlnMicrographName', micName)
            imgRow.setValue('rlnParticleName', particleName)
            
        for label, _ in imgRow:
            if not label in XMIPP_RELION_LABELS:
                imgRow.removeLabel(label)
        self._rowCount += 1
        
        if self._postprocessImageRow:
            self._postprocessImageRow(img, imgRow)


def objectToRow(obj, row, attrDict, extraLabels={}):
    """ This function will convert an EMObject into a XmippMdRow.
    Params:
        obj: the EMObject instance (input)
        row: the XmippMdRow instance (output)
        attrDict: dictionary with the map between obj attributes(keys) and 
            row MDLabels in Xmipp (values).        
        extraLabels: a list with extra labels that could be included
            as _xmipp_labelName
    """
    if obj.isEnabled():
        enabled = 1
    else:
        enabled = -1
    row.setValue(xmipp.RLN_IMAGE_ENABLED, enabled)
    
    for attr, label in attrDict.iteritems():
        if hasattr(obj, attr):
            valueType = getLabelPythonType(label)
            row.setValue(label, valueType(getattr(obj, attr).get()))

    attrLabels = attrDict.values()
    
    for label in extraLabels:
        attrName = '_rln_%s' % xmipp.label2Str(label)
        if label not in attrLabels and hasattr(obj, attrName):
            value = obj.getAttributeValue(attrName) 
            row.setValue(label, value)


def rowToObject(row, obj, attrDict, extraLabels={}):
    """ This function will convert from a XmippMdRow to an EMObject.
    Params:
        row: the XmippMdRow instance (input)
        obj: the EMObject instance (output)
        attrDict: dictionary with the map between obj attributes(keys) and 
            row MDLabels in Xmipp (values).
        extraLabels: a list with extra labels that could be included
            as _xmipp_labelName
    """
    obj.setEnabled(row.getValue(xmipp.RLN_IMAGE_ENABLED, 1) > 0)
    
    for attr, label in attrDict.iteritems():
        value = row.getValue(label)
        if not hasattr(obj, attr):
            setattr(obj, attr, ObjectWrap(value))
        else:
            getattr(obj, attr).set(value)
        
    attrLabels = attrDict.values()
    
    for label in extraLabels:
        if label not in attrLabels and row.hasLabel(label):
            labelStr = xmipp.label2Str(label)
            setattr(obj, '_xmipp_%s' % labelStr, row.getValueAsObject(label))
    
    
def rowFromMd(md, objId):
    row = XmippMdRow()
    row.readFromMd(md, objId)
    return row


def _containsAll(row, labels):
    """ Check if the labels (values) in labelsDict
    are present in the row.
    """
    values = labels.values() if isinstance(labels, dict) else labels
    return all(row.containsLabel(l) for l in values)


def _containsAny(row, labels):
    """ Check if the labels (values) in labelsDict
    are present in the row.
    """
    values = labels.values() if isinstance(labels, dict) else labels
    return any(row.containsLabel(l) for l in values)

    
def _rowToObjectFunc(obj):
    """ From a given object, return the function rowTo{OBJECT_CLASS_NAME}. """
    return globals()['rowTo' + obj.getClassName()]


def _readSetFunc(obj):
    """ From a given object, return the function read{OBJECT_CLASS_NAME}. """
    return globals()['read' + obj.getClassName()]


def setObjId(obj, mdRow, label=xmipp.RLN_IMAGE_ID):
    if mdRow.containsLabel(label):
        obj.setObjId(mdRow.getValue(label))
    else:
        obj.setObjId(None)


def setRowId(mdRow, obj, label=xmipp.RLN_IMAGE_ID):
    mdRow.setValue(label, long(obj.getObjId()))


def acquisitionToRow(acquisition, ctfRow):
    """ Set labels values from acquisition to md row. """
    objectToRow(acquisition, ctfRow, ACQUISITION_DICT)


def rowToAcquisition(acquisitionRow):
    """ Create an acquisition from a row of a metadata. """
    if _containsAll(acquisitionRow, ACQUISITION_DICT):
        acquisition = Acquisition()
        rowToObject(acquisitionRow, acquisition, ACQUISITION_DICT) 
    else:                
        acquisition = None
    
    return acquisition


def setPsdFiles(ctfModel, ctfRow):
    """ Set the PSD files of CTF estimation related
    to this ctfModel. The values will be read from
    the ctfRow if present.
    """
    for attr, label in CTF_PSD_DICT.iteritems():
        if ctfRow.containsLabel(label):
            setattr(ctfModel, attr, String(ctfRow.getValue(label)))


def rowToCtfModel(ctfRow):
    """ Create a CTFModel from a row of a metadata. """
    if _containsAll(ctfRow, CTF_DICT):
        ctfModel = CTFModel()
        rowToObject(ctfRow, ctfModel, CTF_DICT, extraLabels=CTF_EXTRA_LABELS)
        ctfModel.standardize()
        setPsdFiles(ctfModel, ctfRow)
    else:
        ctfModel = None
        
    return ctfModel


def geometryFromMatrix(matrix, inverseTransform):
    from pyworkflow.em.transformations import translation_from_matrix, euler_from_matrix
    from numpy import rad2deg
    
    if inverseTransform:
        from numpy.linalg import inv
        matrix = inv(matrix)
        shifts = -translation_from_matrix(matrix)
    else:
        shifts = translation_from_matrix(matrix)
    angles = -rad2deg(euler_from_matrix(matrix, axes='szyz'))
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


def rowToAlignment(alignmentRow, alignType):
    """
    is2D == True-> matrix is 2D (2D images alignment)
            otherwise matrix is 3D (3D volume alignment or projection)
    invTransform == True  -> for xmipp implies projection
    """
    is2D = alignType == ALIGN_2D
    inverseTransform = alignType == ALIGN_PROJ
    
    if _containsAny(alignmentRow, ALIGNMENT_DICT):
        alignment = Transform()
        angles = numpy.zeros(3)
        shifts = numpy.zeros(3)
        angles[2] = alignmentRow.getValue(xmipp.RLN_ORIENT_PSI, 0.)
        shifts[0] = alignmentRow.getValue(xmipp.RLN_ORIENT_ORIGIN_X, 0.)
        shifts[1] = alignmentRow.getValue(xmipp.RLN_ORIENT_ORIGIN_Y, 0.)
        if not is2D:
            angles[0] = alignmentRow.getValue(xmipp.RLN_ORIENT_ROT, 0.)
            angles[1] = alignmentRow.getValue(xmipp.RLN_ORIENT_TILT, 0.)
            shifts[2] = alignmentRow.getValue(xmipp.RLN_ORIENT_ORIGIN_Z, 0.)

        M = matrixFromGeometry(shifts, angles, inverseTransform)
        alignment.setMatrix(M)
        
        #FIXME: now are also storing the alignment parameters since
        # the conversions to the Transform matrix have not been extensively tested.
        # After this, we should only keep the matrix 
        #for paramName, label in ALIGNMENT_DICT.iteritems():
        #    if alignmentRow.hasLabel(label):
        #        setattr(alignment, paramName, alignmentRow.getValueAsObject(label))    
    else:
        alignment = None
    
    return alignment


def alignmentToRow(alignment, alignmentRow, alignType):
    """
    is2D == True-> matrix is 2D (2D images alignment)
            otherwise matrix is 3D (3D volume alignment or projection)
    invTransform == True  -> for xmipp implies projection
                          -> for xmipp implies alignment
    """
    is2D = alignType == ALIGN_2D
    inverseTransform = alignType == ALIGN_PROJ
    
    shifts, angles = geometryFromMatrix(alignment.getMatrix(), inverseTransform)

    alignmentRow.setValue(xmipp.RLN_ORIENT_ORIGIN_X, shifts[0])
    alignmentRow.setValue(xmipp.RLN_ORIENT_ORIGIN_Y, shifts[1])
    
    if is2D:
        angle = angles[0] + angles[2]
        alignmentRow.setValue(xmipp.RLN_ORIENT_PSI,  angle)
    else:
        alignmentRow.setValue(xmipp.RLN_ORIENT_ORIGIN_Z, shifts[2])
        alignmentRow.setValue(xmipp.RLN_ORIENT_ROT,  angles[0])
        alignmentRow.setValue(xmipp.RLN_ORIENT_TILT, angles[1])
        alignmentRow.setValue(xmipp.RLN_ORIENT_PSI,  angles[2])


def coordinateToRow(coord, coordRow, copyId=True):
    """ Set labels values from Coordinate coord to md row. """
    if copyId:
        setRowId(coordRow, coord)
    objectToRow(coord, coordRow, COOR_DICT)
    if coord.getMicId():
        coordRow.setValue(xmipp.RLN_MICROGRAPH_NAME, str(coord.getMicId()))


def rowToCoordinate(coordRow):
    """ Create a Coordinate from a row of a metadata. """
    # Check that all required labels are present in the row
    if _containsAll(coordRow, COOR_DICT):
        coord = Coordinate()
        rowToObject(coordRow, coord, COOR_DICT)
            
        # Setup the micId if is integer value
        try:
            coord.setMicId(int(coordRow.getValue(xmipp.RLN_MICROGRAPH_NAME)))
        except Exception:
            pass
    else:
        coord = None
        
    return coord


def rowToParticle(partRow, **kwargs):
    """ Create a Particle from a row of a metadata. """
    img = Particle()
    
    # Provide a hook to be used if something is needed to be 
    # done for special cases before converting image to row
    preprocessImageRow = kwargs.get('preprocessImageRow', None)
    if preprocessImageRow:
        preprocessImageRow(img, partRow)
    
    # Decompose Relion filename
    index, filename = relionToLocation(partRow.getValue(xmipp.RLN_IMAGE_NAME))
    img.setLocation(index, filename)
    
#     if partRow.containsLabel(xmipp.MDL_REF):
#         img.setClassId(partRow.getValue(xmipp.MDL_REF))
#     elif partRow.containsLabel(xmipp.MDL_REF3D):
#         img.setClassId(partRow.getValue(xmipp.MDL_REF3D))
    
    if kwargs.get('readCtf', True):
        img.setCTF(rowToCtfModel(partRow))
        
    # alignment is mandatory at this point, it shoud be check
    # and detected defaults if not passed at readSetOf.. level
    alignType = kwargs.get('alignType') 
    
    if alignType != ALIGN_NONE:
        img.setTransform(rowToAlignment(partRow, alignType))
        
    if kwargs.get('readAcquisition', True):
        img.setAcquisition(rowToAcquisition(partRow))
        
    if kwargs.get('magnification', None):
        img.getAcquisition().setMagnification(kwargs.get("magnification"))
    
    setObjId(img, partRow)
    # Read some extra labels
    rowToObject(partRow, img, {}, 
                extraLabels=IMAGE_EXTRA_LABELS + kwargs.get('extraLabels', []))

    # Provide a hook to be used if something is needed to be 
    # done for special cases before converting image to row
    postprocessImageRow = kwargs.get('postprocessImageRow', None)
    if postprocessImageRow:
        postprocessImageRow(img, partRow)
    
    img.setCoordinate(rowToCoordinate(partRow))
    # copy micId if available
    # if not copy micrograph name if available
    try:
        if partRow.hasLabel(xmipp.RLN_MICROGRAPH_ID):
            img.setMicId(partRow.getValue(xmipp.RLN_MICROGRAPH_ID))
    except Exception as e:
        print "Warning:", e.message
    return img


def readSetOfParticles(filename, partSet, **kwargs):
    """read from Relion image metadata.
        filename: The metadata filename where the image are.
        imgSet: the SetOfParticles that will be populated.
        rowToParticle: this function will be used to convert the row to Object
    """    
    imgMd = xmipp.MetaData(filename)
    # By default remove disabled items from metadata
    # be careful if you need to preserve the original number of items
    if kwargs.get('removeDisabled', True):
        imgMd.removeDisabled()
    
    # If the type of alignment is not sent throught the kwargs
#     # try to deduced from the metadata labels
#     if 'alignType' not in kwargs:
#         imgRow = rowFromMd(imgMd, imgMd.firstObject())
#         if _containsAny(imgRow, ALIGNMENT_DICT):
#             if imgRow.getValue(xmipp.RLN_ORIENT_TILT) != 0:
#                 kwargs['alignType'] = ALIGN_PROJ
#             else:
#                 kwargs['alignType'] = ALIGN_2D
#         else:
#             kwargs['alignType'] = ALIGN_NONE
    
    for objId in imgMd:
        imgRow = rowFromMd(imgMd, objId)
        img = rowToParticle(imgRow, **kwargs)
        partSet.append(img)
        
    partSet.setHasCTF(img.hasCTF())
    partSet.setAlignment(kwargs['alignType'])
    
##### DE AQUI EN ADELANTE HAY QUE PASAR DESDE XMIPP.

def writeSetOfParticles(imgSet, starFile,
                        outputDir,
                        originalSet=None, **kwargs):
    """ This function will write a SetOfImages as Relion metadata.
    Params:
        imgSet: the SetOfImages instance.
        starFile: the filename where to write the metadata.
        filesMapping: this dict will help when there is need to replace images names
    """
    import pyworkflow.em.packages.xmipp3 as xmipp3
    filesMapping = convertBinaryFiles(imgSet, outputDir)
    pa = ParticleAdaptor(imgSet, filesMapping, originalSet,
                         preprocessImageRow=kwargs.get('preprocessImageRow', None),
                         postprocessImageRow=kwargs.get('postprocessImageRow', None))
    # Intercept the pre and post processImageRow
    kwargs['preprocessImageRow'] = pa.preprocessImageRow
    kwargs['postprocessImageRow'] = pa.postprocessImageRow
    xmipp3.writeSetOfParticles(imgSet, starFile, **kwargs)


def createClassesFromImages(inputImages, inputStar,
                            classesFn, ClassType, classLabel, classFnTemplate,
                            iter, preprocessImageRow=None):
    """ From an intermediate dataXXX.star file produced by relion, create
    the set of classes in which those images are classified.
    Params:
        inputImages: the SetOfImages that were initially classified by relion. 
        inputStar: the filename with relion star format.
        classesFn: filename where to write the classes.
        ClassType: the output type of the classes set ( usually SetOfClass2D or 3D )
        classLabel: label that is the class reference.
        classFnTemplate: the template to get the classes averages filenames
        iter: the iteration number, just used in Class template
    """
    import pyworkflow.em.packages.xmipp3 as xmipp3
    addRelionLabels(replace=True, extended=True)
    # We asume here that the volumes (classes3d) are in the same folder than imgsFn
    # rootDir here is defined to be used expanding locals()
    if ClassType == SetOfClasses2D:
        alignType = ALIGN_2D
    elif ClassType == SetOfClasses3D:
        alignType = ALIGN_PROJ
    else:
        raise Exception('Invalid ClassType %s' % ClassType)
    
    xmipp3.createClassesFromImages(inputImages, inputStar,
                                   classesFn, ClassType, classLabel, classFnTemplate,
                                   alignType, iter, preprocessImageRow)
    restoreXmippLabels()    
    
# def createClassesFromImages(inputImages, inputStar, classesFn, ):
    

def readSetOfClasses2D(classes2DSet, filename, **args):
    clsSet = SetOfClasses2D(filename=filename)
    classes2DSet.appendFromClasses(clsSet)


def readSetOfClasses3D(classes3DSet, filename, **args):
    clsSet = SetOfClasses3D(filename=filename)
    classes3DSet.appendFromClasses(clsSet)
    

def writeSqliteIterData(imgStar, imgSqlite):
    """ Given a Relion images star file (from some iteration)
    create the corresponding SetOfParticles (sqlite file)
    for this iteration. This file can be visualized sorted
    by the LogLikelihood.
    """
    addRelionLabels(replace=True, extended=True)
    cleanPath(imgSqlite)
    imgSet = SetOfParticles(filename=imgSqlite)
    readSetOfParticles(imgStar, imgSet)
    imgSet.write()
    restoreXmippLabels()
    
    
def writeIterAngularDist(self, inDataStar, outAngularDist, numberOfClasses, prefixes):
    """ Write the angular distribution. Should be overriden in subclasses. """
    addRelionLabels(replace=True, extended=True)
    
    md = xmipp.MetaData(inDataStar)
    
    refsList = range(1, numberOfClasses + 1) 
    print "refsList: ", refsList
    for ref3d in refsList:
        for prefix in prefixes:
            mdGroup = xmipp.MetaData()
            mdGroup.importObjects(md, xmipp.MDValueEQ(xmipp.MDL_REF, ref3d))
            mdDist = xmipp.MetaData()
            mdDist.aggregateMdGroupBy(mdGroup, xmipp.AGGR_COUNT,
                                      [xmipp.MDL_ANGLE_ROT, xmipp.MDL_ANGLE_TILT],
                                      xmipp.MDL_ANGLE_ROT, xmipp.MDL_WEIGHT)
            mdDist.setValueCol(xmipp.MDL_ANGLE_PSI, 0.0)
            blockName = '%sclass%06d_angularDist@' % (prefix, ref3d)
            print "Writing angular distribution to: ", blockName + outAngularDist
            mdDist.write(blockName + outAngularDist, xmipp.MD_APPEND) 
            
    restoreXmippLabels()
    
    
def splitInCTFGroups(imgStar, groups):
    """ Add a new colunm in the image star to separate the particles into ctf groups """
    addRelionLabels(replace=True)
    
    md = xmipp.MetaData(imgStar)
    defocus = [md.getValue(xmipp.MDL_CTF_DEFOCUSU, objId) for objId in md]
    
    minDef = min(defocus)
    maxDef = max(defocus)
    
    increment = (maxDef - minDef) / groups
    counter = 1
    part = 1
    md.sort(xmipp.MDL_CTF_DEFOCUSU)
    
    for objId in md:
        partDef = md.getValue(xmipp.MDL_CTF_DEFOCUSU, objId)
        currDef = minDef + increment
        
        if partDef <= currDef:
            part = part + 1
            md.setValue(xmipp.MDL_SERIE, "ctfgroup_%05d" % counter, objId)
        else:
            if part < 100:
                increment = md.getValue(xmipp.MDL_CTF_DEFOCUSU, objId) - minDef
                part = part + 1
                md.setValue(xmipp.MDL_SERIE, "ctfgroup_%05d" % counter, objId)
            else:
                part = 1
                minDef = md.getValue(xmipp.MDL_CTF_DEFOCUSU, objId)
                counter = counter + 1
                increment = (maxDef - minDef) / (groups - counter)
                md.setValue(xmipp.MDL_SERIE, "ctfgroup_%05d" % counter, objId)
    md.write(imgStar)
    restoreXmippLabels()
    
        
def getMdFirstRow(starfile):
    addRelionLabels(replace=True, extended=True)
    from pyworkflow.em.packages.xmipp3.utils import getMdFirstRow
    row = getMdFirstRow(starfile)
    restoreXmippLabels()
    return row


def findImagesPath(starFile):
    """ Find the path of the images relative to some star file. """
    row = getMdFirstRow(starFile)
    _, imgFile = relionToLocation(row.getValue(xmipp.MDL_IMAGE))
    
    return findRootFrom(starFile, imgFile)
    

def prependToFileName(imgRow, prefixPath):
    """ Prepend some root name to imageRow filename. """
    index, imgPath = relionToLocation(imgRow.getValue(xmipp.MDL_IMAGE))
    newLoc = locationToRelion(index, os.path.join(prefixPath, imgPath))
    imgRow.setValue(xmipp.MDL_IMAGE, newLoc)


def relativeFromFileName(imgRow, prefixPath):
    """ Remove some prefix from filename in row. """
    index, imgPath = relionToLocation(imgRow.getValue(xmipp.MDL_IMAGE))
    imgPath = os.path.relpath(imgPath, prefixPath)
    newLoc = locationToRelion(index, imgPath)
    imgRow.setValue(xmipp.MDL_IMAGE, newLoc)
    
    
def copyOrLinkFileName(imgRow, prefixDir, outputDir, copyFiles=False):
    index, imgPath = relionToLocation(imgRow.getValue(xmipp.MDL_IMAGE))
    baseName = os.path.basename(imgPath)
    newName = os.path.join(outputDir, baseName)
    if not os.path.exists(newName):
        if copyFiles:
            copyFile(os.path.join(prefixDir, imgPath), newName)
        else:
            createLink(os.path.join(prefixDir, imgPath), newName)
            
    imgRow.setValue(xmipp.MDL_IMAGE, locationToRelion(index, newName))
    

def setupCTF(imgRow, sampling):
    """ Do some validations and set some values
    for Relion import.
    """
    imgRow.setValue(xmipp.MDL_SAMPLINGRATE, sampling)
    # TODO: check if we want to move this behaviour to setup CTFModel by default
    hasDefocusU = imgRow.containsLabel(xmipp.MDL_CTF_DEFOCUSU)
    hasDefocusV = imgRow.containsLabel(xmipp.MDL_CTF_DEFOCUSV)
    hasDefocusAngle = imgRow.containsLabel(xmipp.MDL_CTF_DEFOCUS_ANGLE)
    
    if hasDefocusU or hasDefocusV:
        if not hasDefocusU:
            imgRow.setValue(xmipp.MDL_CTF_DEFOCUSU, imgRow.getValue(xmipp.MDL_CTF_DEFOCUSV))
        if not hasDefocusV:
            imgRow.setValue(xmipp.MDL_CTF_DEFOCUSV, imgRow.getValue(xmipp.MDL_CTF_DEFOCUSU))
        if not hasDefocusAngle:
            imgRow.setValue(xmipp.MDL_CTF_DEFOCUS_ANGLE, 0.)
            

def convertBinaryFiles(imgSet, outputDir):
    """ Convert binary images files to a format read by Relion.
    Params:
        imgSet: input image set to be converted.
        outputDir: where to put the converted file(s)
    Return:
        A dictionary with old-file as key and new-file as value
        If empty, not conversion was done.
    """
    filesMapping = {}
    # This approach can be extended when
    # converting from a binary file format that
    # is not read from Relion
    ext = imgSet.getFirstItem().getFileName()
    if ext.endswith('.mrc'):
        uniqueFiles = list(imgSet.getFiles())
        for f in uniqueFiles:
            newFile = os.path.join(outputDir, os.path.basename(f) + 's')
            createLink(f, newFile)
            filesMapping[f] = newFile
    elif ext.endswith('.hdf'):  # assume eman
        uniqueFiles = list(imgSet.getFiles())
        from pyworkflow.em.convert import convertStack
        for f in uniqueFiles:
            newFile = os.path.join(outputDir, os.path.basename(f.rsplit(".", 1)[ 0 ]) + '.spi')
            convertStack(f, newFile)
            filesMapping[f] = newFile

    return filesMapping


def createItemMatrix(self, item, row, align):
    from pyworkflow.em.packages.xmipp3.convert import rowToAlignment
    
    item.setTransform(rowToAlignment(row, alignType=align))

def iterRows(md):
    from pyworkflow.em.packages.xmipp3.utils import iterMdRows
    return iterMdRows(md)
