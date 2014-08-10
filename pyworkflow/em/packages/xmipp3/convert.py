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
"""
This module contains converter functions that will serve to:
1. Write from base classes to Xmipp specific files
2. Read from Xmipp files to base classes
"""

import os
from collections import OrderedDict
from itertools import izip
import numpy

import xmipp
from xmipp3 import XmippMdRow, getLabelPythonType, RowMetaData
from pyworkflow.em import *
from pyworkflow.em.constants import NO_INDEX
from pyworkflow.object import String
from pyworkflow.utils.path import join, dirname, replaceBaseExt


# This dictionary will be used to map
# between CTFModel properties and Xmipp labels
ACQUISITION_DICT = OrderedDict([ 
       ("_amplitudeContrast", xmipp.MDL_CTF_Q0),
       ("_sphericalAberration", xmipp.MDL_CTF_CS),
       ("_voltage", xmipp.MDL_CTF_VOLTAGE)
       ])

COOR_DICT = OrderedDict([
             ("_x", xmipp.MDL_XCOOR), 
             ("_y", xmipp.MDL_YCOOR) 
             ])

CTF_DICT = OrderedDict([
       ("_defocusU", xmipp.MDL_CTF_DEFOCUSU),
       ("_defocusV", xmipp.MDL_CTF_DEFOCUSV),
       ("_defocusAngle", xmipp.MDL_CTF_DEFOCUS_ANGLE)
       ])

CTF_EXTRA_LABELS = [   
    xmipp.MDL_CTF_CA,
    xmipp.MDL_CTF_ENERGY_LOSS,
    xmipp.MDL_CTF_LENS_STABILITY,
    xmipp.MDL_CTF_CONVERGENCE_CONE,
    xmipp.MDL_CTF_LONGITUDINAL_DISPLACEMENT,
    xmipp.MDL_CTF_TRANSVERSAL_DISPLACEMENT,
    xmipp.MDL_CTF_K,
    xmipp.MDL_CTF_BG_GAUSSIAN_K,
    xmipp.MDL_CTF_BG_GAUSSIAN_SIGMAU,
    xmipp.MDL_CTF_BG_GAUSSIAN_SIGMAV,
    xmipp.MDL_CTF_BG_GAUSSIAN_CU,
    xmipp.MDL_CTF_BG_GAUSSIAN_CV,
    xmipp.MDL_CTF_BG_SQRT_K,
    xmipp.MDL_CTF_BG_SQRT_U,
    xmipp.MDL_CTF_BG_SQRT_V,
    xmipp.MDL_CTF_BG_SQRT_ANGLE,
    xmipp.MDL_CTF_BG_BASELINE,
    xmipp.MDL_CTF_BG_GAUSSIAN2_K,
    xmipp.MDL_CTF_BG_GAUSSIAN2_SIGMAU,
    xmipp.MDL_CTF_BG_GAUSSIAN2_SIGMAV,
    xmipp.MDL_CTF_BG_GAUSSIAN2_CU,
    xmipp.MDL_CTF_BG_GAUSSIAN2_CV,
    xmipp.MDL_CTF_BG_GAUSSIAN2_ANGLE,
    xmipp.MDL_CTF_CRIT_FITTINGSCORE,
    xmipp.MDL_CTF_CRIT_FITTINGCORR13,
    xmipp.MDL_CTF_DOWNSAMPLE_PERFORMED,
    xmipp.MDL_CTF_CRIT_PSDVARIANCE,
    xmipp.MDL_CTF_CRIT_PSDPCA1VARIANCE,
    xmipp.MDL_CTF_CRIT_PSDPCARUNSTEST,
    xmipp.MDL_CTF_CRIT_FIRSTZEROAVG,
    xmipp.MDL_CTF_CRIT_DAMPING,
    xmipp.MDL_CTF_CRIT_FIRSTZERORATIO,
    xmipp.MDL_CTF_CRIT_PSDCORRELATION90,
    xmipp.MDL_CTF_CRIT_PSDRADIALINTEGRAL,
    xmipp.MDL_CTF_CRIT_NORMALITY,  
    ]

# Some extra labels to take into account the zscore
IMAGE_EXTRA_LABELS = [
    xmipp.MDL_ZSCORE,
    xmipp.MDL_ZSCORE_HISTOGRAM,
    xmipp.MDL_ZSCORE_RESMEAN,
    xmipp.MDL_ZSCORE_RESVAR,
    xmipp.MDL_ZSCORE_RESCOV,
    xmipp.MDL_ZSCORE_SHAPE1,
    xmipp.MDL_ZSCORE_SHAPE2,
    xmipp.MDL_ZSCORE_SNR1,
    xmipp.MDL_ZSCORE_SNR2,
    ]

ANGLES_DICT = OrderedDict([
       ("_angleY", xmipp.MDL_ANGLE_Y),
       ("_angleY2", xmipp.MDL_ANGLE_Y2),
       ("_angleTilt", xmipp.MDL_ANGLE_TILT)
       ])

ALIGNMENT_EXTRA_LABELS = [
    xmipp.MDL_ANGLE_PSI,
    xmipp.MDL_SHIFT_X,
    xmipp.MDL_SHIFT_Y,
    xmipp.MDL_FLIP,
    ]


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
    for attr, label in attrDict.iteritems():
        if hasattr(obj, attr):
            valueType = getLabelPythonType(label)
            row.setValue(label, valueType(getattr(obj, attr).get()))

    attrLabels = attrDict.values()
    
    for label in extraLabels:
        attrName = '_xmipp_%s' % xmipp.label2Str(label)
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
    for attr, label in attrDict.iteritems():
        value = row.getValue(label)
        if not hasattr(obj, attr):
            setattr(obj, attr, ObjectWrap(value))
        else:
            getattr(obj, attr).set(value)
        
    attrLabels = attrDict.values()
    
    for label in extraLabels:
        if label not in attrLabels and row.hasLabel(label):
            value = row.getValue(label)
            labelStr = xmipp.label2Str(label)
            setattr(obj, '_xmipp_%s' % labelStr, ObjectWrap(value))
    
    
def rowFromMd(md, objId):
    row = XmippMdRow()
    row.readFromMd(md, objId)
    return row


def containsLabels(row, labels):
    """ Check if the labels (values) in labelsDict
    are present in the row.
    """
    values = labels.values() if isinstance(labels, dict) else labels
    return all(row.containsLabel(l) for l in values)
    
#     
# def rowToObject(md, objId, obj, attrDict, extraLabels={}):
#     """ Same as rowToObject, but creating the row from md and objId. """
#     _rowToObject(rowFromMd(md, objId), obj, attrDict, extraLabels)
    
def _rowToObjectFunc(obj):
    """ From a given object, return the function rowTo{OBJECT_CLASS_NAME}. """
    return globals()['rowTo' + obj.getClassName()]

def _readSetFunc(obj):
    """ From a given object, return the function read{OBJECT_CLASS_NAME}. """
    return globals()['read' + obj.getClassName()]
    
def locationToXmipp(index, filename):
    """ Convert an index and filename location
    to a string with @ as expected in Xmipp.
    """
    #TODO: Maybe we need to add more logic dependent of the format
    if index != NO_INDEX:
        return "%06d@%s" % (index, filename)
    
    return filename


def fixVolumeFileName(image):
    """ This method will add :mrc to .mrc volumes
    because for mrc format is not possible to distinguish 
    between 3D volumes and 2D stacks. 
    """
    fn = image.getFileName()
    if isinstance(image, Volume):
        if fn.endswith('.mrc'):
            fn += ':mrc'
        
    return fn   
    
def getImageLocation(image):
    xmippFn = locationToXmipp(image.getIndex(),
                              fixVolumeFileName(image))
        
    return xmippFn


def xmippToLocation(xmippFilename):
    """ Return a location (index, filename) given
    a Xmipp filename with the index@filename structure. """
    if '@' in xmippFilename:
        return xmipp.FileName(xmippFilename).decompose()
    else:
        return NO_INDEX, str(xmippFilename)


def setObjId(obj, mdRow, label=xmipp.MDL_ITEM_ID):
    if mdRow.containsLabel(label):
        obj.setObjId(mdRow.getValue(label))
    else:
        obj.setObjId(None)
    
    
def setRowId(mdRow, obj, label=xmipp.MDL_ITEM_ID):
    mdRow.setValue(label, long(obj.getObjId()))


def micrographToCTFParam(mic, ctfparam):
    """ This function is used to convert a Micrograph object
    to the .ctfparam metadata needed by some Xmipp protograms.
    If the micrograph already comes from xmipp, the ctfparam
    will be returned, if not, the new file. 
    """
    ctf = mic.getCTF()
    
    if hasattr(ctf, '_xmippMd'):
        return ctf._xmippMd.get()
    
    md = xmipp.MetaData()
    md.setColumnFormat(False)
    row = XmippMdRow()
    ctfModelToRow(ctf, row)
    acquisitionToRow(mic.getAcquisition(), row)
    row.writeToMd(md, md.addObject())
    md.write(ctfparam)
    
    return ctfparam
    
    
def imageToRow(img, imgRow, imgLabel, **kwargs):
    setRowId(imgRow, img) # Set the id in the metadata as MDL_ITEM_ID
    index, filename = img.getLocation()
    fn = locationToXmipp(index, filename)
    imgRow.setValue(imgLabel, fn)
    
    if kwargs.get('writeCtf', True) and img.hasCTF():
        ctfModelToRow(img.getCTF(), imgRow)
        
    if kwargs.get('writeAlignment', True) and img.hasAlignment():
        alignmentToRow(img.getAlignment(), imgRow)
        
    if kwargs.get('writeAcquisition', True) and img.hasAcquisition():
        acquisitionToRow(img.getAcquisition(), imgRow)
    
    
    # Write all extra labels to the row    
    objectToRow(img, imgRow, {}, extraLabels=IMAGE_EXTRA_LABELS)
    
        
def rowToImage(imgRow, imgLabel, imgClass, **kwargs):
    """ Create an Image from a row of a metadata. """
    img = imgClass()
    
    preprocessRow = kwargs.get('preprocessRow', None)
    if preprocessRow:
        preprocessRow(imgRow)
    
    # Decompose Xmipp filename
    index, filename = xmippToLocation(imgRow.getValue(imgLabel))
    img.setLocation(index, filename)
    
    if kwargs.get('readCtf', True):
        img.setCTF(rowToCtfModel(imgRow))
        
    if kwargs.get('readAlignment', True):
        img.setAlignment(rowToAlignment(imgRow))
        
    if kwargs.get('readAcquisition', True):
        img.setAcquisition(rowToAcquisition(imgRow))
    
    setObjId(img, imgRow)
    # Read some extra labels
    rowToObject(imgRow, img, {}, extraLabels=IMAGE_EXTRA_LABELS)
    
    return img
    
    
def micrographToRow(mic, micRow, **kwargs):
    """ Set labels values from Micrograph mic to md row. """
    imageToRow(mic, micRow, imgLabel=xmipp.MDL_MICROGRAPH, **kwargs)
    
    
def rowToMicrograph(micRow, **kwargs):
    """ Create a Micrograph object from a row of Xmipp metadata. """
    return rowToImage(micRow, xmipp.MDL_MICROGRAPH, Micrograph, **kwargs)


def volumeToRow(vol, volRow, **kwargs):
    """ Set labels values from Micrograph mic to md row. """
    imageToRow(vol, volRow, imgLabel=xmipp.MDL_IMAGE, **kwargs)


def rowToVolume(volRow, **kwargs):
    """ Create a Volume object from a row of Xmipp metadata. """
    return rowToImage(volRow, xmipp.MDL_IMAGE, Volume, **kwargs)


def coordinateToRow(coord, coordRow, copyId=True):
    """ Set labels values from Coordinate coord to md row. """
    if copyId:
        setRowId(coordRow, coord)
    objectToRow(coord, coordRow, COOR_DICT)
    if coord.getMicId():
        coordRow.setValue(xmipp.MDL_MICROGRAPH, str(coord.getMicId()))


def rowToCoordinate(coordRow):
    """ Create a Coordinate from a row of a metadata. """
    # Check that all required labels are present in the row
    if containsLabels(coordRow, COOR_DICT):
        coord = Coordinate()
        rowToObject(coordRow, coord, COOR_DICT)
            
        # Setup the micId if is integer value
        try:
            coord.setMicId(int(coordRow.getValue(xmipp.MDL_MICROGRAPH)))
        except Exception:
            pass
    else:
        coord = None
        
    return coord


def rowToParticle(partRow, **kwargs):
    """ Create a Particle from a row of a metadata. """
    img = rowToImage(partRow, xmipp.MDL_IMAGE, Particle, **kwargs)
    img.setCoordinate(rowToCoordinate(partRow))
    # Setup the micId if is integer value
    try:
        img.setMicId(int(partRow.getValue(xmipp.MDL_MICROGRAPH)))
    except Exception:
        pass    
    return img
    
    
def particleToRow(part, partRow, **kwargs):
    """ Set labels values from Particle to md row. """
    imageToRow(part, partRow, xmipp.MDL_IMAGE, **kwargs)
    coord = part.getCoordinate()
    if coord is not None:
        coordinateToRow(coord, partRow, copyId=False)


def rowToClass(classRow, classItem):
    """ Method base to create a class2D, class3D or classVol from
    a row of a metadata
    """
    setObjId(classItem, classRow, label=xmipp.MDL_REF)

    if classRow.containsLabel(xmipp.MDL_IMAGE):
        index, filename = xmippToLocation(classRow.getValue(xmipp.MDL_IMAGE))
        img = classItem.REP_TYPE()
        classItem.setObjId(classRow.getObjId())
#         img.copyObjId(classItem)
        img.setLocation(index, filename)
        img.setSamplingRate(classItem.getSamplingRate())
        classItem.setRepresentative(img)
    
    return classItem


def class2DToRow(class2D, classRow):
    """ Set labels values from Class2D to md row. """

    if class2D.hasRepresentative():
        index, filename = class2D.getRepresentative().getLocation()
        fn = locationToXmipp(index, filename)
        classRow.setValue(xmipp.MDL_IMAGE, fn)
    n = long(len(class2D))
    classRow.setValue(xmipp.MDL_CLASS_COUNT, n)
    classRow.setValue(xmipp.MDL_REF, int(class2D.getObjId()))
    
        
def ctfModelToRow(ctfModel, ctfRow):
    """ Set labels values from ctfModel to md row. """
    objectToRow(ctfModel, ctfRow, CTF_DICT)


def defocusGroupSetToRow(defocusGroup, defocusGroupRow):
    """ Set labels values from ctfModel to md row. """
    objectToRow(defocusGroup, defocusGroupRow, CTF_DICT)


def rowToCtfModel(ctfRow):
    """ Create a CTFModel from a row of a metadata. """
    if containsLabels(ctfRow, CTF_DICT):
        ctfModel = CTFModel()
        rowToObject(ctfRow, ctfModel, CTF_DICT, CTF_EXTRA_LABELS)
        ctfModel.standardize()
    else:
        ctfModel = None
        
    return ctfModel


def acquisitionToRow(acquisition, ctfRow):
    """ Set labels values from acquisition to md row. """
    objectToRow(acquisition, ctfRow, ACQUISITION_DICT)


def rowToAcquisition(acquisitionRow):
    """ Create an acquisition from a row of a metadata. """
    if containsLabels(acquisitionRow, ACQUISITION_DICT):
        acquisition = Acquisition()
        rowToObject(acquisitionRow, acquisition, ACQUISITION_DICT) 
    else:                
        acquisition = None
    
    return acquisition
    
    
def readSetOfMicrographs(filename, micSet, **kwargs):    
    readSetOfImages(filename, micSet, rowToMicrograph, **kwargs)


def writeSetOfMicrographs(micSet, filename, rowFunc=None, blockName='Micrographs', **kwargs):
    writeSetOfImages(micSet, filename, micrographToRow, rowFunc, blockName, **kwargs)
    
    
def readSetOfVolumes(filename, volSet, **kwargs):    
    readSetOfImages(filename, volSet, rowToVolume, **kwargs)


def writeSetOfVolumes(volSet, filename, rowFunc=None, blockName='Volumes', **kwargs):
    writeSetOfImages(volSet, filename, volumeToRow, rowFunc, blockName, **kwargs)    
    
    
def readCTFModel(filename, mic):
    """ Read from Xmipp .ctfparam and create a CTFModel object. """
    md = xmipp.MetaData(filename)
    ctfRow = rowFromMd(md, md.firstObject())
    ctfObj = rowToCtfModel(ctfRow)
    ctfObj.setObjId(mic.getObjId())
    ctfObj.setMicrograph(mic)
    
    return ctfObj

        
def writeSetOfCoordinates(posDir, coordSet):
    """ Write a pos file on metadata format for each micrograph 
    on the coordSet. 
    Params:
        posDir: the directory where the .pos files will be written.
        coordSet: the SetOfCoordinates that will be read.
    """
    posFiles = []
    boxSize = coordSet.getBoxSize() or 100  
    
    print "writing coordinates to folder: ", posDir
    # Write pos metadatas (one per micrograph)    
    for mic in coordSet.iterMicrographs():
        micName = mic.getFileName()
        posFn = join(posDir, replaceBaseExt(micName, "pos"))
        print " posFn: ", posFn
        
        md = xmipp.MetaData()
        for coord in coordSet.iterCoordinates(micrograph=mic):
            objId = md.addObject()
            coordRow = XmippMdRow()
            coordinateToRow(coord, coordRow)
            coordRow.writeToMd(md, objId)
            
        if not md.isEmpty():
            md2 = xmipp.MetaData()    
            objId = md2.addObject()
            md2.setValue(xmipp.MDL_PICKING_MICROGRAPH_STATE, 'Manual', objId)
            # Write header block
            md2.write('header@%s' % posFn)
            # Write particles block
            md.write('particles@%s' % posFn, xmipp.MD_APPEND)
            posFiles.append(posFn)
            
    # Write config.xmd metadata
    configFn = join(posDir, 'config.xmd')
    md = xmipp.MetaData()
    # Write properties block
    objId = md.addObject()
    micName = removeBaseExt(micName)
    md.setValue(xmipp.MDL_MICROGRAPH, str(micName), objId)
    #md.setValue(xmipp.MDL_COLOR, int(-16776961), objId)
    md.setValue(xmipp.MDL_PICKING_PARTICLE_SIZE, int(boxSize), objId)
    md.setValue(xmipp.MDL_PICKING_STATE, 'Manual', objId)
    
    md.write('properties@%s' % configFn)

    # Write filters block
    md = xmipp.MetaData()    
    objId = md.addObject()
    md.setValue(xmipp.MDL_MACRO_CMD, 'Gaussian_Blur...', objId)
    md.setValue(xmipp.MDL_MACRO_CMD_ARGS, 'sigma=2', objId)
    md.write('filters@%s' % configFn, xmipp.MD_APPEND)
    
    return posFiles


def readSetOfCoordinates(outputDir, micSet, coordSet):
    """ Read from Xmipp .pos files.
    Params:
        outputDir: the directory where the .pos files are.
            It is also expected a file named: config.xmd
            in this directory where the box size can be read.
        micSet: the SetOfMicrographs to associate the .pos, which 
            name should be the same of the micrographs.
        coordSet: the SetOfCoordinates that will be populated.
    """
    # Read the boxSize from the config.xmd metadata
    configfile = join(outputDir, 'config.xmd')
    if exists(configfile):
        md = xmipp.MetaData('properties@' + join(outputDir, 'config.xmd'))
        boxSize = md.getValue(xmipp.MDL_PICKING_PARTICLE_SIZE, md.firstObject())
        coordSet.setBoxSize(boxSize)
    for mic in micSet:
        posFile = join(outputDir, replaceBaseExt(mic.getFileName(), 'pos'))
        scipionPosFile = join(outputDir, "scipion_" + replaceBaseExt(mic.getFileName(), 'pos'))
        posMd = readPosCoordinates(posFile)
        posMd.addLabel(xmipp.MDL_ITEM_ID)
        
        for objId in posMd:
            coord = rowToCoordinate(rowFromMd(posMd, objId))
            coord.setMicrograph(mic)
            coordSet.append(coord)      
            # Add an unique ID that will be propagated to particles
            posMd.setValue(xmipp.MDL_ITEM_ID, long(coord.getObjId()), objId)
        if not posMd.isEmpty():
            posMd.write("particles@%s"  % scipionPosFile)
    coordSet._xmippMd = String(outputDir)
    

def readPosCoordinates(posFile):
    """ Read the coordinates in .pos file and return corresponding metadata.
    There are two possible blocks with particles:
    particles: with manual/supervised particles
    particles_auto: with automatically picked particles.
    If posFile doesn't exist, the metadata will be empty 
    """
    md = xmipp.MetaData()
    
    if exists(posFile):
        blocks = xmipp.getBlocksInMetaDataFile(posFile)
        
        for b in ['particles', 'particles_auto']:
            if b in blocks:
                mdAux = xmipp.MetaData('%(b)s@%(posFile)s' % locals())
                md.unionAll(mdAux)
        md.removeDisabled()
    
    return md


def readSetOfImages(filename, imgSet, rowToFunc, **kwargs):
    """read from Xmipp image metadata.
        filename: The metadata filename where the image are.
        imgSet: the SetOfParticles that will be populated.
        rowToFunc: this function will be used to convert the row to Object
    """    
    imgMd = xmipp.MetaData(filename)
    imgMd.removeDisabled()
    # Read the sampling rate from the acquisition info file if exists
    
    acqFile = join(dirname(filename), 'acquisition_info.xmd')
    if exists(acqFile):
        md = RowMetaData(acqFile)
        samplingRate = md.getValue(xmipp.MDL_SAMPLING_RATE)
        imgSet.setSamplingRate(samplingRate)
    else:
        pass # TODO: what to do if not exists
    
    hasCtf = False
    hasAlignment = False
    
    for objId in imgMd:
        imgRow = rowFromMd(imgMd, objId)
        img = rowToFunc(imgRow, **kwargs)
        imgSet.append(img)
        hasCtf = img.hasCTF()
        hasAlignment = img.hasAlignment()
        
    imgSet.setHasCTF(hasCtf)
    imgSet.setHasAlignment(hasAlignment)
        

def setOfImagesToMd(imgSet, md, imgToFunc, rowFunc, **kwargs):
    """ This function will fill Xmipp metadata from a SetOfMicrographs
    Params:
        imgSet: the set of images to be converted to metadata
        md: metadata to be filled
        rowFunc: this function can be used to setup the row before 
            adding to metadata.
    """
    for img in imgSet:
        objId = md.addObject()
        imgRow = XmippMdRow()
        imgToFunc(img, imgRow, **kwargs)#drop ctfFn
        
        if rowFunc:
            rowFunc(img, imgRow)
        imgRow.writeToMd(md, objId)


def readAnglesFromMicrographs(micFile, anglesSet):
    """ Read the angles from a micrographs Metadata.
    """
    micMd = xmipp.MetaData(micFile)
#    micMd.removeDisabled()    
    
    for objId in micMd:
        angles = Angles()
        row = rowFromMd(micMd, objId)
        rowToObject(row, angles, ANGLES_DICT)
        angles.setObjId(micMd.getValue(xmipp.MDL_ITEM_ID, objId)) 
        anglesSet.append(angles)
    

def writeSetOfImages(imgSet, filename, imgToFunc, rowFunc, 
                     blockName='Images', **kwargs):
    """ This function will write a SetOfMicrographs as Xmipp metadata.
    Params:
        imgSet: the set of images to be written (particles, micrographs or volumes)
        filename: the filename where to write the metadata.
        rowFunc: this function can be used to setup the row before 
            adding to metadata.
    """
    md = xmipp.MetaData()
    setOfImagesToMd(imgSet, md, imgToFunc, rowFunc, **kwargs)
    md.write('%s@%s' % (blockName, filename))
        
        
def readSetOfParticles(filename, partSet, **kwargs):
    readSetOfImages(filename, partSet, rowToParticle, **kwargs)

def setOfParticlesToMd(imgSet, md, rowFunc=None, **kwargs):
    setOfImagesToMd(imgSet, md, particleToRow, rowFunc, **kwargs)


def setOfMicrographsToMd(imgSet, md, rowFunc=None, **kwargs):
    setOfImagesToMd(imgSet, md, micrographToRow, rowFunc, **kwargs)


def writeSetOfParticles(imgSet, filename, rowFunc=None, blockName='Particles', **kwargs):
    writeSetOfImages(imgSet, filename, particleToRow, rowFunc, blockName, **kwargs)


def writeSetOfCTFs(ctfSet, mdCTF):
    """ Write a ctfSet on metadata format. 
    Params:
        ctfSet: the SetOfCTF that will be read.
        mdCTF: The file where metadata should be written.
    """
    md = xmipp.MetaData()
            
    for ctfModel in ctfSet:
        objId = md.addObject()
        ctfRow = XmippMdRow()
        ctfRow.setValue(xmipp.MDL_MICROGRAPH, ctfModel.getMicFile())
        if ctfModel.getPsdFile():
            ctfRow.setValue(xmipp.MDL_PSD, ctfModel.getPsdFile())
        ctfModelToRow(ctfModel, ctfRow)
        ctfRow.writeToMd(md, objId)
        
    md.write(mdCTF)
    ctfSet._xmippMd = String(mdCTF)


def writeSetOfDefocusGroups(defocusGroupSet, fnDefocusGroup): # also metadata
    """ Write a defocuGroupSet on metadata format. 
    Params:
        defocusGroupSet: the SetOfDefocus that will be read.
        fnDefocusGroup: The file where defocusGroup should be written.
    """
    md = xmipp.MetaData()
            
    for defocusGroup in defocusGroupSet:
        objId = md.addObject()
        defocusGroupRow = XmippMdRow()
        defocusGroupSetToRow(defocusGroup, defocusGroupRow)
        defocusGroupRow.setValue(xmipp.MDL_CTF_GROUP, defocusGroup.getObjId())
        defocusGroupRow.setValue(xmipp.MDL_MIN, defocusGroup.getDefocusMin())
        defocusGroupRow.setValue(xmipp.MDL_MAX, defocusGroup.getDefocusMax())
        defocusGroupRow.setValue(xmipp.MDL_AVG, defocusGroup.getDefocusAvg())
        defocusGroupRow.writeToMd(md, objId)
        
    md.write(fnDefocusGroup)
    defocusGroupSet._xmippMd = String(fnDefocusGroup)


def writeSetOfClasses2D(classes2DSet, filename, classesBlock='classes'):    
    """ This function will write a SetOfClasses2D as Xmipp metadata.
    Params:
        classes2DSet: the SetOfClasses2D instance.
        filename: the filename where to write the metadata.
    """
    classFn = '%s@%s' % (classesBlock, filename)
    classMd = xmipp.MetaData()
    classMd.write(classFn) # Empty write to ensure the classes is the first block
    
    classRow = XmippMdRow()
    for class2D in classes2DSet:        
        class2DToRow(class2D, classRow)
        classRow.writeToMd(classMd, classMd.addObject())
        ref = class2D.getObjId()
        imagesFn = 'class%06d_images@%s' % (ref, filename)
        imagesMd = xmipp.MetaData()
        imgRow = XmippMdRow()
        
        for img in class2D:
            particleToRow(img, imgRow)
            imgRow.writeToMd(imagesMd, imagesMd.addObject())
        imagesMd.write(imagesFn, xmipp.MD_APPEND)
    
    classMd.write(classFn, xmipp.MD_APPEND) # Empty write to ensure the classes is the first block


def writeSetOfMicrographsPairs(uSet, tSet, filename):
    """ This function will write a MicrographsTiltPair as Xmipp metadata.
    Params:
        uSet: the untilted set of micrographs to be written
        tSet: the tilted set of micrographs to be written
        filename: the filename where to write the metadata.
    """
    md = xmipp.MetaData()

    for micU, micT in izip(uSet, tSet):
        objId = md.addObject()
        pairRow = XmippMdRow()
        pairRow.setValue(xmipp.MDL_ITEM_ID, long(micU.getObjId()))
        pairRow.setValue(xmipp.MDL_MICROGRAPH, micU.getFileName())
        pairRow.setValue(xmipp.MDL_MICROGRAPH_TILTED, micT.getFileName())
        pairRow.writeToMd(md, objId)
        
    md.write(filename)   

    
def readSetOfClasses(classesSet, filename, classesBlock='classes', **kwargs):
    """ Read a set of classes from an Xmipp metadata with the given
    convention of a block for classes and another block for each set of 
    images assigned to each class.
    Params:
        classesSet: the SetOfClasses object that will be populated.
        filename: the file path where to read the Xmipp metadata.
        classesBlock (by default 'classes'):
            the block name of the classes group in the metadata.
    """
    blocks = xmipp.getBlocksInMetaDataFile(filename)
    
    classesMd = xmipp.MetaData('%s@%s' % (classesBlock, filename))
    
    for objId in classesMd:
        classItem = classesSet.ITEM_TYPE()
        classItem.setObjId(objId)
        classItem = rowToClass(rowFromMd(classesMd, objId), classItem)
        # FIXME: the following is only valid for SetOfParticles
        SetOfParticles.copyInfo(classItem, classesSet.getImages())
        #classItem.copyInfo(classesSet.getImages())
        classesSet.append(classItem)
        ref = classItem.getObjId()
        b = 'class%06d_images' % ref
        
        if b in blocks:
            #FIXME: we need to adapt the following line
            # when we face classes of volumes and not just particles
            readSetOfParticles('%s@%s' % (b, filename), classItem, **kwargs) 
        # Update with new properties of classItem such as _size
        classesSet.update(classItem)
        
        
def readSetOfClasses2D(classes2DSet, filename, classesBlock='classes', **kwargs):
    """ Just a wrapper to readSetOfClasses. """
    readSetOfClasses(classes2DSet, filename, classesBlock, **kwargs)


def readSetOfClasses3D(classes3DSet, filename, classesBlock='classes', **kwargs):
    """ Just a wrapper to readSetOfClasses. """
    readSetOfClasses(classes3DSet, filename, classesBlock, **kwargs)


def writeSetOfClassesVol(classesVolSet, filename, classesBlock='classes'):    
    """ This function will write a SetOfClassesVol as Xmipp metadata.
    Params:
        classesVolSet: the SetOfClassesVol instance.
        filename: the filename where to write the metadata.
    """
    classFn = '%s@%s' % (classesBlock, filename)
    classMd = xmipp.MetaData()
    classMd.write(classFn) # Empty write to ensure the classes is the first block
    
    classRow = XmippMdRow()
    for classVol in classesVolSet:        
        classVolToRow(classVol, classRow)
        classRow.writeToMd(classMd, classMd.addObject())
        ref = class3D.getObjId()
        imagesFn = 'class%06d_images@%s' % (ref, filename)
        imagesMd = xmipp.MetaData()
        imgRow = XmippMdRow()
        
        for vol in classVol:
            volumeToRow(vol, imgRow)
            imgRow.writeToMd(imagesMd, imagesMd.addObject())
        imagesMd.write(imagesFn, xmipp.MD_APPEND)
  
    classMd.write(classFn, xmipp.MD_APPEND) # Empty write to ensure the classes is the first block


def readSetOfClassesVol(classesVolSet, filename, classesBlock='classes', **args):
    """read from Xmipp image metadata.
        fnImages: The metadata filename where the particles properties are.
        imgSet: the SetOfParticles that will be populated.
        hasCtf: is True if the ctf information exists.
    """
    blocks = xmipp.getBlocksInMetaDataFile(filename)
    classesMd = xmipp.MetaData('%s@%s' % (classesBlock, filename))
    samplingRate = classesVolSet.getImages().getSamplingRate()
    
    for objId in classesMd:
        classVol = ClassVol()
        classVol = rowToClassVol(classesMd, objId, classVol)
        classesVolSet.append(classVol)
        ref = classVol.getObjId()
        b = 'class%06d_images' % ref
        
        if b in blocks:
            classImagesMd = xmipp.MetaData('%s@%s' % (b, filename))
            
            for imgId in classImagesMd:
                img = rowToParticle(classImagesMd, imgId, hasCtf=False)
                img.setSamplingRate(samplingRate)
                classVol.append(img)
        # Update with new properties of class2D such as _size
        classesVolSet.update(classVol)                
        # Check if write function is necessary
        classVol.write()


def writeSetOfMovies(moviesSet, filename, moviesBlock='movies'):    
    """ This function will write a SetOfMovies as Xmipp metadata.
    Params:
        moviesSet: the SetOfMovies instance.
        filename: the filename where to write the metadata.
    """
       
    for movie in moviesSet:        
        
        ref = movie.getObjId()
        micrographsFn = 'movie%06d_micrographs@%s' % (ref, filename)
        micrographsMd = xmipp.MetaData()
        micRow = XmippMdRow()
        
        for mic in movie:
            micrographToRow(mic, micRow)
            micRow.writeToMd(micrographsMd, micrographsMd.addObject())
        micrographsMd.write(micrographsFn, xmipp.MD_APPEND)


def createXmippInputImages(prot, imgSet, rowFunc=None, imagesFn=None):  
    if prot is not None:
        imgsFn = prot._getPath(imagesFn or 'input_images.xmd')
    
    writeSetOfParticles(imgSet, imgsFn, rowFunc)
    return imgsFn


def createXmippInputMicrographs(prot, micSet, rowFunc=None, micsFn=None):    
    if prot is not None:
        micsFn = prot._getPath('input_micrographs.xmd')

    writeSetOfMicrographs(micSet, micsFn, rowFunc)
    return micsFn


def createXmippInputVolumes(prot, volSet, volsFn=None):    
    if volsFn is None:
        volsFn = prot._getPath('input_volumes.xmd')
    
    writeSetOfVolumes(volSet, volsFn)
    return volsFn


def createXmippInputClasses2D(prot, classSet, classFn=None):
    if prot is not None and classFn is None:
        classFn = prot._getPath('input_classes.xmd')
    
    writeSetOfClasses2D(classSet, classFn)
    return classFn
    

def createXmippInputCTF(prot, ctfSet, ctfFn=None):
    ctfMd = getattr(ctfSet, '_xmippMd', None)
    if ctfMd is None:
        if prot is not None:
            ctfFn = prot._getPath('input_ctfs.xmd')
        
        writeSetOfCTFs(ctfSet, ctfFn)
    else:
        ctfFn = ctfMd.get()
    return ctfFn

def geometryFromMatrix(matrix):
    from pyworkflow.em.transformations import translation_from_matrix, euler_from_matrix
    from numpy import rad2deg
    shifts = -translation_from_matrix(matrix)
    angles = rad2deg(euler_from_matrix(matrix, axes='szyz'))
    
    return shifts, angles


def matrixFromGeometry(shifts, angles):
    """ Create the transformation matrix from a given
    2D shifts in X and Y...and the 3 euler angles.
    """
    #TODO: CHECK THIS IS CORRECT
    from pyworkflow.em.transformations import translation_matrix, euler_matrix, concatenate_matrices
    from numpy import deg2rad
    radAngles = deg2rad(angles)
    
    R = euler_matrix(radAngles[0], radAngles[1], radAngles[2], 'szyz')
    T = translation_matrix(shifts)
    
    M = concatenate_matrices(R, T)
    
    return M


def rowToAlignment(alignmentRow):
    """ Fill the alignment matrix reading
    the geometry values from Xmipp metadata row. """
    if containsLabels(alignmentRow, ALIGNMENT_EXTRA_LABELS):
        alignment = Alignment()
        angles = numpy.zeros(3)
        shifts = numpy.zeros(3)
        angles[2] = alignmentRow.getValue(xmipp.MDL_ANGLE_PSI)
        shifts[0] = alignmentRow.getValue(xmipp.MDL_SHIFT_X)
        shifts[1] = alignmentRow.getValue(xmipp.MDL_SHIFT_Y)
        flip = alignmentRow.getValue(xmipp.MDL_FLIP)
        
        M = matrixFromGeometry(shifts, angles)
        alignment.setMatrix(M)
        #FIXME: remove this after the conversions from Transform matrix
        # is completed tested
        rowToObject(alignmentRow, alignment, 
                    {}, ALIGNMENT_EXTRA_LABELS)
    else:
        alignment = None
    
    return alignment


def alignmentToRow(alignment, mdRow):
    #FIXME: we should use the transformation matrix
    #shifts, angles = geometryFromMatrix(alignment.getMatrix())
    #mdRow.setValue(xmipp.MDL_ANGLE_PSI, angles[2])
    #mdRow.setValue(xmipp.MDL_SHIFT_X, shifts[0])
    #mdRow.setValue(xmipp.MDL_SHIFT_Y, shifts[1])
    mdRow.setValue(xmipp.MDL_ANGLE_PSI, alignment._xmipp_anglePsi.get())
    mdRow.setValue(xmipp.MDL_SHIFT_X, alignment._xmipp_shiftX.get())
    mdRow.setValue(xmipp.MDL_SHIFT_Y, alignment._xmipp_shiftY.get())
    mdRow.setValue(xmipp.MDL_FLIP, bool(alignment.getAttributeValue('_xmipp_flip', False)))     


def createClassesFromImages(inputImages, inputMd, classesFn, ClassType, 
                            classLabel, classFnTemplate, iter, preprocessRow=None):
    """ From an intermediate X.xmd file produced by xmipp, create
    the set of classes in which those images are classified.
    Params:
        inputImages: the SetOfImages that were initially classified by relion. 
        inputMd: the filename metadata.
        classesFn: filename where to write the classes.
        ClassType: the output type of the classes set ( usually SetOfClass2D or 3D )
        classLabel: label that is the class reference.
        classFnTemplate: the template to get the classes averages filenames
        iter: the iteration number, just used in Class template
    """
    # We asume here that the volumes (classes3d) are in the same folder than imgsFn
    # rootDir here is defined to be used expanding locals()
    if "@" in inputMd:
        inputFn = inputMd.split('@')[1]
        tmpDir = os.path.dirname(inputFn)
    else:
        tmpDir = os.path.dirname(inputMd)
    rootDir = tmpDir
    md = xmipp.MetaData(inputMd)
    clsDict = {} # Dictionary to store the (classId, classSet) pairs
    clsSet = ClassType(filename=classesFn)
    clsSet.setImages(inputImages)
    hasCtf = inputImages.hasCTF()
    
    for objId in md:
        ref = md.getValue(classLabel, objId)
        if not ref in clsDict: # Register a new class set if the ref was not found.
            cls = clsSet.ITEM_TYPE(objId=ref)
            refFn = classFnTemplate % locals()
            refLocation = xmippToLocation(refFn)
            rep = clsSet.REP_TYPE()
            rep.setLocation(refLocation)
            cls.setRepresentative(rep)
            
            clsDict[ref] = cls
            cls.copyInfo(inputImages)
            cls.setAcquisition(inputImages.getAcquisition())
            clsSet.append(cls)
        classItem = clsDict[ref] # Try to get the class set given its ref number
        # Set images attributes from the md row values
        imgRow = rowFromMd(md, objId)
        if preprocessRow:
            preprocessRow(imgRow)
        img = rowToParticle(imgRow)
        
        classItem.append(img)
        
    for classItem in clsDict.values():
        clsSet.update(classItem)
        
    clsSet.write()
    
