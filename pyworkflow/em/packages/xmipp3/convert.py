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

import xmipp
from xmipp3 import XmippMdRow
from pyworkflow.em import *
from pyworkflow.em.constants import NO_INDEX
from pyworkflow.object import String
from pyworkflow.utils.path import join, dirname, replaceBaseExt

LABEL_TYPES = { 
               xmipp.LABEL_SIZET: long,
               xmipp.LABEL_DOUBLE: float,
               xmipp.LABEL_INT: int,
               xmipp.LABEL_BOOL: bool              
               }

def objectToRow(obj, row, attrDict):
    """ This function will convert an EMObject into a XmippMdRow.
    Params:
        obj: the EMObject instance (input)
        row: the XmippMdRow instance (output)
        attrDict: dictionary with the map between obj attributes(keys) and 
            row MDLabels in Xmipp (values).
    """
    for attr, label in attrDict.iteritems():
        if hasattr(obj, attr):
            labelType = xmipp.labelType(label)
            valueType = LABEL_TYPES.get(labelType, str)
            row.setValue(label, valueType(getattr(obj, attr).get()))

def _rowToObject(row, obj, attrDict):
    """ This function will convert from a XmippMdRow to an EMObject.
    Params:
        row: the XmippMdRow instance (input)
        obj: the EMObject instance (output)
        attrDict: dictionary with the map between obj attributes(keys) and 
            row MDLabels in Xmipp (values).
    """
    for attr, label in attrDict.iteritems():
        if not hasattr(obj, attr):
            setattr(obj, attr, String()) #TODO: change string for the type of label
        getattr(obj, attr).set(row.getValue(label))
    
def rowToObject(md, objId, obj, attrDict):
    """ Same as rowToObject, but creating the row from md and objId. """
    row = XmippMdRow()
    row.readFromMd(md, objId)
    _rowToObject(row, obj, attrDict)
    
def readCTFModel(filename):
    """ Read from Xmipp .ctfparam and create a CTFModel object. """
    md = xmipp.MetaData(filename)
    ctfObj = CTFModel()    
    ctfDict = { 
               "defocusU": xmipp.MDL_CTF_DEFOCUSU,
               "defocusV": xmipp.MDL_CTF_DEFOCUSV,
               "defocusAngle": xmipp.MDL_CTF_DEFOCUS_ANGLE,
               "sphericalAberration": xmipp.MDL_CTF_CS
               }    
    rowToObject(md, md.firstObject(), ctfObj, ctfDict)  
    ctfObj._xmippMd = String(filename) 
    
    return ctfObj


def writeCTFModel(ctfObj, filename):
    """ Write a CTFModel object as Xmipp .ctfparam"""
    ctfDict = { 
           "defocusU": xmipp.MDL_CTF_DEFOCUSU,
           "defocusV": xmipp.MDL_CTF_DEFOCUSV,
           "defocusAngle": xmipp.MDL_CTF_DEFOCUS_ANGLE,
           "sphericalAberration": xmipp.MDL_CTF_CS
           }   
    
    md = xmipp.MetaData()
    md.setColumnFormat(False)
    objId = md.addObject()
    ctfRow = XmippMdRow() 
    objectToRow(ctfObj, ctfRow, ctfDict)
    ctfRow.writeToMd(md, objId)
    md.write(filename)


def locationToXmipp(index, filename):
    """ Convert an index and filename location
    to a string with @ as expected in Xmipp.
    """
    #TODO: Maybe we need to add more logic dependent of the format
    if index != NO_INDEX:
        return "%d@%s" % (index, filename)
    
    return filename

def xmippToLocation(xmippFilename):
    """ Return a location (index, filename) given
    a Xmipp filename with the index@filename structure. """
    if '@' in xmippFilename:
        return xmipp.FileName(xmippFilename).decompose()
    else:
        return NO_INDEX, str(xmippFilename)

def imageToRow(img, imgRow, ctfDir, hasCtf, imgLabel):
    imgDict = { 
               "_id": xmipp.MDL_ITEM_ID,
               }
    index, filename = img.getLocation()
    fn = locationToXmipp(index, filename)
    if hasCtf:
        rootFn = "%06d_%s" % (img.getId(), replaceBaseExt(img.getFileName(), "ctfparam"))
        ctfFn = join(ctfDir, rootFn)
        writeCTFModel(img.getCTF(), ctfFn)
        imgRow.setValue(xmipp.MDL_CTF_MODEL, ctfFn)
    imgRow.setValue(imgLabel, fn)
    objectToRow(img, imgRow, imgDict)
    
def rowToImage(md, objId, imgLabel, imgClass, hasCtf):
    """ Create a Particle from a row of a metadata. """
    imgDict = { 
               "_id": xmipp.MDL_ITEM_ID,
               }
    img = imgClass()
    # Decompose Xmipp filename
    index, filename = xmippToLocation(md.getValue(imgLabel, objId))
    img.setLocation(index, filename)
    if hasCtf:
        ctFilename = md.getValue(xmipp.MDL_CTF_MODEL, objId)
        img.setCTF(readCTFModel(ctFilename))
    rowToObject(md, objId, img, imgDict)
    
    return img
    
def micrographToRow(mic, micRow, ctfDir, hasCtf):
    """ Set labels values from Micrograph mic to md row. """
    imageToRow(mic, micRow, ctfDir, hasCtf, imgLabel=xmipp.MDL_MICROGRAPH)
    
def rowToMicrograph(md, objId, hasCtf):
    """ Create a Micrograph object from a row of Xmipp metadata. """
    return rowToImage(md, objId, xmipp.MDL_MICROGRAPH, Micrograph, hasCtf)

def volumeToRow(vol, volRow):
    """ Set labels values from Micrograph mic to md row. """
    imageToRow(vol, volRow, ctfDir=None, hasCtf=False, imgLabel=xmipp.MDL_IMAGE)
    
def rowToVolume(md, objId):
    """ Create a Micrograph object from a row of Xmipp metadata. """
    return rowToImage(md, objId, xmipp.MDL_IMAGE, Volume, False)

def coordinateToRow(coord, coordRow):
    """ Set labels values from Coordinate coord to md row. """
    coordDict = { 
               "_id": xmipp.MDL_ITEM_ID,
               "_x": xmipp.MDL_XCOOR,
               "_y": xmipp.MDL_YCOOR
               }
#    x, y = coord.getPosition()
#    coordRow.setValue(xmipp.MDL_XCOOR, x)
#    coordRow.setValue(xmipp.MDL_YCOOR, y)
      
    objectToRow(coord, coordRow, coordDict)

def rowToCoordinate(md, objId):
    """ Create a Coordinate from a row of a metadata. """
    coordDict = { 
               "_id": xmipp.MDL_ITEM_ID,
               "_x": xmipp.MDL_XCOOR,
               "_y": xmipp.MDL_YCOOR,
               }
    coord = Coordinate()
    rowToObject(md, objId, coord, coordDict)
    
    return coord

def rowToParticle(md, objId, hasCtf):
    """ Create a Particle from a row of a metadata. """
    return rowToImage(md, objId, xmipp.MDL_IMAGE, Particle, hasCtf)
    
def particleToRow(part, partRow, ctfDir, hasCtf):
    """ Set labels values from Particle to md row. """
    imageToRow(part, partRow, ctfDir, hasCtf, imgLabel=xmipp.MDL_IMAGE)

def rowToClass2D(md, objId, samplingRate):
    """ Create a Class2D from a row of a metadata. """
    classDict = { 
               "_id": xmipp.MDL_REF,
               }
    class2D = Class2D()

    if md.containsLabel(xmipp.MDL_IMAGE):
        index, filename = xmippToLocation(md.getValue(xmipp.MDL_IMAGE, objId))
        img = Image()
        img.setLocation(index, filename)
        img.setSamplingRate(samplingRate)
        class2D.setHasRepresentativeImage(True)
        class2D.setRepresentativeImage(img)
    
    rowToObject(md, objId, class2D, classDict) 
    
    return class2D
    
def rowToImageClassAssignment(md, objId):
    """ Create a ImageClassAssignment from a row of a metadata. """
    imgCADict = { 
               "_imgId": xmipp.MDL_ITEM_ID,
#               "_anglePsi": xmipp.MDL_ANGLE_PSI,
#               "_shiftX": xmipp.MDL_SHIFT_X,
#               "_shiftY": xmipp.MDL_SHIFT_Y,
#               "_flip": xmipp.MDL_FLIP
               }
    
    imageCA = ImageClassAssignment()
    rowToObject(md, objId, imageCA, imgCADict) 
    return imageCA
    
def class2DToRow(class2D, classRow):
    """ Set labels values from Class2D to md row. """
    classDict = { 
               "_id": xmipp.MDL_REF,
               }
    if class2D.getHasRepresentativeImage():
        index, filename = class2D.getRepresentativeImage().getLocation()
        fn = locationToXmipp(index, filename)
        classRow.setValue(xmipp.MDL_IMAGE, fn)
    n = long(len(class2D.getImageAssignments()))
    classRow.setValue(xmipp.MDL_CLASS_COUNT, n)
    objectToRow(class2D, classRow, classDict)
    
def imageClassAssignmentToRow(imgCA, imgCARow, img, ctfDir, ref):
    """ Set label values from ImageClassAssignment to md row. """
    imgCADict = { 
               "_imgId": xmipp.MDL_ITEM_ID,
               #"_ref": xmipp.MDL_REF,
#               "_anglePsi": xmipp.MDL_ANGLE_PSI,
#               "_shiftX": xmipp.MDL_SHIFT_X,
#               "_shiftY": xmipp.MDL_SHIFT_Y,
#               "_flip": xmipp.MDL_FLIP
               }
    
    index, filename = img.getLocation()
    fn = locationToXmipp(index, filename)
    imgCARow.setValue(xmipp.MDL_IMAGE, fn)
    imgCARow.setValue(xmipp.MDL_REF, ref)
    
    if img.hasCTF():
        rootFn = "%06d_%s" % (img.getId(), replaceBaseExt(img.getFileName(), "ctfparam"))
        ctfFn = join(ctfDir, rootFn)
        writeCTFModel(img.getCTF(), ctfFn)
        imgCARow.setValue(xmipp.MDL_CTF_MODEL, ctfFn)
        
    objectToRow(imgCA, imgCARow, imgCADict)
    
    

def readSetOfMicrographs(filename, micSet, hasCtf=False):    
    readSetOfImages(filename, micSet, rowToMicrograph, hasCtf)

def writeSetOfMicrographs(micSet, filename, ctfDir=None, rowFunc=None):
    writeSetOfImages(micSet, filename, micrographToRow, ctfDir, rowFunc)
    
    
def readSetOfVolumes(filename, volSet, hasCtf=False):    
    readSetOfImages(filename, volSet, rowToVolume, False)

def writeSetOfVolumes(volSet, filename):
    writeSetOfImages(volSet, filename, volumeToRow)    
    
    
def readPosCoordinates(posFile):
    """ Read the coordinates in .pos file. 
    and return corresponding metadata. 
    """
    md = xmipp.MetaData()
    blocks = xmipp.getBlocksInMetaDataFile(posFile)
    
    for b in ['particles', 'particles_auto']:
        if b in blocks:
            mdAux = xmipp.MetaData('%(b)s@%(posFile)s' % locals())
            md.unionAll(mdAux)
    
    return md

def writePosCoordinates(posDir, coordSet):
    """ Write a pos file on metadata format for each micrograph 
    on the coordSet. 
    Params:
        posDir: the directory where the .pos files will be written.
        coordSet: the SetOfCoordinates that will be read.
    """
    posFiles = []
    for mic in coordSet.iterMicrographs():
        posFn = join(posDir, replaceBaseExt(mic.getFileName(), "pos"))
        md = xmipp.MetaData()
        for coord in coordSet.iterCoordinates(micrograph=mic):
            objId = md.addObject()
            coordRow = XmippMdRow()
            coordinateToRow(coord, coordRow)
            coordRow.writeToMd(md, objId)
        md.write(posFn)
        posFiles.append(posFn)
        
    return posFiles
    
def readSetOfCoordinates(posDir, micSet, coordSet):
    """ Read from Xmipp .pos files.
    Params:
        posDir: the directory where the .pos files are.
            It is also expected a file named: config.xmd
            in this directory where the box size can be read.
        micSet: the SetOfMicrographs to associate the .pos, which 
            name should be the same of the micrographs.
        coordSet: the SetOfCoordinates that will be populated.
    """
    # Read the boxSize from the config.xmd metadata
    md = xmipp.MetaData('properties@' + join(posDir, 'config.xmd'))
    boxSize = md.getValue(xmipp.MDL_PICKING_PARTICLE_SIZE, md.firstObject())
    
    for mic in micSet:
        posFile = join(posDir, replaceBaseExt(mic.getFileName(), 'pos'))
        posMd = readPosCoordinates(posFile)
        
        for objId in posMd:
            coord = rowToCoordinate(posMd, objId)
            coord.setMicrograph(mic)
            coordSet.append(coord)

    coordSet.setBoxSize(boxSize)


def readSetOfImages(filename, imgSet, rowToFunc, hasCtf):
    """read from Xmipp image metadata.
        filename: The metadata filename where the image are.
        imgSet: the SetOfParticles that will be populated.
        rowToFunc: this function will be used to convert the row to Object
        hasCtf: is True if the ctf information exists.
    """    
    imgMd = xmipp.MetaData(filename)
    for objId in imgMd:
        img = rowToFunc(imgMd, objId, hasCtf)
        imgSet.append(img)   
         
def writeSetOfImages(imgSet, filename, imgToFunc, ctfDir, rowFunc):
    """ This function will write a SetOfMicrographs as Xmipp metadata.
    Params:
        imgSet: the SetOfMicrograph instance.
        filename: the filename where to write the metadata.
        ctfDir: where to write ctfparam files if necessary
        rowFunc: this function can be used to setup the row before 
            adding to metadata.
    """
    md = xmipp.MetaData()
    hasCtf = imgSet.hasCTF()
    if ctfDir is None:
        ctfDir = dirname(filename)
    
    for img in imgSet:
        objId = md.addObject()
        imgRow = XmippMdRow()
        imgToFunc(img, imgRow, ctfDir, hasCtf)
        if rowFunc:
            rowFunc(img, imgRow)
        imgRow.writeToMd(md, objId)
        
    md.write(filename)
    imgSet._xmippMd = String(filename)
        

def readSetOfParticles(filename, partSet, hasCtf=False):
    readSetOfImages(filename, partSet, rowToParticle, hasCtf)
        
def writeSetOfParticles(imgSet, filename, ctfDir=None, rowFunc=None):
    writeSetOfImages(imgSet, filename, particleToRow, ctfDir, rowFunc)
        
        
def writeSetOfClasses2D(classes2DSet, filename, ctfDir=None, classesBlock='classes'):    
    """ This function will write a SetOfClasses2D as Xmipp metadata.
    Params:
        classes2DSet: the SetOfClasses2D instance.
        filename: the filename where to write the metadata.
    """
    imgSet = classes2DSet.getImages()
    imgSet.loadIfEmpty()
    
    classFn = '%s@%s' % (classesBlock, filename)
    classMd = xmipp.MetaData()
    classMd.write(classFn) # Empty write to ensure the classes is the first block
    if ctfDir is None:
        ctfDir = dirname(filename)
    
    for class2D in classes2DSet:        
        classRow = XmippMdRow()
        class2DToRow(class2D, classRow)
        classRow.writeToMd(classMd, classMd.addObject())
        ref = class2D.getId()
        imagesFn = 'class%06d_images@%s' % (ref, filename)
        imagesMd = xmipp.MetaData()
        
        for imgCA in class2D.getImageAssignments():
            imgCARow = XmippMdRow()
            img = imgSet[imgCA.getImageId()]
            imageClassAssignmentToRow(imgCA, imgCARow, img, ctfDir, ref)
            imgCARow.writeToMd(imagesMd, imagesMd.addObject())
        imagesMd.write(imagesFn, xmipp.MD_APPEND)
     
    classMd.write(classFn, xmipp.MD_APPEND) # Empty write to ensure the classes is the first block
    classes2DSet._xmippMd = String(filename)


def readSetOfClasses2D(classes2DSet, filename, classesBlock='classes', **args):
    """read from Xmipp image metadata.
        fnImages: The metadata filename where the particles properties are.
        imgSet: the SetOfParticles that will be populated.
        hasCtf: is True if the ctf information exists.
    """
    blocks = xmipp.getBlocksInMetaDataFile(filename)
    classesMd = xmipp.MetaData('%s@%s' %(classesBlock, filename))
    
    samplingRate = classes2DSet.getImages().getSamplingRate()
    
    for objId in classesMd:           
        class2D = rowToClass2D(classesMd, objId, samplingRate)
        ref = classesMd.getValue(xmipp.MDL_REF, objId)
        b = 'class%06d_images' % ref
        
        if b in blocks:
            imgAssignmentMd = xmipp.MetaData('%s@%s' % (b, filename))
            
            for objCAId in imgAssignmentMd:
                imgCA = ImageClassAssignment()
                imgCA = rowToImageClassAssignment(imgAssignmentMd, objCAId)
                class2D.addImageClassAssignment(imgCA)
        classes2DSet.append(class2D)
            

def createXmippInputImages(self, imgSet, rowFunc=None):    
    imgsFn = self._getPath('input_images.xmd')
    ctfDir = self._getExtraPath()
    writeSetOfParticles(imgSet, imgsFn, ctfDir, rowFunc)
    return imgsFn


def createXmippInputMicrographs(self, micSet, rowFunc=None):    
    micsFn = self._getPath('input_micrographs.xmd')
    ctfDir = self._getExtraPath()
    writeSetOfMicrographs(micSet, micsFn, ctfDir, rowFunc)
    return micsFn
