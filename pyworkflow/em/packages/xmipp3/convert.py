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
from data import *
from xmipp3 import XmippMdRow
from pyworkflow.em.constants import NO_INDEX

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
               "_id": xmipp.MDL_ITEM_ID,
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
    objectToRow(class2D, classRow, classDict)
    
def imageClassAssignmentToRow(imgCA, imgCARow, img,  ctfDir):
    """ Set label values from ImageClassAssignment to md row. """
    imgCADict = { 
               "_id": xmipp.MDL_ITEM_ID,
#               "_anglePsi": xmipp.MDL_ANGLE_PSI,
#               "_shiftX": xmipp.MDL_SHIFT_X,
#               "_shiftY": xmipp.MDL_SHIFT_Y,
#               "_flip": xmipp.MDL_FLIP
               }
    index, filename = img.getLocation()
    fn = locationToXmipp(index, filename)
    imgCARow.setValue(xmipp.MDL_IMAGE, fn)
    
    if img.hasCTF():
        rootFn = "%06d_%s" % (img.getId(), replaceBaseExt(img.getFileName(), "ctfparam"))
        ctfFn = join(ctfDir, rootFn)
        writeCTFModel(img.getCTF(), ctfFn)
        imgCARow.setValue(xmipp.MDL_CTF_MODEL, ctfFn)
    
    objectToRow(img, imgCARow, imgCADict)
    
def readSetOfMicrographs(filename, micSet, hasCtf=False):
    
    micMd = xmipp.MetaData(filename)
    for objId in micMd:
        mic = rowToMicrograph(micMd, objId, hasCtf)
        micSet.append(mic)


def writeSetOfMicrographs(micSet, filename, ctfDir=None, rowFunc=None):
    """ This function will write a SetOfMicrographs as Xmipp metadata.
    Params:
        micSet: the SetOfMicrograph instance.
        filename: the filename where to write the metadata.
        rowFunc: this function can be used to setup the row before 
            adding to metadata.
    """
    md = xmipp.MetaData()
    hasCtf = micSet.hasCTF()
    
    for mic in micSet:
        objId = md.addObject()
        micRow = XmippMdRow()
        micrographToRow(mic, micRow, ctfDir, hasCtf)
        if rowFunc:
            rowFunc(mic, micRow)
        micRow.writeToMd(md, objId)
        
    md.write(filename)
    micSet._xmippMd = String(filename)
    
    
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

def writeSetOfParticles(imgSet, filename, ctfDir=None, rowFunc=None):
    
    """ This function will write a SetOfParticles as Xmipp metadata.
    Params:
        imgSet: the SetOfParticles instance.
        filename: the filename where to write the metadata.
        rowFunc: this function can be used to setup the row before 
            adding to metadata.
    """
    
#     mdFn = getattr(obj, '_xmippMd', None)
#     if mdFn:
#         fn = mdFn.get()
#     else:
#         fn = self._getTmpPath(obj.getName() + '_micrographs.xmd')
#         writeSetOfMicrographs(obj, fn)
    
    md = xmipp.MetaData()
    hasCtf = imgSet.hasCTF()
    
    for img in imgSet:
        objId = md.addObject()
        imgRow = XmippMdRow()
        particleToRow(img, imgRow, ctfDir, hasCtf)
        if rowFunc:
            rowFunc(img, imgRow)
        imgRow.writeToMd(md, objId)
        
    md.write(filename)
    imgSet._xmippMd = String(filename)


def readSetOfParticles(fnImages, imgSet, hasCtf=False):
    """read from Xmipp image metadata.
        fnImages: The metadata filename where the particles properties are.
        imgSet: the SetOfParticles that will be populated.
        hasCtf: is True if the ctf information exists.
    """
    
    imgMd = xmipp.MetaData(fnImages)
    for objId in imgMd:
        part = rowToParticle(imgMd, objId, hasCtf)
        imgSet.append(part)
        
def writeSetOfClasses2D(classes2DSet, filename, ctfDir=None, classesBlock='classes'):
    
    """ This function will write a SetOfClasses2D as Xmipp metadata.
    Params:
        classes2DSet: the SetOfClasses2D instance.
        filename: the filename where to write the metadata.
    """
    md = xmipp.MetaData()
    
    images = classes2DSet.getImages()
        
    for class2D in classes2DSet.iterClasses():
        classRow = XmippMdRow()
        class2DToRow(class2D, classRow)
        objId = md.addObject()
        classRow.writeToMd(md, objId)
        md.write('%s@%s' %(classesBlock, filename))
        ref = class2D.getId() 
        for imgCA in class2D:
            #FIXME: This doesnt work, _idMap is empty (and not even defined!)
            #img = images[imgCA.getImageId()]
            for img1 in images:
                if img1.getId() == imgCA.getImageId():
                    img = img1
                    break
                
            imgCARow = XmippMdRow()
            imageClassAssignmentToRow(imgCA, imgCARow, img, ctfDir)
            objCAId = md.addObject()
            assignmentBlock = 'class%06d_images@%s' % (ref, filename)
            classRow.writeToMd(md, objCAId)  
            md.write(assignmentBlock, xmipp.MD_APPEND)
            
    classes2DSet._xmippMd = String(filename)


def readSetOfClasses2D(classes2DSet, filename, classesBlock='classes', **args):
    """read from Xmipp image metadata.
        fnImages: The metadata filename where the particles properties are.
        imgSet: the SetOfParticles that will be populated.
        hasCtf: is True if the ctf information exists.
    """
    classesMd = xmipp.MetaData('%s@%s' %(classesBlock, filename))
    
    samplingRate = classes2DSet.getImages().getSamplingRate()
    
    for objId in classesMd:           
        class2D = rowToClass2D(classesMd, objId, samplingRate)
        ref = classesMd.getValue(xmipp.MDL_REF, objId)
        assignmentBlock = 'class%06d_images@%s' % (ref, filename)
        imgAssignmentMd = xmipp.MetaData(assignmentBlock)
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
