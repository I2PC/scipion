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
from xmipp3 import XmippMdRow, getLabelPythonType
from pyworkflow.em import *
from pyworkflow.em.constants import NO_INDEX
from pyworkflow.object import String
from pyworkflow.utils.path import join, dirname, replaceBaseExt



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
            valueType = getLabelPythonType(label)
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
    
def rowFromMd(md, objId):
    row = XmippMdRow()
    row.readFromMd(md, objId)
    return row
    
def rowToObject(md, objId, obj, attrDict):
    """ Same as rowToObject, but creating the row from md and objId. """
    _rowToObject(rowFromMd(md, objId), obj, attrDict)
    
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

def setObjId(obj, mdRow, label=xmipp.MDL_ITEM_ID):
    obj.setObjId(mdRow.getValue(label))
    
def setRowId(mdRow, obj, label=xmipp.MDL_ITEM_ID):
    mdRow.setValue(label, long(obj.getObjId()))

def imageToRow(img, imgRow, ctfFn, imgLabel):
    index, filename = img.getLocation()
    fn = locationToXmipp(index, filename)
    if ctfFn is not None:
        imgRow.setValue(xmipp.MDL_CTF_MODEL, ctfFn)
        #writeCTFModel(img.getCTF(), filename)
    imgRow.setValue(imgLabel, fn)
    setRowId(imgRow, img)
        
def rowToImage(md, objId, imgLabel, imgClass, hasCtf):
    """ Create a Particle from a row of a metadata. """
    img = imgClass()
    # Decompose Xmipp filename
    index, filename = xmippToLocation(md.getValue(imgLabel, objId))
    img.setLocation(index, filename)
    if hasCtf:
        ctFilename = md.getValue(xmipp.MDL_CTF_MODEL, objId)
        ctfModel = readCTFModel(ctFilename)
        ctfModel.micFile.set(md.getValue(xmipp.MDL_MICROGRAPH, objId))
        img.setCTF(ctfModel)
    setObjId(img, rowFromMd(md, objId))
    
    return img
    
def micrographToRow(mic, micRow, ctfFn):
    """ Set labels values from Micrograph mic to md row. """
    imageToRow(mic, micRow, ctfFn, imgLabel=xmipp.MDL_MICROGRAPH)
    
def rowToMicrograph(md, objId, hasCtf):
    """ Create a Micrograph object from a row of Xmipp metadata. """
    return rowToImage(md, objId, xmipp.MDL_MICROGRAPH, Micrograph, hasCtf)

def volumeToRow(vol, volRow, ctfFn):
    """ Set labels values from Micrograph mic to md row. """
    imageToRow(vol, volRow, ctfFn=None, imgLabel=xmipp.MDL_IMAGE)
    
def rowToVolume(md, objId, hasCtf):
    """ Create a Volume object from a row of Xmipp metadata. """
    return rowToImage(md, objId, xmipp.MDL_IMAGE, Volume, False)

def coordinateToRow(coord, coordRow):
    """ Set labels values from Coordinate coord to md row. """
    setRowId(coordRow, coord)
    objectToRow(coord, coordRow, {"_x": xmipp.MDL_XCOOR, "_y": xmipp.MDL_YCOOR })

def rowToCoordinate(md, objId):
    """ Create a Coordinate from a row of a metadata. """
    coord = Coordinate()
    #setObjId(coord, rowFromMd(md, objId))
    rowToObject(md, objId, coord, {"_x": xmipp.MDL_XCOOR, "_y": xmipp.MDL_YCOOR })
    
    return coord

def rowToParticle(md, objId, hasCtf):
    """ Create a Particle from a row of a metadata. """
    return rowToImage(md, objId, xmipp.MDL_IMAGE, Particle, hasCtf)
    
def particleToRow(part, partRow, ctfFn):
    """ Set labels values from Particle to md row. """
    imageToRow(part, partRow, ctfFn, imgLabel=xmipp.MDL_IMAGE)

def rowToClass2D(md, objId, class2D):
    """ Create a Class2D from a row of a metadata. """
    setObjId(class2D, rowFromMd(md, objId), label=xmipp.MDL_REF)

    if md.containsLabel(xmipp.MDL_IMAGE):
        index, filename = xmippToLocation(md.getValue(xmipp.MDL_IMAGE, objId))
        img = Particle()
        img.copyObjId(class2D)
        img.setLocation(index, filename)
        class2D.setAverage(img)
    
    return class2D
    
    
def class2DToRow(class2D, classRow):
    """ Set labels values from Class2D to md row. """

    if class2D.hasAverage():
        index, filename = class2D.getAverage().getLocation()
        fn = locationToXmipp(index, filename)
        classRow.setValue(xmipp.MDL_IMAGE, fn)
    n = long(len(class2D.getImages()))
    classRow.setValue(xmipp.MDL_CLASS_COUNT, n)
    setRowId(classRow, class2D, label=xmipp.MDL_REF)
        
def ctfModelToRow(ctfModel, ctfRow):
    """ Set labels values from ctfModel to md row. """    
    ctfDict = { 
           "defocusU": xmipp.MDL_CTF_DEFOCUSU,
           "defocusV": xmipp.MDL_CTF_DEFOCUSV,
           "defocusAngle": xmipp.MDL_CTF_DEFOCUS_ANGLE,
           "sphericalAberration": xmipp.MDL_CTF_CS
           }   
    
    objectToRow(ctfModel, ctfRow, ctfDict)

def rowToCtfModel(md, objId):
    """ Create a CTFModel from a row of a metadata. """
    ctfDict = { 
           "defocusU": xmipp.MDL_CTF_DEFOCUSU,
           "defocusV": xmipp.MDL_CTF_DEFOCUSV,
           "defocusAngle": xmipp.MDL_CTF_DEFOCUS_ANGLE,
           "sphericalAberration": xmipp.MDL_CTF_CS
           }   
    
    ctfModel = CTFModel()
    rowToObject(md, objId, ctfModel, ctfDict) 
    return ctfModel
    
def readSetOfMicrographs(filename, micSet, hasCtf=False):    
    readSetOfImages(filename, micSet, rowToMicrograph, hasCtf)

def writeSetOfMicrographs(micSet, filename, ctfDir=None, rowFunc=None):
    writeSetOfImages(micSet, filename, micrographToRow, ctfDir, rowFunc)
    
    
def readSetOfVolumes(filename, volSet, hasCtf=False):    
    readSetOfImages(filename, volSet, rowToVolume, False)

def writeSetOfVolumes(volSet, filename):
    writeSetOfImages(volSet, filename, volumeToRow, None, None)    
    
def readCTFModel(filename):
    """ Read from Xmipp .ctfparam and create a CTFModel object. """
    md = xmipp.MetaData(filename)
    ctfObj = rowToCtfModel(md, md.firstObject())  
    ctfObj._xmippMd = String(filename) 
    
    return ctfObj

def writeCTFModel(ctfObj, filename):
    """ Write a CTFModel object as Xmipp .ctfparam"""
    md = xmipp.MetaData()
    md.setColumnFormat(False)
    objId = md.addObject()
    ctfRow = XmippMdRow() 
    ctfModelToRow(ctfObj, ctfRow)
    ctfRow.writeToMd(md, objId)
    md.write(filename)
        
def writeSetOfCoordinates(posDir, coordSet):
    """ Write a pos file on metadata format for each micrograph 
    on the coordSet. 
    Params:
        posDir: the directory where the .pos files will be written.
        coordSet: the SetOfCoordinates that will be read.
    """
    posFiles = []
    boxSize = coordSet.getBoxSize()   
    
    # Write pos metadatas (one per micrograph)    
    for mic in coordSet.iterMicrographs():
        micName = mic.getFileName()
        posFn = join(posDir, replaceBaseExt(micName, "pos"))
        
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
        scipionPosFile = join(posDir, "scipion_" + replaceBaseExt(mic.getFileName(), 'pos'))
        posMd = readPosCoordinates(posFile)
        
        for objId in posMd:
            coord = rowToCoordinate(posMd, objId)
            coord.setMicrograph(mic)
            coordSet.append(coord)      
            # Add an unique ID that will be propagated to particles
            posMd.setValue(xmipp.MDL_ITEM_ID, long(coord.getObjId()), objId)
        if not posMd.isEmpty():
            posMd.write("particles@%s"  % scipionPosFile)
    coordSet._xmippMd = String(scipionPosFile)
    coordSet.setBoxSize(boxSize)

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
        
    imgSet._xmippMd = String(filename)
         
def writeSetOfImages(imgSet, filename, imgToFunc, ctfDir, rowFunc):
    """ This function will write a SetOfMicrographs as Xmipp metadata.
    Params:
        imgSet: the SetOfMicrograph instance.
        filename: the filename where to write the metadata.
        ctfDir: where to write ctfparam files if necessary
        rowFunc: this function can be used to setup the row before 
            adding to metadata.
    """
    particlesCtfFn = []
    md = xmipp.MetaData()
    hasCtf = imgSet.hasCTF()
    if ctfDir is None:
        ctfDir = dirname(filename)
    
    for img in imgSet:
        objId = md.addObject()
        imgRow = XmippMdRow()
        ctfFn = None
        if img.hasCTF():
            ctfModel = img.getCTF()
            rootFn = replaceBaseExt(img.getFileName(), "ctfparam")
            #rootFn = "%06d_%s" % (img.getObjId(), replaceBaseExt(img.getFileName(), "ctfparam"))
            ctfFn = join(ctfDir, rootFn)
            ctfModel.sphericalAberration = Float(imgSet.getAcquisition().sphericalAberration.get())
            if ctfModel.getObjId() not in particlesCtfFn:
                particlesCtfFn.append(ctfModel.getObjId())
                writeCTFModel(ctfModel, ctfFn)
        imgToFunc(img, imgRow, ctfFn)
        
        if rowFunc:
            rowFunc(img, imgRow)
        imgRow.writeToMd(md, objId)
        
    md.write(filename)
    imgSet._xmippMd = String(filename)
        
def writeImgToMetadata(md, img, hasCtf, ctfDir, imgToFunc, rowFunc ):
    objId = md.addObject()
    imgRow = XmippMdRow()
    imgToFunc(img, imgRow, ctfDir, hasCtf)
    if rowFunc:
        rowFunc(img, imgRow)
    imgRow.writeToMd(md, objId)        
    

def readSetOfParticles(filename, partSet, hasCtf=False):
    readSetOfImages(filename, partSet, rowToParticle, hasCtf)
        
def writeSetOfParticles(imgSet, filename, ctfDir=None, rowFunc=None):
    writeSetOfImages(imgSet, filename, particleToRow, ctfDir, rowFunc)

def writeSetOfCTFs(ctfSet, mdCTF):
    """ Write a ctfSet on metadata format. 
    Params:
        ctfSet: the SetOfCTF that will be read.
        mdCTF: The file where metadata should be written.
        ctfDir: where to write ctfparam files if necessary
    """
    md = xmipp.MetaData()
            
    for ctfModel in ctfSet:
        objId = md.addObject()
        ctfRow = XmippMdRow()
        ctfModelToRow(ctfModel, ctfRow)
        ctfRow.setValue(xmipp.MDL_MICROGRAPH, ctfModel.micFile.get())
        ctfRow.writeToMd(md, objId)
        
    md.write(mdCTF)
    ctfSet._xmippMd = String(mdCTF)
    
def writeSetOfClasses2D(classes2DSet, filename, ctfDir=None, classesBlock='classes'):    
    """ This function will write a SetOfClasses2D as Xmipp metadata.
    Params:
        classes2DSet: the SetOfClasses2D instance.
        filename: the filename where to write the metadata.
    """
    classFn = '%s@%s' % (classesBlock, filename)
    classMd = xmipp.MetaData()
    classMd.write(classFn) # Empty write to ensure the classes is the first block
    if ctfDir is None:
        ctfDir = dirname(filename)
    
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
    classes2DSet._xmippMd = String(filename)


def readSetOfClasses2D(classes2DSet, filename, classesBlock='classes', **args):
    """read from Xmipp image metadata.
        fnImages: The metadata filename where the particles properties are.
        imgSet: the SetOfParticles that will be populated.
        hasCtf: is True if the ctf information exists.
    """
    blocks = xmipp.getBlocksInMetaDataFile(filename)
    classesMd = xmipp.MetaData('%s@%s' % (classesBlock, filename))
    samplingRate = classes2DSet.getImages().getSamplingRate()
    averages = None
    
    if classesMd.containsLabel(xmipp.MDL_IMAGE):
        averages = classes2DSet.createAverages()
    
    for objId in classesMd:
        class2D = Class2D()
        class2D = rowToClass2D(classesMd, objId, class2D)
        classes2DSet.append(class2D)
        ref = class2D.getObjId()
        b = 'class%06d_images' % ref
        
        if b in blocks:
            classImagesMd = xmipp.MetaData('%s@%s' % (b, filename))
            
            for imgId in classImagesMd:
                img = rowToParticle(classImagesMd, imgId, hasCtf=False)
                img.setSamplingRate(samplingRate)
                class2D.append(img)
                
        class2D.write()
    
        if averages is not None:
            index, avgFn = xmippToLocation(classesMd.getValue(xmipp.MDL_IMAGE, objId))
            avg = Particle()
            avg.setLocation(index, avgFn)
            avg.copyObjId(class2D)
            averages.append(avg)
            
    classes2DSet._xmippMd = String(filename)
    print "end of readSetOfClasses2D"
         

def createXmippInputImages(self, imgSet, rowFunc=None, imagesFn=None):  
    imgsMd = getattr(imgSet, '_xmippMd', None)
    if imgsMd is None:
        imgsFn = imagesFn    
        ctfDir = None
        if self is not None:
            imgsFn = self._getPath(imagesFn or 'input_images.xmd')
            ctfDir = self._getExtraPath()

        writeSetOfParticles(imgSet, imgsFn, ctfDir, rowFunc)
    else:
        imgsFn = imgsMd.get()
    return imgsFn


def createXmippInputMicrographs(self, micSet, rowFunc=None, micsFn=None):    
    micsMd = getattr(micSet, '_xmippMd', None)
    if micsMd is None:
        ctfDir = None
        if self is not None:
            micsFn = self._getPath('input_micrographs.xmd')
            ctfDir = self._getExtraPath()
        writeSetOfMicrographs(micSet, micsFn, ctfDir, rowFunc)
    else:
        micsFn = micsMd.get()
    return micsFn

def createXmippInputVolumes(self, volSet, volsFn=None):    
    volsMd = getattr(volSet, '_xmippMd', None)
    if volsMd is None:
        if self is not None:
            volsFn = self._getPath('input_volumes.xmd')
        
        writeSetOfVolumes(volSet, volsFn)
    else:
        volsFn = volsMd.get()
    return volsFn

def createXmippInputClasses2D(self, classSet, classFn=None):
    classMd = getattr(classSet, '_xmippMd', None)
    if classMd is None:
        if self is not None:
            classFn = self._getPath('input_classes.xmd')
        
        writeSetOfClasses2D(classSet, classFn)
    else:
        classFn = classMd.get()
    return classFn
    
def createXmippInputCTF(self, ctfSet, ctfFn=None):
    ctfMd = getattr(ctfSet, '_xmippMd', None)
    if ctfMd is None:
        if self is not None:
            ctfFn = self._getPath('input_ctfs.xmd')
        
        writeSetOfCTFs(ctfSet, ctfFn)
    else:
        ctfFn = ctfMd.get()
    return ctfFn
