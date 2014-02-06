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
from protlib_xmipp import RowMetaData

# This dictionary will be used to map
# between CTFModel properties and Xmipp labels
ACQUISITION_DICT = { 
       "_amplitudeContrast": xmipp.MDL_CTF_Q0,
       "_sphericalAberration": xmipp.MDL_CTF_CS,
       "_voltage": xmipp.MDL_CTF_VOLTAGE,
       }

COOR_DICT = {"_x": xmipp.MDL_XCOOR, 
             "_y": xmipp.MDL_YCOOR 
             }

CTF_DICT = { 
       "_defocusAngle": xmipp.MDL_CTF_DEFOCUS_ANGLE,
       "_defocusU": xmipp.MDL_CTF_DEFOCUSU,
       "_defocusV": xmipp.MDL_CTF_DEFOCUSV,
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


def micrographToCTFParam(mic, ctfparam):
    """ This function is used to convert a Micrograph object
    to the .ctfparam metadata needed by some Xmipp protograms.
    If the micrograph already comes from xmipp, the ctfparam
    will be returned, if not, the new file. 
    """
#    if hasattr(mic, '_xmippMd'):
#        return mic._xmippMd.get()
    
    md = xmipp.MetaData()
    row = XmippMdRow()
    ctfModelToRow(mic.getCTF(), row)
    acquisitionToRow(mic.getAcquisition(), row)
    row.writeToMd(md, md.addObject())
    md.write(ctfparam)
    print 'MD------'
    print md
    print 'mic-----'
    mic.printAll()
    
    return ctfparam
    
    
def imageToRow(img, imgRow, imgLabel):
    setRowId(imgRow, img) # Set the id in the metadata as MDL_ITEM_ID
    index, filename = img.getLocation()
    fn = locationToXmipp(index, filename)
    imgRow.setValue(imgLabel, fn)
       
    if img.hasCTF():
        ctfModelToRow(img.getCTF(), imgRow)
    acquisitionToRow(img.getAcquisition(), imgRow)
        
        
def rowToImage(md, objId, imgLabel, imgClass, hasCtf):
    """ Create a Particle from a row of a metadata. """
    img = imgClass()
    # Decompose Xmipp filename
    index, filename = xmippToLocation(md.getValue(imgLabel, objId))
    img.setLocation(index, filename)
    if hasCtf:
        ctfModel = CTFModel()
        rowToCtfModel(md, objId)
        #TODO: CHECK NEXT LINE
        #ctfModel.setMicFile(md.getValue(xmipp.MDL_MICROGRAPH, objId))
        img.setCTF(ctfModel)
    setObjId(img, rowFromMd(md, objId))
    
    return img
    
    
def micrographToRow(mic, micRow):
    """ Set labels values from Micrograph mic to md row. """
    imageToRow(mic, micRow, imgLabel=xmipp.MDL_MICROGRAPH)
    
    
def rowToMicrograph(md, objId, hasCtf):
    """ Create a Micrograph object from a row of Xmipp metadata. """
    return rowToImage(md, objId, xmipp.MDL_MICROGRAPH, Micrograph, hasCtf)

def volumeToRow(vol, volRow):
    """ Set labels values from Micrograph mic to md row. """
    imageToRow(vol, volRow, imgLabel=xmipp.MDL_IMAGE)

def rowToVolume(md, objId, hasCtf):
    """ Create a Volume object from a row of Xmipp metadata. """
    return rowToImage(md, objId, xmipp.MDL_IMAGE, Volume, False)


def coordinateToRow(coord, coordRow):
    """ Set labels values from Coordinate coord to md row. """
    setRowId(coordRow, coord)
    objectToRow(coord, coordRow, COOR_DICT)


def rowToCoordinate(md, objId):
    """ Create a Coordinate from a row of a metadata. """
    coord = Coordinate()
    rowToObject(md, objId, coord, COOR_DICT)
        
    return coord


def rowToParticle(md, objId, hasCtf):
    """ Create a Particle from a row of a metadata. """
    return rowToImage(md, objId, xmipp.MDL_IMAGE, Particle, hasCtf)
    
    
def particleToRow(part, partRow):
    """ Set labels values from Particle to md row. """
    imageToRow(part, partRow, imgLabel=xmipp.MDL_IMAGE)


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
    
def rowToClass3D(md, objId, class3D):
    """ Create a Class2D from a row of a metadata. """
    setObjId(class3D, rowFromMd(md, objId), label=xmipp.MDL_REF)

    if md.containsLabel(xmipp.MDL_IMAGE):
        index, filename = xmippToLocation(md.getValue(xmipp.MDL_IMAGE, objId))
        img = Volume()
        img.copyObjId(class3D)
        img.setLocation(index, filename)
        class3D.setAverage(img)
    
    return class3D
    
def class2DToRow(class2D, classRow):
    """ Set labels values from Class2D to md row. """

    if class2D.hasAverage():
        index, filename = class2D.getAverage().getLocation()
        fn = locationToXmipp(index, filename)
        classRow.setValue(xmipp.MDL_IMAGE, fn)
    n = long(len(class2D))
    classRow.setValue(xmipp.MDL_CLASS_COUNT, n)
    setRowId(classRow, class2D, label=xmipp.MDL_REF)
        
        
def ctfModelToRow(ctfModel, ctfRow):
    """ Set labels values from ctfModel to md row. """
    objectToRow(ctfModel, ctfRow, CTF_DICT)


def rowToCtfModel(md, objId):
    """ Create a CTFModel from a row of a metadata. """
    ctfModel = CTFModel()
    rowToObject(md, objId, ctfModel, CTF_DICT) 
    return ctfModel
    

def acquisitionToRow(acquisition, ctfRow):
    """ Set labels values from acquisition to md row. """
    objectToRow(acquisition, ctfRow, ACQUISITION_DICT)

def rowToAcquisition(md, objId):
    """ Create an acquisition from a row of a metadata. """
    acquisition = Acquisition()
    rowToObject(md, objId, acquisition, ACQUISITION_DICT) 
    return acquisition
    
    
def readSetOfMicrographs(filename, micSet, hasCtf=False):    
    readSetOfImages(filename, micSet, rowToMicrograph, hasCtf)


def writeSetOfMicrographs(micSet, filename, rowFunc=None):
    writeSetOfImages(micSet, filename, micrographToRow, rowFunc)
    
    
def readSetOfVolumes(filename, volSet, hasCtf=False):    
    readSetOfImages(filename, volSet, rowToVolume, False)


def writeSetOfVolumes(volSet, filename, rowFunc=None):
    writeSetOfImages(volSet, filename, volumeToRow, rowFunc)    
    
    
def readCTFModel(filename):
    """ Read from Xmipp .ctfparam and create a CTFModel object. """
    md = xmipp.MetaData(filename)
    ctfObj = rowToCtfModel(md, md.firstObject())  
    ctfObj._xmippMd = String(filename) 
    
    return ctfObj

        
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
    # Read the sampling rate from the acquisition info file if exists
    
    acqFile = join(dirname(filename), 'acquisition_info.xmd')
    print acqFile
    if exists(acqFile):
        md = RowMetaData(acqFile)
        samplingRate = md.getValue(xmipp.MDL_SAMPLING_RATE)
        imgSet.setSamplingRate(samplingRate)
    else:
        pass # TODO: what to do if not exists
        
    imgSet._xmippMd = String(filename)
         
 
def setOfImagesToMd(imgSet, md, imgToFunc, rowFunc):
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
        imgToFunc(img, imgRow)#drop ctfFn
        
        if rowFunc:
            rowFunc(img, imgRow)
        imgRow.writeToMd(md, objId)
            
def writeSetOfImages(imgSet, filename, imgToFunc, rowFunc):
    """ This function will write a SetOfMicrographs as Xmipp metadata.
    Params:
        imgSet: the set of images to be written (particles, micrographs or volumes)
        filename: the filename where to write the metadata.
        rowFunc: this function can be used to setup the row before 
            adding to metadata.
    """
    md = xmipp.MetaData()
    setOfImagesToMd(imgSet, md, imgToFunc, rowFunc)
    md.write(filename)
    imgSet._xmippMd = String(filename)
        
        
def writeImgToMetadata(md, img, hasCtf, imgToFunc, rowFunc ):
    objId = md.addObject()
    imgRow = XmippMdRow()
    imgToFunc(img, imgRow, hasCtf)
    if rowFunc:
        rowFunc(img, imgRow)
    imgRow.writeToMd(md, objId)        
    

def readSetOfParticles(filename, partSet, hasCtf=False):
    readSetOfImages(filename, partSet, rowToParticle, hasCtf)
     
def setOfParticlesToMd(imgSet, md, rowFunc=None):
    setOfImagesToMd(imgSet, md, particleToRow, rowFunc)
      
def setOfMicrographsToMd(imgSet, md, rowFunc=None):
    setOfImagesToMd(imgSet, md, micrographToRow, rowFunc)
      
def writeSetOfParticles(imgSet, filename, rowFunc=None):
    writeSetOfImages(imgSet, filename, particleToRow, rowFunc)

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
        ctfModelToRow(ctfModel, ctfRow)
        ctfRow.setValue(xmipp.MDL_MICROGRAPH, ctfModel.getMicFile())
        ctfRow.writeToMd(md, objId)
        
    md.write(mdCTF)
    ctfSet._xmippMd = String(mdCTF)
    
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

def writeSetOfClasses3D(classes3DSet, filename, classesBlock='classes'):    
    """ This function will write a SetOfClasses3D as Xmipp metadata.
    Params:
        classes3DSet: the SetOfClasses3D instance.
        filename: the filename where to write the metadata.
    """
    classFn = '%s@%s' % (classesBlock, filename)
    classMd = xmipp.MetaData()
    classMd.write(classFn) # Empty write to ensure the classes is the first block
    
    classRow = XmippMdRow()
    for class3D in classes3DSet:        
        class3DToRow(class3D, classRow)
        classRow.writeToMd(classMd, classMd.addObject())
        ref = class3D.getObjId()
        imagesFn = 'class%06d_images@%s' % (ref, filename)
        imagesMd = xmipp.MetaData()
        imgRow = XmippMdRow()
        
        for vol in class3D:
            volumeToRow(vol, imgRow)
            imgRow.writeToMd(imagesMd, imagesMd.addObject())
        imagesMd.write(imagesFn, xmipp.MD_APPEND)
    
    classMd.write(classFn, xmipp.MD_APPEND) # Empty write to ensure the classes is the first block
    classes3DSet._xmippMd = String(filename)

def readSetOfClasses3D(classes3DSet, filename, classesBlock='classes', **args):
    """read from Xmipp image metadata.
        fnImages: The metadata filename where the particles properties are.
        imgSet: the SetOfParticles that will be populated.
        hasCtf: is True if the ctf information exists.
    """
    blocks = xmipp.getBlocksInMetaDataFile(filename)
    classesMd = xmipp.MetaData('%s@%s' % (classesBlock, filename))
    samplingRate = classes3DSet.getVolumes().getSamplingRate()
    averages = None
    
    if classesMd.containsLabel(xmipp.MDL_IMAGE):
        averages = classes3DSet.createAverages()
    
    for objId in classesMd:
        class3D = Class3D()
        class3D = rowToClass3D(classesMd, objId, class3D)
        classes3DSet.append(class3D)
        ref = class3D.getObjId()
        b = 'class%06d_images' % ref
        
        if b in blocks:
            classImagesMd = xmipp.MetaData('%s@%s' % (b, filename))
            
            for imgId in classImagesMd:
                img = rowToParticle(classImagesMd, imgId, hasCtf=False)
                img.setSamplingRate(samplingRate)
                class3D.append(img)
                
        class3D.write()
    
        if averages is not None:
            index, avgFn = xmippToLocation(classesMd.getValue(xmipp.MDL_IMAGE, objId))
            avg = Volume()
            avg.setLocation(index, avgFn)
            avg.copyObjId(class3D)
            averages.append(avg)
    classes3DSet._xmippMd = String(filename)

def createXmippInputImages(self, imgSet, rowFunc=None, imagesFn=None):  
    imgsMd = getattr(imgSet, '_xmippMd', None)
    if imgsMd is None:
        imgsFn = imagesFn    
        if self is not None:
            imgsFn = self._getPath(imagesFn or 'input_images.xmd')

        writeSetOfParticles(imgSet, imgsFn, rowFunc)
    else:
        imgsFn = imgsMd.get()
    return imgsFn


def createXmippInputMicrographs(self, micSet, rowFunc=None, micsFn=None):    
    micsMd = getattr(micSet, '_xmippMd', None)
    if micsMd is None:
        if self is not None:
            micsFn = self._getPath('input_micrographs.xmd')
        writeSetOfMicrographs(micSet, micsFn, rowFunc)
    else:
        micsFn = micsMd.get()
    return micsFn

def createXmippInputVolumes(prot, volSet, volsFn=None):    
    volsMd = getattr(volSet, '_xmippMd', None)
    if volsMd is None:
        if volsFn is None:
            volsFn = prot._getPath('input_volumes.xmd')
        
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
