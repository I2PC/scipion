# **************************************************************************
# *
# * Authors:     Mohsen Kazemi  (mkazemi@cnb.csic.es)
# *              Joaquin Oton   (joton@cnb.csic.es)
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
This modules contains basic hierarchy
for ET data objects like: Tilt Series, Focal Series and others
"""

import os
import json

from pyworkflow.object import *
#import pyworkflow.utils as pwutils

from constants import *
from pyworkflow.em.convert import *
import numpy as np
# import xmipp





class ETObject(OrderedObject):
    """Base object for all ET classes"""
    def __init__(self, **args):
        OrderedObject.__init__(self, **args)
        
    def __str__(self):
        return self.getClassName()
    
    def getFiles(self):
        """ Get all filePaths """
        return None




    
class TiltSeries(ETObject):
    """Represents an ET Tilt-series object"""
    def __init__(self, **args):
        ETObject.__init__(self, **args)
        # Image location is composed by an index and a filename
        self._index = Integer(0)
        self._filename = String()
        self._samplingRate = Float()
        self._ctfModel = None
        self._acquisition = None
        # _transform property will store the transformation matrix
        # this matrix can be used for 2D/3D alignment or
        # to represent projection directions
        self._transform = None

    def getSamplingRate(self):
        """ Return image sampling rate. (A/pix) """
        return self._samplingRate.get()
    
    def setSamplingRate(self, sampling):
        self._samplingRate.set(sampling)

    def getFormat(self):
        pass
    
    def getDataType(self):
        pass

    def getDimensions(self):
        """getDimensions is redundant here but not in setOfVolumes
         create it makes easier to crete protocols for both images
         and sets of images
        """
        self.getDim()

    def getDim(self):
        """Return image dimensions as tuple: (Xdim, Ydim, Zdim)"""
        fn = self.getFileName()
        if fn is not None and os.path.exists(fn.replace(':mrc', '')):
            x, y, z, n = ImageHandler().getDimensions(self)
            return x, y, z
        return None

    def getXDim(self):
        return self.getDim()[0] if self.getDim() is not None else 0
    
    def getIndex(self):
        return self._index.get()
    
    def setIndex(self, index):
        self._index.set(index)

    def getFileName(self):
        """ Use the _objValue attribute to store filename. """
        return self._filename.get()
    
    def setFileName(self, filename):
        """ Use the _objValue attribute to store filename. """
        self._filename.set(filename)
        
    def getLocation(self):
        """ This function return the image index and filename.
        It will only differs from getFileName, when the image
        is contained in a stack and the index make sense. 
        """
        return (self.getIndex(), self.getFileName())
    
    def setLocation(self, *args):
        """ Set the image location, see getLocation. 
        Params:
            First argument can be:
             1. a tuple with (index, filename)
             2. a index, this implies a second argument with filename
             3. a filename, this implies index=NO_INDEX
        """
        first = args[0]
        t = type(first)
        if t == tuple:
            index, filename = first
        elif t == int or t == long:
            index , filename = first, args[1]
        elif t == str or t == unicode:
            index, filename = NO_INDEX, first
        else:
            raise Exception('setLocation: unsupported type %s as input.' % t)
            
        self.setIndex(index)
        self.setFileName(filename)
        
    def copyInfo(self, other):
        """ Copy basic information """
        self.copyAttributes(other, '_samplingRate')
        
    def copyLocation(self, other):
        """ Copy location index and filename from other image. """
        self.setIndex(other.getIndex())
        self.setFileName(other.getFileName())
        
    #def hasCTF(self):
    #    return self._ctfModel is not None
    
    #def getCTF(self):
    #    """ Return the CTF model """
    #    return self._ctfModel
    
    #def setCTF(self, newCTF):
    #    self._ctfModel = newCTF
        
    def hasAcquisition(self):
        return (self._acquisition is not None and
                self._acquisition.getVoltage() is not None and
                self._acquisition.getMagnification() is not None
                )
        
    def getAcquisition(self):
        return self._acquisition

    def setAcquisition(self, acquisition):
        self._acquisition = acquisition
        
    def hasTransform(self):
        return self._transform is not None
    
    def getTransform(self):
        return self._transform
    
    def setTransform(self, newTransform):
        self._transform = newTransform

    def __str__(self):
        """ String representation of an Image. """
        dim = self.getDim()
        if dim:
            dimStr = str(ImageDim(*dim))
        else:
            dimStr = 'No-Dim'
        return "%s (%s, %0.2f A/px)" % (self.getClassName(), dimStr, self.getSamplingRate() or 99999.)
    
    def getFiles(self):
        filePaths = set()
        filePaths.add(self.getFileName())
        return filePaths


class ETSet(Set, ETObject):
    
    def _loadClassesDict(self):
        import pyworkflow.em as em
        classDict = em.getObjects()
        classDict.update(globals())
        
        return classDict
    
    def copyInfo(self, other):
        """ Define a dummy copyInfo function to be used
        for some generic operations on sets.
        """
        pass
    
    def clone(self):
        """ Override the clone defined in Object
        to avoid copying _mapperPath property
        """
        pass
    
    def copyItems(self, otherSet, 
                  updateItemCallback=None, 
                  itemDataIterator=None,
                  copyDisabled=False):
        """ Copy items from another set.
        If the updateItemCallback is passed, it will be 
        called with each item (and optionally with a data row).
        This is a place where items can be updated while copying.
        This is useful to set new attributes or update values
        for each item.
        """
        for item in otherSet:
            # copy items if enabled or copyDisabled=True
            if copyDisabled or item.isEnabled():
                newItem = item.clone()
                if updateItemCallback:
                    row = None if itemDataIterator is None else next(itemDataIterator)
                    updateItemCallback(newItem, row)
                # If updateCallBack function returns attribute
                # _appendItem to False do not append the item
                if getattr(newItem, "_appendItem", True):
                    self.append(newItem)
            else:
                if itemDataIterator is not None:
                    next(itemDataIterator) # just skip disabled data row
                    
    def getFiles(self):
        return Set.getFiles(self)


  
class FocalSeries(ETSet):
    """ Represents a set of Images """
    ITEM_TYPE = TiltSeries
    
    def __init__(self, **kwargs):
        ETSet.__init__(self, **kwargs)
        self._samplingRate = Float()
        self._hasCtf = Boolean(kwargs.get('ctf', False))
        self._alignment = String(ALIGN_NONE)
        self._isPhaseFlipped = Boolean(False)
        self._isAmplitudeCorrected = Boolean(False)
        self._acquisition = Acquisition()
        self._firstDim = ImageDim() # Dimensions of the first image

    def getAcquisition(self):
        return self._acquisition
        
    def setAcquisition(self, acquisition):
        self._acquisition = acquisition
        
    def hasAcquisition(self):
        return self._acquisition.getMagnification() is not None
        
    #def hasCTF(self):
    #    """Return True if the SetOfImages has associated a CTF model"""
    #    return self._hasCtf.get()  
    
    #def setHasCTF(self, value):
    #    self._hasCtf.set(value)
        
    def hasAlignment(self):
        return self._alignment != ALIGN_NONE
    
    def hasAlignment2D(self):
        return self._alignment == ALIGN_2D

    def hasAlignment3D(self):
        return self._alignment == ALIGN_3D

    def hasAlignmentProj(self):
        return self._alignment == ALIGN_PROJ

    def getAlignment(self):
        return self._alignment.get()

    def setAlignment(self, value):
        if not value in ALIGNMENTS:
            raise Exception('Invalid alignment value: "%s"' % value)
        self._alignment.set(value)
        
    def setAlignment2D(self):
        self.setAlignment(ALIGN_2D)
        
    def setAlignment3D(self):
        self.setAlignment(ALIGN_3D)

    def setAlignmentProj(self):
        self.setAlignment(ALIGN_PROJ)
    
    #def isPhaseFlipped(self):
    #    return self._isPhaseFlipped.get()
    
    #def setIsPhaseFlipped(self, value):
    #    self._isPhaseFlipped.set(value)
        
    #def isAmplitudeCorrected(self):
    #    return self._isAmplitudeCorrected.get()
    
    #def setIsAmplitudeCorrected(self, value):
    #    self._isAmplitudeCorrected.set(value)

    def append(self, image):
        """ Add a image to the set. """
        # If the sampling rate was set before, the same value
        # will be set for each image added to the set
        if self.getSamplingRate() or not image.getSamplingRate():
            image.setSamplingRate(self.getSamplingRate())
        # Copy the acquistion from the set to images
        if self.hasAcquisition(): # only override image acquisition if setofImages acquisition is not none
            #TODO: image acquisition should not be overwritten
            if not image.hasAcquisition():
                image.setAcquisition(self.getAcquisition())
        # Store the dimensions of the first image, just to 
        # avoid reading image files for further queries to dimensions
        if self.getSize() == 0: # only check this for first time append is called
            if self._firstDim.isEmpty():
                self._firstDim.set(image.getDim())
        ETSet.append(self, image)

    def copyInfo(self, other):
        """ Copy basic information (sampling rate and ctf)
        from other set of images to current one"""
        self.copyAttributes(other, '_samplingRate', '_isPhaseFlipped', '_isAmplitudeCorrected', '_alignment')
        self._acquisition.copyInfo(other._acquisition)
        
    def getFiles(self):
        filePaths = set()
        uniqueFiles = self.aggregate(['count'],'_filename',['_filename'])

        for row in uniqueFiles:
            filePaths.add(row['_filename'])
        return filePaths
    
    def setDownsample(self, downFactor):
        """ Update the values of samplingRate and scannedPixelSize
        after applying a downsampling factor of downFactor.
        """
        self.setSamplingRate(self.getSamplingRate() * downFactor)        
        
    def setSamplingRate(self, samplingRate):
        """ Set the sampling rate and adjust the scannedPixelSize. """
        self._samplingRate.set(samplingRate)
        
    def getSamplingRate(self):
        return self._samplingRate.get()
    
    def writeStack(self, fnStack, orderBy='id', direction='ASC',
                   applyTransform=False):
        # TODO create empty file to improve efficiency
        ih = ImageHandler()
        applyTransform = applyTransform and self.hasAlignment2D()

        for i, img in enumerate(self.iterItems(orderBy=orderBy,
                                               direction=direction)):
            transform = img.getTransform() if applyTransform else None
            ih.convert(img, (i+1, fnStack), transform=transform)
    
    # TODO: Check whether this function can be used.
    # for example: protocol_apply_mask
    def readStack(self, fnStack, postprocessImage=None):
        """ Populate the set with the images in the stack """
        _,_,_, ndim = ImageHandler().getDimensions(fnStack)
        img = self.ITEM_TYPE()
        for i in range(1, ndim+1):
            img.setObjId(None)
            img.setLocation(i, fnStack)
            if postprocessImage is not None:
                postprocessImage(img)
            self.append(img)
    
    def getDim(self):
        """ Return the dimensions of the first image in the set. """
        if self._firstDim.isEmpty():
            return None
        x, y, z = self._firstDim
        return x, y, z
    
    def getXDim(self):
        return self.getDim()[0] if self.getDim() is not None else 0
    
    def isOddX(self):
        """ Return True if the first item x dimension is odd. """
        return self.getXDim() % 2 == 1
    
    def getDimensions(self):
        """Return first image dimensions as a tuple: (xdim, ydim, zdim)"""
        return self.getFirstItem().getDim()
    
    def __str__(self):
        """ String representation of a set of images. """
        sampling = self.getSamplingRate()
        
        if not sampling:
            print "FATAL ERROR: Object %s has no sampling rate!!!" % self.getName()
            sampling = -999.0

        s = "%s (%d items, %s, %0.2f A/px)" % (self.getClassName(),
                                               self.getSize(),
                                               self._firstDim, sampling)
        return s

    def iterItems(self, orderBy='id', direction='ASC'):
        """ Redefine iteration to set the acquisition to images. """
        for img in Set.iterItems(self, orderBy=orderBy, direction=direction):
            # Sometimes the images items in the set could
            # have the acquisition info per data row and we
            # don't want to override with the set acquisition for this case
            if not img.hasAcquisition():
                img.setAcquisition(self.getAcquisition())
            yield img

    def appendFromImages(self, imagesSet):
        """ Iterate over the images and append 
        every image that is enabled. 
        """
        for img in imagesSet:
            if img.isEnabled():
                self.append(img)
                                    
    def appendFromClasses(self, classesSet):
        """ Iterate over the classes and the element inside each
        class and append to the set all that are enabled. 
        """
        for cls in classesSet:
            if cls.isEnabled() and cls.getSize() > 0:
                for img in cls:
                    if img.isEnabled():               
                        self.append(img)




class ETFile(ETObject):
    """ Class to link usually to text files. """
    def __init__(self, filename=None, **args):
        ETObject.__init__(self, **args)
        self._filename = String(filename)
        
    def getFileName(self):
        """ Use the _objValue attribute to store filename. """
        return self._filename.get()
    
    def setFileName(self, filename):
        """ Use the _objValue attribute to store filename. """
        self._filename.set(filename)    

    





