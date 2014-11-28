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
This modules contains basic hierarchy
for EM data objects like: Image, SetOfImage and others
"""

import json

from pyworkflow.mapper.sqlite import SqliteMapper, SqliteFlatMapper
from pyworkflow.object import *
from pyworkflow.utils.path import cleanPath, dirname, join, replaceExt, exists

from constants import *
from convert import ImageHandler
import numpy as np
import xmipp


class EMObject(OrderedObject):
    """Base object for all EM classes"""
    def __init__(self, **args):
        OrderedObject.__init__(self, **args)
        
    def __str__(self):
        return self.getClassName()
    
    def getFiles(self):
        """ Get all filePaths """
        return None


class Acquisition(EMObject):
    """Acquisition information"""
    def __init__(self, **args):
        EMObject.__init__(self, **args)
        self._magnification = Float(args.get('magnification', None)) 
        # Microscope voltage in kV
        self._voltage = Float(args.get('voltage', None))
        # Spherical aberration in mm
        self._sphericalAberration = Float(args.get('sphericalAberration', None)) 
        self._amplitudeContrast = Float(args.get('amplitudeContrast', None))
        
    def copyInfo(self, other):
        self.copyAttributes(other, '_magnification', '_voltage', 
                            '_sphericalAberration', '_amplitudeContrast')
        
    def getMagnification(self):
        return self._magnification.get()
        
    def setMagnification(self, value):
        self._magnification.set(value)
        
    def getVoltage(self):
        return self._voltage.get()
        
    def setVoltage(self, value):
        self._voltage.set(value)
        
    def getSphericalAberration(self):
        return self._sphericalAberration.get()
    
    def setSphericalAberration(self, value):
        self._sphericalAberration.set(value)
        
    def getAmplitudeContrast(self):
        return self._amplitudeContrast.get()
    
    def setAmplitudeContrast(self, value):
        self._amplitudeContrast.set(value)        

    def __str__(self):
        return "\n    mag=%f\n    volt= %f\n    Cs=%f\n    Q0=%f\n\n"%(self._magnification.get(),
                                                                     self._voltage.get(),
                                                                     self._sphericalAberration.get(),
                                                                     self._amplitudeContrast.get())


    
class CTFModel(EMObject):
    """ Represents a generic CTF model. """
    def __init__(self, **args):
        EMObject.__init__(self, **args)
        self._defocusU = Float(args.get('defocusU', None))
        self._defocusV = Float(args.get('defocusV', None))
        self._defocusAngle = Float(args.get('defocusAngle', None))
        self._defocusRatio = Float()
        self._psdFile = String()
#         self._micFile = String()
        self._micObj  = None

    def __str__(self):
        ctfStr = "defocus(U,V,a) = (%0.2f,%0.2f,%0.2f)" % (self._defocusU.get(),
                                                           self._defocusV.get(), 
                                                           self._defocusAngle.get())
        if self._micObj:
            ctfStr + " mic=%s" % self._micObj
        return ctfStr

    def getDefocusU(self):
        return self._defocusU.get()
        
    def setDefocusU(self, value):
        self._defocusU.set(value)
        
    def getDefocusV(self):
        return self._defocusV.get()
        
    def setDefocusV(self, value):
        self._defocusV.set(value)
        
    def getDefocusAngle(self):
        return self._defocusAngle.get()
        
    def setDefocusAngle(self, value):
        self._defocusAngle.set(value)
        
    def getDefocusRatio(self):
        return self._defocusRatio.get()

    def copyInfo(self, other):
        self.copyAttributes(other, '_defocusU', '_defocusV','_defocusAngle',
                            '_defocusRatio', '_psdFile', '_micFile')
        
    def getPsdFile(self):
        return self._psdFile.get()
    
    def setPsdFile(self, value):
        self._psdFile.set(value)
        
    def getMicrograph(self):
        self._micObj.copyObjId(self)
        return self._micObj
    
    def setMicrograph(self, mic):
        self._micObj = mic

    def getDefocus(self):
        """ Returns defocusU, defocusV and defocusAngle. """
        return (self._defocusU.get(), 
                self._defocusV.get(), 
                self._defocusAngle.get())

    def setStandardDefocus(self, defocusU, defocusV, defocusAngle):
        """ Set defocus values following emx conventions. 
        See _standardize function."""
        self._defocusU.set(defocusU)
        self._defocusV.set(defocusV)
        self._defocusAngle.set(defocusAngle)
        self.standardize()

    def standardize(self):
        """ Modify defocusU, defocusV and defocusAngle to conform 
        the EMX standard: defocusU > defocusV, 0 <= defocusAngle < 180
        and the defocusAnges is between x-axis and defocusU. Also
        determine the defocusRatio(defocusU/defocusV).
        For more details see:
        http://i2pc.cnb.csic.es/emx/LoadDictionaryFormat.htm?type=Convention#ctf
        """
        if self._defocusV > self._defocusU:
            self._defocusV.swap(self._defocusU) # exchange defocuU by defocuV
            self._defocusAngle.sum(90.)
        if self._defocusAngle >= 180.:
            self._defocusAngle.sum(-180.)
        elif self._defocusAngle < 0.:
            self._defocusAngle.sum(180.)
        # At this point defocusU is always greater than defocusV
        # following the EMX standard
        self._defocusRatio.set(self.getDefocusU()/self.getDefocusV())
        
    def equalAttributes(self, other, ignore=[], verbose=False):
        """ Override default behaviour to compare two
        CTF objects, now ignoring the psdFile.
        """
        return (self._defocusU == other._defocusU and
                self._defocusV == other._defocusV and
                self._defocusAngle == other._defocusAngle
                )


class DefocusGroup(EMObject):
    """ Groups CTFs by defocus"""
    def __init__(self, **args):
        EMObject.__init__(self, **args)
        self._defocusMin = Float()
        self._defocusMax = Float()
        self._defocusAvg = Float()
        self._size = Integer()
        
    def getDefocusMin(self):
        return self._defocusMin.get()
        
    def setDefocusMin(self, value):
        self._defocusMin.set(value)

    def getDefocusMax(self):
        return self._defocusMax.get()
        
    def setDefocusMax(self, value):
        self._defocusMax.set(value)
        
    def getDefocusAvg(self):
        return self._defocusAvg.get()
        
    def setDefocusAvg(self, value):
        self._defocusAvg.set(value)
        
    def setSize(self, value):
        self._size.set(value)
        
    def getSize(self):
        return self._size.get()


class ImageDim(CsvList):
    """ Just a wrapper to a CsvList to store image dimensions
    as X, Y and Z. 
    """
    def __init__(self, x=None, y=None, z=None):
        CsvList.__init__(self, pType=int)
        if x is not None and y is not None:
            self.append(x)
            self.append(y)
            if z is not None:
                self.append(z)
        
    def getX(self):
        return self[0]
    
    def getY(self):
        return self[1]
    
    def getZ(self):
        return self[2]
    
    def __str__(self):
        if self.isEmpty():
            s = 'No-Dim'
        else:
            s = '%d x %d' % (self.getX(), self.getY())
            if self.getZ() > 1:
                s += ' x %d' % self.getZ()
        return s

    
class Image(EMObject):
    """Represents an EM Image object"""
    def __init__(self, **args):
        EMObject.__init__(self, **args)
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
    
    def getDim(self):
        """Return image dimensions as tuple: (Xdim, Ydim, Zdim)"""
        if exists(self.getFileName().replace(':mrc', '')):
            x, y, z, n = ImageHandler().getDimensions(self)
            return x, y, z
        return None
    
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
        
    def hasCTF(self):
        return self._ctfModel is not None
    
    def getCTF(self):
        """ Return the CTF model """
        return self._ctfModel
    
    def setCTF(self, newCTF):
        self._ctfModel = newCTF
        
    def hasAcquisition(self):
        # This doesn't work
        return self._acquisition is not None
        #FIXME: check this later, very very IMPORTANT!!!
        #return (self._acquisition is not None and 
        #        self._acquisition.getMagnification() is not None)
        
    def getAcquisition(self):
        return self._acquisition
    
    def setAcquisition(self, acquisition):
        self._acquisition = acquisition
        
    def hasTransfrom(self):
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
        return "%s (%s, %0.2f A/px)" % (self.getClassName(), dimStr, self.getSamplingRate())
    
    def getFiles(self):
        filePaths = set()
        filePaths.add(self.getFileName())
        return filePaths

class Micrograph(Image):
    """ Represents an EM Micrograph object """
    pass


class Particle(Image):
    """ Represents an EM Particle object """
    def __init__(self, **args):
        Image.__init__(self, **args)
        # This may be redundant, but make the Particle
        # object more indenpent for tracking coordinates
        self._coordinate = None
        self._micId = Integer()
        
    def hasCoordinate(self):
        return self._coordinate is not None
    
    def setCoordinate(self, coordinate):
        self._coordinate = coordinate
        
    def getCoordinate(self):
        return self._coordinate
    
    def getMicId(self):
        """ Return the micrograph id if the coordinate is not None.
        or have set the _micId property.
        """
        if self._micId.hasValue():
            return self._micId.get()
        if self.hasCoordinate():
            return self.getCoordinate().getMicId()
        
        return None
    
    def setMicId(self, micId):
        self._micId.set(micId)
        
    def hasMicId(self):
        return self.getMicId() is not None


class Mask(Particle):
    """ Represent a mask. """
    pass


class Volume(Image):
    """ Represents an EM Volume object """
    def __init__(self, **args):
        Image.__init__(self, **args)

        
class VolumeMask(Volume):
    """ A 3D mask to be used with volumes. """
    pass


class EMFile(EMObject):
    """ Class to link usually to text files. """
    def __init__(self, filename=None, **args):
        EMObject.__init__(self, **args)
        self._filename = String(filename)
        
    def getFileName(self):
        """ Use the _objValue attribute to store filename. """
        return self._filename.get()
    
    def setFileName(self, filename):
        """ Use the _objValue attribute to store filename. """
        self._filename.set(filename)    

    
class PdbFile(EMFile):
    """Represents an PDB file. """
    def __init__(self, filename=None, pseudoatoms=False, **args):
        EMFile.__init__(self, filename, **args)
        self._pseudoatoms = Boolean(pseudoatoms)
        
    def getPseudoAtoms(self):
        return self._pseudoatoms.get()
    
    def setPseudoAtoms(self, value):
        self._pseudoatoms.set(value)
        
    def __str__(self):
        return "%s (pseudoatoms=%s)" % (self.getClassName(), self.getPseudoAtoms())
    
    
class EMXObject(EMObject):
    """Represents EMX data object, mainly comprising two files:
    1- XML file specifiying metadata information
    2- A binary data file of either Micrographs or Particles.
    """
    def __init__(self, xmlFile=None, binaryFile=None, **args):
        EMObject.__init__(self, **args)
        self._xmlFile = String(xmlFile)
        self._binaryFile = String(binaryFile)
        
    def getXmlFile(self):
        return self._xmlFile.get()
    
    def getBinaryFile(self):
        return self._binaryFile.get()        
                
#       
# class EMSet(EMObject, Set):
#     def __init__(self, *args, **kwargs):
#         Set.__init__(self, *args, **kwargs)
#         EMObject.__init__(self, *args, **kwargs)

class EMSet(Set, EMObject):
    def _loadClassesDict(self):
        return globals()
    
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
                self.append(newItem)
            else:
                if itemDataIterator is not None:
                    next(itemDataIterator) # just skip disabled data row
  
  
class SetOfImages(EMSet):
    """ Represents a set of Images """
    ITEM_TYPE = Image

    
    def __init__(self, **args):
        EMSet.__init__(self, **args)
        self._samplingRate = Float()
        self._hasCtf = Boolean(args.get('ctf', False))
        #self._hasAlignment = Boolean(args.get('alignmet', False))
        #self._hasProjection = Boolean(False)
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
        
    def hasCTF(self):
        """Return True if the SetOfImages has associated a CTF model"""
        return self._hasCtf.get()  
    
    def setHasCTF(self, value):
        self._hasCtf.set(value)
        
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
    
    def isPhaseFlipped(self):
        return self._isPhaseFlipped.get()
    
    def setIsPhaseFlipped(self, value):
        self._isPhaseFlipped.set(value)
        
    def isAmplitudeCorrected(self):
        return self._isAmplitudeCorrected.get()
    
    def setIsAmplitudeCorrected(self, value):
        self._isAmplitudeCorrected.set(value)
        
    def append(self, image):
        """ Add a image to the set. """
        # If the sampling rate was set before, the same value
        # will be set for each image added to the set
        if self.getSamplingRate() or not image.getSamplingRate():
            image.setSamplingRate(self.getSamplingRate())
        # Copy the acquistion from the set to images
        if self.hasAcquisition(): # only override acquisition if not None
            image.setAcquisition(self.getAcquisition())
        # Store the dimensions of the first image, just to 
        # avoid reading image files for further queries to dimensions
        if self._firstDim.isEmpty():
            self._firstDim.set(image.getDim())
        EMSet.append(self, image)
    
    def copyInfo(self, other):
        """ Copy basic information (sampling rate, scannedPixelSize and ctf)
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
    
    def writeStack(self, fnStack):
        # TODO creaty empty file to improve efficiency
        ih = ImageHandler()
        for i, img in enumerate(self):
            ih.convert(img, (i+1, fnStack))
    
    # TODO: Check whether this function can be used.
    # for example: protocol_apply_mask
    def readStack(self, fnStack):
        """ Populate the set with the images in the stack """
        _,_,_, ndim = ImageHandler().getDimensions(fnStack)
        img = self.ITEM_TYPE()
        for i in range(1, ndim+1):
            img.setLocation(i, fnStack)
            self.append(img)
            img.cleanObjId()
    
    def getDim(self):
        """ Return the dimensions of the first image in the set. """
        if self._firstDim.isEmpty():
            return None
        x, y, z = self._firstDim
        return x, y, z
    
    def isOddX(self):
        """ Return True if the first item x dimension is odd. """
        x, _, _ = self.getDim()
        return x % 2 == 1
    
    def getDimensions(self):
        """Return first image dimensions as a tuple: (xdim, ydim, zdim)"""
        return self.getFirstItem().getDim()
    
    def __str__(self):
        """ String representation of a set of images. """
        sampling = self.getSamplingRate()
        
        if not sampling:
            print "FATAL ERROR: Object %s has no sampling rate!!!" % self.getName()
            sampling = -999.0
        if self._firstDim.isEmpty():
            try:
                self._firstDim.set(self.getFirstItem().getDim())
            except Exception, ex:
                print "Error reading dimension: ", ex
                import traceback
                traceback.print_exc()
        dimStr = str(self._firstDim)
        s = "%s (%d items, %s, %0.2f A/px)" % (self.getClassName(), self.getSize(), dimStr, sampling)
        return s

    def __iter__(self, random=False):
        """ Redefine iteration to set the acquisition to images. """
        for img in self._iterItems(random=random):
            # Sometimes the images items in the set could
            # have the acquisition info per data row and we
            # dont want to override with the set acquistion for this case
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
            if cls.isEnabled():
                for img in cls:
                    if img.isEnabled():                
                        self.append(img)


class SetOfMicrographs(SetOfImages):
    """Represents a set of Micrographs"""
    ITEM_TYPE = Micrograph
    
    def __init__(self, **args):
        SetOfImages.__init__(self, **args)
        self._scannedPixelSize = Float()
    
    def copyInfo(self, other):
        """ Copy basic information (voltage, spherical aberration and sampling rate)
        from other set of micrographs to current one.
        """
        SetOfImages.copyInfo(self, other)
        self._scannedPixelSize.set(other.getScannedPixelSize())
    
    def setSamplingRate(self, samplingRate):
        """ Set the sampling rate and adjust the scannedPixelSize. """
        self._samplingRate.set(samplingRate)
        mag = self._acquisition.getMagnification()
        if mag is None:
            self._scannedPixelSize.set(None)
        else:
            self._scannedPixelSize.set(1e-4 * samplingRate * mag)
    
    def getScannedPixelSize(self):
        return self._scannedPixelSize.get()
    
    def setScannedPixelSize(self, scannedPixelSize):
        """ Set scannedPixelSize and update samplingRate. """
        mag = self._acquisition.getMagnification()
        if mag is None:
            raise Exception("SetOfMicrographs: cannot set scanned pixel size if Magnification is not set.")
        self._scannedPixelSize.set(scannedPixelSize)
        self._samplingRate.set((1e+4 * scannedPixelSize) / mag)


class SetOfParticles(SetOfImages):
    """ Represents a set of Particles.
    The purpose of this class is to separate the
    concepts of Micrographs and Particles, even if
    both are considered Images
    """
    ITEM_TYPE = Particle
    REP_TYPE = Particle
    
    def __init__(self, **args):
        SetOfImages.__init__(self, **args)
        self._coordsPointer = Pointer()
        
    def hasCoordinates(self):
        return self._coordsPointer.hasValue()
    
    def getCoordinates(self):
        """ Returns the SetOfCoordinates associated with 
        this SetOfParticles"""
        return self._coordsPointer.get()
    
    def setCoordinates(self, coordinates):
        """ Set the SetOfCoordinates associates with
        this set of particles.
         """
        self._coordsPointer.set(coordinates)    
        
    def copyInfo(self, other):
        """ Copy basic information (voltage, spherical aberration and sampling rate)
        from other set of micrographs to current one.
        """
        SetOfImages.copyInfo(self, other)
        self.setHasCTF(other.hasCTF())    


class SetOfAverages(SetOfParticles):
    """Represents a set of Averages.
    It is a SetOfParticles but it is useful to differentiate outputs."""
    def __init__(self, **args):
        SetOfParticles.__init__(self, **args)


class SetOfVolumes(SetOfImages):
    """Represents a set of Volumes"""
    ITEM_TYPE = Volume
    
    def __init__(self, **args):
        SetOfImages.__init__(self, **args)


class SetOfCTF(EMSet):
    """ Contains a set of CTF models estimated for a set of images."""
    ITEM_TYPE = CTFModel
    
    def __init__(self, **args):
        EMSet.__init__(self, **args)
        self._micrographsPointer = Pointer()
        
    def getMicrographs(self):
        """ Return the SetOFImages used to create the SetOfClasses. """
        return self._micrographsPointer.get()
    
    def setMicrographs(self, micrographs):
        self._micrographsPointer.set(micrographs)


class SetOfDefocusGroup(EMSet):
    """ Contains a set of DefocusGroup.
        if min/max/avg exists the corresponding flag must be
        set to true.
    """
    ITEM_TYPE = DefocusGroup
        
    def __init__(self, **args):
        EMSet.__init__(self, **args) 
        self._minSet=False
        self._maxSet=False
        self._avgSet=False
        
    def getMinSet(self):
        return self._minSet.get()
        
    def setMinSet(self, value):
        self._minSet.set(value)
        
    def getMaxSet(self):
        return self._maxSet.get()
        
    def setMaxSet(self, value):
        self._maxSet.set(value)
        
    def getAvgSet(self):
        return self._avgSet.get()
        
    def setAvgSet(self, value):
        self._avgSet.set(value)


class Coordinate(EMObject):
    """This class holds the (x,y) position and other information
    associated with a coordinate"""
    def __init__(self, **kwargs):
        EMObject.__init__(self, **kwargs)
        self._micrographPointer = Pointer(objDoStore=False)
        self._x = Integer(kwargs.get('x', None))
        self._y = Integer(kwargs.get('y', None))
        self._micId = Integer()
        
    def getX(self):
        return self._x.get()
    
    def setX(self, x):
        self._x.set(x)
        
    def shiftX(self, shiftX):
        self._x.sum(shiftX)
    
    def getY(self):
        return self._y.get()
    
    def setY(self, y):
        self._y.set(y)
        
    def shiftY(self, shiftY):
        self._y.sum(shiftY)    
    
    def getPosition(self):
        """ Return the position of the coordinate as a (x, y) tuple.
        mode: select if the position is the center of the box
        or in the top left corner.
        """
        return (self.getX(), self.getY())

    def setPosition(self, x, y):
        self.setX(x)
        self.setY(y)
    
    def getMicrograph(self):
        """ Return the micrograph object to which
        this coordinate is associated.
        """
        return self._micrographPointer.get()
    
    def setMicrograph(self, micrograph):
        """ Set the micrograph to which this coordinate belongs. """
        self._micrographPointer.set(micrograph)
        self._micId.set(micrograph.getObjId())
    
    def copyInfo(self, coord):
        """ Copy information from other coordinate. """
        self.setPosition(*coord.getPosition())
        self.setObjId(coord.getObjId())
        self.setBoxSize(coord.getBoxSize())
        
    def getMicId(self):
        return self._micId.get()
    
    def setMicId(self, micId):
        self._micId.set(micId)
    
    def invertY(self):
        if not self.getMicrograph() is None:
            dims = self.getMicrograph().getDim()
            height = dims[1]
            self.setY(height - self.getY())
        #else: error TODO


class SetOfCoordinates(EMSet):
    """ Encapsulate the logic of a set of particles coordinates.
    Each coordinate has a (x,y) position and is related to a Micrograph
    The SetOfCoordinates can also have information about TiltPairs.
    """
    ITEM_TYPE = Coordinate
    
    def __init__(self, **args):
        EMSet.__init__(self, **args)
        self._micrographsPointer = Pointer()
        self._boxSize = Integer()

    def getBoxSize(self):
        """ Return the box size of the particles.
        """
        return self._boxSize.get()
    
    def setBoxSize(self, boxSize):
        """ Set the box size of the particles. """
        self._boxSize.set(boxSize)
    
    def iterMicrographs(self):
        """ Iterate over the micrographs set associated with this
        set of coordinates.
        """
        return self.getMicrographs()
    
    def iterMicrographCoordinates(self, micrograph):
        """ Iterates over the set of coordinates belonging to that micrograph. """
        pass
    
    def iterCoordinates(self, micrograph=None):
        """ Iterate over the coordinates associated with a micrograph.
        If micrograph=None, the iteration is performed over the whole set of coordinates.
        """
        for coord in self:
            if micrograph is None:
                yield coord 
            else:
                if coord.getMicId() == micrograph.getObjId():
                    yield coord
    
    def getMicrographs(self):
        """ Returns the SetOfMicrographs associated with 
        this SetOfCoordinates"""
        return self._micrographsPointer.get()
    
    def setMicrographs(self, micrographs):
        """ Set the SetOfMicrograph associates with 
        this set of coordinates.
         """
        self._micrographsPointer.set(micrographs)
        
    def getFiles(self):
        filePaths = set()
        filePaths.add(self.getFileName())
        return filePaths
    
    def __str__(self):
        """ String representation of a set of coordinates. """
        if self._boxSize.hasValue():
            boxSize = self._boxSize.get()
            boxStr = ' %d x %d' % (boxSize, boxSize)
        else:
            boxStr = 'No-Box'
        s = "%s (%d items, %s)" % (self.getClassName(), self.getSize(), boxStr)
        
        return s
    

class Matrix(Scalar):
    def __init__(self, **args):
        Scalar.__init__(self, **args)
        self._matrix = np.eye(4)
        
    def _convertValue(self, value):
        """Value should be a str with comman separated values
        or a list.
        """
        self._matrix = np.array(json.loads(value))
            
    def getObjValue(self):
        self._objValue = json.dumps(self._matrix.tolist())
        return self._objValue
    
    def setValue(self, i, j, value):
        self._matrix[i, j] = value
        
    def getMatrix(self):
        """ Return internal numpy matrix. """
        return self._matrix
    
    def setMatrix(self, matrix):
        """ Override internal numpy matrix. """
        self._matrix = matrix
        
    def __str__(self):
        return np.array_str(self._matrix)
    
        
class Transform(EMObject):
    """ This class will contain a transformation matrix
    that can be applied to 2D/3D objects like images and volumes.
    It should contain information about euler angles, translation(or shift)
    and mirroring.
    """

    def __init__(self, matrix=None, **args):
        EMObject.__init__(self, **args)
        self._matrix = Matrix()
        if matrix is not None:
            self.setMatrix(matrix)

    def getMatrix(self):
        return self._matrix.getMatrix()
    
    def setMatrix(self, matrix):
        self._matrix.setMatrix(matrix)

    def __str__(self):
        return str(self._matrix)

    def scale(self, factor):
        m = self.getMatrix()
        m *= factor
        m[3, 3] = 1.


class Class2D(SetOfParticles):
    """ Represent a Class that groups Particles objects.
    Usually the representative of the class is another Particle 
    (some kind of average particle from the particles assigned
    to the class) 
    """
    def copyInfo(self, other):
        """ Copy basic information (id and other properties) but not _mapperPath or _size
        from other set of micrographs to current one.
        """
        self.copy(other, copyId=False, ignoreAttrs=['_mapperPath', '_size'])
        
    def clone(self):
        clone = self.getClass()()
        clone.copy(self, ignoreAttrs=['_mapperPath', '_size'])
        return clone
        
    
class Class3D(SetOfParticles):
    """ Represent a Class that groups Particles objects.
    Usually the representative of the class is a Volume 
    reconstructed from the particles assigned to the class.
    """
    def copyInfo(self, other):
        """ Copy basic information (id and other properties) but not _mapperPath or _size
        from other set of micrographs to current one.
        """
        self.copy(other, copyId=False, ignoreAttrs=['_mapperPath', '_size'])
        
    def clone(self):
        clone = self.getClass()()
        clone.copy(self, ignoreAttrs=['_mapperPath', '_size'])
        return clone
        

class ClassVol(SetOfVolumes):
    """ Represent a Class that groups Volume objects.
    Usually the representative of the class is another Volume. 
    """
    pass


class SetOfClasses(EMSet):
    """ Store results from a classification. """
    ITEM_TYPE = None # type of classes stored in the set
    REP_TYPE = None # type of the representatives of each class
    
    def __init__(self, **args):
        EMSet.__init__(self, **args)
        self._representatives = Boolean(False) # Store the average images of each class(SetOfParticles)
        self._imagesPointer = Pointer()

    def iterClassImages(self):
        """ Iterate over the images of a class. """
        pass
    
    def hasRepresentatives(self):
        return self._representatives.get()
    
    def getImages(self):
        """ Return the SetOFImages used to create the SetOfClasses. """
        return self._imagesPointer.get()
    
    def setImages(self, images):
        self._imagesPointer.set(images)
    
    def getDimensions(self):
        """Return first image dimensions as a tuple: (xdim, ydim, zdim)"""
        if self.hasRepresentatives():
            return self.getFirstItem().getRepresentative().getDim()
        return None
    
    def _setItemMapperPath(self, classItem):
        """ Set the mapper path of this class according to the mapper
        path of the SetOfClasses and also the prefix acording to class id
        """
        classPrefix = 'Class%03d' % classItem.getObjId()
        classItem._mapperPath.set('%s,%s' % (self.getFileName(), classPrefix))
        classItem._mapperPath.setStore(False)
        
    def _insertItem(self, classItem):
        """ Create the SetOfImages assigned to a class.
        If the file exists, it will load the Set.
        """
        if classItem.hasRepresentative():
            self._representatives.set(True)
                        
        self._setItemMapperPath(classItem)
        EMSet._insertItem(self, classItem)
        classItem.write()#Set.write(self)
        
    def __getitem__(self, itemId):
        """ Setup the mapper classes before returning the item. """
        classItem = self._mapper.selectById(itemId)
        self._setItemMapperPath(classItem)
        return classItem

    def _iterItems(self):    
        for classItem in self._mapper.selectAll():
            self._setItemMapperPath(classItem)
            yield classItem
            
    def getSamplingRate(self):
        return self.getImages().getSamplingRate()

    def appendFromClasses(self, classesSet):
        """ Iterate over the classes and the elements inside each
        class and append classes and items that are enabled.
        """
        for cls in classesSet:
            if cls.isEnabled():
                newCls = self.ITEM_TYPE()
                newCls.copyInfo(cls)
                newCls.setObjId(cls.getObjId())
                self.append(newCls)
                for img in cls:
                    if img.isEnabled():                
                        newCls.append(img)
                self.update(newCls)
                
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
                self.append(newItem)
                # copy items inside the class
                newItem.copyItems(item, copyDisabled=copyDisabled)
                self.update(newItem)
            else:
                if itemDataIterator is not None:
                    next(itemDataIterator) # just skip disabled data row
                

class SetOfClasses2D(SetOfClasses):
    """ Store results from a 2D classification of Particles. """
    ITEM_TYPE = Class2D
    REP_TYPE = Particle

    pass


class SetOfClasses3D(SetOfClasses):
    """ Store results from a 3D classification of Particles. """
    ITEM_TYPE = Class3D
    REP_TYPE = Volume
    
    pass
       

class SetOfClassesVol(SetOfClasses3D):
    """ Store results from a classification of Volumes. """
    ITEM_TYPE = ClassVol

    pass
    

class NormalMode(EMObject):
    """ Store normal mode information. """
    def __init__(self, **kwargs):
        EMObject.__init__(self, **kwargs)
        self._modeFile = String(kwargs.get('modeFile', None))
        self._collectivity = Float(kwargs.get('collectivity', None))
        self._score = Float(kwargs.get('score', None))
        
    def getModeFile(self):
        return self._modeFile.get()
    
    def setModeFile(self, value):
        self._modeFile.set(value)
        
    def getCollectivity(self):
        return self._collectivity.get()
    
    def setCollectivity(self, value):
        self._collectivity.set(value)
        
    def getScore(self):
        return self._score.get()
    
    def setScore(self, value):
        self._score.set(value)
        
    
class SetOfNormalModes(EMSet):
    """ Set containing NormalMode items. """
    ITEM_TYPE = NormalMode
    
    def __init__(self, **kwargs):
        EMSet.__init__(self, **kwargs)
        # Store a pointer to the PdbFile object
        # from which this normal modes where computed.
        self._pdbPointer = Pointer()
        
    def setPdb(self, pdbFile):
        self._pdbPointer.set(pdbFile)
        
    def getPdb(self):
        return self._pdbPointer.get()
    
    def copyInfo(self, other):
        self._pdbPointer.copy(other._pdbPointer, copyId=False)


class Movie(Micrograph):
    """ Represent a set of frames of micrographs.
    """
    def __init__(self, filename=None, **kwargs):
        Micrograph.__init__(self, filename=filename, **kwargs)
        self._alignment = None

    def isCompressed(self):
        return self.getFileName().endswith('bz2') 
        
    def getDim(self):
        """Return image dimensions as tuple: (Xdim, Ydim, Zdim)
        Consider compressed Movie files"""
        if not self.isCompressed():
            return Image.getDim(self)
        return None
    
    def hasAlignment(self):
        return self._alignment is not None
    
    def getAlignment(self):
        return self._alignment
    
    def setAlignment(self, alignment):
        self._alignment = alignment
    
    
class MovieAlignment(EMObject):
    """ Store the alignment between the different Movie frames.
    Also store the first and last frames used for alignment.
    """
    def __init__(self, first=0, last=0, shifts=[], **kwargs):
        EMObject.__init__(self, **kwargs)
        self._first = Integer(first)
        self._last = Integer(last)
        self._shifts = CsvList(pType=float)
        self._shifts.set(shifts)
        
    def getRange(self):
        """ Return the first and last frames used for alignment.
        The first frame in a movie stack is 0.
        """
        return self._firstFrame.get(), self._lastFrame.get()
    
    def getShifts(self):
        """ Return the list of alignment between one frame
        to another, from first to last frame used.
        """
        return self._shifts


class SetOfMovies(SetOfMicrographs):
    """ Represents a set of Movies. """
    ITEM_TYPE = Movie
    
    def __init__(self, **kwargs):
        SetOfMicrographs.__init__(self, **kwargs)
        self._gainFile = String()
        
    def setGain(self, gain):
        self._gainFile.set(gain)
        
    def getGain(self):
        return self._gainFile.get()
    
    
class SetOfMovieParticles(SetOfParticles):
    """ This is just to distinguish the special case
    when the particles have been extracted from a set of movies.
    """
    pass
    
    
