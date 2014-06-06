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

import numpy as np
import json

from constants import *
from convert import ImageHandler
from pyworkflow.object import *
from pyworkflow.mapper.sqlite import SqliteMapper, SqliteFlatMapper
from pyworkflow.utils.path import cleanPath, dirname, join, replaceExt, exists
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
       
    
class CTFModel(EMObject):
    """ Represents a generic CTF model. """
    def __init__(self, **args):
        EMObject.__init__(self, **args)
        self._defocusU = Float(args.get('defocusU', None))
        self._defocusV = Float(args.get('defocusV', None))
        self._defocusAngle = Float(args.get('defocusAngle', None))
        self._psdFile = String()
        self._micFile = String()
        
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
        
    def copyInfo(self, other):
        self.copyAttributes(other, '_defocusU', '_defocusV',
                            '_defocusAngle', '_psdFile', '_micFile')
        
    def getPsdFile(self):
        return self._psdFile.get()
    
    def setPsdFile(self, value):
        self._psdFile.set(value)
        
    def getMicFile(self):
        return self._micFile.get()
    
    def setMicFile(self, value):
        self._micFile.set(value)


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
        """Return image dimensions as tuple: (Xdim, Ydim, Zdim, N)"""
        i, fn = self.getLocation()
        if exists(fn):
            x, y, z, n = ImageHandler().getDimensions(self.getLocation())
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
        
    def getAcquisition(self):
        return self._acquisition
    
    def setAcquisition(self, acquisition):
        self._acquisition = acquisition
        
    def hasAcquisition(self):
        return (self._acquisition is not None and 
                self._acquisition.getMagnification() is not None)
        
    def __str__(self):
        """ String representation of an Image. """
        return "%s (%s, %0.2f A/px)" % (self.getClassName(), ImageDim(*self.getDim()), self.getSamplingRate())


class Micrograph(Image):
    """ Represents an EM Micrograph object """
    def __init__(self, **args):
        Image.__init__(self, **args)


class Particle(Image):
    """ Represents an EM Particle object """
    def __init__(self, **args):
        Image.__init__(self, **args)
        # This may be redundant, but make the Particle
        # object more indenpent for tracking coordinates
        self._coordinate = None
        
    def hasCoordinate(self):
        return self._coordinate is not None
    
    def setCoordinate(self, coordinate):
        self._coordinate = coordinate
        
    def getCoordinate(self):
        return self._coordinate


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
                
      
class EMSet(Set, EMObject):
    def _loadClassesDict(self):
        return globals()
  
  
class SetOfImages(EMSet):
    """ Represents a set of Images """
    ITEM_TYPE = Image
    
    def __init__(self, **args):
        EMSet.__init__(self, **args)
        self._samplingRate = Float()
        self._hasCtf = Boolean(args.get('ctf', False))
        self._hasAlignment = Boolean(args.get('alignmet', False))
        self._hasProjectionMatrix = Boolean(False)
        self._isPhaseFlippled = Boolean(False)
        self._isAmplitudeCorrected = Boolean(False)
        self._acquisition = Acquisition()
        self._firstDim = ImageDim() # Dimensions of the first image
           
    def getAcquisition(self):
        return self._acquisition
        
    def setAcquisition(self, acquisition):
        self._acquisition = acquisition
        
    def hasCTF(self):
        """Return True if the SetOfImages has associated a CTF model"""
        return self._hasCtf.get()  
    
    def setHasCTF(self, value):
        self._hasCtf.set(value)
        
    def hasAlignment(self):
        return self._hasAlignment.get()
    
    def setHasAlignment(self, value):
        self._hasAlignment.set(value)
        
    def hasProjectionMatrix(self):
        return self._hasProjectionMatrix.get()
    
    def setHasProjectionMatrix(self, value):
        self._hasProjectionMatrix.set(value)
        
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
        if not image.getSamplingRate():
            image.setSamplingRate(self.getSamplingRate())
        if self._firstDim.isEmpty():
            self._firstDim.set(image.getDim())
        EMSet.append(self, image)
    
    def copyInfo(self, other):
        """ Copy basic information (sampling rate, scannedPixelSize and ctf)
        from other set of images to current one"""
        self.copyAttributes(other, '_samplingRate')
        self._acquisition.copyInfo(other._acquisition)
        
    def getFiles(self):
        filePaths = set()
        filePaths.add(self.getFileName())
        for item in self:
            # item is an XmippImage or an Image
            filePaths.add(item.getFileName())
            # If it has CTF we must include ctf file
#            if item.hasCTF():
#                # ctf is a XMippCTFModel
#                filePaths.update(item.getCTF().getFiles())
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
            ih.convert(img.getLocation(), (i+1, fnStack))
    
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
    
    def getDimensions(self):
        """Return first image dimensions as a tuple: (xdim, ydim, zdim)"""
        return self.getFirstItem().getDim()
    
    def __str__(self):
        """ String representation of a set of images. """
        sampling = self.getSamplingRate()
        dimStr = "No Dim"
        
        if not sampling:
            print "FATAL ERROR: Object %s has no sampling rate!!!" % self.getName()
            sampling = -999.0
        if self._firstDim.isEmpty():
            try:
                self._firstDim.set(self.getFirstItem().getDim())
                dimStr = str(self._firstDim)
            except Exception, ex:
                print "Error reading dimension: ", ex
                import traceback
                traceback.print_exc()
        s = "%s (%d items, %s, %0.2f A/px)" % (self.getClassName(), self.getSize(), dimStr, sampling)
        return s

    def __iter__(self):
        """ Redefine iteration to set the acquisition to images. """
        for img in self._iterItems():
            img.setAcquisition(self.getAcquisition())
            
            yield img


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
        
        
class SetOfDefocusGroup(EMSet):
    """ Contains a set of DefocusGroup.
        id min/max/avg exists the corresponding flaf must be
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
    def __init__(self, **args):
        EMObject.__init__(self, **args)
        self._micrographPointer = Pointer(objDoStore=False)
        self._x = Integer()
        self._y = Integer()
        self._micId = Integer()
        
    def getX(self):
        return self._x.get()
    
    def setX(self, x):
        self._x.set(x)
    
    def getY(self):
        return self._y.get()
    
    def setY(self, y):
        self._y.set(y)        
    
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
    def __init__(self, **args):
        EMObject.__init__(self, **args)
        self._matrix = Matrix()
        
    def getMatrix(self):
        return self._matrix.getMatrix()
    
    def setMatrix(self, matrix):
        self._matrix.setMatrix(matrix)
        
    def __str__(self):
        return str(self._matrix)
        

class SetOfAlignment(EMSet):
    """ An Aligment is an particular type of Transform.
    A set of transform is usually the result of alignment or multi-reference
    alignment of a SetOfPartices. Each Transformation modifies the original
    image to be the same of a given reference.
    """
    ITEM_TYPE = Transform
    
    def __init__(self, **args):
        EMSet.__init__(self, **args)
        self._particlesPointer = Pointer()

    def getParticles(self):
        """ Return the SetOfParticles from which the SetOfAligment was obtained. """
        return self._particlesPointer.get()
    
    def setParticles(self, particles):
        """ Set the SetOfParticles associated with this SetOfAlignment..
         """
        self._particlesPointer.set(particles)
        

class TransformParams(object):
    """ Class to store transform parameters in the way
    expected by Xmipp/Spider.
    """
    def __init__(self, **args):
        defaults = {'shiftX': 0., 'shiftY': 0., 'shiftZ': 0.,
                    'angleRot': 0., 'angleTilt': 0., 'anglePsi': 0.,
                    'scale': 1., 'mirror': False}.update(args)
        for k, v in defaults.iteritems():
            setattr(self, k, v)


class Class2D(SetOfParticles):
    """ Represent a Class that groups Particles objects.
    Usually the representative of the class is another Particle 
    (some kind of average particle from the particles assigned
    to the class) 
    """
    pass
        
    
class Class3D(SetOfParticles):
    """ Represent a Class that groups Particles objects.
    Usually the representative of the class is a Volume 
    reconstructed from the particles assigned to the class.
    """
    pass


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
    
    def _insertItem(self, classItem):
        """ Create the SetOfImages assigned to a class.
        If the file exists, it will load the Set.
        """
        if classItem.hasRepresentative():
            self._representatives.set(True)
            
        if classItem.getFileName() is None:
            classPrefix = 'Class%03d' % classItem.getObjId()
            classItem._mapperPath.set('%s,%s' % (self.getFileName(), classPrefix))
        EMSet._insertItem(self, classItem)
        classItem.write()#Set.write(self)
           
    
    def write(self):
        """ Override super method to also write the representatives. """
        EMSet.write(self)
#         if self._representatives: # Write if not None
#             self._representatives.write()
            
    def getSamplingRate(self):
        return self.getImages().getSamplingRate()


class SetOfClasses2D(SetOfClasses):
    """ Store results from a 2D classification of Particles. """
    ITEM_TYPE = Class2D
    REP_TYPE = SetOfParticles

    pass


class SetOfClasses3D(SetOfClasses):
    """ Store results from a 3D classification of Particles. """
    ITEM_TYPE = Class3D
    REP_TYPE = SetOfVolumes
    
    pass
       

class SetOfClassesVol(SetOfClasses3D):
    """ Store results from a classification of Volumes. """
    ITEM_TYPE = ClassVol

    pass
    

class NormalModes(EMObject):
    """ Store results from a 2D classification. """
    def __init__(self, filename=None, **args):
        EMObject.__init__(self, **args)
        self._filename = String(filename)
        
    def getFileName(self):
        return self._filename.get()


class Movie(SetOfMicrographs):
    """ Represent a set of frames of micrographs.
    """
    def __init__(self, **args):
        SetOfMicrographs.__init__(self, **args)


class SetOfMovies(EMSet):
    """ Represents a set of Movies. """
    ITEM_TYPE = Movie
    
    def __init__(self, **args):
        EMSet.__init__(self, **args)
        self._acquisition = Acquisition()
        self._samplingRate = Float()
        self._scannedPixelSize = Float()
        self._representatives = Boolean(False)
        self._imagesPointer = Pointer()
    
    def setSamplingRate(self, samplingRate):
        """ Set the sampling rate and adjust the scannedPixelSize. """
        self._samplingRate.set(samplingRate)
        self._scannedPixelSize.set(1e-4 * samplingRate * self._acquisition.getMagnification())
    
    def getScannedPixelSize(self):
        return self._scannedPixelSize.get()
    
    def setScannedPixelSize(self, scannedPixelSize):
        """ Set scannedPixelSize and update samplingRate. """
        self._scannedPixelSize.set(scannedPixelSize)
        self._samplingRate.set((1e+4 * scannedPixelSize) / self._acquisition.getMagnification())
    
    def getSamplingRate(self):
        return self._samplingRate.get()
    
    def getAcquisition(self):
        return self._acquisition
    
    def setAcquisition(self, acquisition):
        self._acquisition = acquisition
    
    def copyInfo(self, other):
        """ Copy basic information (sampling rate, scannedPixelSize and ctf)
        from other set of movies to current one"""
        self.copyAttributes(other, '_samplingRate')
        self._acquisition.copyInfo(other._acquisition)
    
    def getFiles(self):
        filePaths = set()
        filePaths.add(self.getFileName())
        for item in self:
            filePaths.add(item.getFileName())
        return filePaths
    
    def iterMovieFrames(self):
        """ Iterate over the frames of a movie. """
        pass
    
    def hasRepresentatives(self):
        return self._representatives.get()
    
    def getMicrographs(self):
        """ Return the SetOfMicrographs used to create the SetOfMovies. """
        return self._imagesPointer.get()
    
    def setMicrographs(self, micrographs):
        self._imagesPointer.set(micrographs)
    
    def getDimensions(self):
        """Return first micrograph dimensions as a tuple: (xdim, ydim, zdim)"""
        if self.hasRepresentatives():
            return self.getRepresentatives().getDimensions()
    
    def _insertItem(self, movie):
        """ Create the SetOfMicrographs assigned to a Movie.
        If the file exists, it will load the Set.
        """
        if movie.getFileName() is None:
            moviePrefix = 'Movie%03d' % movie.getObjId()
            movie._mapperPath.set('%s,%s' % (self.getFileName(), moviePrefix))
        EMSet._insertItem(self, movie)
        movie.write()#Set.write(self)
