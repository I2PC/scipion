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

from constants import *
from pyworkflow.object import *
from pyworkflow.mapper.sqlite import SqliteMapper
from posixpath import join
from pyworkflow.utils.utils import getUniqueItems
from pyworkflow.utils.path import exists
import xmipp


class EMObject(OrderedObject):
    """Base object for all EM classes"""
    def __init__(self, **args):
        OrderedObject.__init__(self, **args)
        
    def __str__(self):
        return str(self.get())
    
    def getFiles(self):
        """ Get all filePaths """
        return None


class Microscope(EMObject):
    """Microscope information"""
    def __init__(self, **args):
        EMObject.__init__(self, **args)
        self.magnification = Float(60000)
        self.voltage = Float(300)
        self.sphericalAberration = Float(2.0)
        
    def copyInfo(self, other):
        self.magnification.set(other.magnification.get())
        self.voltage.set(other.voltage.get())
        self.sphericalAberration.set(other.sphericalAberration.get())
    

# TODO: Move this class and Set to a separated base module
class Item(Object):
    """ This class should be subclasses to be used as elements of Set.
    """
    def __init__(self, **args):
        Object.__init__(self, **args)
        self._id =  Integer()
        
    #TODO: replace this id with objId
    def getId(self):
        return self._id.get()
        
    def setId(self, imgId):
        """ This id identifies the element inside a set """
        self._id.set(imgId)
        
    def hasId(self):
        return self._id.hasValue()
        

class CTFModel(EMObject, Item):
    """ Represents a generic CTF model. """
    def __init__(self, **args):
        EMObject.__init__(self, **args)
        Item.__init__(self, **args)
        self.defocusU = Float()
        self.defocusV = Float()
        self.defocusAngle = Float()
        self.psdFile = String()
        self.micFile = String()
        
    def copyInfo(self, other):
        self.copyAttributes(other, 'defocusU', 'defocusV',
                            'defocusAngle', 'psdFile')


class Image(EMObject, Item):
    """Represents an EM Image object"""
    def __init__(self, **args):
        EMObject.__init__(self, **args)
        Item.__init__(self, **args)
        #TODO: replace this id with objId
        # Image location is composed by an index and a filename
        self._index = Integer(0)
        self._filename = String()
        self._samplingRate = Float()
        self._ctfModel = None
        
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
        """Return image dimensions as tuple: (Ydim, Xdim)"""
        imagePath = '%s@%s' % (self._index, self._filename) 
            
        imgXmipp = xmipp.Image(imagePath)
        
        return imgXmipp.getDimensions()
    
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
    
    def setLocation(self, index, filename):
        """ Set the image location, see getLocation. """
        self.setIndex(index)
        self.setFileName(filename)
        
    def hasCTF(self):
        return self._ctfModel is not None
    
    def getCTF(self):
        """ Return the CTF model """
        return self._ctfModel
    
    def setCTF(self, newCTF):
        self._ctfModel = newCTF
        
        
class Micrograph(Image):
    """ Represents an EM Micrograph object """
    def __init__(self, **args):
        Image.__init__(self, **args)
        
        
class Particle(Image):
    """ Represents an EM Particle object """
    def __init__(self, **args):
        Image.__init__(self, **args)


class Volume(Image):
    """ Represents an EM Volume object """
    def __init__(self, **args):
        Image.__init__(self, **args)


class Set(EMObject):
    """ This class will be a container implementation for elements.
    It will use an extra sqlite file to store the elements.
    All items will have an unique id that identifies each element in the set.
    """
    def __init__(self, **args):
        # Use the object value to store the filename
        EMObject.__init__(self, **args)
        self._filename = String() # sqlite filename
        self._mapper = None
        self._idCount = 0
        self._size = Integer(0) # cached value of the number of images    
       
    def getSize(self):
        """Return the number of images"""
        return self._size.get()
    
    def getFileName(self):
        return self._filename.get()
    
    def setFileName(self, filename):
        self._filename.set(filename)
    
    def write(self):
        """This method will be used to persist in a file the
        list of images path contained in this Set
        path: output file path
        images: list with the images path to be stored
        """
        #TODO: If mapper is in memory, do commit and dump to disk
        self._mapper.commit()
    
    def load(self):
        """ Load extra data from files. """
        if self.getFileName() is None:
            raise Exception("Set filename before calling load()")
        self._mapper = SqliteMapper(self.getFileName(), globals())
        
    def loadIfEmpty(self):
        """ Load data only if the main set is empty. """
        if self._mapper is None:
            self.load()
            
    def append(self, item):
        """ Add a image to the set. """
        if not item.hasId():
            self._idCount += 1
            item.setId(self._idCount)
        self.loadIfEmpty()
        self._mapper.insert(item)
        self._size.set(self._size.get() + 1)
                
    
class SetOfImages(Set):
    """ Represents a set of Images """
    def __init__(self, **args):
        Set.__init__(self, **args)
        self._samplingRate = Float()        
        self._hasCtf = Boolean(args.get('ctf', False))
        self._hasAlignment = Boolean(args.get('alignmet', False))
        self._hasProjectionMatrix = Boolean(False)
        self._isPhaseFlippled = Boolean(False)
        self._isAmplitudeCorrected = Boolean(False)
        self._mapper = None
        self._idCount = 0
        self._size = Integer(0) # cached value of the number of images
        self._idMap = {}#FIXME, remove this after id is the one in mapper
    
    def load(self):
        """ Load data only if the main set is empty. """
        Set.load(self)
        self._idMap = {} #FIXME, remove this after id is the one in mapper
        for img in self._mapper.selectByClass("Image", iterate=True):
            self._idMap[img.getId()] = img
            
    def __getitem__(self, imgId):
        """ Get the image with the given id. """
        self.loadIfEmpty() 
            
        return self._idMap.get(imgId, None)

    def __iter__(self):
        """ Iterate over the set of images. """
        self.loadIfEmpty()
        for img in self._mapper.selectByClass("Image", iterate=True):
            yield img            
        
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
        Set.append(self, image)
        self._idMap[image.getId()] = image

    def copyInfo(self, other):
        """ Copy basic information (sampling rate, scannedPixelSize and ctf)
        from other set of images to current one"""
        self.copyAttributes(other, '_samplingRate', '_hasCtf')
        
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
    
    
class SetOfMicrographs(SetOfImages):
    """Represents a set of Micrographs"""
    def __init__(self, **args):
        SetOfImages.__init__(self, **args)
        self._microscope = Microscope()
        self._scannedPixelSize = Float()
        
    def getMicroscope(self, index=0):
        return self._microscope
        
    def copyInfo(self, other):
        """ Copy basic information (voltage, spherical aberration and sampling rate)
        from other set of micrographs to current one.
        """
        SetOfImages.copyInfo(self, other)
        self._microscope.copyInfo(other._microscope)
        self._scannedPixelSize.set(other.getScannedPixelSize())
        
    def setSamplingRate(self, samplingRate):
        """ Set the sampling rate and adjust the scannedPixelSize. """
        self._samplingRate.set(samplingRate)
        self._scannedPixelSize.set(1e-4 * samplingRate * self._microscope.magnification.get())
               
    def getScannedPixelSize(self):
        return self._scannedPixelSize.get()
                       
    def setScannedPixelSize(self, scannedPixelSize):
        """ Set scannedPixelSize and update samplingRate. """
        self._scannedPixelSize.set(scannedPixelSize)
        self._samplingRate.set((1e+4 * scannedPixelSize) / self._microscope.magnification.get())


class SetOfParticles(SetOfImages):
    """ Represents a set of Particles.
    The purpose of this class is to separate the
    concepts of Micrographs and Particles, even if
    both are considered Images
    """
    def __init__(self, **args):
        SetOfImages.__init__(self, **args)
        self._coordsPointer = Pointer()
        
    def getCoordinates(self):
        """ Returns the SetOfCoordinates associated with 
        this SetOfParticles"""
        return self._coordsPointer.get()
    
    def setCoordinates(self, coordinates):
        """ Set the SetOfCoordinates associates with 
        this set of particles.
         """
        self._coordsPointer.set(coordinates)    


class SetOfVolumes(SetOfImages):
    """Represents a set of Volumes"""
    def __init__(self, **args):
        SetOfImages.__init__(self, **args)
        
        
class SetOfCTF(Set):
    """ Contains a set of CTF models estimated for a set of images."""
    def __init__(self, **args):
        Set.__init__(self, **args)      
        
    def __iter__(self):
        """ Iterate over the set of ctfs. """
        self.loadIfEmpty()
        for ctfModel in self._mapper.selectByClass("CTFModel", iterate=True):
            yield ctfModel   
        

class Coordinate(EMObject, Item):
    """This class holds the (x,y) position and other information
    associated with a coordinate"""
    def __init__(self, **args):
        EMObject.__init__(self, **args)
        Item.__init__(self, **args)
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
        self._micId.set(micrograph.getId())
    
    def copyInfo(self, coord):
        """ Copy information from other coordinate. """
        self.setPosition(*coord.getPosition())
        self.setId(coord.getId())
        self.setBoxSize(coord.getBoxSize())
        
    def getMicId(self):
        return self._micId.get()
        
    
class SetOfCoordinates(Set):
    """ Encapsulate the logic of a set of particles coordinates.
    Each coordinate has a (x,y) position and is related to a Micrograph
    The SetOfCoordinates can also have information about TiltPairs.
    """
    
    def __init__(self, **args):
        Set.__init__(self, **args)
        self._micrographsPointer = Pointer()
        self._boxSize = Integer()

    def getBoxSize(self):
        """ Return the box size of the particles.
        """
        return self._boxSize.get()
    
    def setBoxSize(self, boxSize):
        """ Set the box size of the particles.
        """
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
        self.loadIfEmpty()
        for coord in self._mapper.selectByClass("Coordinate", iterate=True):
            if micrograph is None:
                yield coord 
            else:
                if coord.getMicId() == micrograph.getId():
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

    
class ImageClassAssignment(EMObject):
    """ This class represents the relation of
    an image assigned to a class. It serve to
    store additional information like weight, transformation
    or others. 
    """
    def __init__(self, **args):
        EMObject.__init__(self, **args)
        #self._imagePointer = Pointer() # Pointer to image
        # This parameters will dissappear when transformation matrix is used
#         self._anglePsi = Float()
#         self._shiftX = Float()
#         self._shiftY = Float()
#         self._flip = Boolean()
        self._imgId = Integer()
        
#    def setImage(self, image):
#        """ Set associated image. """
#        self._imagePointer.set(image)
#        
#    def getImage(self):
#        """ Get associated image. """
#        return self._imagePointer.get()

    def setImageId(self, imgId):
        """ Set associated image Id. """
        self._imgId.set(imgId)
        
    def getImageId(self):
        """ Get associated image Id. """
        return self._imgId.get()
    
    def setAnglePsi(self, anglePsi):
        self._anglePsi.set(anglePsi)
        
#     def getAnglePsi(self):
#         return self._anglePsi.get()
#     
#     def setShiftX(self, shiftX):
#         self._shiftX.set(shiftX)
#         
#     def getShiftX(self):
#         return self.shiftX.get()
# 
#     def setShiftY(self, shiftY):
#         self._shiftY.set(shiftY)
#         
#     def getShiftY(self):
#         return self.shiftY.get()   
# 
#     def setFlip(self, flip):
#         self._flip.set(flip)
#         
#     def getFlip(self):
#         return self.flip.get()       
     
    
class Class2D(EMObject):
    """ Represent a Class that group some elements 
    from a classification. 
    """
    def __init__(self, **args):
        EMObject.__init__(self, **args)
        self._id =  Integer()
        self._hasRepresentativeImage = Boolean(False)
        self._representativeImage = Image()
        self._imageAssignments = List()
    
    #TODO: replace this id with objId
    def getId(self):
        return self._id.get()
        
    def setId(self, imgId):
        """ This id identifies the element inside a set """
        self._id.set(imgId)
        
    def hasId(self):
        return self._id.hasValue()
    
    def __iter__(self):
        """ Iterate over the assigments of images
        to this particular class.
        """
        for imgCA in self._imageAssignments:
            yield imgCA
            
    def getImageAssignments(self):
        return self._imageAssignments
    
    def addImageClassAssignment(self, imgCA):
        self._imageAssignments.append(imgCA)
    
    def setRepresentativeImage(self, representativeImage):
        self._representativeImage = representativeImage
    
    def getRepresentativeImage(self):
        """ Usually the representative is an average of 
        the images assigned to that class.
        """
        return self._representativeImage
    
    def getHasRepresentativeImage(self):
        """ Return true if have an average image. """
        return self._hasRepresentativeImage.get()
    
    def setHasRepresentativeImage(self, hasRepresentativeImage):
        self._hasRepresentativeImage.set(hasRepresentativeImage)
        
        
class SetOfClasses2D(Set):
    """ Store results from a 2D classification. """
    def __init__(self, **args):
        Set.__init__(self, **args)
        self._hasRepresentativeImages = True # True if the classes have associated average image
        self._imagesPointer = Pointer()
        
    def __iter__(self):
        """ Iterate over all classes. """
        self.loadIfEmpty()
        self._idMap = {} #FIXME, remove this after id is the one in mapper
        for class2D in self._mapper.selectByClass("Class2D", iterate=True):
            self._idMap[class2D.getId()] = class2D
            yield class2D  
            
    def iterImagesClass(self):
        """ Iterate over the images of a class. """
        pass
    
    def hasRepresentativeImages(self):
        return self._hasRepresentativeImages
    
    def getRepresentativeImages(self):
        """ Return a SetOfImages composed by all the average images 
        of the 2D classes. """
        imgs = SetOfImages()
        
        for class2D in self:
            imgs.append(class2D.getHasRepresentativeImage())
            
        return imgs
    
    def getImages(self):
        """ Return the SetOFImages used to create the SetOfClasses2D. """
        return self._imagesPointer.get()
    
    def setImages(self, images):
        self._imagesPointer.set(images)
    
    def __getitem__(self, classId):
        """ Get the class with the given id. """
        #FIXME, remove this after id is the one in mapper
        return self._idMap.get(classId, None)
     
     
