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
from convert import ImageHandler
from pyworkflow.object import *
from pyworkflow.mapper.sqlite import SqliteMapper, SqliteFlatMapper
from pyworkflow.utils.path import cleanPath, dirname, join, replaceExt
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
        return ImageHandler().getDimensions(self.getLocation())
    
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
        
    def __str__(self):
        """ String representation of an Image. """
        return "%s (index=%d, filename=%s)" % (self.getClassName(), self.getIndex(), self.getFileName())
        
        
class Micrograph(Image):
    """ Represents an EM Micrograph object """
    def __init__(self, **args):
        Image.__init__(self, **args)
        
        
class Particle(Image):
    """ Represents an EM Particle object """
    def __init__(self, **args):
        Image.__init__(self, **args)


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


class PdbFile(EMObject):
    """Represents an EM Image object"""
    def __init__(self, filename=None, pseudoatoms=False, **args):
        EMObject.__init__(self, **args)
        self._filename = String(filename)
        self._pseudoatoms = Boolean(pseudoatoms)
        
    def getFileName(self):
        """ Use the _objValue attribute to store filename. """
        return self._filename.get()
    
    def setFileName(self, filename):
        """ Use the _objValue attribute to store filename. """
        self._filename.set(filename)
        
    def getPseudoAtoms(self):
        return self._pseudoatoms.get()
    
    def setPseudoAtoms(self, value):
        self._pseudoatoms.set(value)
        
    def __str__(self):
        return "%s (pseudoatoms=%s)" % (self.getClassName(), self.getPseudoAtoms())
        
        
class Set(EMObject):
    """ This class will be a container implementation for elements.
    It will use an extra sqlite file to store the elements.
    All items will have an unique id that identifies each element in the set.
    """
    def __init__(self, filename=None, prefix='', mapperClass=SqliteFlatMapper, **args):
        # Use the object value to store the filename
        EMObject.__init__(self, **args)
        self._mapper = None
        self._idCount = 0
        self._size = Integer(0) # cached value of the number of images  
        #self._idMap = {}#FIXME, remove this after id is the one in mapper
        self.setMapperClass(mapperClass)
        self._mapperPath = CsvList() # sqlite filename
        self._mapperPath.trace(self.load) # Load the mapper whenever the filename is changed

#TODO
#        self.itemType
# ESta propierdad indicara que tipo de objeto guarda el set. Tiene que ser inicializada cuando el primer objeto es insertado/cargado.
# Esta variable sera guardada en la DB y la utilizaremos para preguntar. ie si preguntamos por volume tb buscaremos todos los set que tenga
# esta propiedad inicialida a volume
# En caso de que el usuario seleccione uno de los volume del set hay que ver como guardarmos ese value. ie. si el set tiene id = 3 y el volumen 
# que cogemos tiene id=100, se podria guardar como 3, 100


        # If filename is passed in the constructor, it means that
        # we want to create a new object, so we need to delete it if
        # the file exists
        if filename:
            self._mapperPath.set('%s, %s' % (filename, prefix)) # This will cause the creation of the mapper           
        
    def setMapperClass(self, MapperClass):
        """ Set the mapper to be used for storage. """
        Object.__setattr__(self, '_MapperClass', MapperClass)
        
    def __getitem__(self, itemId):
        """ Get the image with the given id. """
        return self._mapper.selectById(itemId)

    def __iterItems(self):
        return self._mapper.selectAll(iterate=True)
    
    def getFirstItem(self):
        """ Return the first item in the Set. """
        return self._mapper.selectFirst()
    
    def __iter__(self):
        """ Iterate over the set of images. """
        return self.__iterItems()
       
    def __len__(self):
        return self._size.get()
    
    def getSize(self):
        """Return the number of images"""
        return self._size.get()
    
    def getFileName(self):
        if len(self._mapperPath):
            return self._mapperPath[0]
        return None
    
    def getPrefix(self):
        if len(self._mapperPath) > 1:
            return self._mapperPath[1]
        return None
    
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
        if self._mapperPath.isEmpty():
            raise Exception("Set.load:  mapper path and prefix not set.")
        fn, prefix = self._mapperPath
        self._mapper = self._MapperClass(fn, globals(), prefix)
            
    def append(self, item):
        """ Add a image to the set. """
        if not item.hasObjId():
            self._idCount += 1
            item.setObjId(self._idCount)
        self._insertItem(item)
        self._size.increment()
#        self._idMap[item.getObjId()] = item
        
    def _insertItem(self, item):
        self._mapper.insert(item)
        
    def update(self, item):
        """ Update an existing item. """
        self._mapper.update(item)
                
    def __str__(self):
        return "%-20s (%d items)" % (self.getClassName(), self.getSize())
    
    def getDimensions(self):
        """Return first image dimensions as a tuple: (xdim, ydim, zdim, n)"""
        return self.getFirstItem().getDim()
    
    def getSubset(self, n):
        """ Return a subset of n element, making a clone of each. """
        subset = []
        for i, item in enumerate(self):
            subset.append(item.clone())
            if i == n:
                break
        return subset
            
                
    
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
        self._acquisition = Acquisition()
           
    def getAcquisition(self, index=0):
        return self._acquisition
        
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
    
# TODO: Create this general purpose function. See example in protocol_apply_mask
#    def readFromStack(self, fnStack):
#        """ Populate the set with the images in the stack """
#        (_,_,_,ndim)=ImageHandler().getDimensions(fnStack)
#        for i in range(1,ndim+1):
#            image=ImageHandler().convert(img.getLocation(), (i, fnStack))
    
    def __str__(self):
        """ String representation of a set of images. """
        if self.getSamplingRate() is None:
            raise Exception("FATAL ERROR: Object %s has no sampling rate!!!" % self.getName())
        s = "%s (%d items, %0.2f A/px)" % (self.getClassName(), self.getSize(), self.getSamplingRate())
        return s
    
    
class SetOfMicrographs(SetOfImages):
    """Represents a set of Micrographs"""
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
        self._scannedPixelSize.set(1e-4 * samplingRate * self._acquisition.getMagnification())
               
    def getScannedPixelSize(self):
        return self._scannedPixelSize.get()
                       
    def setScannedPixelSize(self, scannedPixelSize):
        """ Set scannedPixelSize and update samplingRate. """
        self._scannedPixelSize.set(scannedPixelSize)
        self._samplingRate.set((1e+4 * scannedPixelSize) / self._acquisition.getMagnification())


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
        
    def copyInfo(self, other):
        """ Copy basic information (voltage, spherical aberration and sampling rate)
        from other set of micrographs to current one.
        """
        SetOfImages.copyInfo(self, other)
        self.setHasCTF(other.hasCTF())    


class SetOfVolumes(SetOfImages):
    """Represents a set of Volumes"""
    def __init__(self, **args):
        SetOfImages.__init__(self, **args)
        
        
class SetOfCTF(Set):
    """ Contains a set of CTF models estimated for a set of images."""
    def __init__(self, **args):
        Set.__init__(self, **args)    
        

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


class Transform(EMObject):
    """ This class will contain a transformation matrix
    that can be applied to 2D/3D objects like images and volumes.
    It should contain information about euler angles, translation(or shift)
    and mirroring.
    """
    def __init__(self, **args):
        EMObject.__init__(self, **args)
        from numpy import eye
        self._matrix = eye(4)
        self._matrix[3, 3] = 0.
      
      
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

        
#class ImageClassAssignment(EMObject):
#    """ This class represents the relation of
#    an image assigned to a class. It serve to
#    store additional information like weight, transformation
#    or others. 
#    """
#    def __init__(self, **args):
#        EMObject.__init__(self, **args)
#        #self._imagePointer = Pointer() # Pointer to image
#        # This parameters will dissappear when transformation matrix is used
##         self._anglePsi = Float()
##         self._shiftX = Float()
##         self._shiftY = Float()
##         self._flip = Boolean()
#        self._imgId = Integer()
#        
##    def setImage(self, image):
##        """ Set associated image. """
##        self._imagePointer.set(image)
##        
##    def getImage(self):
##        """ Get associated image. """
##        return self._imagePointer.get()
#
#    def setImageId(self, imgId):
#        """ Set associated image Id. """
#        self._imgId.set(imgId)
#        
#    def getImageId(self):
#        """ Get associated image Id. """
#        return self._imgId.get()
#    
#    def setAnglePsi(self, anglePsi):
#        self._anglePsi.set(anglePsi)
        
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
     
    
class Class2D(SetOfParticles):
    """ Represent a Class that group some elements 
    from a classification. 
    """
    def __init__(self, **args):
        SetOfParticles.__init__(self, **args)
        # This properties should be set when retrieving from the SetOfClasses2D
        self._average = None
    
    def setAverage(self, avgImage):
        self._average = avgImage
    
    def getAverage(self):
        """ Usually the representative is an average of 
        the images assigned to that class.
        """
        return self._average
    
    def hasAverage(self):
        """ Return true if have an average image. """
        return self._average is not None


class SetOfClasses2D(Set):
    """ Store results from a 2D classification. """
    def __init__(self, **args):
        Set.__init__(self, **args)
        self._averages = None # Store the averages images of each class(SetOfParticles)
        self._imagesPointer = Pointer()

    def iterClassImages(self):
        """ Iterate over the images of a class. """
        pass
    
    def hasAverages(self):
        return self._averages is not None
    
    def getAverages(self):
        """ Return a SetOfImages composed by all the average images 
        of the 2D classes. """
        return self._averages
    
    def createAverages(self):
        self._averages = SetOfParticles(filename=self.getFileName(), prefix='Averages')
        if not self.getImages().hasValue():
            raise Exception("SetOfClasses2D.createAverages: you must set the images before creating the averages!!!")
        self._averages.copyInfo(self.getImages())
        return self._averages
    
    def getImages(self):
        """ Return the SetOFImages used to create the SetOfClasses2D. """
        return self._imagesPointer.get()
    
    def setImages(self, images):
        self._imagesPointer.set(images)
        
    def getDimensions(self):
        """Return first image dimensions as a tuple: (xdim, ydim, zdim, n)"""
        if self.hasAverages():
            return self.getAverages().getDimensions()
        
    def _insertItem(self, class2D):
        """ Create the SetOfImages assigned to a class.
        If the file exists, it will load the Set.
        """
        if class2D.getFileName() is None:
            classPrefix = 'Class%03d' % class2D.getObjId()
            class2D._mapperPath.set('%s,%s' % (self.getFileName(), classPrefix))
        Set._insertItem(self, class2D)
        class2D.write()#Set.write(self)
        
    def write(self):
        """ Override super method to also write the averages. """
        Set.write(self)
        if self._averages: # Write if not None
            self._averages.write()
            

class Class3D(SetOfVolumes):
    """ Represent a Class that group some elements 
    from a classification. 
    """
    def __init__(self, **args):
        SetOfVolumes.__init__(self, **args)
        # This properties should be set when retrieving from the SetOfClasses2D
        self._average = None
    
    def setAverage(self, avgVolume):
        self._average = avgVolume
    
    def getAverage(self):
        """ Usually the representative is an average of 
        the volumes assigned to that class.
        """
        return self._average
    
    def hasAverage(self):
        """ Return true if have an average volume. """
        return self._average is not None


class SetOfClasses3D(Set):
    """ Store results from a 3D classification. """
    def __init__(self, **args):
        Set.__init__(self, **args)
        self._averages = None # Store the averages images of each class(SetOfParticles)
        self._imagesPointer = Pointer()

    def iterClassImages(self):
        """ Iterate over the images of a class. """
        pass
    
    def hasAverages(self):
        return self._averages is not None
    
    def getAverages(self):
        """ Return a SetOfVolumes composed by all the average volumes 
        of the 3D classes. """
        return self._averages
    
    def createAverages(self):
        self._averages = SetOfVolumes(filename=self.getFileName(), prefix='Averages')
        if not self.getVolumes().hasValue():
            raise Exception("SetOfClasses3D.createAverages: you must set the volumes before creating the averages!!!")
        self._averages.copyInfo(self.getVolumes())
        
        return self._averages
    
    def getVolumes(self):
        """ Return the SetOfVolumes used to create the SetOfClasses3D. """
        return self._imagesPointer.get()
    
    def setVolumes(self, volumes):
        self._imagesPointer.set(volumes)
        
    def getDimensions(self):
        """Return first volume dimensions as a tuple: (xdim, ydim, zdim, n)"""
        if self.hasAverages():
            return self.getAverages().getDimensions()
        
    def _insertItem(self, class3D):
        """ Create the SetOfVolumes assigned to a class.
        If the file exists, it will load the Set.
        """
        if class3D.getFileName() is None:
            classPrefix = 'Class%03d' % class3D.getObjId()
            class3D._mapperPath.set('%s,%s' % (self.getFileName(), classPrefix))
        Set._insertItem(self, class3D)
        class3D.write()#Set.write(self)
        
    def write(self):
        """ Override super method to also write the averages. """
        Set.write(self)
        if self._averages: # Write if not None
            self._averages.write()


class NormalModes(EMObject):
    """ Store results from a 2D classification. """
    def __init__(self, filename=None, **args):
        EMObject.__init__(self, **args)
        self._filename = String(filename)
        
    def getFileName(self):
        return self._filename.get()
    
