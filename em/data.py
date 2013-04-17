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

from pyworkflow.object import *


class EMObject(Object):
    """Base object for all EM classes"""
    def __init__(self, **args):
        Object.__init__(self, **args)


class Microscope(EMObject):
    """Microscope information"""
    def __init__(self, **args):
        EMObject.__init__(self, **args)
        self.magnification = Float(60000)
        self.voltage = Float(300)
        self.sphericalAberration = Float(1.2)
        

class Image(EMObject):
    """Represents an EM Image object"""
    def __init__(self, **args):
        EMObject.__init__(self, **args)
        self.samplingRate = Float()
        
    def getFormat(self):
        pass
    
    def getDataType(self):
        pass
        
    def getSamplingRate(self):
        pass
    
    def getDim(self):
        """Return image dimensions as tuple: (Ydim, Xdim)"""
        pass
    
    def getFileName(self):
        return self._objValue
    
    def setFileName(self, newFileName):
        self._objValue = newFileName
        
        
class Micrograph(Image):
    """Represents an EM Image object"""
    def __init__(self, **args):
        Image.__init__(self, **args)
        self.ctfModel = None
        
    def getMicroscope(self):
        pass
    
    def hasCTF(self):
        return self.ctfModel is not None

class SetOfImages(EMObject):
    """Represents a set of Images"""
    def __init__(self, **args):
        if 'filename' in args:
            args['value'] = args['filename']
        EMObject.__init__(self, **args)
        self.samplingRate = Float()
        
    def getSize(self):
        """Return the number of images"""
        pass
    
    def getFileName(self):
        return self._objValue
    
    def setFileName(self, newFileName):
        self._objValue = newFileName
    
    
class SetOfMicrographs(SetOfImages):
    """Represents a set of Micrographs"""
    def __init__(self, **args):
        SetOfImages.__init__(self, **args)
        self.microscope = Microscope()
        self._tiltPairs = Boolean(args.get('tiltPairs', False))
        self._ctf = Boolean(args.get('ctf', False))
        self._files = []
        
    def getMicroscope(self, index=0):
        return self.microscope
    
    def hasTiltPairs(self):
        return self._tiltPairs.get()
    
    def hasCTF(self):
        """Return True if the SetOfMicrographs has associated a CTF model"""
        return self._ctf.get()
    
    def append(self, path):
        """Add simply a micrograph path to the set"""
        self._files.append(path)        
    
    def __loadFiles(self):
        """Read files from text files"""
        f = open(self.getFileName())
        self._files = [l.strip() for l in f]
        f.close()
        
    def __iter__(self):
        """Iterate over the set of micrographs in a .txt file"""
        self.__loadFiles()
        for f in self._files:
            m = Micrograph()
            m.setFileName(f)       
            yield m
        
    def iterPairs(self):
        """Iterate over the tilt pairs if is the case"""
        #FIXME, we need to store the real relation 
        # between micrographs
        # TODO: Validate number of micrographs is even for Tilt Pairs
        n = len(self._files) / 2
        for i in range(n):
            mU = Micrograph()
            mU.setFileName(self._files[2*i])
            mT = Micrograph()
            mT.setFileName(self._files[2*i+1])
            yield (mU, mT)
        
    def writeToFile(self, path):
        """This method will be used to persist in a file the
        list of micrographs path contained in this Set
        path: output file path
        micrographs: list with the micrographs path to be stored
        """
        micFile = open(path, 'w+')
        for f in self._files:            
            print >> micFile, f
        micFile.close()
        
    def copyInfo(self, other):
        """ Copy basic information (voltage, spherical aberration and sampling rate)
        from other set of micrographs to current one"""
        self.microscope.voltage.set(other.microscope.voltage.get())
        self.microscope.sphericalAberration.set(other.microscope.sphericalAberration.get())
        self.samplingRate.set(other.samplingRate.get())
        self._tiltPairs.set(other._tiltPairs.get())
        self._ctf.set(other._ctf.get())
    

class Coordinate(EMObject):
    """This class holds the (x,y) position and other information
    associated with a coordinate"""
    POS_CENTER = 0
    POS_TOPLEFT = 1
    
    def __init__(self, **args):
        self._micrograph = None
    
    def getPosition(self, mode=POS_CENTER):
        """Return the position of the coordinate.
        mode: select if the position is the center of the box
          or in the top left corner."""
        pass
    
    def getMicrograph(self):
        """Return the micrograph object to which
        this coordinate is associated"""
        return self._micrograph
    
    def setMicrograph(self, micrograph):
        """Set the micrograph to which this coordinate belongs"""
        self._micrograph = micrograph
    
    def getPair(self):
        """It should return the paired coordinate associate to self.
        If self is an untilted coordinate, getPaired will return the 
        tilted one and viceversa"""
        pass 
    
    
class SetOfCoordinates(EMObject):
    """Encapsulate the logic of a set of particles coordinates.
    Each coordinate has a (x,y) position and is related to a Micrograph
    The SetOfCoordinates can also have information about TiltPairs"""
    
    def __init__(self, **args):
        self._micrographsPointer = Pointer()
    
    def getBoxSize(self):
        """Return the box size of the future particles.
        This can be None, since when the POS_CENTER mode is used,
        the box size is only relevant when extraction"""
        return None
    
    def iterCoordinates(self):
        """Itearates over the whole set of coordinates.
        If the SetOfMicrographs has tilted pairs, the coordinates
        should have the information related to its paired coordinate."""
        pass
    
    def hasTiltPairs(self):
        """Returns True if the SetOfMicrographs has tilted pairs"""
        return self.getMicrographs().hasTiltPairs()
    
    def getMicrographs(self):
        """Returns the SetOfMicrographs associated with 
        this SetOfCoordinates"""
        return self._micrographsPointer.get()
    
    def setMicrographs(self, micrographs):
        """ Set the SetOfMicrograph associates with 
        this set of coordinates """
        self._micrographsPointer.set(micrographs)
        
    
class CTFModel(EMObject):
    """Represents a generic CTF model"""
#TODO: See how this can be generic (no pointing to Xmipp labels
#    ctfParams = {
#                 "ctfSamplingRate":MDL_CTF_SAMPLING_RATE,
#                 "ctfVoltage":MDL_CTF_VOLTAGE,
#                 "ctfDefocusU":MDL_CTF_DEFOCUSU,
#                 "ctfDefocusV":MDL_CTF_DEFOCUSV,
#                 "ctfDefocusAngle":MDL_CTF_DEFOCUS_ANGLE,
#                 "ctfSphericalAberration":MDL_CTF_CS,
#                 "ctfQ0":MDL_CTF_Q0,
#                 "ctfK":MDL_CTF_K
#                }
    def __init__(self, **args):
        EMObject.__init__(self, **args)
        
        #Aqui meter los comunes
        samplingRate = Float()
