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
    pass


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
        
    def getMicroscope(self):
        pass
    

class SetOfImages(EMObject):
    """Represents a set of Images"""
    def __init__(self, **args):
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
    """Represents a set of Images"""
    def __init__(self, **args):
        SetOfImages.__init__(self, **args)
        self.microscope = Microscope()
        self.tiltPairs = Boolean(False)
        self._files = []
        
    def getMicroscope(self, index=0):
        return self.microscope
    
    def hasTiltPairs(self):
        return self.tiltPairs.get()
    
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
        
    