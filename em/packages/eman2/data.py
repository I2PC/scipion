# **************************************************************************
# *
# * Authors:     Antonio Poza Ballesteros (apoza@cnb.csic.es)
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
# *  e-mail address 'apoza@cnb.csic.es'
# *
# **************************************************************************
"""
This modules contains basic hierarchy
for specific Xmipp3 EM data objects
"""

from pyworkflow.em import *  
from pyworkflow.utils.path import removeBaseExt
    
class EmanSetOfMicrographs(SetOfMicrographs):
    
    def __iter__(self):
        """Iterate over the set of micrographs in the MetaData"""
        metaDataFile=open(self.getFileName(),"r")
        for line in metaDataFile:
            micrograph = Micrograph()
            micrograph.setFileName(line)       
            yield micrograph
        metaDataFile.close()
            
class EmanCoordinate(Coordinate):
    """This class holds the (x,y) position and other information
    associated with a EMAN coordinate (Eman coordinates are POS_TOPLEFT mode)"""    
    
    def getPosition(self, mode=Coordinate.POS_TOPLEFT):
        """Return the position of the coordinate.
        mode: select if the position is the center of the box
          or in the top left corner."""
        if mode == Coordinate.POS_TOPLEFT:
            return self.x, self.y
        elif mode == Coordinate.POS_CENTER: 
            return (self.x + self.boxSize/2, self.y + self.boxSize/2)
        else:
            raise Exception("No coordinate mode registered for : " + str(mode)) 
    
    def setPosition(self, x, y):
        self.x = x
        self.y = y
    
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
    
    
class EmanSetOfCoordinates(SetOfCoordinates):
    """Encapsulate the logic of a set of particles coordinates for EMAN.
    Each coordinate has a (x,y) position and is related to a Micrograph
    The SetOfCoordinates can also have information about TiltPairs.
    EMAN coordinates are taken from top left"""
    
    def iterCoordinates(self):
        """Iterates over the whole set of coordinates.
        If the SetOfMicrographs has tilted pairs, the coordinates
        should have the information related to its paired coordinate."""
        
        for micrograph in self._micrographsPointer.get():
            micrographFileNameNE = removeBaseExt(micrograph.getFileName())
            boxFileName = micrographFileNameNE + ".box"
            boxFile =  open (boxFileName,"r")
            for line in boxFile:
                if len(line.strip()) > 0:
                    parts = line.strip().split()
                    x = parts[0]
                    y = parts[1]
                    coordinate = EmanCoordinate()
                    coordinate.setMicrograph(micrograph)
                    coordinate.setPosition(x, y)
                    coordinate.setBoxSize(self.boxSize)
                    yield coordinate
                else:
                    pass
            boxFile.close()
    
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
        