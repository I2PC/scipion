# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Laura del Cano         (ldelcano@cnb.csic.es)
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
This modules contains data classes related to Random Connical Tilt workflow.
"""

#NOTE: Some of this importS are needed by the mapper,
# not directly in the code
from pyworkflow.em.data import (EMObject, EMSet, Micrograph, 
                                Acquisition, Particle, Coordinate)
from pyworkflow.object import Float, Pointer, Integer, String, Object 



class TiltPair(EMObject):
    def __init__(self, untilted=None, tilted=None, **kwargs):
        EMObject.__init__(self, **kwargs)
        self._untilted = untilted
        self._tilted = tilted
        
        
class MicrographsTiltPair(EMSet):
    """Represents a Micrographs Tilt Pair"""
    ITEM_TYPE = TiltPair
    
    def __init__(self, **args):
        EMSet.__init__(self, **args)
        self._tilted = None#SetOfMicrographs()
        self._untilted = None#SetOfMicrographs()
        
    def getUntilted(self):
        return self._untilted
    
    def getTilted(self):
        return self._tilted
    
    def setUntilted(self, untilted):
        self._untilted = untilted
        
    def setTilted(self, tilted):
        self._tilted = tilted
        
    def getFiles(self):
        filePaths = set()
        filePaths.add(self.getTilted().getFiles)
        filePaths.add(self.getUntilted().getFiles)
        
        return filePaths 
    
    def _loadClassesDict(self):
        return globals() 
    
    
class CoordinatesTiltPair(EMObject):
    """Represents a Coordinates Tilt Pair"""
    
    def __init__(self, **args):
        EMObject.__init__(self, **args)
        self._tilted = None#SetOfMicrographs()
        self._untilted = None#SetOfMicrographs()
        self._angles = SetOfAngles()
        self._micsPair = Pointer()
        
    def getUntilted(self):
        return self._untilted
    
    def getTilted(self):
        return self._tilted
    
    def getAngles(self):
        return self._angles
    
    def getMicsPair(self):
        return self._micsPair.get()
    
    def setUntilted(self, untilted):
        self._untilted = untilted
        
    def setTilted(self, tilted):
        self._tilted = tilted
        
    def setAngles(self, setAngles):
        self._angles = setAngles
        
    def setMicsPair(self, micsPair):
        self._micsPair.set(micsPair)
        
    def getFiles(self):
        filePaths = set()
        filePaths.add(self.getTilted().getFiles)
        filePaths.add(self.getUntilted().getFiles)
        
        return filePaths
         
class Angles(EMObject):
    """Represents a triplet of angles"""

    def __init__(self, **args):
        EMObject.__init__(self, **args)
        self._angleY = Float()    
        self._angleY2 = Float()
        self._angleTilt = Float()
        
    def setAngles(self, angleY, angleY2, angleTilt):
        self._angleY.set(angleY)
        self._angleY2.set(angleY2)
        self._angleTilt.set(angleTilt)
        
    def getAngles(self):
        return (self._angleY.get(), self._angleY2.get(), self._angleTilt.get())
    
    
class SetOfAngles(EMSet):
    """ Represents a set of Images """
    ITEM_TYPE = Angles
    
    def __init__(self, **args):
        EMSet.__init__(self, **args)
        
    def _loadClassesDict(self):
        return globals()  
    
    
class ParticlesTiltPair(EMSet):
    """Represents a Particles Tilt Pair"""
    ITEM_TYPE = TiltPair
    
    def __init__(self, **args):
        EMSet.__init__(self, **args)
        self._tilted = None#SetOfImages()
        self._untilted = None#SetOfImages()
        self._coordsPair = Pointer()
        
    def getUntilted(self):
        return self._untilted
    
    def getTilted(self):
        return self._tilted

    def getCoordsPair(self):
        return self._coordsPair.get()
        
    def setUntilted(self, untilted):
        self._untilted = untilted
        
    def setTilted(self, tilted):
        self._tilted = tilted

    def setCoordsPair(self, coordsPair):
        self._coordsPair.set(coordsPair)
                
    def getFiles(self):
        filePaths = set()
        filePaths.add(self.getTilted().getFiles)
        filePaths.add(self.getUntilted().getFiles)
        
        return filePaths    
    
    def _loadClassesDict(self):
        return globals() 