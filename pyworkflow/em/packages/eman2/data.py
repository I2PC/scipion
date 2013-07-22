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
for specific EMAN2 EM data objects
"""

from pyworkflow.em import * 
from pyworkflow.utils.path import replaceBaseExt, exists
from glob import glob
import os
            
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
            return (int(self.x) + self._boxSize/2, int(self.y) + self._boxSize/2)
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
    def __init__(self, filename=None, **args):
        # Use object value to store filename
        # Here filename is the path where pos filePaths can be found
        SetOfCoordinates.__init__(self, value=filename, **args)
        
    def getFileName(self):
        return self.get()   
    
    def getSize(self):
        """ Return the number of coordinates on the set """
        from glob import glob
        size = 0
        boxDir = glob(join(self.getFileName(), '*.box'))
        for boxFn in boxDir:
            size += sum(1 for line in open(boxFn))
        return size      
        
    def iterMicrographCoordinates(self, micrograph):
        """ Iterates over the set of coordinates belonging to that micrograph. """
        path = self.getFileName()
        pathBox = join(path, replaceBaseExt(micrograph.getFileName(), 'box'))
        if exists(pathBox):
            boxFile =  open (pathBox,"r")
            for line in boxFile:
                if len(line.strip()) > 0:
                    parts = line.strip().split()
                    x = parts[0]
                    y = parts[1]
                    coordinate = EmanCoordinate()
                    coordinate.setMicrograph(micrograph)
                    coordinate.setPosition(x, y)
                    coordinate.setBoxSize(self.boxSize.get())
                    yield coordinate
                else:
                    pass
            boxFile.close()
        
    def iterCoordinates(self):
        """ Iterates over the whole set of coordinates.
        If the SetOfMicrographs has tilted pairs, the coordinates
        should have the information related to its paired coordinate.
        """
        for mic in self.getMicrographs():
            for coord in self.iterMicrographCoordinates(mic):
                yield coord

    def getFiles(self):
        filePaths = set()
        path = self.getFileName()
        for mic in self.getMicrographs():            
            filePath = join(path, replaceBaseExt(mic.getFileName(), 'box'))
            filePaths.add(filePath)
        return filePaths
    
    def hasTiltPairs(self):
        """Returns True if the SetOfMicrographs has tilted pairs"""
        return self.getMicrographs().hasTiltPairs()

class EmanImage(Image):
    """Eman implementation for Image"""
    
    def __init__(self, filename=None, **args):
        Image.__init__(self, filename, **args)
    
class EmanSetOfImages(SetOfImages):
    """Represents a set of Images for Eman"""
    BDB_FORMAT = 'bdb'
    HDF_FORMAT = 'hdf'
    SPI_FORMAT = 'spi'
    IMG_FORMAT = 'img'
    
    _format = BDB_FORMAT

    def load(self):
        """ Load extra data from bdb. """
        if self.getFileName() is None:
            raise Exception("Set filename before calling load()")
        
        os.chdir(os.path.join(os.getcwd(), self.getFileName()))
        self._format = str(EmanDbd.getEmanParamValue('format'))
        os.chdir(os.getcwd())
        
    def loadIfEmpty(self):
        """ Load format only if None. """
        if self._format is None:
            self.load()
                
    def __init__(self, filename=None, **args):
        SetOfImages.__init__(self, filename, **args)
        self._format = args.get('format')
  
    def __iter__(self):
        """ Iterate over the set of images. """
        self.loadIfEmpty()
        if self._format == self.BDB_FORMAT:
            # TODO: Find out how to iterate over images on bdb
            pass
        else:
            imgPaths = glob('*_ptcls.*')
            if len(imgPaths) == 0:
                raise Exception('There are not particles.')
            for imgFn in imgPaths:
                yield EmanImage(imgFn)


class EmanSetOfMicrographs(EmanSetOfImages, SetOfMicrographs):
    """Represents a set of particles for eman2"""
    def __init__(self, filename=None, **args):
        EmanSetOfImages.__init__(self, filename, **args)
        SetOfMicrographs.__init__(self, filename, **args)

                
class EmanSetOfParticles(EmanSetOfImages, SetOfParticles):
    """Represents a set of particles for eman2"""
    def __init__(self, filename=None, **args):
        EmanSetOfImages.__init__(self, filename, **args)
        SetOfParticles.__init__(self, filename, **args)


    
class EmanDbd():
    """ Utility class to access the Eman dbd database """

    @staticmethod    
    def getEmanParamValue(paramName):
        """ Recover a parameter value from EMAN Berkeley data base. """        
        command = "e2bdb.py -D bdb:emboxerbase"
        pipe = os.popen(command)
        stOutput = pipe.readlines()
        pipe.close()
        auxValue = None
        for line in stOutput:
            if (paramName in line):
                auxValue = line.split(" : ")[1]
        if auxValue is None:
            raise Exception("Error getting the stored paramter with command: " + command) 
        return auxValue


class EmanVolume(Image):
    """Eman implementation for Volume"""
    def __init__(self, filename=None, **args):
        Image.__init__(self, filename, **args)
#        Volume.__init__(self, filename, **args)    


class EmanSetOfVolumes(EmanSetOfImages, SetOfVolumes):
    """Represents a set of Volumes for Xmipp"""
    def __init__(self, filename=None, **args):
        SetOfVolumes.__init__(self, filename, **args)
        EmanSetOfImages.__init__(self, filename, **args)

    
def getBoxSize():
    """Method to read boxsize from EMANDB"""
    from EMAN2db import db_open_dict,db_check_dict,db_close_dict

    EMBOXERBASE_DB = "bdb:emboxerbase"
    db = db_open_dict(EMBOXERBASE_DB)
    boxsize = db.get("box_size",dfl=128)
    return boxsize