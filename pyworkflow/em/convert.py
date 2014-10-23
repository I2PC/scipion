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
This module contains several conversion utilities
"""

import os
from constants import NO_INDEX
import xmipp
from constants import *

# TODO: remove dependency from Xmipp


class ImageHandler(object):
    """ Class to provide several Image manipulation utilities. """
    def __init__(self):
        # Now it will use Xmipp image library
        # to read and write most of formats, in the future
        # if we want to be indepent of Xmipp, we should have
        # our own image library
        from packages.xmipp3 import fixVolumeFileName
        
        self._img = xmipp.Image()
        self._imgClass = xmipp.Image
        self._fixVolumeFileName = fixVolumeFileName

    def convertStack(self, inputFn, outputFn, inFormat=None, outFormat=None):
        """ convert format of stack file. Output/input format is
        specified by outFormat/inFormat. If outFormat/inFomat=None then
        there will be inferred from extension
        """
        #get input dim
        (x,y,z,n) = xmipp.getImageSize(inputFn)
        #Create empty output stack for efficiency
        xmipp.createEmptyFile(outputFn,x,y,z,n)
        # handle image formats
        for i in range(1, n+1):
            self.convert((i, inputFn), (i, outputFn))

    def _convertToLocation(self, location):
        """ Get a location in a tuple format (index, filename).
        location could be:
            tuple -> (index, filename)
            string -> (NO_INDEX, filename)
            image -> (image.getIndex(), image.getFileName())
        """
        if isinstance(location, tuple):
            outLocation = location
        
        elif isinstance(location, str):
            outLocation = (NO_INDEX, location)
            
        elif hasattr(location, 'getLocation'): #this include Image and subclasses
            # In this case inputLoc should be a subclass of Image
            outLocation = (location.getIndex(), self._fixVolumeFileName(location))
            
        else:
            raise Exception('Can not convert object %s to (index, location)' % type(location))
        
        return outLocation
        
    def convert(self, inputObj, outputObj):
        """ Convert from one image to another.
        inputObj and outputObj can be: tuple, string, or Image subclass 
        (see self._convertToLocation)
        """
        # Read from input
        self._img.read(self._convertToLocation(inputObj))
        # Write to output
        self._img.write(self._convertToLocation(outputObj))
        
    def getDimensions(self, locationObj):
        """ It will return a tuple with the images dimensions.
        The tuple will contains:
            (x, y, z, n) where x, y, z are image dimensions (z=1 for 2D) and 
            n is the number of elements if stack.
        """
        location = self._convertToLocation(locationObj)
        self._img.read(location, xmipp.HEADER)
        
        return self._img.getDimensions()
    
    def read(self, inputObj):
        """ Create a new Image class from inputObj 
        (inputObj can be tuple, str or Image subclass). """
        location = self._convertToLocation(inputObj)
        
        return self._imgClass(location)
    
    def write(self, image, outputObj):
        """ Write to disk an image from outputObj 
        (outputObj can be tuple, str or Image subclass). """
        location = self._convertToLocation(outputObj)
        image.write(location)

