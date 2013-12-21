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


class ImageHandler(object):
    """ Class to provide several Image manipulation utilities. """
    def __init__(self):
        # Now it will use Xmipp image library
        # to read and write most of formats, in the future
        # if we want to be indepent of Xmipp, we should have
        # our own image library
        import xmipp
        from packages.xmipp3 import locationToXmipp
        
        self._img = xmipp.Image()
        self._imgClass = xmipp.Image
        self._locationToStr = locationToXmipp
        
    def convert(self, inputLoc, outputLoc):
        """ Convert from one image to another.
        Params:
            inputLoc: input location (index and filename)
            outputLoc: output location (index and filename)
        """
        if isinstance(inputLoc, str):
            inputLoc = (NO_INDEX, inputLoc)
        if isinstance(outputLoc, str):
            outputLoc = (NO_INDEX, outputLoc)
        # Read from input
        inputStr = self._locationToStr(*inputLoc)
        self._img.read(inputStr)
        # Write to output
        outputStr = self._locationToStr(*outputLoc)
        self._img.write(outputStr)
        
    def getDimensions(self, location):
        """ It will return a tuple with the images dimensions.
        The tuple will contains:
            (x, y, z, n) where x, y, z are image dimensions (z=1 for 2D) and 
            n is the number of elements if stack.
        """
        if isinstance(location, str):
            location = (NO_INDEX, location)
        img = self.read(location)
        return img.getDimensions()
    
    def read(self, location):
        """ Read an image from a location. """
        return self._imgClass(self._locationToStr(*location))
    
    def write(self, image, location):
        """ Write an image to disk. """
        image.write(self._locationToStr(*location))
        
        
def paramsToTransform(params):
    """ Convert from geometry params (with Spider/Xmipp standard) to
    the Transform object (transformation matrix)
    """
    # TODO: Remove dependency to Xmipp
    import xmipp
    from data import Transform
    t = Transform()
    e = xmipp.Euler_angles2matrix(params.angleRot, params.angleTilt, params.anglePsi)
    m = t._matrix
    m[:3, :3] = e
    m[:, 3] = [params.shiftX, params.shiftY, params.shiftZ, 0.]
    # TODO: Consider mirror and scale
    return t

       
def transformToParams(transform):
    """ Inverse convertion from a Transform to param. """
    import xmipp
    from data import TransformParams
    p = TransformParams()    
    m = transform._matrix
    e = m[:3, :3]
    p.angleRot, p.angleTilt, p.anglePsi = xmipp.Euler_matrix2angles(e)
    p.shiftX, p.shiftY, p.shiftZ = m[:, 3]
    # TODO: Consider mirror and scale
    return p

