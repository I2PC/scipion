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
import sys
import PIL

from constants import NO_INDEX
import xmipp
from constants import *
from pyworkflow.utils import runJob, getExt

# TODO: remove dependency from Xmipp
DT_FLOAT = xmipp.DT_FLOAT


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

    def _convertToLocation(self, location):
        """ Get a location in a tuple format (index, filename).
        location could be:
            tuple -> (index, filename)
            string -> (NO_INDEX, filename)
            image -> (image.getIndex(), image.getFileName())
        """
        if isinstance(location, tuple):
            outLocation = location
        
        elif isinstance(location, basestring):
            outLocation = (NO_INDEX, location)
            
        elif hasattr(location, 'getLocation'): #this include Image and subclasses
            # In this case inputLoc should be a subclass of Image
            outLocation = (location.getIndex(), self._fixVolumeFileName(location))
            
        else:
            raise Exception('Can not convert object %s to (index, location)' % type(location))
        
        return outLocation
    
    def _existsLocation(self, location):
        """ Return True if a given location exists. 
        Location have the same meaning than in _convertToLocation.
        """
        if isinstance(location, tuple):
            fn = location[1]
        elif isinstance(location, basestring):
            fn = location
        elif hasattr(location, 'getLocation'): #this include Image and subclasses
            # In this case inputLoc should be a subclass of Image
            fn = location.getLocation()[1]
        else:
            raise Exception('Can not match object %s to (index, location)' % type(location))

        return os.path.exists(fn.replace(':mrc', ''))
        
    def convert(self, inputObj, outputObj, dataType=None):
        """ Convert from one image to another.
        inputObj and outputObj can be: tuple, string, or Image subclass 
        (see self._convertToLocation)
        """
        # Read from input
        self._img.read(self._convertToLocation(inputObj))
        
        if dataType is not None:
            self._img.convert2DataType(dataType)
        # Write to output
        self._img.write(self._convertToLocation(outputObj))
        
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
        
    def getDimensions(self, locationObj):
        """ It will return a tuple with the images dimensions.
        The tuple will contains:
            (x, y, z, n) where x, y, z are image dimensions (z=1 for 2D) and 
            n is the number of elements if stack.
        """
        
        if self._existsLocation(locationObj):
            
            location = self._convertToLocation(locationObj)
            fn = location[1]
            ext = getExt(fn).lower()
            
#             print "Extension %s" % ext
            
            if ext == '.png' or ext == '.jpg':
#                 print "Reading with PIL"
                im = PIL.Image.open(fn)
                x, y = im.size # (width,height) tuple
                return x, y, 1, 1
            
            else:
#                 print "Reading with Xmipp"
                self._img.read(location, xmipp.HEADER)
                x, y, z, n = self._img.getDimensions()
                return x, y, z, n
        
        else: 
            return None, None, None, None
        
    def read(self, inputObj):
        """ Create a new Image class from inputObj 
        (inputObj can be tuple, str or Image subclass). """
        location = self._convertToLocation(inputObj)
        
        return self._imgClass(location)
    
    def createImage(self):
        return self._imgClass()
    
    def write(self, image, outputObj):
        """ Write to disk an image from outputObj 
        (outputObj can be tuple, str or Image subclass). """
        location = self._convertToLocation(outputObj)
        image.write(location)
        
    def compareData(self, locationObj1, locationObj2, tolerance=0.0001):
        """ Compare if two locations have the same binary data.
        """
        loc1 = self._convertToLocation(locationObj1)
        loc2 = self._convertToLocation(locationObj2)
        
        return xmipp.compareTwoImageTolerance(loc1, loc2, tolerance)

    def computeAverage(self, inputSet):
        n = inputSet.getSize()
        if n > 0:
            imageIter = iter(inputSet)
            img = imageIter.next()
            avgImage = self.read(img)

            for img in imageIter:
                self._img.read(self._convertToLocation(img))
                avgImage.inplaceAdd(self._img)

            avgImage.inplaceDivide(n)
            return avgImage
        else:
            return None
        
    def invertStack(self, inputFn, outputFn):
        #get input dim
        (x,y,z,n) = xmipp.getImageSize(inputFn)
        #Create empty output stack for efficiency
        xmipp.createEmptyFile(outputFn,x,y,z,n)
        # handle image formats
        for i in range(1, n+1):
            self.invert((i, inputFn), (i, outputFn))
    
    def invert(self, inputObj, outputObj):
        """ invert the pixels.
        inputObj and outputObj can be: tuple, string, or Image subclass 
        (see self._convertToLocation)
        """
        # Read from input
        self._img.read(self._convertToLocation(inputObj))
        self._img.inplaceMultiply(-1)
        # Write to output
        self._img.write(self._convertToLocation(outputObj))
        
    def createCircularMask(self, radius, refImage, outputFile):
        """ Create a circular mask with the given radius (pixels)
        and with the same dimensions of the refImage.
        The radius should be less or equal dim(refImage)/2
        The mask will be stored in 'outputFile'
        """
        #TODO: right now we need to call an xmipp program to create 
        # the spherical mask, it would be nicer to have such utility in the binding
        import pyworkflow.em.packages.xmipp3 as xmipp3
        xmippEnv = xmipp3.getEnviron()
        inputRef = xmipp3.getImageLocation(refImage)
        runJob(None, 'xmipp_transform_mask', 
                    '-i %s --create_mask  %s --mask circular -%d' % (inputRef, outputFile, radius),
                    env=xmippEnv)
        
    def isImageFile(self, imgFn):
        """ Check if imgFn has an image extension. The function
        is implemented in the xmipp binding."""
        return xmipp.FileName(imgFn).isImage()


def downloadPdb(pdbId, pdbFile, log=None):
    pdbGz = pdbFile + ".gz"
    result = (__downloadPdb(pdbId, pdbGz, log) and 
              __unzipPdb(pdbGz, pdbFile, log))
    return result
    
def __downloadPdb(pdbId, pdbGz, log):
    import ftplib
    """Download a pdb file given its id. """
    if log:
        log.info("File to download and unzip: %s" % pdbGz)
    
    pdborgHostname = "ftp.wwpdb.org"
    pdborgDirectory = "/pub/pdb/data/structures/all/pdb/"
    prefix = "pdb"
    suffix = ".ent.gz"
    success = True
    # Log into serverhttp://www.rcsb.org/pdb/files/2MP1.pdb.gz
    ftp = ftplib.FTP()
    try:
        ftp.connect(pdborgHostname)
        ftp.login()
    except ftplib.error_temp:
        if log:
            log.error("ERROR! Timeout reached!")
        success = False
    
    if success:
        # Download  file
        _fileIn = "%s/%s%s%s" % (pdborgDirectory, prefix, pdbId, suffix) 
        _fileOut = pdbGz
        try:
            ftp.retrbinary("RETR %s" % _fileIn, open(_fileOut, "wb").write)
        except ftplib.error_perm:
            os.remove(_fileOut)
            if log:
                log.error("ERROR!  %s could not be retrieved!" % _fileIn)
            success = False
        # Log out
        ftp.quit()
        
    return success

# TODO unzip may go to utilities
def __unzipPdb(pdbGz, pdbFile, log, cleanFile=True):
    """
    Unzip a pdb file.
    Params:
        pdbGz: zipped pdb file.
        pdbFile: output pdb file.
        cleanFile: remove the zipped file.
    """
    import gzip
    success = True
    try:
        f = gzip.open(pdbGz, 'r')
        g = open(pdbFile, 'w')
        g.writelines(f.readlines())
        f.close()
        g.close()
    except:
        e = sys.exc_info()[0]
        if log:
            log.error('ERROR opening gzipped file %s: %s' % (pdbGz, e))
        success = False
    
    try:
        if success:
            os.remove(pdbGz)
    except:
        e = sys.exc_info()[0]
        if log:
            log.error('ERROR deleting gzipped file: %s' % e)
        success = False
        
    return success
