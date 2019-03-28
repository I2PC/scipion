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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
from itertools import izip
import PIL

try:
  import xmippLib
except Exception as e:
  print("\nXmipp installation have failed or xmipp's C++ binnding have some problem:")
  print(e)
  print("\n  > Please, remove the Xmipp installation (usually 'rm software/em/xmipp') "
        "and re-install it\n")
  sys.exit(-1)

import pyworkflow.utils as pwutils
from pyworkflow.em.constants import *


class ImageHandler(object):
    """ Class to provide several Image manipulation utilities. """
    # TODO: remove dependency from Xmipp
    DT_DEFAULT = xmippLib.DT_DEFAULT
    DT_UNKNOWN = xmippLib.DT_UNKNOWN
    DT_UCHAR = xmippLib.DT_UCHAR
    DT_SCHAR = xmippLib.DT_SCHAR
    DT_USHORT = xmippLib.DT_USHORT
    DT_SHORT = xmippLib.DT_SHORT
    DT_UINT = xmippLib.DT_UINT
    DT_INT = xmippLib.DT_INT
    DT_LONG = xmippLib.DT_LONG
    DT_FLOAT = xmippLib.DT_FLOAT
    DT_DOUBLE = xmippLib.DT_DOUBLE
    DT_COMPLEXSHORT = xmippLib.DT_COMPLEXSHORT
    DT_COMPLEXINT = xmippLib.DT_COMPLEXINT
    DT_COMPLEXFLOAT = xmippLib.DT_COMPLEXFLOAT
    DT_COMPLEXDOUBLE = xmippLib.DT_COMPLEXDOUBLE
    DT_BOOL = xmippLib.DT_BOOL
    DT_LASTENTRY = xmippLib.DT_LASTENTRY
    
    def __init__(self):
        # Now it will use Xmipp image library
        # to read and write most of formats, in the future
        # if we want to be independent of Xmipp, we should have
        # our own image library
        self._img = xmippLib.Image()
        self._imgClass = xmippLib.Image
    
    @classmethod
    def fixXmippVolumeFileName(cls, image):
        """ This method will add :mrc to .mrc volumes
        because for mrc format is not possible to distinguish
        between 3D volumes and 2D stacks.
        """
        # We can not import Volume from top level since
        # data depends on this module
        from pyworkflow.em.data import Volume, Movie
        fn = image.getFileName()
        if isinstance(image, Volume):
            if fn.endswith('.mrc') or fn.endswith('.map'):
                fn += ':mrc'
        elif isinstance(image, Movie):
            if fn.endswith('.mrc'):
                fn += ':mrcs'
            elif fn.endswith('.em'):
                fn += ':ems'
        
        return fn
    
    @classmethod
    def locationToXmipp(cls, location):
        """ Convert an index and filename location
        to a string with @ as expected in Xmipp.
        """
        index, filename = cls._convertToLocation(location)
        if index != NO_INDEX:
            return "%06d@%s" % (index, filename)
        return filename
    
    @classmethod
    def _convertToLocation(cls, location):
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
        
        elif hasattr(location, 'getLocation'):
            # This case includes Image and its subclasses
            outLocation = (location.getIndex(),
                           cls.fixXmippVolumeFileName(location))
        else:
            raise Exception('Can not convert object %s to (index, location)'
                            % type(location))
        
        return outLocation
    
    @classmethod
    def existsLocation(cls, locationObj):
        """ Return True if a given location exists. 
        Location have the same meaning than in _convertToLocation.
        """
        if locationObj is None:
            fn = None
        elif isinstance(locationObj, tuple):
            fn = locationObj[1]
        elif isinstance(locationObj, basestring):
            fn = locationObj
        elif hasattr(locationObj, 'getLocation'):
            # This case includes Image and its subclasses
            fn = locationObj.getLocation()[1]
        else:
            raise Exception('Can not match object %s to '
                            '(index, location)' % type(locationObj))
        
        # If either the location is None or location
        if fn is None:
            return False
        
        # Remove filename format specification such as :mrc, :mrcs or :ems
        if ':' in fn:
            fn = fn.split(':')[0]
        
        return os.path.exists(fn)
    
    @classmethod
    def getSupportedDataType(cls, inDataType, outputFilename):
        """ Returns the most similar data type supported by the
        output format"""
        outDataType = inDataType

        if outputFilename.endswith(".mrc") or outputFilename.endswith(".mrcs"):
            if inDataType == cls.DT_SCHAR:
                outDataType = cls.DT_USHORT

        return outDataType

    def convert(self, inputObj, outputObj, dataType=None, transform=None):
        """ Convert from one image to another.
        inputObj and outputObj can be: tuple, string, or Image subclass 
        (see self._convertToLocation)
        transform: if not None, apply this transformation
        """
        inputLoc = self._convertToLocation(inputObj)
        outputLoc = self._convertToLocation(outputObj)
        
        if outputLoc[1].lower().endswith('.img'):
            # FIXME Since now we can not read dm4 format in Scipion natively
            # we are opening an Eman2 process to read the dm4 file
            convertImage = pwutils.importFromPlugin('eman2.convert',
                                                    'convertImage')
            convertImage(inputLoc, outputLoc)
        else:
            # Read from input
            self._img.read(inputLoc)

            if dataType is not None:
                self._img.convert2DataType(dataType)
            if transform is not None:
                self._img.applyTransforMatScipion(transform.getMatrixAsList())
            # Write to output
            self._img.write(outputLoc)
    
    def convertStack(self, inputFn, outputFn, firstImg=None, lastImg=None,
                     inFormat=None, outFormat=None):
        """ Convert an input stack file into another.
        It is possible to only use a subset of frames to be written in the
            output stack.
        If outFormat/inFomat=None then there will be
        inferred from extension.If firstFrame/lastFrame are not None, the output
        stack will be a subset of input stack. If it are none, the conversion is
        over the whole stack. If the input format is ".dm4" or  ".img" only is
        allowed the conversion of the whole stack.
        """
        # inputLower = inputFn.lower()
        outputLower = outputFn.lower()
        if outputLower.endswith('.img'):
            if (firstImg and lastImg) is None:
                # FIXME Since now we can not read dm4 format in Scipion natively
                # or writing recent .img format
                # we are opening an Eman2 process to read the dm4 file
                convertImage = pwutils.importFromPlugin('eman2.convert',
                                                        'convertImage')
                convertImage(inputFn, outputFn)
            else:
                ext = os.path.splitext(outputFn)[1]
                raise Exception("if convert from %s, firstImg and lastImg "
                                "must be None" % ext)
        # elif inputLower.endswith('.tif'):
        #     # FIXME: It seems that we have some flip problem with compressed
        #     # tif files, we need to check that
        #     if outputLower.endswith('.mrc'):
        #         self.runJob('tif2mrc', '%s %s' % (inputFn, outputFn))
        #     else:
        #         raise Exception("Conversion from tif to %s is not "
        #                        "implemented yet. " % pwutils.getExt(outputFn))
        else:
            # get input dim
            (x, y, z, n) = xmippLib.getImageSize(inputFn)
            
            location = self._convertToLocation(inputFn)
            self._img.read(location, xmippLib.HEADER)
            dataType = self.getSupportedDataType(self._img.getDataType(),
                                                 outputLower)

            if (firstImg and lastImg) is None:
                n = max(z, n)
                firstImg = 1
                lastImg = n
            else:
                n = lastImg - firstImg + 1
            
            # Create empty output stack file to reserve desired space
            xmippLib.createEmptyFile(outputFn, x, y, 1, n, dataType)
            for i, j in izip(range(firstImg, lastImg + 1), range(1, n+1)):
                self.convert((i, inputFn), (j, outputFn))
    
    def getDimensions(self, locationObj):
        """ It will return a tuple with the images dimensions.
        The tuple will contains:
            (x, y, z, n) where x, y, z are image dimensions (z=1 for 2D) and 
            n is the number of elements if stack.
        """
        if self.existsLocation(locationObj):
            location = self._convertToLocation(locationObj)
            fn = location[1]
            ext = pwutils.getExt(fn).lower()
            
            if ext == '.png' or ext == '.jpg':
                im = PIL.Image.open(fn)
                x, y = im.size # (width,height) tuple
                return x, y, 1, 1
            elif ext == '.img':
                # FIXME Since now we can not read dm4 format in Scipion natively
                # or recent .img format
                # we are opening an Eman2 process to read the dm4 file
                getImageDimensions = pwutils.importFromPlugin(
                                        'eman2.convert', 'getImageDimensions')
                return getImageDimensions(fn) # we are ignoring index here
            else:
                self._img.read(location, xmippLib.HEADER)
                return self._img.getDimensions()
        else:
            return None, None, None, None
    
    def getDataType(self, locationObj):
        if self.existsLocation(locationObj):
            location = self._convertToLocation(locationObj)
            self._img.read(location, xmippLib.HEADER)
            return self._img.getDataType()
        else:
            return None
    
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
        
        return xmippLib.compareTwoImageTolerance(loc1, loc2, tolerance)
    
    def computeAverage(self, inputSet):
        """ Compute the average image either from filename or set.
        If inputSet is a filename, we will read the whole stack
        and compute the average from all images.
        If inputSet is a SetOfImages subclass, we will iterate
        and compute the average from all images.
        """
        if isinstance(inputSet, basestring):
            _, _, _, n = self.getDimensions(inputSet)
            if n:
                avgImage = self.read((1, inputSet))
                
                for i in range(2, n+1):
                    self._img.read((i, inputSet))
                    avgImage.inplaceAdd(self._img)
                
                avgImage.inplaceDivide(n)
                return avgImage
        else:
            n = inputSet.getSize()
            if n:
                imageIter = iter(inputSet)
                img = imageIter.next()
                avgImage = self.read(img)
                
                for img in imageIter:
                    self._img.read(self._convertToLocation(img))
                    avgImage.inplaceAdd(self._img)
                
                avgImage.inplaceDivide(n)
                return avgImage
        
        return None
    
    def invertStack(self, inputFn, outputFn):
        # get input dim
        (x, y, z, n) = xmippLib.getImageSize(inputFn)
        # Create empty output stack for efficiency
        xmippLib.createEmptyFile(outputFn, x, y, z, n)
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

    @classmethod
    def __runXmippProgram(cls, program, args):
        """ Internal shortcut function to launch a Xmipp program. """
        xmipp3 = pwutils.importFromPlugin('xmipp3')
        xmipp3.Plugin.runXmippProgram(program, args)

    @classmethod
    def __runEman2Program(cls, program, args):
        """ Internal workaround to launch an EMAN2 program. """
        eman2 = pwutils.importFromPlugin('eman2')
        from pyworkflow.utils.process import runJob
        runJob(None, eman2.Plugin.getProgram(program), args,
               env=eman2.Plugin.getEnviron())
    
    def createCircularMask(self, radius, refImage, outputFile):
        """ Create a circular mask with the given radius (pixels)
        and with the same dimensions of the refImage.
        The radius should be less or equal dim(refImage)/2
        The mask will be stored in 'outputFile'
        """
        inputRef = self.locationToXmipp(refImage)
        self.__runXmippProgram('xmipp_transform_mask',
                               '-i %s --create_mask  %s --mask circular -%d'
                               % (inputRef, outputFile, radius))
    
    def addNoise(self, inputFile, outputFile, std=1., avg=0.):
        """ Add Gaussian noise to an input image (or stack)
        to produce noisy images.
        Params:
            inputFile: the filename of the input images
            outputFile: the filename of the output noisy images
            noiseStd: standard deviation for the Gaussian noise.
        """
        self.__runXmippProgram('xmipp_transform_add_noise',
                               '-i %s -o %s --type gaussian %f %f'
                               % (inputFile, outputFile, std, avg))

    def truncateMask(self, inputFile, outputFile, newDim=None):
        """ Forces the values of a mask to be between 0 and 1.
        Additionally, the output mask can be scaled to a new dimension.

        Params:
            inputFile: the filename of the input either image or volume
            outputFile: the filename of the output either image or volume
            newDim: scale the output Mask to a new dimension if not None
        """
        if inputFile.endswith('.mrc'):
            inputFile += ':mrc'

        ioStr = '-i %s -o %s' % (inputFile, outputFile)

        if newDim:
            self.__runXmippProgram('xmipp_image_resize',
                                   '%s --dim %d' % (ioStr, newDim))
            ioStr = '-i %s' % outputFile
        self.__runXmippProgram('xmipp_transform_threshold',
                               '%s --select below 0 --substitute '
                               'value 0' % ioStr)

        self.__runXmippProgram('xmipp_transform_threshold',
                               '-i %s --select above 1 --substitute '
                               'value 1' % outputFile)

    @classmethod
    def createEmptyImage(cls, fnOut, xDim=1, yDim=1, zDim=1, nDim=1,
                         dataType=None):
        dt = dataType or cls.DT_FLOAT
        xmippLib.createEmptyFile(fnOut, xDim, yDim, zDim, nDim, dt)

    @classmethod
    def isImageFile(cls, imgFn):
        """ Check if imgFn has an image extension. The function
        is implemented in the Xmipp binding.
        """
        return xmippLib.FileName(imgFn).isImage()

    def computeThumbnail(self, inputFn, outputFn, scaleFactor=6, flipOnY=False,
                         flipOnX=False):
        """ Compute a thumbnail of inputFn, save to ouptutFn.
        Optionally choose a scale factor eg scaleFactor=6 will make
        a thumbnail 6 times smaller.
        """
        outputFn = outputFn or self.getThumbnailFn(inputFn)
        args = '"%s" "%s" ' % (inputFn, outputFn)

        process = "--process normalize"
        process += '' if not flipOnY else " --process=xform.flip:axis=y"
        process += '' if not flipOnX else " --process=xform.flip:axis=x"

        args += "--fouriershrink %s %s" % (scaleFactor, process)

        self.__runEman2Program('e2proc2d.py', args)

        return outputFn

    @staticmethod
    def getThumbnailFn(inputFn):
        """Replace the extension in inputFn with thumb.png"""
        return pwutils.replaceExt(inputFn, "thumb.png")

    @classmethod
    def getVolFileName(cls, location):
        if isinstance(location, tuple):
            fn = location[1]
        elif isinstance(location, basestring):
            fn = location
        elif hasattr(location, 'getLocation'):
            fn = location.getLocation()[1]
        else:
            raise Exception('Can not match object %s to (index, location)'
                            % type(location))
        
        if fn.endswith('.mrc') or fn.endswith('.map'):
            fn += ':mrc'
        
        return fn

    @classmethod
    def removeFileType(cls, fileName):
        # Remove filename format specification such as :mrc, :mrcs or :ems
        if ':' in fileName:
            fileName = fileName.split(':')[0]
        return fileName

    @classmethod
    def scaleFourier(cls, inputFn, outputFn, scaleFactor):
        """ Scale an image by cropping in Fourier space. """
        # TODO: Avoid using xmipp program for this
        cls.__runXmippProgram("xmipp_transform_downsample",
                              "-i %s -o %s --step %f --method fourier"
                              % (inputFn, outputFn, scaleFactor))


DT_FLOAT = ImageHandler.DT_FLOAT
