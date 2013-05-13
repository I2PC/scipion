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
for specific Xmipp3 EM data objects
"""

from pyworkflow.em import *
from xmipp3 import XmippMdRow, XmippSet

from pyworkflow.utils.path import replaceBaseExt, exists
import xmipp
    
      
    
class XmippImage(XmippMdRow, Image):
    """Xmipp implementation for Image"""
    def __init__(self, filename=None, label=xmipp.MDL_IMAGE, **args):
        self._label = label
        XmippMdRow.__init__(self)
        Image.__init__(self, filename, **args)
        
    def setFileName(self, filename):
        self.setValue(self._label, filename)
        
    def getFileName(self):
        return self.getValue(self._label)
        
    @staticmethod
    def convert(img):
        """ Create a XmippImage from a general Image instance. """
        if isinstance(img, XmippImage):
            return img
        imgXmipp = XmippImage(img.getFileName())
        # TODO: copyInfo??
        return imgXmipp
        
        
class XmippSetOfImages(XmippSet, SetOfImages):
    """Represents a set of Images for Xmipp"""
    def __init__(self, filename=None, label=xmipp.MDL_IMAGE, **args):
        self._label = label
        SetOfImages.__init__(self, filename, **args)
        XmippSet.__init__(self)
       
    def _getItemLabel(self):
        return xmipp.MDL_IMAGE
    
    def _getItemClass(self):
        return XmippImage
            
    @staticmethod
    def convert(setOfImgs, filename):
        return XmippSet.convert(setOfImgs, XmippSetOfImages, filename)

                
class XmippMicrograph(XmippImage, Micrograph):
    """Xmipp implementation for Micrograph"""
    def __init__(self, filename=None, **args):
        XmippImage.__init__(self, filename, label=xmipp.MDL_MICROGRAPH, **args)
        Micrograph.__init__(self, filename, **args)
            
    @staticmethod
    def convert(mic):
        """Convert from Micrograph to XmippMicrograph"""
        if isinstance(mic, XmippMicrograph):
            return mic
        micXmipp = XmippMicrograph(mic.getFileName())
        if mic.hasCTF():
            ctf = mic.ctfModel
            micXmipp.setValue(xmipp.MDL_CTF_DEFOCUSU, ctf.defocusU.get())
            micXmipp.setValue(xmipp.MDL_CTF_DEFOCUSV, ctf.defocusV.get())
            print "psd: ", mic.ctfModel.getPsdFile()
            print type(mic.ctfModel.getPsdFile())
            micXmipp.setValue(xmipp.MDL_IMAGE1, mic.ctfModel.getPsdFile())
        # TODO: copyInfo??
        # from mic to micXmipp??  
        return micXmipp
    
    
class XmippSetOfMicrographs(XmippSetOfImages, SetOfMicrographs):
    """Represents a set of Micrographs for Xmipp"""
    def __init__(self, filename=None, **args):
        SetOfMicrographs.__init__(self, filename, **args)
        XmippSetOfImages.__init__(self, filename, **args)
        
    def _getItemLabel(self):
        return xmipp.MDL_MICROGRAPH
    
    def _getItemClass(self):
        return XmippMicrograph
            
    @staticmethod
    def convert(setOfMics, filename):
        if isinstance(setOfMics, XmippSetOfMicrographs):
            return setOfMics
        xmippMics = XmippSet.convert(setOfMics, XmippSetOfMicrographs, filename)
        
        # If there are tilt pairs add tilt pairs block the the micrograph metadata
        if setOfMics.hasTiltPairs():
            mdPairs = 'TiltedPairs@' + filename
            md = xmipp.MetaData()
            for micU, micT in setOfMics.iterTiltPairs():
                objId = md.addObject()
                md.setValue(xmipp.MDL_MICROGRAPH, str(micU.getFileName()), objId)
                md.setValue(xmipp.MDL_MICROGRAPH_TILTED, str(micT.getFileName()), objId)
            md.write(mdPairs, xmipp.MD_APPEND)

        return xmippMics
                    
class XmippCoordinate(Coordinate):
    """This class holds the (x,y) position and other information
    associated with a Xmipp coordinate (Xmipp coordinates are POS_CENTER mode)"""
    
    def getPosition(self, mode=Coordinate.POS_CENTER):
        """Return the position of the coordinate.
        mode: select if the position is the center of the box
          or in the top left corner."""
        if mode == Coordinate.POS_CENTER:
            return self.x, self.y
        elif mode == Coordinate.POS_TOPLEFT: 
            return (int(self.x) - self._boxSize / 2, int(self.y) - self._boxSize / 2)
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
    

class XmippCTFModel(CTFModel, XmippMdRow):
    # Mapping to general CTFModel attributes
    ctfParams = {
                 "samplingRate":xmipp.MDL_CTF_SAMPLING_RATE,
                 "voltage":xmipp.MDL_CTF_VOLTAGE,
                 "defocusU":xmipp.MDL_CTF_DEFOCUSU,
                 "defocusV":xmipp.MDL_CTF_DEFOCUSV,
                 "defocusAngle":xmipp.MDL_CTF_DEFOCUS_ANGLE,
                 "sphericalAberration":xmipp.MDL_CTF_CS,
                 "ampContrast":xmipp.MDL_CTF_Q0
                 }

    def __init__(self, filename = None, **args):
        CTFModel.__init__(self, **args)
        XmippMdRow.__init__(self)
        self.getFileName = self.get
        self.setFileName = self.set
        self.setFileName(filename)  
          
        if filename is not None:
            # Use the object value to store the filename

            md = xmipp.MetaData(filename)
            objId = md.firstObject()
            self.getFromMd(md, objId)
            
            for key, label in  self.ctfParams.iteritems():
                mdVal = md.getValue(label, objId)
                getattr(self, key).set(mdVal)
      
    
    def write(self, fn):

        for key, label in  self.ctfParams.iteritems():
            self.setValue(label, getattr(self, key).get())
        
        md = xmipp.MetaData()
        md.setColumnFormat(False)
        objId = md.addObject()
        self.setToMd(md, objId)
        md.write(fn)
        
        self.setFileName(fn)
          
    
class XmippSetOfCoordinates(SetOfCoordinates):
    """Implementation of SetOfCoordinates for Xmipp"""
    def __init__(self, filename=None, **args):
        # Use object value to store filename
        # Here filename is the path where pos files can be found
        SetOfCoordinates.__init__(self, value=filename, **args)
        self.family = String()
        
    def getFileName(self):
        return self.get()       
    
    def iterMicrographCoordinates(self, micrograph):
        """ Iterates over the set of coordinates belonging to that micrograph. """
        path = self.getFileName()
        template = self.family.get() + '@%s'
        
        pathPos = join(path, replaceBaseExt(micrograph.getFileName(), 'pos'))
            
        if exists(pathPos):
            mdPos = xmipp.MetaData(template % pathPos)
                            
            for objId in mdPos:
                x = mdPos.getValue(xmipp.MDL_XCOOR, objId)
                y = mdPos.getValue(xmipp.MDL_YCOOR, objId)
                coordinate = XmippCoordinate()
                coordinate.setPosition(x, y)
                coordinate.setMicrograph(micrograph)
                coordinate.setBoxSize(self.boxSize.get())
                yield coordinate
                
    def iterCoordinates(self):
        """Iterates over the whole set of coordinates.
        If the SetOfMicrographs has tilted pairs, the coordinates
        should have the information related to its paired coordinate."""
        
        for mic in self.getMicrographs():
            self.iterMicrographCoordinates(mic)

    def iterTiltPairs(self):
        """Iterate over the tilt pairs if is the case"""
        #TODO: Implement this method
        pass
            
class XmippImageClassAssignment(ImageClassAssignment, XmippMdRow):
    """ Image-class assignments in Xmipp are stored in a MetaData row. """
    def __init__(self, **args):
        XmippMdRow.__init__(self, **args)
        ImageClassAssignment.__init__(self, **args)
        
    def setImage(self, image):
        """ Set associated image. """
        self.setValue(xmipp.MDL_IMAGE, image.getFileName())
        
    def getImage(self):
        """ Get associated image. """
        return XmippImage(self.getValue(xmipp.MDL_IMAGE))
    
    
class XmippClass2D(Class2D):
    """ Xmipp classes are stored in blocks of the MetaData. """
    def __init__(self, classNumber, filename, representative=None, **args):
        Class2D.__init__(self, **args)
        self._number = classNumber
        self._filename = filename
        self._representative = representative
        
    def iterImageAssignemts(self):
        md = xmipp.MetaData('class%06d_images@%s' % 
                            (self._number, self._filename))
        for objId in md:
            imgCA = XmippImageClassAssignment()
            imgCA.readFromMd(md, objId)
            yield imgCA
    
    def getClassRepresentative(self):
        return self._representative
        
        
class XmippClassification2D(Classification2D):
    """ Store results from a 2D classification. """
    def __init__(self, filename=None, classesBlock='classes', **args):
        Classification2D.__init__(self, **args)
        self.getFileName = self.get
        self.setFileName = self.set
        self.setFileName(filename)
        self.classesBlock = String(classesBlock)
        
    def iterClasses(self):
        fn = self.getFileName()
        block = self.classesBlock.get()
        md = xmipp.MetaData('%(block)s@%(fn)s' % locals())
        for objId in md:
            ref = md.getValue(xmipp.MDL_REF, objId)
            img = md.getValue(xmipp.MDL_IMAGE, objId)
            yield XmippClass2D(ref, fn, img)
            
    def getClassesMdFileName(self):
        """ Return the filename with block pointing
        to the classes. 
        """
        return "%s@%s" % (self.classesBlock.get(),
                          self.getFileName())
  
