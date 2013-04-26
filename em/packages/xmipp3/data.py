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
from pyworkflow.utils.path import replaceBaseExt, exists
import xmipp
    
class _XmippMdRow():
    """ Store label and value pairs corresponding to a Metadata row """
    def __init__(self):
        self._labelDict = {} # Dictionary containing labels and values
    
    def setValue(self, *args):
        """args: this list should contains tuples with 
        MetaData Label and the desired value"""
        for label, value in args:
            self._labelDict[label] = value
            
    def iterItems(self):
        return self._labelDict.iteritems()
            
            
class XmippMicrograph(Micrograph, _XmippMdRow):
    """Xmipp implementation for Micrograph"""
    def __init__(self, filename=None, **args):
        _XmippMdRow.__init__(self)
        Micrograph.__init__(self, filename, **args)
        self.setValue((xmipp.MDL_MICROGRAPH, filename))
            
        
class XmippSetOfMicrographs(SetOfMicrographs):
    """Represents a set of Micrographs for Xmipp"""
    def __init__(self, filename=None, **args):
        SetOfMicrographs.__init__(self, filename, **args)
        self._md = xmipp.MetaData()
        
    def append(self, micrograph):
        """Add a micrograph to the set"""
        objId = self._md.addObject()
        # Convert to xmipp micrograph if necessary
        micXmipp = convertMicrograph(micrograph)
        for label, value in micXmipp.iterItems():
            # TODO: Check how to handle correctly unicode type
            # in Xmipp and Scipion
            if type(value) is unicode:
                value = str(value)
#            print "setting label: %s with value: %s" % (label2Str(label), str(value))
#            print " type(value): ", type(value)
            self._md.setValue(label, value, objId)
            
    def sort(self):
        """Sort the set according to MDL_MICROGRAPH"""
        self._md.sort(xmipp.MDL_MICROGRAPH)
        
    def write(self):
        self._md.write(self.getFileName())
        
    def __iter__(self):
        """Iterate over the set of micrographs in the MetaData"""
        md = xmipp.MetaData(self.getFileName())
        
        for objId in md:    
            m = Micrograph(md.getValue(xmipp.MDL_MICROGRAPH, objId))
            if self.hasCTF():
                m.ctfModel = XmippCTFModel(md.getValue(xmipp.MDL_CTF_MODEL, objId)) 
            yield m


class XmippImage(Image, _XmippMdRow):
    """Xmipp implementation for Image"""
    def __init__(self, filename=None, **args):
        _XmippMdRow.__init__(self)
        Image.__init__(self, filename, **args)
        self.setValue((xmipp.MDL_IMAGE, filename))
        
        
class XmippSetOfImages(SetOfImages):
    """Represents a set of Images for Xmipp"""
    def __init__(self, filename=None, **args):
        SetOfImages.__init__(self, filename, **args)
        self._md = xmipp.MetaData()
                
    def __iter__(self):
        """Iterate over the set of images in the MetaData"""
        md = xmipp.MetaData(self.getFileName())
        
        for objId in md:    
            m = Image(md.getValue(xmipp.MDL_IMAGE, objId))
            if self.hasCTF():
                m.ctfModel = XmippCTFModel(md.getValue(xmipp.MDL_CTF_MODEL, objId)) 
            yield m
            
    def setMd(self, md):
        self._md = md
        
    def sort(self):
        """Sort the set according to MDL_IMAGE"""
        self._md.sort(xmipp.MDL_IMAGE)
        
    def write(self):
        self._md.write(self.getFileName())
                    
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
            return (self.x - self.boxSize / 2, self.y - self.boxSize / 2)
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
    

class XmippCTFModel(CTFModel):
    
    ctfParams = {
                 "samplingRate":xmipp.MDL_CTF_SAMPLING_RATE,
                 "voltage":xmipp.MDL_CTF_VOLTAGE,
                 "defocusU":xmipp.MDL_CTF_DEFOCUSU,
                 "defocusV":xmipp.MDL_CTF_DEFOCUSV,
                 "defocusAngle":xmipp.MDL_CTF_DEFOCUS_ANGLE,
                 "sphericalAberration":xmipp.MDL_CTF_CS,
                 "chromaticAberration":xmipp.MDL_CTF_CA,
                 "energyLoss":xmipp.MDL_CTF_ENERGY_LOSS,
                 "lensStability":xmipp.MDL_CTF_LENS_STABILITY,
                 "convergenceCone":xmipp.MDL_CTF_CONVERGENCE_CONE,
                 "longitudinalDisplacement":xmipp.MDL_CTF_LONGITUDINAL_DISPLACEMENT,
                 "transversalDisplacement":xmipp.MDL_CTF_TRANSVERSAL_DISPLACEMENT,
                 "q0":xmipp.MDL_CTF_Q0,
                 "k":xmipp.MDL_CTF_K,
                 "bgGaussianK":xmipp.MDL_CTF_BG_GAUSSIAN_K,
                 "bgGaussianSigmaU":xmipp.MDL_CTF_BG_GAUSSIAN_SIGMAU,
                 "bgGaussianSigmaV":xmipp.MDL_CTF_BG_GAUSSIAN_SIGMAV,
                 "bgGaussianCU":xmipp.MDL_CTF_BG_GAUSSIAN_CU,
                 "bgGaussianCV":xmipp.MDL_CTF_BG_GAUSSIAN_CV,
                 "bgGaussianAngle":xmipp.MDL_CTF_BG_GAUSSIAN_ANGLE,
                 "bgSqrtK":xmipp.MDL_CTF_BG_SQRT_K,
                 "bgSqrtU":xmipp.MDL_CTF_BG_SQRT_U,
                 "bgSqrtV":xmipp.MDL_CTF_BG_SQRT_V,
                 "bgSqrtAngle":xmipp.MDL_CTF_BG_SQRT_ANGLE,
                 "bgBaseline":xmipp.MDL_CTF_BG_BASELINE,
                 "bgGaussian2K":xmipp.MDL_CTF_BG_GAUSSIAN2_K,
                 "bgGaussian2SigmaU":xmipp.MDL_CTF_BG_GAUSSIAN2_SIGMAU,
                 "bgGaussian2SigmaV":xmipp.MDL_CTF_BG_GAUSSIAN2_SIGMAV,
                 "bgGaussian2CU":xmipp.MDL_CTF_BG_GAUSSIAN2_CU,
                 "bgGaussian2CV":xmipp.MDL_CTF_BG_GAUSSIAN2_CV,
                 "bgGaussian2Angle":xmipp.MDL_CTF_BG_GAUSSIAN2_ANGLE,
#                 "X0":xmipp.MDL_CTF_X0,
#                 "XF":xmipp.MDL_CTF_XF,
#                 "Y0":xmipp.MDL_CTF_Y0,
#                 "YF":xmipp.MDL_CTF_YF,
                 "critFitting":xmipp.MDL_CTF_CRIT_FITTINGSCORE,
                 "critCorr13":xmipp.MDL_CTF_CRIT_FITTINGCORR13,
#                 "downsampleFactor":xmipp.MDL_CTF_DOWNSAMPLE_PERFORMED,
                 "critPsdStdQ":xmipp.MDL_CTF_CRIT_PSDVARIANCE,
                 "critPsdPCA1":xmipp.MDL_CTF_CRIT_PSDPCA1VARIANCE,
                 "critPsdPCARuns":xmipp.MDL_CTF_CRIT_PSDPCARUNSTEST
                 }

    def __init__(self, filename = None, **args):

        args['value'] = filename
                    
        CTFModel.__init__(self, **args)       
        if filename is not None:
            # Use the object value to store the filename

            md = xmipp.MetaData(filename)
            objId = md.firstObject()
            
            for key, val in  self.ctfParams.iteritems():
                mdVal = md.getValue(val, objId)
                if not hasattr(self, key):
                    setattr(self, key, Float(mdVal))
                else:
                    getattr(self, key).set(mdVal)
                
    def getFileName(self):
        return self.get()          
    
    def write(self, fn):

        md = xmipp.MetaData()
        objId = md.addObject()
        for key, val in  self.ctfParams.iteritems():
            if hasattr(self, key):
                md.setValue(val, getattr(self, key).get(), objId)
                
        md.write(fn)
        
        self.set(fn)
          
    
class XmippSetOfCoordinates(SetOfCoordinates):
    """Implementation of SetOfCoordinates for Xmipp"""
    def __init__(self, filename=None, **args):
        # Use object value to store filename
        SetOfCoordinates.__init__(self, value=filename, **args)
        self.family = String()
        
    def getFileName(self):
        return self.get()
    
    def iterMicrographs(self):
        """Iterate over the micrographs set associated with this
        set of coordinates
        """
        return self.getMicrographs()
        
    
    def iterCoordinates(self, micrograph=None):
        """Iterates over the whole set of coordinates.
        If the SetOfMicrographs has tilted pairs, the coordinates
        should have the information related to its paired coordinate."""
        
        path = self.getFileName()
        template = self.family.get() + '@%s'
        
        if micrograph is None:
            for mic in self.getMicrographs():
                pathPos = join(path, replaceBaseExt(mic.getFileName(), 'pos'))
                
                if exists(pathPos):
                    mdPos = xmipp.MetaData(template % pathPos)
                                
                    for objId in mdPos:
                        x = mdPos.getValue(xmipp.MDL_XCOOR, objId)
                        y = mdPos.getValue(xmipp.MDL_YCOOR, objId)
                        coordinate = XmippCoordinate()
                        coordinate.setPosition(x, y)
                        coordinate.setMicrograph(mic)
                        
                        yield coordinate
        else:
            pathPos = join(path, replaceBaseExt(micrograph.getFileName(), 'pos'))
                
            if exists(pathPos):
                mdPos = xmipp.MetaData(template % pathPos)
                                
                for objId in mdPos:
                    x = mdPos.getValue(xmipp.MDL_XCOOR, objId)
                    y = mdPos.getValue(xmipp.MDL_YCOOR, objId)
                    coordinate = XmippCoordinate()
                    coordinate.setPosition(x, y)
                    coordinate.setMicrograph(micrograph)
                        
                    yield coordinate
        
    
