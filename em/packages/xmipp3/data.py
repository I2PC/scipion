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
from xmipp import *
    
    
class XmippMicrograph(Micrograph):
    """Xmipp implementation for Micrograph"""
    def __init__(self, filename=None, **args):
        Micrograph.__init__(self, filename, **args)
        self._labelDict = {MDL_MICROGRAPH: filename}
        
    def setValue(self, *args):
        """args: this list should contains tuples with 
        MetaData Label and the desired value"""
        for label, value in args:
            self._labelDict[label] = value
            
        
class XmippSetOfMicrographs(SetOfMicrographs):
    """Represents a set of Micrographs for Xmipp"""
    def __init__(self, filename=None, **args):
        SetOfMicrographs.__init__(self, filename, **args)
        self._md = MetaData()
        
    def append(self, micrograph):
        """Add a micrograph to the set"""
        objId = self._md.addObject()
        # Convert to xmipp micrograph if necessary
        micXmipp = convertMicrograph(micrograph)
        for label, value in micXmipp._labelDict.iteritems():
            # TODO: Check how to handle correctly unicode type
            # in Xmipp and Scipion
            if type(value) is unicode:
                value = str(value)
#            print "setting label: %s with value: %s" % (label2Str(label), str(value))
#            print " type(value): ", type(value)
            self._md.setValue(label, value, objId)
            
    def sort(self):
        """Sort the set according to MDL_MICROGRAPH"""
        self._md.sort(MDL_MICROGRAPH)
        
    def write(self):
        self._md.write(self.getFileName())
        
    def __iter__(self):
        """Iterate over the set of micrographs in the MetaData"""
        md = MetaData(self.getFileName())
        
        for objId in md:    
            m = Micrograph()
            m.setFileName(md.getValue(MDL_MICROGRAPH, objId))
            if self.hasCTF():
                m.ctfModel = XmippCTFModel(md.getValue(MDL_CTF_MODEL, objId)) 
            yield m

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
                 "samplingRate":MDL_CTF_SAMPLING_RATE,
                 "voltage":MDL_CTF_VOLTAGE,
                 "defocusU":MDL_CTF_DEFOCUSU,
                 "defocusV":MDL_CTF_DEFOCUSV,
                 "defocusAngle":MDL_CTF_DEFOCUS_ANGLE,
                 "sphericalAberration":MDL_CTF_CS,
                 "chromaticAberration":MDL_CTF_CA,
                 "energyLoss":MDL_CTF_ENERGY_LOSS,
                 "lensStability":MDL_CTF_LENS_STABILITY,
                 "convergenceCone":MDL_CTF_CONVERGENCE_CONE,
                 "longitudinalDisplacement":MDL_CTF_LONGITUDINAL_DISPLACEMENT,
                 "transversalDisplacement":MDL_CTF_TRANSVERSAL_DISPLACEMENT,
                 "q0":MDL_CTF_Q0,
                 "k":MDL_CTF_K,
                 "bgGaussianK":MDL_CTF_BG_GAUSSIAN_K,
                 "bgGaussianSigmaU":MDL_CTF_BG_GAUSSIAN_SIGMAU,
                 "bgGaussianSigmaV":MDL_CTF_BG_GAUSSIAN_SIGMAV,
                 "bgGaussianCU":MDL_CTF_BG_GAUSSIAN_CU,
                 "bgGaussianCV":MDL_CTF_BG_GAUSSIAN_CV,
                 "bgGaussianAngle":MDL_CTF_BG_GAUSSIAN_ANGLE,
                 "bgSqrtK":MDL_CTF_BG_SQRT_K,
                 "bgSqrtU":MDL_CTF_BG_SQRT_U,
                 "bgSqrtV":MDL_CTF_BG_SQRT_V,
                 "bgSqrtAngle":MDL_CTF_BG_SQRT_ANGLE,
                 "bgBaseline":MDL_CTF_BG_BASELINE,
                 "bgGaussian2K":MDL_CTF_BG_GAUSSIAN2_K,
                 "bgGaussian2SigmaU":MDL_CTF_BG_GAUSSIAN2_SIGMAU,
                 "bgGaussian2SigmaV":MDL_CTF_BG_GAUSSIAN2_SIGMAV,
                 "bgGaussian2CU":MDL_CTF_BG_GAUSSIAN2_CU,
                 "bgGaussian2CV":MDL_CTF_BG_GAUSSIAN2_CV,
                 "bgGaussian2Angle":MDL_CTF_BG_GAUSSIAN2_ANGLE,
#                 "X0":MDL_CTF_X0,
#                 "XF":MDL_CTF_XF,
#                 "Y0":MDL_CTF_Y0,
#                 "YF":MDL_CTF_YF,
                 "critFitting":MDL_CTF_CRIT_FITTINGSCORE,
                 "critCorr13":MDL_CTF_CRIT_FITTINGCORR13,
#                 "downsampleFactor":MDL_CTF_DOWNSAMPLE_PERFORMED,
                 "critPsdStdQ":MDL_CTF_CRIT_PSDVARIANCE,
                 "critPsdPCA1":MDL_CTF_CRIT_PSDPCA1VARIANCE,
                 "critPsdPCARuns":MDL_CTF_CRIT_PSDPCARUNSTEST
                 }
    
    # Implementar el constructor para crear las variables del modelo usando el params de arriba
    # y leyendo del metadata. No hace falta el __getattr__
    
    def __init__(self, filename):
        md = MetaData(filename)
        objId = md.firstObject()
        
        for key, val in  self.ctfParams.iteritems():
            mdVal = md.getValue(val, objId)
            if not hasattr(self, key):
                setattr(self, key, Float(mdVal))
            else:
                getattr(self, key).set(mdVal)
                
                
class XmippSetOfCoordinates(SetOfCoordinates):
    """Implementation of SetOfCoordinates for Xmipp"""
    def __init__(self, filename=None, **args):
        # Use object value to store filename
        SetOfCoordinates.__init__(self, value=filename, **args)
        self.family = String()
        
        
    def iterCoordinates(self):
        """Iterates over the whole set of coordinates.
        If the SetOfMicrographs has tilted pairs, the coordinates
        should have the information related to its paired coordinate."""
        
        path = self.getFileName()
        
        for mic in self.getMicrographs():
            pathPos = join(path, replaceBaseExt(mic.getFilename(), 'pos'))
            
            mdPos = MetaData('%s@%s' % (family.get(), pathPos))
                        
            for objId in mdPos:
                x = mdPos.getValue(MDL_XCOOR, objId)
                y = mdPos.getValue(MDL_YCOOR, objId)
                coordinate = XmippCoordinate()
                coordinate.setPosition(x, y)
                coordinate.setMicrograph(mic)
                
                yield coordinate

        
# Group of converter fuctions
def convertMicrograph(mic):
    """Convert from Micrograph to XmippMicrograph"""
    if type(mic) is XmippMicrograph:
        return mic
    
    micXmipp = XmippMicrograph(mic.getFileName())
    # TODO: copyInfo??
    # from mic to micXmipp??  
    return micXmipp
       
def convertSetOfMicrographs(setOfMics, filename):
    """Method to convert from a general SetOfMicrographs to XmippSetOfMicrographs"""
    if type(setOfMics) is XmippSetOfMicrographs:
        return setOfMics
        
    micsOut = XmippSetOfMicrographs(filename)
    micsOut.copyInfo(setOfMics)

    for mic in setOfMics:
        micsOut.append(mic)

    micsOut.write()
        
    return micsOut
    
