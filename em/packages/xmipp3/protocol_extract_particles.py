# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Laura del Cano (ldelcano@cnb.csic.es)
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
This sub-package contains the XmippProtExtractParticles protocol
"""


from pyworkflow.em import *  
from pyworkflow.em.packages.xmipp3.data import *
from pyworkflow.utils.path import makePath, removeBaseExt, join, exists
from pyworkflow.utils import runJob
from protlib_particles import runNormalize
from glob import glob
import xmipp


class XmippDefExtractParticles(Form):
    """Create the definition of parameters for
    the XmippProtExtractParticles protocol"""
    def __init__(self):
        Form.__init__(self)
    
        self.addSection(label='Input')
        self.addParam('inputCoordinates', PointerParam, label="Coordinates", 
                      pointerClass='SetOfCoordinates',
                      help='Select the SetOfCoordinates ')
        self.addParam('downsampleType', EnumParam, choices=['original', 'same as picking', 'other'], 
                      default=1, important=True, label='Downsampling type', 
                      help='Select the downsampling type.')
        self.addParam('downFactor', FloatParam, default=2, condition='downsampleType==2',
                      label='Downsampling factor',
                      help='This factor is always referred to the original sampling rate. '
                      'You may use independent downsampling factors for extracting the '
                      'particles, picking them and estimating the CTF. All downsampling '
                      'factors are always referred to the original sampling rate, and '
                      'the differences are correctly handled by Xmipp.')        
        self.addParam('inputMicrographs', PointerParam, label="Micrographs", 
                      condition='downsampleType in (1,2)',
                      pointerClass='SetOfMicrographs',
                      help='Select the original SetOfMicrographs')
        
        self.addParam('family', StringParam, default='DefaultFamily',
                      label='Family',
                      help='For Xmipp coordinates specify the family.')

        self.addParam('boxSize', IntParam, default=0,
                      label='Particle box size',
                      help='In pixels. The box size is the size of the boxed particles, ' 
                      'actual particles may be smaller than this.')
        
        self.addSection(label='Preprocess')
        self.addParam('DoRemoveDust', BooleanParam, default=True, important=True,
                      label='Dust removal (Recommended)', 
                      help='Sets pixels with unusually large values to random values from a Gaussian '
                      'with zero-mean and unity-standard deviation.')
        self.addParam('thresholdDust', FloatParam, default=3.5, condition='DoRemoveDust',
                      label='Threshold for dust removal',
                      help='Pixels with a signal higher or lower than this value times the standard '
                      'deviation of the image will be affected. For cryo, 3.5 is a good value.'
                      'For high-contrast negative stain, the signal itself may be affected so '
                      'that a higher value may be preferable.',
                      expertLevel=LEVEL_ADVANCED)
        self.addParam('DoFlip', BooleanParam, default=True, important=True,
                      label='Phase flipping (Recommended)', 
                      help='Use the information from the CTF to compensate for phase reversals.')
        self.addParam('doInvert', BooleanParam, default=False, important=True,
                      label='Phase flipping (Recommended)', 
                      help='Invert the contrast if your particles are black over a white background.')
        self.addParam('doNormalize', BooleanParam, default=True, important=True,
                      label='Normalize (Recommended)', 
                      help='It subtract a ramp in the gray values and normalizes so that in the '
                      'background there is 0 mean and standard deviation 1.')
        self.addParam('normaType', EnumParam, choices=['OldXmipp','NewXmipp','Ramp'], 
                      default=0, important=True, condition='doNormalize',
                      label='Normalization type', 
                      help='OldXmipp (mean(Image)=0, stddev(Image)=1). '
                      'NewXmipp (mean(background)=0, stddev(background)=1)'
                      'Ramp (subtract background+NewXmipp).',
                      expertLevel=LEVEL_ADVANCED)
        self.addParam('backRadius', IntParam, default=-1, important=True, condition='doNormalize',
                      label='Background radius',
                      help='Pixels outside this circle are assumed to be noise and their stddev '
                      'is set to 1. Radius for background circle definition (in pix.). '
                      'If this value is 0, then half the box size is used.', 
                      expertLevel=LEVEL_ADVANCED)
                
class XmippProtExtractParticles(ProtExtractParticles):
    """Protocol to extract particles from a set of coordinates in the project"""
    
    _definition = XmippDefExtractParticles()
    
#    def __init__(self, **args):
#        
#        ProtExtractParticles.__init__(self, **args)
        
    def _defineSteps(self):
        '''for each micrograph insert the steps to preprocess it
        '''
        print "En defineSteps"
        
        # Set sampling rate and inputMics according to downsample type
        self.inputCoords = self.inputCoordinates.get() 
        self.inputMics = self.inputMicrographs.get()
        
        self.samplingOriginal = self.inputMics.samplingRate.get()

        if self.downsampleType == 1:
            # If 'same as picking' get samplingRate from input micrographs  
            self.inputMics = self.inputCoords.getMicrographs()
            self.samplingFinal = self.inputMics.samplingRate.get()
        elif self.downsampleType == 0:
            # If 'original' get sampling rate from original micrographs
            self.samplingFinal = self.inputMics.samplingRate.get()
        else:
            # IF 'other' multiply the original sampling rate by the factor provided
            self.samplingFinal = self.inputMics.samplingRate.get()*self.downFactor

        # For each micrograph insert the steps
        for mic in self.inputMics:
            fn = mic.getFileName()
            micrographToExtract = fn
        
            # If downsample type is 'other' perform a downsample
            if self.downsampleType == 2:
                fnDownsampled = self._getTmpPath(removeBaseExt(fn)+"_downsampled.xmp")
                downFactor = self.downFactor.get()
                args = "-i %(micrographToExtract)s -o %(fnDownsampled)s --step %(downFactor)f --method fourier" % locals()
                self._insertRunJobStep("xmipp_transform_downsample",args)
                micographToExtract = fnDownsampled
            # If remove dust 
            if self.DoRemoveDust:
                fnNoDust = self._getTmpPath(removeBaseExt(fn)+"_noDust.xmp")
                
                thresholdDust = self.thresholdDust.get() #TODO: remove this extra variable
                args=" -i %(micrographToExtract)s -o %(fnNoDust)s --bad_pixels outliers %(thresholdDust)f" % locals()
                self._insertRunJobStep("xmipp_transform_filter", args)
                micographToExtract = fnNoDust
                
            # Flipping?
            if self.DoFlip:
                fnFlipped = self._getTmpPath(removeBaseExt(fn)+"_flipped.xmp")
                #TODO: Y si el ctfmodel no es de xmipp? como saco un fichero ctf.param?
                if mic.hasCTF():
                    ctf = mic.ctfModel.filename
                    args=" -i %(micrographToExtract)s --ctf %(ctf)s -o %(fnFlipped)s"
                    # If some downsampling has been performed (either before picking or now) pass the downsampling factor 
                    if self.downsampleType != 1:
                        args+=" --downsampling "+str(self.samplingFinal/self.samplingOriginal)
                        self._insertRunJobStep("xmipp_ctf_phase_flip", args)
                        micographToExtract = fnFlipped
                else:
                    #TODO: Raise some type of error!!
                    pass
                
            # Actually extract
            self._insertFunctionStep('extractParticles', micrographToExtract)
                
        # Delete temporary files
        
        
                
        # Insert step to create output objects       
        self._insertFunctionStep('createOutput')
 
    def extractParticles(self, micrographToExtract):
        ''' Extract particles from one micrograph '''
        # Extract 
        outputRoot = self._getExtraPath(removeBaseExt(micrographToExtract))
        
        
        fnExtractList = self._getExtraPath('extract_list.xmd')
        self.createExtractList(fnExtractList)
        
        posMetadata = removeBaseExt(micrographToExtract) + "@" + fnExtractList
        
        print 'posMetadata' , posMetadata
        
        boxSize = self.boxSize.get()
        args = "-i %(micrographToExtract)s --pos %(posMetadata)s -o %(outputRoot)s --Xdim %(boxSize)d" % locals()
        if self.downsampleType != 1:
            args += " --downsampling " + str(self.samplingFinal/self.samplingOriginal)
        if self.doInvert:
            args += " --invert"
    
        runJob(None,"xmipp_micrograph_scissor", args)
        # Normalize 
        if self.doNormalize:
            runNormalize(None, outputRoot + '.stk',self.normaType.get(), self.backRadius.get(), 1)          
                                   
    def createExtractList(self, fnExtractList):
        ''' Create xmipp metadata extract_list with the coordinates for all micrographs '''
        
        #Iterate over the micrographs on the SetOfCoordinates
        for mic in self.inputCoords.iterMicrographs():
            micName = removeBaseExt(mic.getFileName())
            
            # Create micrograph datablock
            coorMD = xmipp.MetaData()
            #Iterate over the coordinates on that micrograph
            for coords in self.inputCoords.iterCoordinates(mic):
                print "%d, %d" % (coords.x, coords.y)        
                coorId = coorMD.addObject()
                coorMD.setValue(xmipp.MDL_XCOOR, int(coords.x), coorId)
                coorMD.setValue(xmipp.MDL_YCOOR, int(coords.y), coorId)
                                
            coorMD.write(micName + "@%s" % fnExtractList, xmipp.MD_APPEND)        
        
    def createOutput(self):
        # Create the SetOfImages object on the database
        mdOut = xmipp.FileName(self._getPath('images.xmd'))        
        imgSet = XmippSetOfImages(str(mdOut))
        imgSet.copyInfo(self.inputMics)
        
        stackFiles = glob(join(self._getExtraPath(),"*.stk"))
        stackFiles.sort()
        imagesMd = xmipp.MetaData()
        for stack in stackFiles:
            fn = stack.replace(".stk",".xmd")
            md = xmipp.MetaData(fn)
            imagesMd.unionAll(md)

        imgSet.setMd(imagesMd)
        imgSet.sort()
        imgSet.write()
    
        self._defineOutputs(outputImages=imgSet)
        

