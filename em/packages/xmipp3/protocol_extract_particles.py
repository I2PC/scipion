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
from convert import convertCTFModel 
from pyworkflow.em.packages.xmipp3.data import *
from pyworkflow.utils.path import makePath, removeBaseExt, join, exists
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
                      default=1, important=True, label='Downsampling type', display=EnumParam.DISPLAY_COMBO, 
                      help='Select the downsampling type.')
        self.addParam('downFactor', FloatParam, default=2, condition='downsampleType==2',
                      label='Downsampling factor',
                      help='This factor is always referred to the original sampling rate. '
                      'You may use independent downsampling factors for extracting the '
                      'particles, picking them and estimating the CTF. All downsampling '
                      'factors are always referred to the original sampling rate, and '
                      'the differences are correctly handled by Xmipp.')        
        self.addParam('inputMicrographs', PointerParam, label="Micrographs", 
                      condition='downsampleType != 1',
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
        self.addParam('doRemoveDust', BooleanParam, default=True, important=True,
                      label='Dust removal (Recommended)', 
                      help='Sets pixels with unusually large values to random values from a Gaussian '
                      'with zero-mean and unity-standard deviation.')
        self.addParam('thresholdDust', FloatParam, default=3.5, condition='doRemoveDust',
                      label='Threshold for dust removal',
                      help='Pixels with a signal higher or lower than this value times the standard '
                      'deviation of the image will be affected. For cryo, 3.5 is a good value.'
                      'For high-contrast negative stain, the signal itself may be affected so '
                      'that a higher value may be preferable.',
                      expertLevel=LEVEL_ADVANCED)
        self.addParam('doFlip', BooleanParam, default=True, important=True,
                      label='Phase flipping (Recommended)', 
                      help='Use the information from the CTF to compensate for phase reversals.')
        self.addParam('doInvert', BooleanParam, default=False, important=True,
                      label='Invert contrast', 
                      help='Invert the contrast if your particles are black over a white background.')
        self.addParam('doNormalize', BooleanParam, default=True, important=True,
                      label='Normalize (Recommended)', 
                      help='It subtract a ramp in the gray values and normalizes so that in the '
                      'background there is 0 mean and standard deviation 1.')
        self.addParam('normType', EnumParam, choices=['OldXmipp','NewXmipp','Ramp'], 
                      default=2, condition='doNormalize', display=EnumParam.DISPLAY_COMBO,
                      label='Normalization type', 
                      help='OldXmipp (mean(Image)=0, stddev(Image)=1).\n'
                           'NewXmipp (mean(background)=0, stddev(background)=1)\n'
                           'Ramp (subtract background+NewXmipp).\n',
                      expertLevel=LEVEL_ADVANCED)
        self.addParam('backRadius', IntParam, default=-1, condition='doNormalize',
                      label='Background radius',
                      help='Pixels outside this circle are assumed to be noise and their stddev '
                      'is set to 1. Radius for background circle definition (in pix.). '
                      'If this value is 0, then half the box size is used.', 
                      expertLevel=LEVEL_ADVANCED)
                
class XmippProtExtractParticles(ProtExtractParticles):
    """Protocol to extract particles from a set of coordinates in the project"""
    
    _definition = XmippDefExtractParticles()
    ORIGINAL = 0
    SAME_AS_PICKING = 1
    OTHER = 2
#    def __init__(self, **args):
#        
#        ProtExtractParticles.__init__(self, **args)
        
    def _defineSteps(self):
        """for each micrograph insert the steps to preprocess it
        """       
        # Set sampling rate and inputMics according to downsample type
        self.inputCoords = self.inputCoordinates.get() 

        self.samplingInput = self.inputCoords.getMicrographs().samplingRate.get()
        
        if self.downsampleType.get() == self.SAME_AS_PICKING:
            # If 'same as picking' get samplingRate from input micrographs  
            self.inputMics = self.inputCoords.getMicrographs()
            self.samplingFinal = self.samplingInput
        else:
            self.samplingOriginal = self.inputMicrographs.get().samplingRate.get()
            self.inputMics = self.inputMicrographs.get()
            if self.downsampleType.get() == self.ORIGINAL:
                # If 'original' get sampling rate from original micrographs
                self.samplingFinal = self.samplingOriginal
            else:
                # IF 'other' multiply the original sampling rate by the factor provided
                self.samplingFinal = self.samplingOriginal*self.downFactor.get()
                
        # For each micrograph insert the steps
        for mic in self.inputMics:
            fn = mic.getFileName()
            micrographToExtract = fn
        
            # If downsample type is 'other' perform a downsample
            if self.downsampleType == self.OTHER:
                fnDownsampled = self._getTmpPath(removeBaseExt(fn)+"_downsampled.xmp")
                downFactor = self.downFactor.get()
                args = "-i %(micrographToExtract)s -o %(fnDownsampled)s --step %(downFactor)f --method fourier"
                self._insertRunJobStep("xmipp_transform_downsample", args % locals())
                micographToExtract = fnDownsampled
            # If remove dust 
            if self.doRemoveDust:
                fnNoDust = self._getTmpPath(removeBaseExt(fn)+"_noDust.xmp")
                
                thresholdDust = self.thresholdDust.get() #TODO: remove this extra variable
                args=" -i %(micrographToExtract)s -o %(fnNoDust)s --bad_pixels outliers %(thresholdDust)f"
                self._insertRunJobStep("xmipp_transform_filter", args % locals())
                micographToExtract = fnNoDust
                
            # Flipping if micrograph has CTF model
            if self.doFlip:
                if mic.hasCTF():
                    fnFlipped = self._getTmpPath(removeBaseExt(fn)+"_flipped.xmp")
                    fnCTFTmp = self._getTmpPath("tmp.ctfParam")
                    xmippCTF = convertCTFModel(mic.ctfModel, fnCTFTmp)
                    fnCTF = xmippCTF.getFileName()
                    args = " -i %(micrographToExtract)s --ctf %(fnCTF)s -o %(fnFlipped)s"
                    # If some downsampling has been performed (either before picking or now) pass the downsampling factor 
                    if self.downsampleType != self.ORIGINAL:
                        downFactor = self.samplingFinal/xmippCTF.samplingRate.get()
                        args += " --downsampling %(downFactor)f"
                        self._insertRunJobStep("xmipp_ctf_phase_flip", args % locals())
                        micographToExtract = fnFlipped
                else:
                    #TODO: Raise some type of error!!
                    pass
                
            # Actually extract
            self._insertFunctionStep('extractParticles', micrographToExtract, mic)
                
        # TODO: Delete temporary files
                
        # Insert step to create output objects       
        self._insertFunctionStep('createOutput')
           
 
    def extractParticles(self, micrographToExtract, mic):
        """ Extract particles from one micrograph """
        # Extract 
        outputRoot = str(self._getExtraPath(removeBaseExt(micrographToExtract)))
        
        fnPosFile = self._getExtraPath(removeBaseExt(mic.getFileName()) + '.pos')
        # If it has coordinates extract the particles
        if self._createPosFile(mic, fnPosFile):
            posMetadata = removeBaseExt(micrographToExtract) + "@" + fnPosFile
           
            boxSize = self.boxSize.get()
            args = "-i %(micrographToExtract)s --pos %(posMetadata)s -o %(outputRoot)s --Xdim %(boxSize)d" % locals()
            if self.downsampleType.get() != self.SAME_AS_PICKING:
                args += " --downsampling " + str(self.samplingFinal/self.samplingInput)
            if self.doInvert:
                args += " --invert"
        
            self.runJob(None,"xmipp_micrograph_scissor", args)
            # Normalize 
            if self.doNormalize:
                runNormalize(None, outputRoot + '.stk',self.normaType.get(), self.backRadius.get(), 1)          
                           
            #WHY THIS?????        
            if self.downsampleType.get() == self.OTHER:
                selfile = outputRoot + ".xmd"
                md = xmipp.MetaData(selfile)
                downsamplingFactor=self.samplingFinal/self.samplingInput
                md.operate("Xcoor=Xcoor*%f"%downsamplingFactor)
                md.operate("Ycoor=Ycoor*%f"%downsamplingFactor)
                md.write(selfile)
            
    def _createPosFile(self, mic, fnPosFile):
        """ Create xmipp metadata extract_list with the coordinates for a micrograph """
        micName = removeBaseExt(mic.getFileName())
        mdExtractList = xmipp.MetaData()
        #Iterate over the coordinates on that micrograph
        hasCoords = False
        for coord in self.inputCoords.iterMicrographCoordinates(mic):
            x, y = coord.getPosition(Coordinate.POS_CENTER)    
            coorId = mdExtractList.addObject()
            mdExtractList.setValue(xmipp.MDL_XCOOR, int(x), coorId)
            mdExtractList.setValue(xmipp.MDL_YCOOR, int(y), coorId)
            hasCoords = True
                                
        # Write block only if there are coordinates for this micrograph       
        if hasCoords:
            mdExtractList.write(micName + "@%s" % fnPosFile, xmipp.MD_OVERWRITE)        
        
        return hasCoords
    
    def createOutput(self):
        # Create the SetOfImages object on the database
        imgSet = XmippSetOfImages(self._getPath('images.xmd'))
        imgSet.copyInfo(self.inputMics)
                
        if self.downsampleType == self.OTHER:
            imgSet.samplingRate.set(self.inputMics.samplingRate.get()*self.downFactor.get())
        
        stackFiles = glob(join(self._getExtraPath(),"*.stk"))
        stackFiles.sort()

        for stack in stackFiles:
            fn = stack.replace(".stk",".xmd")
            imgSet.appendFromMd(xmipp.MetaData(fn))
            
        # TODO: If input coordinates have tilt pairs copy the relation
#        if self.inputCoords.hasTiltPairs():
#            imgSet.copyTiltPairs(self.inputCoords)

        #imgSet.setMd(imagesMd)
        imgSet.sort() # We need sort?
        imgSet.write()
    
        self._defineOutputs(outputImages=imgSet)
        

