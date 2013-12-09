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

from convert import writeCTFModel, writeSetOfCoordinates, readSetOfParticles
from pyworkflow.utils.path import makePath, removeBaseExt, join, exists
from xmipp3 import XmippProtocol
from glob import glob
import xmipp
   
                
class XmippProtExtractParticles(ProtExtractParticles, XmippProtocol):
    """Protocol to extract particles from a set of coordinates in the project"""
    
    # Normalization type constants
    ORIGINAL = 0
    SAME_AS_PICKING = 1
    OTHER = 2
    
    # Rejection method constants
    NONE = 0
    MAXZSCORE = 1
    PERCENTAGE = 2
        
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputCoordinates', PointerParam, label="Coordinates", 
                      pointerClass='SetOfCoordinates',
                      help='Select the SetOfCoordinates ')
        form.addParam('downsampleType', EnumParam, choices=['original', 'same as picking', 'other'], 
                      default=1, important=True, label='Downsampling type', display=EnumParam.DISPLAY_COMBO, 
                      help='Select the downsampling type.')
        form.addParam('downFactor', FloatParam, default=2, condition='downsampleType==2',
                      label='Downsampling factor',
                      help='This factor is always referred to the original sampling rate. '
                      'You may use independent downsampling factors for extracting the '
                      'particles, picking them and estimating the CTF. All downsampling '
                      'factors are always referred to the original sampling rate, and '
                      'the differences are correctly handled by Xmipp.')        
        form.addParam('inputMicrographs', PointerParam, label="Micrographs", 
                      condition='downsampleType != 1',
                      pointerClass='SetOfMicrographs',
                      help='Select the original SetOfMicrographs')

        form.addParam('ctfRelations', RelationParam, allowNull=True,
                      label='CTF relations', relationName=RELATION_CTF, relationParent='getInputMicrographs', 
                      relationReverse=True, help='Choose the CTF.  \n  ')     

        form.addParam('boxSize', IntParam, default=0,
                      label='Particle box size', validators=[Positive],
                      help='In pixels. The box size is the size of the boxed particles, ' 
                      'actual particles may be smaller than this.')
        
        form.addParam('rejectionMethod', EnumParam, choices=['None','MaxZscore', 'Percentage'], 
                      default=0, display=EnumParam.DISPLAY_COMBO,
                      label='Automatic particle rejection',
                      help='How to automatically reject particles. It can be none (no rejection),'
                      ' maxZscore (reject a particle if its Zscore is larger than this value), '
                      'Percentage (reject a given percentage in each one of the screening criteria).', 
                      expertLevel=LEVEL_ADVANCED)
        
        form.addParam('maxZscore', IntParam, default=3, condition='rejectionMethod==1',
                      label='Maximum Zscore',
                      help='Maximum Zscore above which particles are rejected.', 
                      expertLevel=LEVEL_ADVANCED)
        
        form.addParam('percentage', IntParam, default=5, condition='rejectionMethod==2',
                      label='Percentage (%)',
                      help='Percentage of particles to reject', expertLevel=LEVEL_ADVANCED)
        
        form.addSection(label='Preprocess')
        form.addParam('doRemoveDust', BooleanParam, default=True, important=True,
                      label='Dust removal (Recommended)', 
                      help='Sets pixels with unusually large values to random values from a Gaussian '
                      'with zero-mean and unity-standard deviation.')
        form.addParam('thresholdDust', FloatParam, default=3.5, condition='doRemoveDust',
                      label='Threshold for dust removal',
                      help='Pixels with a signal higher or lower than this value times the standard '
                      'deviation of the image will be affected. For cryo, 3.5 is a good value.'
                      'For high-contrast negative stain, the signal itself may be affected so '
                      'that a higher value may be preferable.',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('doFlip', BooleanParam, default=True, important=True,
                      label='Phase flipping (Recommended)', 
                      help='Use the information from the CTF to compensate for phase reversals.')
        form.addParam('doInvert', BooleanParam, default=False, important=True,
                      label='Invert contrast', 
                      help='Invert the contrast if your particles are black over a white background.')
        form.addParam('doNormalize', BooleanParam, default=True, important=True,
                      label='Normalize (Recommended)', 
                      help='It subtract a ramp in the gray values and normalizes so that in the '
                      'background there is 0 mean and standard deviation 1.')
        form.addParam('normType', EnumParam, choices=['OldXmipp','NewXmipp','Ramp'], 
                      default=2, condition='doNormalize', display=EnumParam.DISPLAY_COMBO,
                      label='Normalization type', 
                      help='OldXmipp (mean(Image)=0, stddev(Image)=1).  \n  '
                           'NewXmipp (mean(background)=0, stddev(background)=1)  \n  '
                           'Ramp (subtract background+NewXmipp).  \n  ',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('backRadius', IntParam, default=-1, condition='doNormalize',
                      label='Background radius',
                      help='Pixels outside this circle are assumed to be noise and their stddev '
                      'is set to 1. Radius for background circle definition (in pix.). '
                      'If this value is 0, then half the box size is used.', 
                      expertLevel=LEVEL_ADVANCED)
        
    def getInputMicrographs(self):
        """ Return the micrographs associated to the SetOfCoordinates. """
        if self.inputCoordinates.get() is None:
            return None
        return self.inputCoordinates.get().getMicrographs()
        
    def _defineSteps(self):
        """for each micrograph insert the steps to preprocess it
        """       
        # Set sampling rate and inputMics according to downsample type
        self.inputCoords = self.inputCoordinates.get() 

        self.samplingInput = self.inputCoords.getMicrographs().getSamplingRate()
        
        if self.downsampleType.get() == self.SAME_AS_PICKING:
            # If 'same as picking' get samplingRate from input micrographs  
            self.inputMics = self.inputCoords.getMicrographs()
            self.samplingFinal = self.samplingInput
        else:
            self.inputMics = self.inputMicrographs.get()
            self.samplingOriginal = self.inputMics.getSamplingRate()
            if self.downsampleType.get() == self.ORIGINAL:
                # If 'original' get sampling rate from original micrographs
                self.samplingFinal = self.samplingOriginal
            else:
                # IF 'other' multiply the original sampling rate by the factor provided
                self.samplingFinal = self.samplingOriginal*self.downFactor.get()
                
        # Convert input SetOfCoordinates to Xmipp if needed
#        self._insertConvertStep('inputCoords', XmippSetOfCoordinates, 
#                                 self._getExtraPath('scipion_micrographs_coordinates.xmd'))
                
        # Write pos files for each micrograph
        self._insertFunctionStep('writePosFiles')
                
        #if self.doFlip.get():
        ctfSet = self.ctfRelations.get()
           
        # For each micrograph insert the steps
        for mic in self.inputMics:
            micrographToExtract = mic.getFileName()
            micName = removeBaseExt(mic.getFileName())
            micId = mic.getId()
                                            
            # If downsample type is 'other' perform a downsample
            if self.downsampleType == self.OTHER:
                fnDownsampled = self._getTmpPath(micName+"_downsampled.xmp")
                downFactor = self.downFactor.get()
                args = "-i %(micrographToExtract)s -o %(fnDownsampled)s --step %(downFactor)f --method fourier"
                self._insertRunJobStep("xmipp_transform_downsample", args % locals())
                micrographToExtract = fnDownsampled
            # If remove dust 
            if self.doRemoveDust:
                fnNoDust = self._getTmpPath(micName+"_noDust.xmp")
                
                thresholdDust = self.thresholdDust.get() #TODO: remove this extra variable
                args=" -i %(micrographToExtract)s -o %(fnNoDust)s --bad_pixels outliers %(thresholdDust)f"
                self._insertRunJobStep("xmipp_transform_filter", args % locals())
                micrographToExtract = fnNoDust
                
            #if self.doFlip.get():
            if self.ctfRelations.hasValue():
                mic.ctfModel = ctfSet[micId]       
                        
            #self._insertFunctionStep('getCTF', micId, micName, micrographToExtract)
            micName = removeBaseExt(mic.getFileName())
            fnCTF = None
            #FIXME: Check only if mic has CTF when implemented ok
            #if self.doFlip or mic.hasCTF():
            if self.ctfRelations.hasValue():
                # Write CTF metadata if it does not exist
                mdFn = getattr(mic, '_xmippMd', None)
                if mdFn:
                    fnCTF = mdFn.get()
                else:
                    fnCTF = self._getTmpPath("%s.ctfParam" % micName)
                mic.ctfModel.sphericalAberration = Float(self.inputMics.getAcquisition().sphericalAberration.get())
                writeCTFModel(mic.ctfModel, fnCTF)  
                 
                # Insert step to flip micrograph
                if self.doFlip:
                    self._insertFunctionStep('flipMicrograph', micName, fnCTF, micrographToExtract)
                    micrographToExtract = self._getTmpPath(micName +"_flipped.xmp")
                           
            # Actually extract
            self._insertFunctionStep('extractParticles', micId, micName, fnCTF, micrographToExtract)
                
        # TODO: Delete temporary files
                        
        # Insert step to create output objects       
        self._insertFunctionStep('createOutput')
        
        
    def writePosFiles(self):
        """ Write the pos file for each micrograph on metadata format. """
        #self.posFiles = writeSetOfCoordinates(self._getExtraPath(), self.inputCoords)
        writeSetOfCoordinates(self._getExtraPath(), self.inputCoords)

    def flipMicrograph(self, micName, fnCTF, micrographToExtract):
        """ Flip micrograph. """           
        fnFlipped = self._getTmpPath(micName +"_flipped.xmp")

        args = " -i %(micrographToExtract)s --ctf %(fnCTF)s -o %(fnFlipped)s"
        # If some downsampling has been performed (either before picking or now) pass the downsampling factor 
        if self.downsampleType != self.ORIGINAL:
            downFactor = self.samplingFinal/self.samplingInput
            args += " --downsampling %(downFactor)f"
            self.runJob(None, "xmipp_ctf_phase_flip", args % locals())
        
        
    def extractParticles(self, micId, micName, fnCTF, micrographToExtract):
        """ Extract particles from one micrograph """
        #If flip selected and exists CTF model use the flip output
#        if self.doFlip and self.fnCTF:
#            micrographToExtract = self._getTmpPath(micName +"_flipped.xmp")
                
        outputRoot = str(self._getExtraPath(micName))

        #fnPosFile = self.getConvertedInput('inputCoords').getMicrographCoordFile(micId)
        fnPosFile =  self._getExtraPath(micName + ".pos")

        # If it has coordinates extract the particles      
        particlesMd = 'particles@%s' % fnPosFile

        #if fnPosFile is not None and xmipp.existsBlockInMetaDataFile(particlesMd):
        if exists(fnPosFile):
            boxSize = self.boxSize.get()
            args = "-i %(micrographToExtract)s --pos %(particlesMd)s -o %(outputRoot)s --Xdim %(boxSize)d" % locals()
            if self.downsampleType.get() != self.SAME_AS_PICKING:
                args += " --downsampling " + str(self.samplingFinal/self.samplingInput)
            if self.doInvert:
                args += " --invert"
        
            self.runJob(None,"xmipp_micrograph_scissor", args)
            # Normalize 
            if self.doNormalize:
                self.runNormalize(outputRoot + '.stk', self.normType.get(), self.backRadius.get())          
                               
            if (self.downsampleType.get() == self.OTHER) or (fnCTF is not None):
                selfile = outputRoot + ".xmd"
                md = xmipp.MetaData(selfile)
                if self.downsampleType.get() == self.OTHER:
                    downsamplingFactor = self.samplingFinal/self.samplingInput
                    md.operate("Xcoor=Xcoor*%f" % downsamplingFactor)
                    md.operate("Ycoor=Ycoor*%f" % downsamplingFactor)
                # If micrograph has CTF info copy it
                if fnCTF is not None:
                    md.setValueCol(xmipp.MDL_CTF_MODEL, fnCTF)
                md.write(selfile)
                
    def runNormalize(self, stack, normType, bgRadius):
        program = "xmipp_transform_normalize"
        args = "-i %(stack)s "
        
        if bgRadius <= 0:
            particleSize = xmipp.MetaDataInfo(stack)[0]
            bgRadius = int(particleSize/2)
        
        if normType=="OldXmipp":
            args += "--method OldXmipp"
        elif normType=="NewXmipp":
            args += "--method NewXmipp --background circle %(bgRadius)d"
        else:
            args += "--method Ramp --background circle %(bgRadius)d"
        self.runJob(None, program, args % locals())
                    
    def getImgIdFromCoord(self, coordId):
        """ Get the image id from the related coordinate id. """
        '%s:%06d'
        parts = coordId.split(':')
        imgFn = self._getExtraPath(replaceBaseExt(parts[0], "stk")) 
        
        return '%06d@%s' %(int(parts[1]), imgFn)
        
    def createOutput(self):
        # Create the SetOfImages object on the database
        #imgSet = XmippSetOfParticles(self._getPath('images.xmd'))
                  
        #Create images.xmd metadata
        fnImages = self._getPath('images.xmd')
        imgsXmd = xmipp.MetaData() 
        posFiles = glob(self._getExtraPath('*.pos')) 
        for posFn in posFiles:
            xmdFn = self._getExtraPath(replaceBaseExt(posFn, "xmd"))
            md = xmipp.MetaData(xmdFn)
            mdPos = xmipp.MetaData('particles@%s' % posFn)
            mdPos.merge(md) 
            #imgSet.appendFromMd(mdPos)
            imgsXmd.unionAll(mdPos)
   
        imgsXmd.sort(xmipp.MDL_IMAGE)
        imgsXmd.write(fnImages)
        #imgSet.sort() # We need sort?
        #imgSet.write()

        # Run xmipp_image_sort_by_statistics to add zscore info to images.xmd
        args="-i %(fnImages)s --addToInput"
        if self.rejectionMethod==self.MAXZSCORE:
            maxZscore = self.maxZscore.get()
            args+=" --zcut "+str(maxZscore)
        elif self.rejectionMethod==self.PERCENTAGE:
            percentage = self.percentage.get()
            args+=" --percent "+str(percentage)

        self.runJob(None, "xmipp_image_sort_by_statistics", args % locals())
        md = xmipp.MetaData(fnImages) # Should have ZScore label after runJob
        md.sort(xmipp.MDL_ZSCORE)
        md.write(fnImages) 
        
        mdavgZscore = xmipp.MetaData()
        md.read(fnImages)
        mdavgZscore.aggregate(md, xmipp.AGGR_AVG, xmipp.MDL_MICROGRAPH, xmipp.MDL_ZSCORE, xmipp.MDL_ZSCORE)
        
        # Create output SetOfParticles
        imgSet = self._createSetOfParticles()
        imgSet.load()
        imgSet.copyInfo(self.inputMics)
        #imgSet.setHasCTF(self.fnCTF is not None)       
        imgSet.setHasCTF(self.ctfRelations.get() is not None)
                 
        if self.downsampleType == self.OTHER:
            imgSet.setSamplingRate(self.inputMics.getSamplingRate()*self.downFactor.get())
        imgSet.setCoordinates(self.inputCoords)
        readSetOfParticles(fnImages, imgSet, imgSet.hasCTF())
        # For each particle retrieve micId from SetOFCoordinates and set it on the CTFModel
        for img in imgSet:
            coord = self.inputCoords[img.getId()]
            ctfModel = img.getCTF()
            if ctfModel is not None:
                ctfModel.setId(coord.getMicId())
                img.setCTF(ctfModel)
                imgSet.update(img)
        imgSet.write()
        
        self._defineOutputs(outputParticles=imgSet)
        self._defineDataSource(self.inputCoords, imgSet)
        #TODO: pass CTF relation from input micrographs to imgSet
    
    def _summary(self):
        downsampleTypeText = {
                              self.ORIGINAL:'Original micrographs',
                              self.SAME_AS_PICKING:'Same as picking',
                              self.OTHER: 'Other downsampling factor'}
        summary = []
        if not hasattr(self, 'outputParticles'):
            summary.append("Output images not ready yet.") 
        else:
            summary.append("*Input coordinates* http://www.google.com : %s" % self.inputCoordinates.get().getName())
            summary.append("_Downsample type_: %s" % downsampleTypeText.get(self.downsampleType.get()))
            if self.downsampleType.get() == self.OTHER:
                summary.append("Downsampling factor: %d" % self.downFactor.get())
            summary.append("Particle size %d" % self.boxSize.get())
            summary.append("Particles extracted: %d" % (self.outputParticles.getSize()))
        return summary
    
    def _validate(self):
        validateMsgs = []
        # doFlip can only be True if CTF information is available on picked micrographs
        if self.doFlip.get() and not self.ctfRelations.hasValue():
            validateMsgs.append('Phase flipping cannot be performed unless CTF information is provided.')
        return validateMsgs
            
        

