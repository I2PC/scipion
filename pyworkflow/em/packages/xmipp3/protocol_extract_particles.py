# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Laura del Cano (ldelcano@cnb.csic.es)
# *              Adrian Quintana (aquintana@cnb.csic.es)
# *              Javier Vargas (jvargas@cnb.csic.es)
# *              Grigory Sharov (sharov@igbmc.fr)
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

from glob import glob
from itertools import izip
from os.path import exists, basename

import pyworkflow.em.metadata as md
import pyworkflow.utils as pwutils
from pyworkflow.em.packages.xmipp3.constants import SAME_AS_PICKING, OTHER
from pyworkflow.protocol.constants import STEPS_PARALLEL, LEVEL_ADVANCED, STATUS_FINISHED
import pyworkflow.protocol.params as params
from pyworkflow.em.protocol import ProtExtractParticles
from pyworkflow.em.data import SetOfParticles
from pyworkflow.em.constants import RELATION_CTF

from convert import writeSetOfCoordinates, readSetOfParticles, micrographToCTFParam
from xmipp3 import XmippProtocol

# Rejection method constants
REJECT_NONE = 0
REJECT_MAXZSCORE = 1
REJECT_PERCENTAGE = 2


class XmippProtExtractParticles(ProtExtractParticles, XmippProtocol):
    """Protocol to extract particles from a set of coordinates"""
    _label = 'extract particles'
    
    def __init__(self, **args):
        ProtExtractParticles.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL

    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('inputCoordinates', params.PointerParam,
                      pointerClass='SetOfCoordinates',
                      important=True,
                      label="Input coordinates",
                      help='Select the SetOfCoordinates ')

        # The name for the followig param is because historical reasons
        # now it should be named better 'micsSource' rather than
        # 'downsampleType', but this could make inconsistent previous executions
        # of this protocols, we will keep the name
        form.addParam('downsampleType', params.EnumParam,
                      choices=['same as picking', 'other'],
                      default=0, important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='Micrographs source',
                      help='By default the particles will be extracted \n'
                           'from the micrographs used in the picking \n'
                           'step ( _same as picking_ option ). \n'
                           'If you select _other_ option, you must provide \n'
                           'a different set of micrographs to extract from.\n'
                           '*Note*: In the _other_ case, ensure that provided \n'
                           'micrographs and coordinates are related \n'
                           'by micName or by micId. Difference in pixel size will \n'
                           'be handled automatically.')

        form.addParam('inputMicrographs', params.PointerParam,
                      pointerClass='SetOfMicrographs',
                      condition='downsampleType != %s' % SAME_AS_PICKING,
                      important=True, label='Input micrographs',
                      help='Select the SetOfMicrographs from which to extract.')

        form.addParam('ctfRelations', params.RelationParam, allowsNull=True,
                      relationName=RELATION_CTF,
                      attributeName='getInputMicrographs',
                      label='CTF estimation',
                      help='Choose some CTF estimation related to input micrographs. \n'
                           'CTF estimation is needed if you want to do phase flipping or \n'
                           'you want to associate CTF information to the particles.')

        # downFactor should always be 1.0 or greater
        geOne = params.GE(1.0,
                          error='Value should be greater or equal than 1.0')

        form.addParam('downFactor', params.FloatParam, default=1.0,
                      validators=[geOne],
                      label='Downsampling factor',
                      help='Select a value greater than 1.0 to reduce the size '
                           'of micrographs before extracting the particles. '
                           'If 1.0 is used, no downsample is applied. '
                           'Non-integer downsample factors are possible. ')

        form.addParam('boxSize', params.IntParam,
                      label='Particle box size (px)',
                      validators=[params.Positive],
                      help='This is size of the boxed particles (in pixels). '
                           'Note that if you use downsample option, the '
                           'particles are boxed out after downsampling. '
                           'Use the wizard to check boxSize changes after '
                           'downsampling or using a different pixel size. ')

        form.addParam('doSort', params.BooleanParam, default=True,
                      label='Perform sort by statistics',
                      help='Perform sort by statistics to compute a zscore '
                           'for each particle.')

        form.addParam('rejectionMethod', params.EnumParam,
                      choices=['None','MaxZscore', 'Percentage'],
                      default=REJECT_NONE, condition='doSort',
                      display=params.EnumParam.DISPLAY_COMBO,
                      label='Automatic particle rejection',
                      expertLevel=LEVEL_ADVANCED,
                      help='How to automatically reject particles. It can be: '
                           '*None*: no rejection '
                           '*MaxZscore*: reject a particle if its Zscore is '
                           '             larger than this value), '
                           '*Percentage*: reject a given percentage in each '
                           '             one of the screening criteria).')

        form.addParam('maxZscore', params.IntParam, default=3,
                      expertLevel=LEVEL_ADVANCED,
                      condition='doSort and rejectionMethod==%d' % REJECT_MAXZSCORE,
                      label='Maximum Zscore',
                      help='Maximum Zscore above which particles are rejected.')
        
        form.addParam('percentage', params.IntParam, default=5,
                      expertLevel=LEVEL_ADVANCED,
                      condition='rejectionMethod==%d' % REJECT_PERCENTAGE,
                      label='Percentage (%)',
                      help='Percentage of particles to reject')
        
        form.addSection(label='Preprocess')

        form.addParam('doRemoveDust', params.BooleanParam, default=True,
                      label='Dust removal (Recommended)', important=True,
                      help='Sets pixels with unusually large values to random '
                           'values from a Gaussian with zero-mean and '
                           'unity-standard deviation.')

        form.addParam('thresholdDust', params.FloatParam, default=3.5,
                      condition='doRemoveDust', expertLevel=LEVEL_ADVANCED,
                      label='Threshold for dust removal',
                      help='Pixels with a signal higher or lower than this value times the standard '
                           'deviation of the image will be affected. For cryo, 3.5 is a good value.'
                           'For high-contrast negative stain, the signal itself may be affected so '
                           'that a higher value may be preferable.')

        form.addParam('doInvert', params.BooleanParam, default=None,
                      label='Invert contrast', 
                      help='Invert the contrast if your particles are black over a white background.\n'
                           'Xmipp, Spider, Relion and Eman require white particles over a black background\n'
                           'Frealign (up to v9.07) requires black particles over a white background')
        
        form.addParam('doFlip', params.BooleanParam, default=None,
                      label='Phase flipping',
                      help='Use the information from the CTF to compensate for '
                           'phase reversals.\n'
                           'Phase flip is recommended in Xmipp or Eman\n'
                           '(even Wiener filtering and bandpass filter are '
                           'recommended for obtaining better 2D classes)\n'
                           'Otherwise (Frealign, Relion, Spider, ...), '
                           'phase flip is not recommended.')

        form.addParam('doNormalize', params.BooleanParam, default=True,
                      label='Normalize (Recommended)', 
                      help='It subtract a ramp in the gray values and normalizes so that in the '
                           'background there is 0 mean and standard deviation 1.')
        form.addParam('normType', params.EnumParam,
                      choices=['OldXmipp','NewXmipp','Ramp'], default=2,
                      condition='doNormalize', expertLevel=LEVEL_ADVANCED,
                      display=params.EnumParam.DISPLAY_COMBO,
                      label='Normalization type', 
                      help='OldXmipp (mean(Image)=0, stddev(Image)=1).  \n  '
                           'NewXmipp (mean(background)=0, stddev(background)=1)  \n  '
                           'Ramp (subtract background+NewXmipp).  \n  ')
        form.addParam('backRadius', params.IntParam, default=-1,
                      condition='doNormalize',
                      label='Background radius', expertLevel=LEVEL_ADVANCED,
                      help='Pixels outside this circle are assumed to be noise and their stddev '
                           'is set to 1. Radius for background circle definition (in pix.). '
                           'If this value is 0, then half the box size is used.')

        form.addParallelSection(threads=4, mpi=1)

    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        """for each micrograph insert the steps to preprocess it"""
        self._setupBasicProperties()
        # Write pos files for each micrograph
        firstStepId = self._insertFunctionStep('writePosFilesStep')
        
        # For each micrograph insert the steps, run in parallel
        deps = []
        
        for mic in self.inputMics:
            localDeps = [firstStepId]
            fnLast = mic.getFileName()
            baseMicName = pwutils.removeBaseExt(mic.getFileName())

            if self.ctfRelations.hasValue():
                micKey = self.micKey(mic)
                mic.setCTF(self.ctfDict[micKey])

            def getMicTmp(suffix):
                return self._getTmpPath(baseMicName + suffix)

            # Create a list with micrographs operations (programs in xmipp) and
            # the required command line parameters (except input/ouput files)
            micOps = []

            # Check if it is required to downsample your micrographs
            downFactor = self.downFactor.get()

            if self.notOne(downFactor):
                fnDownsampled = getMicTmp("_downsampled.xmp")
                args = "-i %s -o %s --step %f --method fourier"
                micOps.append(('xmipp_transform_downsample',
                               args % (fnLast, fnDownsampled, downFactor)))
                fnLast = fnDownsampled

            if self.doRemoveDust:
                fnNoDust = getMicTmp("_noDust.xmp")
                args = " -i %s -o %s --bad_pixels outliers %f"
                micOps.append(('xmipp_transform_filter',
                               args % (fnLast, fnNoDust, self.thresholdDust)))
                fnLast = fnNoDust

            if self.ctfRelations.hasValue():
                # If the micrograph doesn't come from Xmipp, we need to write
                # a Xmipp ctfparam file to perform the phase flip on the micrograph
                fnCTF = self._getTmpPath("%s.ctfParam" % baseMicName)
                # FIXME: Weird return of fnCTF, could be an internal xmipp one
                fnCTF = micrographToCTFParam(mic, fnCTF)
                # Insert step to flip micrograph
                if self.doFlip:
                    fnFlipped = getMicTmp('_flipped.xmp')
                    args = " -i %s -o %s --ctf %s --sampling %f"
                    micOps.append(('xmipp_ctf_phase_flip',
                                   args % (fnLast, fnFlipped, fnCTF,
                                           self._getNewSampling())))
                    fnLast = fnFlipped
            else:
                fnCTF = None

            # Actually extract
            deps.append(self._insertFunctionStep('extractParticlesStep',
                                                 mic.getObjId(), baseMicName,
                                                 fnCTF, fnLast, micOps,
                                                 self.doInvert.get(),
                                                 self._getNormalizeArgs(),
                                                 prerequisites=localDeps))

        # Insert step to create output objects
        metaDeps = self._insertFunctionStep('createMetadataImageStep',
                                            prerequisites=deps)
        
        if self.doSort:
            screenDep = self._insertFunctionStep('screenParticlesStep',
                                                 prerequisites=[metaDeps])
            finalDeps = [screenDep]
        else:
            finalDeps = [metaDeps]
        
        self._insertFunctionStep('createOutputStep', prerequisites=finalDeps)

    #--------------------------- STEPS functions -------------------------------
    def writePosFilesStep(self):
        """ Write the pos file for each micrograph on metadata format. """
        writeSetOfCoordinates(self._getExtraPath(), self.inputCoords,
                              scale=self.getBoxScale())
        # We need to find the mapping (either by micName (without ext) or micId)
        # between the micrographs in the SetOfCoordinates and
        # the Other micrographs if necessary
        if self._micsOther():
            micDict = {}
            coordMics = self.inputCoords.getMicrographs()
            for mic in coordMics:
                micBase = pwutils.removeBaseExt(mic.getFileName())
                micPos = self._getExtraPath(micBase + ".pos")
                micDict[pwutils.removeBaseExt(mic.getMicName())] = micPos
                micDict[mic.getObjId()] = micPos
                
            if any(pwutils.removeBaseExt(mic.getMicName()) in micDict
                   for mic in self.inputMics):
                micKey = lambda mic: pwutils.removeBaseExt(mic.getMicName())
            elif any(mic.getObjId() in micDict for mic in self.inputMics):
                self.warning('Could not match input micrographs and coordinates'
                             ' micrographs by micName, using micId.')
                micKey = lambda mic: mic.getObjId()
            else:
                raise Exception('Could not match input micrographs and '
                                'coordinates neither by micName or micId.')
            
            for mic in self.inputMics: # micrograph from input (other)
                mk = micKey(mic)
                if mk in micDict:
                    micPosCoord = micDict[mk]
                    if exists(micPosCoord):
                        micBase = pwutils.removeBaseExt(mic.getFileName())
                        micPos = self._getExtraPath(micBase + ".pos")
                        if micPos != micPosCoord:
                            self.info('Moving %s -> %s' % (micPosCoord, micPos))
                            pwutils.moveFile(micPosCoord, micPos)
                        
    def extractParticlesStep(self, micId, baseMicName, fnCTF,
                             micrographToExtract, micOps,
                             doInvert, normalizeArgs):
        """ Extract particles from one micrograph """
        outputRoot = str(self._getExtraPath(baseMicName))
        fnPosFile = self._getExtraPath(baseMicName + ".pos")

        # If it has coordinates extract the particles      
        particlesMd = 'particles@%s' % fnPosFile
        boxSize = self.boxSize.get()
        boxScale = self.getBoxScale()
        print "boxScale: ", boxScale

        if exists(fnPosFile):
            # Apply first all operations required for the micrograph
            for program, args in micOps:
                self.runJob(program, args)

            args = " -i %s --pos %s" % (micrographToExtract, particlesMd)
            args += " -o %s --Xdim %d" % (outputRoot, boxSize)

            if doInvert:
                args += " --invert"

            if fnCTF:
                args += " --ctfparam " + fnCTF

            self.runJob("xmipp_micrograph_scissor", args)

            # Normalize
            if normalizeArgs:
                self.runJob('xmipp_transform_normalize',
                            '-i %s.stk %s' % (outputRoot, normalizeArgs))
        else:
            self.warning("The micrograph %s hasn't coordinate file! "
                         % baseMicName)
            self.warning("Maybe you picked over a subset of micrographs")

        # Let's clean the temporary mrc micrographs
        if not pwutils.envVarOn("SCIPION_DEBUG_NOCLEAN"):
            pwutils.cleanPattern(self._getTmpPath(baseMicName) + '*')

    def _getNormalizeArgs(self):
        if not self.doNormalize:
            return ''

        normType = self.getEnumText("normType")
        args = "--method %s " % normType

        if normType != "OldXmipp":
            bgRadius = self.backRadius.get()
            if bgRadius <= 0:
                bgRadius = int(self.boxSize.get() / 2)
            args += " --background circle %d" % bgRadius

        return args

    def createMetadataImageStep(self):
        #Create images.xmd metadata
        fnImages = self._getOutputImgMd()
        imgsXmd = md.MetaData()
        posFiles = glob(self._getExtraPath('*.pos')) 
        for posFn in posFiles:
            xmdFn = self._getExtraPath(pwutils.replaceBaseExt(posFn, "xmd"))
            if exists(xmdFn):
                mdFn = md.MetaData(xmdFn)
                mdPos = md.MetaData('particles@%s' % posFn)
                mdPos.merge(mdFn) 
                imgsXmd.unionAll(mdPos)
            else:
                self.warning("The coord file %s wasn't used for extraction! "
                             % basename(posFn))
                self.warning("Maybe you are extracting over a subset of "
                             "micrographs")
        imgsXmd.write(fnImages)

    def screenParticlesStep(self):
        # If selected run xmipp_image_sort_by_statistics
        # to add zscore info to images.xmd
        args = " -i %s --addToInput" % self._getOutputImgMd()
        if self.rejectionMethod == REJECT_MAXZSCORE:
            args += " --zcut %s" % self.maxZscore
        elif self.rejectionMethod == REJECT_PERCENTAGE:
            args += " --percent %s" % self.percentage
        
        self.runJob("xmipp_image_sort_by_statistics", args)
    
    def createOutputStep(self):
        # Create the SetOfImages object on the database
        fnImages = self._getOutputImgMd()
        # Create output SetOfParticles
        imgSet = self._createSetOfParticles()
        imgSet.copyInfo(self.inputMics)
        # set coords from the input, will update later if needed
        imgSet.setCoordinates(self.inputCoords)

        if self.doFlip:
            imgSet.setIsPhaseFlipped(not self.inputMics.isPhaseFlipped())

        boxScale = self.getBoxScale()
        doScale = self.notOne(boxScale)

        imgSet.setSamplingRate(self._getNewSampling())

        # Create a temporary set to read from the metadata file
        # and later create the good one with the coordinates 
        # properly set. We need this because the .update is not
        # working in the mapper when new attributes are added.
        imgSet.setHasCTF(self.ctfRelations.hasValue())
        readSetOfParticles(fnImages, imgSet)
        
        # Since the coordinates are properly get in readSetOfParticles (
        # scaled if
        # necessary), is not necessary iterate again over the whole
        # setOfParticles.
        self._storeMethodsInfo(fnImages)
        self._defineOutputs(outputParticles=imgSet)
        self._defineSourceRelation(self.inputCoordinates, imgSet)
        if self.ctfRelations.hasValue():
            self._defineSourceRelation(self.ctfRelations.get(), imgSet)

    #--------------------------- INFO functions --------------------------------
    def _validate(self):
        errors = []
        # doFlip can only be selected if CTF information
        # is available on picked micrographs
        if self.doFlip and not self.ctfRelations.hasValue():
            errors.append('Phase flipping cannot be performed unless '
                          'CTF information is provided.')

        if self.doNormalize:
            if self.backRadius > int(self.boxSize.get() / 2):
                errors.append("Background radius for normalization should be "
                              "equal or less than half of the box size.")

        self._setupCtfProperties() # setup self.micKey among others
        if self.ctfRelations.hasValue() and self.micKey is None:
            errors.append('Some problem occurs matching micrographs and CTF.\n'
                          'There were micrographs for which CTF was not found\n'
                          'either using micName or micId.\n')
        return errors
    
    def _citations(self):
        return ['Vargas2013b']
        
    def _summary(self):
        summary = []
        summary.append("Micrographs source: %s"
                       % self.getEnumText("downsampleType"))
        summary.append("Particle box size: %d" % self.boxSize)
        
        if not hasattr(self, 'outputParticles'):
            summary.append("Output images not ready yet.") 
        else:
            summary.append("Particles extracted: %d" %
                           self.outputParticles.getSize())
            
        return summary
    
    def _methods(self):
        methodsMsgs = []

        if self.getStatus() == STATUS_FINISHED:
            msg = "A total of %d particles of size %d were extracted" % (self.getOutput().getSize(),
                                                                         self.boxSize)
            if self._micsOther():
                msg += " from another set of micrographs: %s" % self.getObjectTag('inputMicrographs')

            msg += " using coordinates %s" % self.getObjectTag('inputCoordinates')
            msg += self.methodsVar.get('')
            methodsMsgs.append(msg)

            if self.doRemoveDust:
                methodsMsgs.append("Removed dust over a threshold of %s."
                                   % self.thresholdDust)
            if self.doInvert:
                methodsMsgs.append("Inverted contrast on images.")
            if self._doDownsample():
                methodsMsgs.append("Particles downsampled by a factor of %0.2f."
                                   % self.downFactor)
            if self.doNormalize:
                methodsMsgs.append("Normalization: %s."
                                   % self.getEnumText('normType'))

        return methodsMsgs

    #--------------------------- UTILS functions -------------------------------
    def _micsOther(self):
        """ Return True if other micrographs are used for extract. """
        return self.downsampleType == OTHER

    def _doDownsample(self):
        return self.downFactor > 1.0

    def notOne(self, value):
        return abs(value - 1) > 0.0001

    def _getNewSampling(self):
        newSampling = self.samplingMics

        if self._doDownsample():
            # Set new sampling, it should be the input sampling of the used
            # micrographs multiplied by the downFactor
            newSampling *= self.downFactor.get()

        return newSampling

    def _setupBasicProperties(self):
        # Set sampling rate (before and after doDownsample) and inputMics
        # according to micsSource type
        self.inputCoords = self.getCoords()
        self.inputMics = self.getInputMicrographs()
        self.samplingInput = self.inputCoords.getMicrographs().getSamplingRate()
        self.samplingMics = self.inputMics.getSamplingRate()
        self.samplingFactor = float(self.samplingMics / self.samplingInput)

        self._setupCtfProperties()
        
    def _setupCtfProperties(self):
        inputMics = self.getInputMicrographs()
        if self.ctfRelations.hasValue():
            # Load CTF dictionary for all micrographs, all CTF should be present
            self.ctfDict = {}
            
            for ctf in self.ctfRelations.get():
                ctfMic = ctf.getMicrograph()
                newCTF = ctf.clone()
                self.ctfDict[ctfMic.getMicName()] = newCTF
                self.ctfDict[ctfMic.getObjId()] = newCTF
            
            if all(mic.getMicName() in self.ctfDict for mic in inputMics):
                self.micKey = lambda mic: mic.getMicName()
            elif all(mic.getObjId() in self.ctfDict for mic in inputMics):
                self.micKey = lambda mic: mic.getObjId()
            else:
                self.micKey = None # some problem matching CTF
            
    def getInputMicrographs(self):
        """ Return the micrographs associated to the SetOfCoordinates or
        Other micrographs. """
        if not self._micsOther():
            return self.inputCoordinates.get().getMicrographs()
        else:
            return self.inputMicrographs.get()

    def _storeMethodsInfo(self, fnImages):
        """ Store some information when the protocol finishes. """
        mdImgs = md.MetaData(fnImages)
        total = mdImgs.size() 
        mdImgs.removeDisabled()
        zScoreMax = mdImgs.getValue(md.MDL_ZSCORE, mdImgs.lastObject())
        numEnabled = mdImgs.size()
        numRejected = total - numEnabled
        msg = ""

        if self.doSort:
            if self.rejectionMethod != REJECT_NONE:
                msg = " %d of them were rejected with Zscore greater than %.2f." % (numRejected, zScoreMax)
        if self.doFlip:
            msg += "\nPhase flipping was performed."

        self.methodsVar.set(msg)

    def getCoords(self):
        return self.inputCoordinates.get()

    def getOutput(self):
        if (self.hasAttribute('outputParticles') and
            self.outputParticles.hasValue()):
            return self.outputParticles
        else:
            return None

    def getCoordSampling(self):
        return self.getCoords().getMicrographs().getSamplingRate()

    def getMicSampling(self):
        return self.getInputMicrographs().getSamplingRate()

    def getBoxScale(self):
        """ Computing the sampling factor between input and output.
        We should take into account the differences in sampling rate between
        micrographs used for picking and the ones used for extraction.
        The downsampling factor could also affect the resulting scale.
        """
        samplingPicking = self.getCoordSampling()
        samplingExtract = self.getMicSampling()
        f = samplingPicking / samplingExtract
        return f / self.downFactor.get() if self._doDownsample() else f

    def getBoxSize(self):
        # This function is needed by the wizard
        return int(self.getCoords().getBoxSize() * self.getBoxScale())

    def _getOutputImgMd(self):
        return self._getPath('images.xmd')
