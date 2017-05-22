# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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

import pyworkflow.utils as pwutils
from pyworkflow.protocol.constants import STATUS_FINISHED
import pyworkflow.protocol.params as params
import pyworkflow.em as em

from convert import writeSetOfCoordinates, writeSetOfMicrographs, rowToParticle
from protocol_base import ProtRelionBase


# Rejection method constants
REJECT_NONE = 0
REJECT_MAXZSCORE = 1
REJECT_PERCENTAGE = 2

# Micrograph type constants for particle extraction
SAME_AS_PICKING = 0
OTHER = 1


class ProtRelionExtractParticles(em.ProtExtractParticles, ProtRelionBase):
    """Protocol to extract particles from a set of coordinates"""
    _label = 'particles extraction'
    
    def __init__(self, **kwargs):
        em.ProtExtractParticles.__init__(self, **kwargs)

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
                      help='By default the particles will be extracted '
                           'from the micrographs used in the picking '
                           'step ( _same as picking_ option ). \n'
                           'If you select _other_ option, you must provide '
                           'a different set of micrographs to extract from. \n'
                           '*Note*: In the _other_ case, ensure that provided '
                           'micrographs and coordinates are related '
                           'by micName or by micId. Difference in pixel size '
                           'will be handled automatically.')

        form.addParam('inputMicrographs', params.PointerParam,
                      pointerClass='SetOfMicrographs',
                      condition='downsampleType != %s' % SAME_AS_PICKING,
                      important=True, label='Input micrographs',
                      help='Select the SetOfMicrographs from which to extract.')

        form.addParam('ctfRelations', params.RelationParam, allowsNull=True,
                      relationName=em.RELATION_CTF,
                      attributeName='getInputMicrographs',
                      label='CTF estimation',
                      help='Choose some CTF estimation related to input '
                           'micrographs. \n CTF estimation is needed if you '
                           'want to do phase flipping or you want to '
                           'associate CTF information to the particles.')

        # downFactor should always be 1.0 or greater
        geOne = params.GE(1.0,
                          error='Value should be greater or equal than 1.0')

        form.addParam('boxSize', params.IntParam,
                      label='Particle box size (px)',
                      validators=[params.Positive],
                      help='This is size of the boxed particles (in pixels). '
                           'Note that if you use downsample option, the '
                           'particles are boxed out after downsampling. '
                           'Use the wizard to check boxSize changes after '
                           'downsampling or using a different pixel size. ')

        form.addParam('doRescale', params.BooleanParam, default=False,
                      label='Rescale particles?',
                      help='If set to Yes, particles will be re-scaled. '
                           'Note that the re-scaled size below will be in '
                           'the down-scaled images.')

        form.addParam('rescaledSize', params.IntParam,
                      validators=[params.Positive],
                      condition='doRescale',
                      label='Re-scaled size (px)',
                      help='Final size in pixels of the extracted particles. '
                           'The provided value should be an even number. ')

        form.addSection(label='Preprocess')

        form.addParam('doInvert', params.BooleanParam, default=None,
                      label='Invert contrast?',
                      help='Invert the contrast if your particles are black '
                           'over a white background.  Xmipp, Spider, Relion '
                           'and Eman require white particles over a black '
                           'background. Frealign (up to v9.07) requires black '
                           'particles over a white background')

        form.addParam('doNormalize', params.BooleanParam, default=True,
                      label='Normalize particles?',
                      help='If set to Yes, particles will be normalized in '
                           'the way RELION prefers it.')

        form.addParam('backDiameter', params.IntParam, default=-1,
                      condition='doNormalize',
                      label='Diameter background circle (px)',
                      help='Particles will be normalized to a mean value of '
                           'zero and a standard-deviation of one for all '
                           'pixels in the background area. The background area '
                           'is defined as all pixels outside a circle with '
                           'this given diameter in pixels (before rescaling). '
                           'When specifying a negative value, a default value '
                           'of 75% of the Particle box size will be used.')

        form.addParam('stddevWhiteDust', params.FloatParam, default=-1,
                      condition='doNormalize',
                      label='Stddev for white dust removal: ',
                      help='Remove very white pixels from the extracted '
                           'particles. Pixels values higher than this many '
                           'times the image stddev will be replaced with '
                           'values from a Gaussian distribution. \n'
                           'Use negative value to switch off dust removal.')

        form.addParam('stddevBlackDust', params.FloatParam, default=-1,
                      condition='doNormalize',
                      label='Stddev for black dust removal: ',
                      help='Remove very black pixels from the extracted '
                           'particles. Pixels values higher than this many '
                           'times the image stddev will be replaced with '
                           'values from a Gaussian distribution. \n'
                           'Use negative value to switch off dust removal.')

        form.addParallelSection(threads=0, mpi=4)
    
    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        firstStepId = self._insertFunctionStep('convertInputStep')

        deps = []

        inputMics = self.getInputMicrographs()

        # Actually extract
        # Since the program will be run in the run working dir
        # we don't need to use the workingDir prefix here
        micStarFile = 'input_micrographs.star'
        partStarFile = os.path.join('extra', 'output_particles.star')

        params = self.getExtractParams(micStarFile, partStarFile)
        deps.append(self._insertFunctionStep('extractParticlesStep',
                                             inputMics.getObjId(), params,
                                             prerequisites=[firstStepId]))

        self._insertFunctionStep('createOutputStep')

    # -------------------------- STEPS functions -------------------------------

    def getExtractParams(self, inputMicStarFile, outputPartStarFile):
        # The following parameters are executing 'relion_preprocess' to
        # extract the particles of a given micrographs
        # The following is assumed:
        # - relion_preproces will be executed from the protocol workingDir
        # - the micrographs (or links) and coordinate files will be in 'extra'
        # - coordinate files have the 'coords.star' suffix
        params = ' --i %s' % inputMicStarFile
        params += ' --coord_dir "."'
        params += ' --coord_suffix .coords.star'
        params += ' --part_star %s' % outputPartStarFile
        params += ' --part_dir "." --extract '
        params += ' --extract_size %d' % self.boxSize

        if self.backDiameter <= 0:
            diameter = self.boxSize.get() * 0.75 / self._getDownFactor()
        else:
            diameter = self.backDiameter.get()

        params += ' --bg_radius %d' % int(diameter/2)

        if self.doInvert:
            params += ' --invert_contrast'

        if self.doNormalize:
            params += ' --norm'

        if self._doDownsample():
            params += ' --scale %d' % self.rescaledSize

        if self.stddevWhiteDust > 0:
            params += ' --white_dust %0.3f' % self.stddevWhiteDust

        if self.stddevBlackDust > 0:
            params += ' --black_dust %0.3f' % self.stddevBlackDust

        return params

    def convertInputStep(self):
        self._ih = em.ImageHandler() # used to convert micrographs
        # Match ctf information against the micrographs
        self.ctfDict = {}
        if self.ctfRelations.get() is not None:
            for ctf in self.ctfRelations.get():
                self.ctfDict[ctf.getMicrograph().getMicName()] = ctf.clone()

        micStar = self._getPath('input_micrographs.star')

        writeSetOfMicrographs(self.getInputMicrographs(), micStar,
                              alignType=em.ALIGN_NONE,
                              preprocessImageRow=self._preprocessMicrographRow)

        # We need to compute a scale factor for the coordinates if extracting
        # from other micrographs with a different pixel size
        micDict = {}
        for mic in self.getInputMicrographs():
            micDict[mic.getMicName()] = mic.getFileName()

        def _getCoordsStarFile(mic):
            micName = mic.getMicName()
            if not micName in micDict:
                return None
            micFn = micDict[micName]
            return self._getExtraPath(pwutils.replaceBaseExt(micFn,
                                                             'coords.star'))

        self.info("Using scale: %s" % self.getScaleFactor())
        writeSetOfCoordinates(self._getExtraPath(), self.getInputCoords(),
                              _getCoordsStarFile, scale=self.getScaleFactor())

    def extractParticlesStep(self, inputId, params):
        """ Extract particles from one micrograph, ignore if the .star
        with the coordinates is not present. """
        self.runJob(self._getProgram('relion_preprocess'), params,
                    cwd=self.getWorkingDir())

    def createOutputStep(self):
        inputMics = self.getInputMicrographs()
        inputCoords = self.getInputCoords()
        coordMics = inputCoords.getMicrographs()

        # Create output SetOfParticles
        partSet = self._createSetOfParticles()
        partSet.copyInfo(inputMics)
        # set coords from the input, will update later if needed
        partSet.setCoordinates(inputCoords)


        hasCTF = self.ctfRelations.hasValue()
        partSet.setSamplingRate(self._getNewSampling())
        partSet.setHasCTF(hasCTF)

        ctfDict = {}

        if hasCTF:
            # Load CTF dictionary for all micrographs, all CTF should be present
            for ctf in self.ctfRelations.get():
                ctfMic = ctf.getMicrograph()
                newCTF = ctf.clone()
                ctfDict[ctfMic.getMicName()] = newCTF

        # Keep a dictionary between the micName and its corresponding
        # particles stack file and CTF model object
        micDict = {}
        for mic in inputMics:
            stackFile = self._getMicStackFile(mic)
            ctfModel = ctfDict[mic.getMicName()] if hasCTF else None
            micDict[mic.getMicName()] = (stackFile, ctfModel)

        scaleFactor = self.getBoxScale()
        doScale = self.notOne(scaleFactor)

        # Track last micName to know when to load particles from another
        # micrograph stack
        lastMicName = None
        count = 0 # Counter for the particles of a given micrograph

        for coord in inputCoords.iterItems(orderBy='_micId'):
            micName = coordMics[coord.getMicId()].getMicName()
            # If Micrograph Source is "other" and extract from a subset
            # of micrographs, micName key should be checked if it exists.
            if micName in micDict.keys():
                # Load the new particles mrcs file and reset counter
                if micName != lastMicName:
                    stackFile, ctfModel = micDict[micName]
                    count = 1
                    lastMicName = micName
    
                p = em.Particle()
                p.setLocation(count, stackFile)
                if hasCTF:
                    p.setCTF(ctfModel)
                p.setCoordinate(coord)
                # Copy objId and micId from the coordinate
                p.copyObjId(coord)
                p.setMicId(coord.getMicId())
    
                if doScale:
                    p.scaleCoordinate(scaleFactor)
    
                partSet.append(p)
                count += 1

        self._defineOutputs(outputParticles=partSet)
        self._defineSourceRelation(self.inputCoordinates, partSet)

        if self.ctfRelations.hasValue():
            self._defineSourceRelation(self.ctfRelations.get(), partSet)

    #--------------------------- INFO functions --------------------------------
    def _validate(self):
        errors = []

        if self.doNormalize and self.backDiameter > self.boxSize:
            errors.append("Background diameter for normalization should "
                          "be equal or less than the box size.")

        self._setupCtfProperties() # setup self.micKey among others
        if self.ctfRelations.hasValue() and self.micKey is None:
            errors.append('Some problem occurs matching micrographs and CTF.\n'
                          'There were micrographs for which CTF was not found '
                          'either using micName or micId.\n')
        return errors
    
    def _citations(self):
        return ['Scheres2012b']
        
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
            msg = ("A total of %d particles of size %d were extracted"
                   % (self.getOutput().getSize(), self.boxSize))
            
            if self._micsOther():
                msg += (" from another set of micrographs: %s"
                        % self.getObjectTag('inputMicrographs'))

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

    def _getDownFactor(self):
        if self.doRescale:
            return self.boxSize.get() / self.rescaledSize.get()
        return 1

    def _doDownsample(self):
        return self.doRescale and self.rescaledSize != self.boxSize

    def notOne(self, value):
        return abs(value - 1) > 0.0001

    def _getNewSampling(self):
        newSampling = self.getInputMicrographs().getSamplingRate()

        if self._doDownsample():
            # Set new sampling, it should be the input sampling of the used
            # micrographs multiplied by the downFactor
            newSampling *= self._getDownFactor()

        return newSampling

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

    def getInputCoords(self):
        return self.inputCoordinates.get()

    def getOutput(self):
        if (self.hasAttribute('outputParticles') and
            self.outputParticles.hasValue()):
            return self.outputParticles
        else:
            return None

    def getCoordSampling(self):
        return

    def getScaleFactor(self):
        """ Returns the scaling factor that needs to be applied to the input
        coordinates to adapt for the input micrographs.
        """
        coordsSampling = self.getInputCoords().getMicrographs().getSamplingRate()
        micsSampling = self.getInputMicrographs().getSamplingRate()
        return coordsSampling / micsSampling

    def getBoxScale(self):
        """ Computing the sampling factor between input and output.
        We should take into account the differences in sampling rate between
        micrographs used for picking and the ones used for extraction.
        The downsampling factor could also affect the resulting scale.
        """
        f = self.getScaleFactor()
        return f / self._getDownFactor() if self._doDownsample() else f

    def getBoxSize(self):
        # This function is needed by the wizard
        return int(self.getInputCoords().getBoxSize() * self.getBoxScale())

    def _getOutputImgMd(self):
        return self._getPath('images.xmd')

    def createParticles(self, item, row):
        particle = rowToParticle(row, readCtf=self.ctfRelations.hasValue())
        coord = particle.getCoordinate()
        item.setY(coord.getY())
        item.setX(coord.getX())
        particle.setCoordinate(item)

        item._appendItem = False

    def _preprocessMicrographRow(self, img, imgRow):
        # Temporarly convert the few micrographs to tmp and make sure
        # they are in 'mrc' format
        # Get basename and replace extension by 'mrc'
        newName = pwutils.replaceBaseExt(img.getFileName(), 'mrc')
        newPath = self._getExtraPath(newName)

        # If the micrographs are in 'mrc' format just create a link
        # if not, convert to 'mrc'
        if img.getFileName().endswith('mrc'):
            pwutils.createLink(img.getFileName(), newPath)
        else:
            self._ih.convert(img, newPath)
        # The command will be launched from the working dir
        # so, let's make the micrograph path relative to that
        img.setFileName(os.path.join('extra', newName))
        if self.ctfRelations.get() is not None:
            img.setCTF(self.ctfDict[img.getMicName()])

    def __getMicFile(self, mic, ext):
        """ Return a filename based on the micrograph.
        The filename will be located in the extra folder and with
        the given extension.
        """
        return self._getExtraPath(pwutils.replaceBaseExt(mic.getFileName(),
                                                         ext))
    def _getMicStarFile(self, mic):
        return self.__getMicFile(mic, 'star')

    def _getMicStackFile(self, mic):
        return self.__getMicFile(mic, 'mrcs')

