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
import pyworkflow.em.metadata as md

from convert import (writeSetOfCoordinates, writeMicCoordinates,
                     relionToLocation, writeSetOfMicrographs, rowToParticle,
                     isVersion2)

from protocol_base import ProtRelionBase


# Micrograph type constants for particle extraction
SAME_AS_PICKING = 0
OTHER = 1


class ProtRelionExtractParticles(em.ProtExtractParticles, ProtRelionBase):
    """Protocol to extract particles from a set of coordinates"""
    _label = 'particles extraction'
    
    @classmethod
    def isDisabled(cls):
        return not isVersion2()

    def __init__(self, **kwargs):
        em.ProtExtractParticles.__init__(self, **kwargs)

    #--------------------------- DEFINE param functions ------------------------
    def _definePreprocessParams(self, form):
        form.addParam('boxSize', params.IntParam,
                      label='Particle box size (px)',
                      validators=[params.Positive],
                      help='This is size of the boxed particles (in pixels).')
    
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
        inputMics = self.getInputMicrographs()
        # self.doFlip is important because in _updateOutputPartSet is used
        self.doFlip = False
        self.streamIsOpen = self._isStreamOpen()

        if self.streamIsOpen:
            # If the input is in streaming, follow the base class policy
            # about inserting new steps and discovery new input/output
            
            # TODO: We have to change the way to process the micrographs. It
            # is now done one by one, rather than by group.
            
            em.ProtExtractParticles._insertAllSteps(self)
        else:
            # If not in streaming, then we will just insert a single step to
            # extract all micrographs at once since it is much faster
            deps = []
            micStarFn, partStarFn = self._getStarFiles()
            firstStepId = self._insertFunctionStep('convertInputStep',
                                                   inputMics.getObjId())

            args = self._getExtractArgs()
            self._insertFunctionStep('extractNoStreamParticlesStep', *args,
                                                 prerequisites=[firstStepId])
            self._insertFunctionStep('createOutputStep')

            # Disable streaming functions:
            self._insertFinalSteps = self._doNothing
            self._stepsCheck = self._doNothing

    def _insertInitialSteps(self):
        self._setupBasicProperties()
    
        return []
    
    def _doNothing(self, *args):
        pass # used to avoid some streaming functions

    # -------------------------- STEPS functions -------------------------------
    def convertInputStep(self, micsId):
        self._ih = em.ImageHandler() # used to convert micrographs
        # Match ctf information against the micrographs
        self.ctfDict = {}
        inputMics = self.getInputMicrographs()
        if self.ctfRelations.get() is not None:
            for ctf in self.ctfRelations.get():
                self.ctfDict[ctf.getMicrograph().getMicName()] = ctf.clone()
        micStar, _ = self._getStarFiles()
        writeSetOfMicrographs(inputMics, self._getPath(micStar),
                              alignType=em.ALIGN_NONE,
                              preprocessImageRow=self._preprocessMicrographRow)

        # We need to compute a scale factor for the coordinates if extracting
        # from other micrographs with a different pixel size
        micDict = {}
        for mic in inputMics:
            micDict[mic.getMicName()] = mic.getFileName()

        def _getCoordsStarFile(mic):
            micName = mic.getMicName()
            if not micName in micDict:
                return None
            micFn = micDict[micName]
            return self._getExtraPath(pwutils.replaceBaseExt(micFn,
                                                             'coords.star'))

        self.info("Using scale: %s" % self.getScaleFactor())
        writeSetOfCoordinates(self._getExtraPath(), self.getCoords(),
                              _getCoordsStarFile, scale=self.getScaleFactor())

    def _convertCoordinates(self, mic, coordList):
        micStar, _ = self._getStarFiles(mic)
        writeSetOfMicrographs([mic], self._getPath(micStar),
                              alignType=em.ALIGN_NONE,
                              preprocessImageRow=self._preprocessMicrographRow)

        writeMicCoordinates(mic, coordList, self._getMicPos(mic),
                            getPosFunc=self._getPos)
    
    def _getMicPos(self, mic):
        """ Return the corresponding .pos file for a given micrograph. """
        micBase = pwutils.removeBaseExt(mic.getFileName())
        return self._getExtraPath(micBase + ".coords.star")

    def extractNoStreamParticlesStep(self, *args):
        self._extractMicrograph(None, *args)
    
    def _extractMicrograph(self, mic, params):
        """ Extract particles from one micrograph, ignore if the .star
        with the coordinates is not present. """
        
        micsStar, partStar = self._getStarFiles(mic)
        print "params: ", params
        args = params % locals()
        self.runJob(self._getProgram('relion_preprocess'), args,
                    cwd=self.getWorkingDir())

    def createOutputStep(self):
        if self.streamIsOpen:
            pass
        else:
            inputMics = self.getInputMicrographs()
            inputCoords = self.getCoords()
            coordMics = inputCoords.getMicrographs()
    
            # Create output SetOfParticles
            partSet = self._createSetOfParticles()
            partSet.copyInfo(inputMics)
            # set coords from the input, will update later if needed
            partSet.setCoordinates(inputCoords)
    
    
            hasCTF = self._useCTF()
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
            # Check missing mics to report error once
            missingMics = set()
 
            for coord in inputCoords.iterItems(orderBy='_micId'):
                micId = coord.getMicId()
                mic = coordMics[micId] 
                if mic is None:
                    if not micId in missingMics:
                        print("Ignoring wrong micId: ", micId)
                        missingMics.add(micId)
                    micName = None
                else:
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
    
            if self._useCTF():
                self._defineSourceRelation(self.ctfRelations.get(), partSet)

    #--------------------------- INFO functions --------------------------------
    def _validate(self):
        errors = []

        if self.doNormalize and self.backDiameter > self.boxSize:
            errors.append("Background diameter for normalization should "
                          "be equal or less than the box size.")

        # We cannot check this if the protocol is in streaming.
        
        # self._setupCtfProperties() # setup self.micKey among others
        # if self._useCTF() and self.micKey is None:
        #     errors.append('Some problem occurs matching micrographs and '
        #                   'CTF.\n There were micrographs for which CTF '
        #                   'was not found either using micName or micId.\n')
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
    def _setupBasicProperties(self):
        # Set sampling rate (before and after doDownsample) and inputMics
        # according to micsSource type
        inputCoords = self.getCoords()
        mics = inputCoords.getMicrographs()
        self.samplingInput = inputCoords.getMicrographs().getSamplingRate()
        self.samplingMics = self.getInputMicrographs().getSamplingRate()
        self.samplingFactor = float(self.samplingMics / self.samplingInput)

        scale = self.getBoxScale()
        self.debug("Scale: %f" % scale)
        if self.notOne(scale):
            # If we need to scale the box, then we need to scale the coordinates
            getPos = lambda coord: (int(coord.getX() * scale),
                                    int(coord.getY() * scale))
        else:
            getPos = lambda coord: coord.getPosition()
        # Store the function to be used for scaling coordinates
        self._getPos = getPos

    def _getExtractArgs(self):
        # The following parameters are executing 'relion_preprocess' to
        # extract the particles of a given micrographs
        # The following is assumed:
        # - relion_preproces will be executed from the protocol workingDir
        # - the micrographs (or links) and coordinate files will be in 'extra'
        # - coordinate files have the 'coords.star' suffix
        params = ' --i %(micsStar)s'
        params += ' --coord_dir "."'
        params += ' --coord_suffix .coords.star'
        params += ' --part_star %(partStar)s'
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

        return [params]

    def readPartsFromMics(self, micList, outputParts):
        """ Read the particles extract for the given list of micrographs
        and update the outputParts set with new items.
        """
        p = em.Particle()
        for mic in micList:
            # We need to make this dict because there is no ID in the
            # coord.star file
            coordDict = {}
            for coord in self.coordDict[mic.getObjId()]:
                coordDict[self._getPos(coord)] = coord
        
            _, partStarFn = self._getStarFiles(mic)
        
            for row in md.iterRows(self._getPath(partStarFn)):
                pos = (row.getValue(md.RLN_IMAGE_COORD_X),
                       row.getValue(md.RLN_IMAGE_COORD_Y))
            
                coord = coordDict.get(pos, None)
                if coord is not None:
                    # scale the coordinates according to particles dimension.
                    coord.scale(self.getBoxScale())
                    p.copyObjId(coord)
                    idx, fn = relionToLocation(row.getValue(md.RLN_IMAGE_NAME))
                    p.setLocation(idx, self._getPath(fn[1:]))
                    p.setCoordinate(coord)
                    p.setMicId(mic.getObjId())
                    p.setCTF(mic.getCTF())
                    outputParts.append(p)
        
            # Release the list of coordinates for this micrograph since it
            # will not be longer needed
            del self.coordDict[mic.getObjId()]

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

    # def _setupCtfProperties(self):
    #     inputMics = self.getInputMicrographs()
    #     if self._useCTF():
    #         # Load CTF dictionary for all micrographs, all CTF should be present
    #         self.ctfDict = {}
    #
    #         for ctf in self.ctfRelations.get():
    #             ctfMic = ctf.getMicrograph()
    #             newCTF = ctf.clone()
    #             self.ctfDict[ctfMic.getMicName()] = newCTF
    #             self.ctfDict[ctfMic.getObjId()] = newCTF
    #
    #         if all(mic.getMicName() in self.ctfDict for mic in inputMics):
    #             self.micKey = lambda mic: mic.getMicName()
    #         elif all(mic.getObjId() in self.ctfDict for mic in inputMics):
    #             self.micKey = lambda mic: mic.getObjId()
    #         else:
    #             self.micKey = None # some problem matching CTF
            
    def getInputMicrographs(self):
        """ Return the micrographs associated to the SetOfCoordinates or
        Other micrographs. """
        if not self._micsOther():
            return self.inputCoordinates.get().getMicrographs()
        else:
            return self.inputMicrographs.get()

    def getCoords(self):
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
        coordsSampling = self.getCoords().getMicrographs().getSamplingRate()
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
        return int(self.getCoords().getBoxSize() * self.getBoxScale())

    def _getOutputImgMd(self):
        return self._getPath('images.xmd')

    def createParticles(self, item, row):
        particle = rowToParticle(row, readCtf=self._useCTF())
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
        if self._useCTF() and not img.hasCTF():
            img.setCTF(self.ctfDict[img.getMicName()])

        if self._useCTF():
            # add phaseShift to micrograph Row
            ctf = img.getCTF()
            if ctf is not None and ctf.getPhaseShift():
                imgRow.setValue(md.RLN_CTF_PHASESHIFT, ctf.getPhaseShift())

    def __getMicFile(self, mic, ext):
        """ Return a filename based on the micrograph.
        The filename will be located in the extra folder and with
        the given extension.
        """
        return self._getExtraPath(pwutils.replaceBaseExt(mic.getFileName(),ext))
    
    def _useCTF(self):
        return self.ctfRelations.hasValue()

    def _getMicStarFile(self, mic):
        return self.__getMicFile(mic, 'star')

    def _getMicStackFile(self, mic):
        return self.__getMicFile(mic, 'mrcs')

    def _getStarFiles(self, mic=None):
        # Actually extract
        # Since the program will be run in the run working dir
        # we don't need to use the workingDir prefix here
        if mic is None:
            micsStar = 'input_micrographs.star'
            partStar = os.path.join('extra', 'output_particles.star')
        else:
            micsStar = 'input_micrographs_%s.star' % mic.getObjId()
            partStar = os.path.join('extra',
                                    'output_particles_%s.star' % mic.getObjId())
        
        return micsStar, partStar

    def _isStreamOpen(self):
        if self._useCTF():
            ctfStreamOpen = self.ctfRelations.get().isStreamOpen()
        else:
            ctfStreamOpen = False
        
        return (self.getInputMicrographs().isStreamOpen() or
                ctfStreamOpen or self.getCoords().isStreamOpen())
