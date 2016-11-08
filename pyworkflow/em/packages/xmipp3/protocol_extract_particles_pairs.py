# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Laura del Cano (ldelcano@cnb.csic.es)
# *              Adrian Quintana (aquintana@cnb.csic.es)
# *              Javier Vargas (jvargas@cnb.csic.es)
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

from itertools import izip
from os.path import exists

import pyworkflow.em.metadata as md
import pyworkflow.utils as pwutils
from pyworkflow.em.packages.xmipp3.protocol_extract_particles import XmippProtExtractParticles
from pyworkflow.em.packages.xmipp3.constants import SAME_AS_PICKING, OTHER
from pyworkflow.protocol.constants import STEPS_PARALLEL, LEVEL_ADVANCED, STATUS_FINISHED
from pyworkflow.protocol.params import (PointerParam, EnumParam, FloatParam, IntParam,
                                        BooleanParam, Positive, GE)
from convert import writeSetOfCoordinates, readSetOfParticles
from pyworkflow.em.data_tiltpairs import ParticlesTiltPair, TiltPair
from pyworkflow.em.data import SetOfMicrographs, SetOfParticles



class XmippProtExtractParticlesPairs(XmippProtExtractParticles):
    """Protocol to extract particles from a set of tilted pairs coordinates"""
    _label = 'extract particle pairs'

    def __init__(self, **args):
        XmippProtExtractParticles.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL

    # --------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputCoordinatesTiltedPairs', PointerParam,
                      important=True, label="Coordinates tilted pairs",
                      pointerClass='CoordinatesTiltPair',
                      help='Select the CoordinatesTiltPairs')
        form.addParam('downsampleType', EnumParam, choices=['same as picking', 'other'],
                      default=0, important=True,
                      label='Micrographs source', display=EnumParam.DISPLAY_HLIST,
                      help='By default the particles will be extracted \n'
                           'from the micrographs used in the picking \n'
                           'step ( _same as picking_ option ). \n'
                           'If you select _other_ option, you must provide \n'
                           'a different set of micrographs to extract from.\n'
                           '*Note*: In the _other_ case, ensure that provided \n'
                           'micrographs and coordinates are related \n'
                           'by micName or by micId. Difference in pixel size will \n'
                           'be handled automatically.')
        form.addParam('inputMicrographsTiltedPair', PointerParam,
                      pointerClass='MicrographsTiltPair',
                      condition='downsampleType != 0',
                      important=True, label='Input tilt pair micrographs',
                      help='Select the tilt pair micrographs from which to extract.')

        # downFactor should always be 1.0 or greater
        geOne = GE(1.0, error='Value should be greater or equal than 1.0')

        form.addParam('downFactor', FloatParam, default=1.0,
                      validators=[geOne],
                      label='Downsampling factor',
                      help='Select a value greater than 1.0 to reduce the size '
                           'of micrographs before extracting the particles. '
                           'If 1.0 is used, no downsample is applied. '
                           'Non-integer downsample factors are possible. ')

        form.addParam('boxSize', IntParam, default=0,
                      label='Particle box size', validators=[Positive],
                      help='In pixels. The box size is the size of the boxed particles, '
                           'actual particles may be smaller than this. If you do downsampling '
                           'after extraction, provide final box size here.')

        form.addSection(label='Preprocess')
        form.addParam('doRemoveDust', BooleanParam, default=True, important=True,
                      label='Dust removal (Recommended)',
                      help='Sets pixels with unusually large values to random values from a Gaussian '
                           'with zero-mean and unity-standard deviation.')
        form.addParam('thresholdDust', FloatParam, default=3.5,
                      condition='doRemoveDust', expertLevel=LEVEL_ADVANCED,
                      label='Threshold for dust removal',
                      help='Pixels with a signal higher or lower than this value times the standard '
                           'deviation of the image will be affected. For cryo, 3.5 is a good value.'
                           'For high-contrast negative stain, the signal itself may be affected so '
                           'that a higher value may be preferable.')
        form.addParam('doInvert', BooleanParam, default=None,
                      label='Invert contrast',
                      help='Invert the contrast if your particles are black over a white background.\n'
                           'Xmipp, Spider, Relion and Eman require white particles over a black background\n'
                           'Frealign (up to v9.07) requires black particles over a white background')
        form.addParam('doNormalize', BooleanParam, default=True,
                      label='Normalize (Recommended)',
                      help='It subtract a ramp in the gray values and normalizes so that in the '
                           'background there is 0 mean and standard deviation 1.')
        form.addParam('normType', EnumParam, choices=['OldXmipp', 'NewXmipp', 'Ramp'], default=2,
                      condition='doNormalize', expertLevel=LEVEL_ADVANCED,
                      display=EnumParam.DISPLAY_COMBO,
                      label='Normalization type',
                      help='OldXmipp (mean(Image)=0, stddev(Image)=1).  \n  '
                           'NewXmipp (mean(background)=0, stddev(background)=1)  \n  '
                           'Ramp (subtract background+NewXmipp).  \n  ')
        form.addParam('backRadius', IntParam, default=-1, condition='doNormalize',
                      label='Background radius', expertLevel=LEVEL_ADVANCED,
                      help='Pixels outside this circle are assumed to be noise and their stddev '
                           'is set to 1. Radius for background circle definition (in pix.). '
                           'If this value is 0, then half the box size is used.')

        form.addParallelSection(threads=4, mpi=1)

    def _insertAllSteps(self):
        self._setupBasicProperties()
        # Write pos files for each micrograph
        firstStepId = self._insertFunctionStep('writePosFilesStep')

        # For each micrograph insert the steps, run in parallel
        deps = []
        for mic in self.inputMics:
            localDeps = [firstStepId]
            fnLast = mic.getFileName()
            micName = pwutils.removeBaseExt(mic.getFileName())

            def getMicTmp(suffix):
                return self._getTmpPath(micName + suffix)

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

            # TODO: implement CTF
            # Actually extract
            deps.append(self._insertFunctionStep('extractParticlesStep',
                                                 mic.getObjId(), micName,
                                                 None, fnLast, micOps,
                                                 self.doInvert.get(),
                                                 self._getNormalizeArgs(),
                                                 prerequisites=localDeps))

        metaDeps = self._insertFunctionStep('createMetadataImageStep', prerequisites=deps)

        # Insert step to create output objects
        self._insertFunctionStep('createOutputStep', prerequisites=[metaDeps])

    # --------------------------- STEPS functions --------------------------------------------
    def writePosFilesStep(self):
        """ Write the pos file for each micrograph in metadata format
        (both untilted and tilted). """
        writeSetOfCoordinates(self._getExtraPath(),
                              self.inputCoords.getUntilted(),
                              scale=self.getBoxScale())
        writeSetOfCoordinates(self._getExtraPath(),
                              self.inputCoords.getTilted(),
                              scale=self.getBoxScale())

        # We need to find the mapping by micName (without ext) between the micrographs in
        # the SetOfCoordinates and the Other micrographs
        if self._micsOther():
            micDict = {}
            # create tmp set with all mics from coords set
            coordMics = SetOfMicrographs(filename=':memory:')
            coordMics.copyInfo(self.inputCoords.getUntilted().getMicrographs())

            for micU, micT in izip(self.inputCoords.getUntilted().getMicrographs(),
                                   self.inputCoords.getTilted().getMicrographs()):
                micU.cleanObjId()
                micT.cleanObjId()
                coordMics.append(micU)
                coordMics.append(micT)

            for mic in coordMics:
                micBase = pwutils.removeBaseExt(mic.getFileName())
                micPos = self._getExtraPath(micBase + ".pos")
                micDict[pwutils.removeExt(mic.getMicName())] = micPos

            # now match micDict and inputMics
            if any(pwutils.removeExt(mic.getMicName()) in micDict for mic in self.inputMics):
                micKey = lambda mic: pwutils.removeExt(mic.getMicName())
            else:
                raise Exception('Could not match input micrographs and coordinates '
                                'by micName.')

            for mic in self.inputMics:  # micrograph from input (other)
                mk = micKey(mic)
                if mk in micDict:
                    micPosCoord = micDict[mk]
                    if exists(micPosCoord):
                        micBase = pwutils.removeBaseExt(mic.getFileName())
                        micPos = self._getExtraPath(micBase + ".pos")
                        if micPos != micPosCoord:
                            self.info('Moving %s -> %s' % (micPosCoord, micPos))
                            pwutils.moveFile(micPosCoord, micPos)

    def createMetadataImageStep(self):
        mdUntilted = md.MetaData()
        mdTilted = md.MetaData()
        # for objId in mdPairs:
        for uMic, tMic in izip(self.uMics, self.tMics):
            umicName = pwutils.removeBaseExt(uMic.getFileName())
            fnMicU = self._getExtraPath(umicName + ".xmd")
            fnPosU = self._getExtraPath(umicName + ".pos")
            # Check if there are picked particles in these micrographs
            if pwutils.exists(fnMicU):
                mdMicU = md.MetaData(fnMicU)
                mdPosU = md.MetaData('particles@%s' % fnPosU)
                mdPosU.merge(mdMicU)
                mdUntilted.unionAll(mdPosU)
                tmicName = pwutils.removeBaseExt(tMic.getFileName())
                fnMicT = self._getExtraPath(tmicName + ".xmd")
                fnPosT = self._getExtraPath(tmicName + ".pos")
                mdMicT = md.MetaData(fnMicT)
                mdPosT = md.MetaData('particles@%s' % fnPosT)
                mdPosT.merge(mdMicT)
                mdTilted.unionAll(mdPosT)

        # Write image metadata (check if it is really necessary)
        fnTilted = self._getExtraPath("images_tilted.xmd")
        fnUntilted = self._getExtraPath("images_untilted.xmd")
        mdUntilted.write(fnUntilted)
        mdTilted.write(fnTilted)

    def createOutputStep(self):
        fnTilted = self._getExtraPath("images_tilted.xmd")
        fnUntilted = self._getExtraPath("images_untilted.xmd")

        # Create outputs SetOfParticles both for tilted and untilted
        imgSetU = self._createSetOfParticles(suffix="Untilted")
        imgSetU.copyInfo(self.uMics)
        imgSetT = self._createSetOfParticles(suffix="Tilted")
        imgSetT.copyInfo(self.tMics)

        sampling = self.samplingMics if self._micsOther() else self.samplingInput
        if self._doDownsample():
            sampling *= self.downFactor.get()
        imgSetU.setSamplingRate(sampling)
        imgSetT.setSamplingRate(sampling)

        # set coords from the input, will update later if needed
        imgSetU.setCoordinates(self.inputCoordinatesTiltedPairs.get().getUntilted())
        imgSetT.setCoordinates(self.inputCoordinatesTiltedPairs.get().getTilted())

        # Read untilted and tilted particles on a temporary object (also disabled particles)
        imgSetAuxU = SetOfParticles(filename=':memory:')
        imgSetAuxU.copyInfo(imgSetU)
        readSetOfParticles(fnUntilted, imgSetAuxU, removeDisabled=False)

        imgSetAuxT = SetOfParticles(filename=':memory:')
        imgSetAuxT.copyInfo(imgSetT)
        readSetOfParticles(fnTilted, imgSetAuxT, removeDisabled=False)

        # calculate factor for coords scaling
        factor = 1 / self.samplingFactor
        if self._doDownsample():
            factor /= self.downFactor.get()

        coordsT = self.inputCoordinatesTiltedPairs.get().getTilted()
        # For each untilted particle retrieve micId from SetOfCoordinates untilted
        for imgU, coordU in izip(imgSetAuxU, self.inputCoordinatesTiltedPairs.get().getUntilted()):
            # FIXME: Remove this check when sure that objIds are equal
            id = imgU.getObjId()
            if id != coordU.getObjId():
                raise Exception('ObjIds in untilted img and coord are not the same!!!!')
            imgT = imgSetAuxT[id]
            coordT = coordsT[id]

            # If both particles are enabled append them
            if imgU.isEnabled() and imgT.isEnabled():
                if self._micsOther() or self._doDownsample():
                    coordU.scale(factor)
                    coordT.scale(factor)
                imgU.setCoordinate(coordU)
                imgSetU.append(imgU)
                imgT.setCoordinate(coordT)
                imgSetT.append(imgT)

        imgSetU.write()
        imgSetT.write()

        # Define output ParticlesTiltPair
        outputset = ParticlesTiltPair(filename=self._getPath('particles_pairs.sqlite'))
        outputset.setTilted(imgSetT)
        outputset.setUntilted(imgSetU)
        for imgU, imgT in izip(imgSetU, imgSetT):
            outputset.append(TiltPair(imgU, imgT))

        outputset.setCoordsPair(self.inputCoordinatesTiltedPairs.get())
        self._defineOutputs(outputParticlesTiltPair=outputset)
        self._defineSourceRelation(self.inputCoordinatesTiltedPairs, outputset)

    # --------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []

        if self.doNormalize:
            if self.backRadius > int(self.boxSize.get() / 2):
                errors.append("Background radius for normalization should be "
                              "equal or less than half of the box size.")
        return errors

    def _citations(self):
        return ['Vargas2013b']

    def _summary(self):
        summary = []
        summary.append("Micrographs source: %s"
                       % self.getEnumText('downsampleType'))
        summary.append("Particle box size: %d" % self.boxSize)

        if not hasattr(self, 'outputParticlesTiltPair'):
            summary.append("Output images not ready yet.")
        else:
            summary.append("Particle pairs extracted: %d" %
                           self.outputParticlesTiltPair.getSize())

        return summary

    def _methods(self):
        methodsMsgs = []

        if self.getStatus() == STATUS_FINISHED:
            msg = "A total of %d particle pairs of size %d were extracted" % (self.getOutput().getSize(),
                                                                              self.boxSize)
            if self._micsOther():
                msg += " from another set of micrographs: %s" % self.getObjectTag('inputMicrographsTiltedPair')

            msg += " using coordinates %s" % self.getObjectTag('inputCoordinatesTiltedPairs')
            msg += self.methodsVar.get('')
            methodsMsgs.append(msg)

            if self.doRemoveDust:
                methodsMsgs.append("Removed dust over a threshold of %s." % (self.thresholdDust))
            if self.doInvert:
                methodsMsgs.append("Inverted contrast on images.")
            if self._doDownsample():
                methodsMsgs.append("Particles downsampled by a factor of %0.2f." % self.downFactor)
            if self.doNormalize:
                methodsMsgs.append("Normalization performed of type %s." % (self.getEnumText('normType')))

        return methodsMsgs

    # --------------------------- UTILS functions --------------------------------------------
    def _setupBasicProperties(self):
        # Get sampling rate and inputMics according to micsSource type
        self.inputCoords = self.getCoords()
        self.uMics, self.tMics = self.getInputMicrographs()
        self.samplingInput = self.inputCoords.getUntilted().getMicrographs().getSamplingRate()
        self.samplingMics = self.uMics.getSamplingRate()
        self.samplingFactor = float(self.samplingMics / self.samplingInput)

        # create tmp set with all mic pairs
        self.inputMics = SetOfMicrographs(filename=':memory:')
        self.inputMics.copyInfo(self.uMics)
        self.inputMics.setStore(False)

        for micU, micT in izip(self.uMics, self.tMics):
            micU.cleanObjId()
            micT.cleanObjId()
            self.inputMics.append(micU)
            self.inputMics.append(micT)

    def getInputMicrographs(self):
        """ Return pairs of micrographs associated to the SetOfCoordinates or
        Other micrographs. """
        if not self._micsOther():
            return self.inputCoordinatesTiltedPairs.get().getUntilted().getMicrographs(), \
                   self.inputCoordinatesTiltedPairs.get().getTilted().getMicrographs()
        else:
            return self.inputMicrographsTiltedPair.get().getUntilted(), \
                   self.inputMicrographsTiltedPair.get().getTilted()

    def getCoords(self):
        return self.inputCoordinatesTiltedPairs.get()

    def getOutput(self):
        if self.hasAttribute('outputParticlesTiltPair') and self.outputParticlesTiltPair.hasValue():
            return self.outputParticlesTiltPair
        else:
            return None

    def getCoordSampling(self):
        return self.getCoords().getUntilted().getMicrographs().getSamplingRate()

    def getMicSampling(self):
        return self.getInputMicrographs()[0].getSamplingRate()

