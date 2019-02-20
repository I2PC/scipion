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
from os.path import relpath

import pyworkflow.protocol.params as params
from pyworkflow.em.protocol import ProtParticlePickingAuto
from pyworkflow.em.constants import RELATION_CTF, ALIGN_NONE
from pyworkflow.em.convert import ImageHandler
from pyworkflow.utils.properties import Message
import pyworkflow.utils as pwutils
import pyworkflow.em.metadata as md
from pyworkflow.em import getSubsetByDefocus

from protocol_base import ProtRelionBase
from convert import (writeSetOfMicrographs, writeReferences,
                     readSetOfCoordinates, isVersion2, micrographToRow)


REF_AVERAGES = 0
REF_BLOBS = 1

RUN_OPTIMIZE = 0 # Run only on several micrographs to optimize parameters
RUN_COMPUTE = 1 # Run the picking for all micrographs after optimize

MICS_AUTO = 0
MICS_SUBSET = 1


class ProtRelion2Autopick(ProtParticlePickingAuto, ProtRelionBase):
    """
    This Relion protocol uses the 'relion_autopick' program to pick particles
    from micrographs, either using templates or gaussian blobs.

    The picking with this protocol is divided in three steps:
    1) Run with 'Optimize' option for several (less than 30) micrographs.
    2) Execute the wizard to refine the picking parameters.
    3) Run with 'Pick all' option to pick particles from all micrographs.

    The first steps will use internally the option '--write-fom-maps' to write
    to disk the FOM maps. The expensive part of this calculation is to calculate
    a probability-based figure-of-merit (related to the cross-correlation
    coefficient between each rotated reference and all positions in the
    micrographs. That's why it is only done in an small subset of the
    micrographs, where one should use representative micrographs for the entire
    data set, e.g. a high and a low-defocus one, and/or with thin or thick ice.

    Step 2 uses a much cheaper peak-detection algorithm that uses the threshold
    and minimum distance parameters.
    """
    _label = 'auto-picking'

    @classmethod
    def isDisabled(cls):
        return not isVersion2()

    # -------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputMicrographs', params.PointerParam,
                      pointerClass='SetOfMicrographs',
                      label='Input micrographs', important=True,
                      help='Select the input micrographs. '
                           'If using the *Optimize* mode, just a subset of '
                           'micrographs are used to compute the FOM maps. '
                           'If in *Compute* mode, all micrographs will be '
                           'auto-picked.')
        form.addParam('ctfRelations', params.RelationParam,
                      relationName=RELATION_CTF,
                      attributeName='getInputMicrographs',
                      label='CTF estimation',
                      help='Choose some CTF estimation related to the '
                           'input micrographs.')

        form.addParam('runType', params.EnumParam, default=RUN_OPTIMIZE,
                      choices=['Optimize params', 'Pick all micrographs'],
                      display=params.EnumParam.DISPLAY_LIST,
                      label='Run type: ',
                      help='Usually, first you should use the *Optimize* mode '
                           'to compute the FOM maps for a few micrographs and '
                           'use them to tune the picking parameters using the '
                           'wizard. After that you can run the job in *Compute*'
                           ' mode and auto-pick all the micrographs. ')

        group = form.addGroup('Micrographs for optimization',
                              condition='runType==%d' % RUN_OPTIMIZE)

        group.addParam('micrographsSelection', params.EnumParam,
                       default=MICS_AUTO,
                       choices=['automatic selection', 'input subset'],
                       display=params.EnumParam.DISPLAY_HLIST,
                       label='Choose micrographs by',
                       help='If you choose "automatic selection", you only '
                            'need to provide the number of microgrphs to use '
                            'and that number will be selected to cover the '
                            'defocus range. ')
        group.addParam('micrographsNumber', params.IntParam, default='10',
                      condition='micrographsSelection==%d' % MICS_AUTO,
                      label='Micrographs for optimization:',
                      help='Select the number of micrographs that you want'
                           'to be used for the parameters optimization. ')
        group.addParam('micrographsSubset', params.PointerParam,
                       condition='micrographsSelection==%d' % MICS_SUBSET,
                       pointerClass='SetOfMicrographs',
                       label='Subset of micrographs',
                       help='Choose as input a subset of micrographs that '
                            'you have previously selected. '
                            '(Probably covering the defocus range).')

        # From Relion 2.+, it can be picked with gaussian blobs, so we
        # need to add these parameters
        refCondition = 'referencesType==%s' % REF_AVERAGES

        group = form.addGroup('References')
        group.addParam('referencesType', params.EnumParam,
                      choices=['References', 'Gaussian blobs'],
                      default=REF_AVERAGES,
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='References type',
                      help='You may select "Gaussian blobs" to be used as '
                           'references. The preferred way to autopick is '
                           'by providing 2D references images that were '
                           'by 2D classification. \n'
                           'The Gaussian blob references may be useful to '
                           'kickstart a new data set.')

        group.addParam('gaussianPeak', params.FloatParam, default=0.1,
                      condition='referencesType==%s' % REF_BLOBS,
                      label='Gaussian peak value',
                      help='The peak value of the Gaussian blob. '
                           'Weaker data will need lower values.')

        group.addParam('inputReferences', params.PointerParam,
                      pointerClass='SetOfAverages',
                      condition=refCondition,
                      label='Input references', important=True,
                      help='Input references (SetOfAverages) for auto-pick. \n\n'
                           'Note that the absolute greyscale needs to be correct, \n'
                           'so only use images with proper normalization.')

        group.addParam('particleDiameter', params.IntParam, default=-1,
                      label='Mask diameter (A)',
                      help='Diameter of the circular mask that will be applied '
                           'around the templates in Angstroms. When set to a '
                           'negative value, this value is estimated '
                           'automatically from the templates themselves.')

        form.addSection('References')

        form.addParam('lowpassFilterRefs', params.IntParam, default=20,
                      condition=refCondition,
                      label='Lowpass filter references (A)',
                      help='Lowpass filter that will be applied to the '
                           'references before template matching. \n'
                           'Do NOT use very high-resolution templates to '
                           'search your micrographs. \n'
                           'The signal will be too weak at high resolution '
                           'anyway, and you may find Einstein from noise...')

        form.addParam('highpassFilterMics', params.IntParam, default=-1,
                      label='Highpass filter (A)',
                      help='Highpass filter that will be applied to the '
                           'micrographs. This may be useful to get rid of '
                           'background ramps due to uneven ice distributions. '
                           'Give a negative value to skip the highpass '
                           'filter.  Useful values are often in the range '
                           'of 200-400 Angstroms.')

        form.addParam('angularSampling', params.IntParam, default=5,
                      label='Angular sampling (deg)',
                      help='Angular sampling in degrees for exhaustive searches '
                           'of the in-plane rotations for all references.')

        form.addParam('refsHaveInvertedContrast', params.BooleanParam,
                      default=True,
                      label='References have inverted contrast?',
                      help='Set to Yes to indicate that the reference have '
                           'inverted contrast with respect to the particles '
                           'in the micrographs.')

        form.addParam('refsCtfCorrected', params.BooleanParam, default=True,
                      condition=refCondition,
                      label='Are References CTF corrected?',
                      help='Set to Yes if the references were created with '
                           'CTF-correction inside RELION.\n'
                           'If set to Yes, the input micrographs should contain '
                           'the CTF information.')

        form.addParam('ignoreCTFUntilFirstPeak', params.BooleanParam,
                      condition=refCondition,
                      default=False, expertLevel=params.LEVEL_ADVANCED,
                      label='Ignore CTFs until first peak?',
                      help='Set this to Yes, only if this option was also used '
                           'to generate the references.')

        form.addSection('Autopicking')

        group = form.addGroup('Autopick')
        group.addParam('pickingThreshold', params.FloatParam, default=0.25,
                      label='Picking threshold:',
                      help='Use lower thresholds to pick more particles '
                           '(and more junk probably)')

        group.addParam('interParticleDistance', params.IntParam, default=-1,
                      label='Minimum inter-particle distance (A):',
                      help='Particles closer together than this distance \n'
                           'will be consider to be a single cluster. \n'
                           'From each cluster, only one particle will be '
                           'picked.')

        group.addParam('maxStddevNoise', params.FloatParam, default=1.1,
                      label='Maximum stddev noise:',
                      help='This is useful to prevent picking in carbon areas, '
                           'or areas with big contamination features. Peaks in '
                           'areas where the background standard deviation in '
                           'the normalized micrographs is higher than this '
                           'value will be ignored. Useful values are probably '
                           'in the range 1.0 to 1.2. Set to -1 to switch off '
                           'the feature to eliminate peaks due to high '
                           'background standard deviations.')

        group = form.addGroup('Computing')
        group.addParam('shrinkFactor', params.FloatParam, default=0,
                      validators=[params.Range(0, 1, "value should be "
                                                     "between 0 and 1. ")],
                      label='Shrink factor',
                      help='This is useful to speed up the calculations, '
                           'and to make them less memory-intensive. The '
                           'micrographs will be downscaled (shrunk) to '
                           'calculate the cross-correlations, and peak '
                           'searching will be done in the downscaled FOM '
                           'maps. When set to 0, the micrographs will de '
                           'downscaled to the lowpass filter of the '
                           'references, a value between 0 and 1 will '
                           'downscale the micrographs by that factor. '
                           'Note that the results will not be exactly '
                           'the same when you shrink micrographs!')

        group.addParam('doGpu', params.BooleanParam, default=True,
                      label='Use GPU acceleration?',
                      help='If set to Yes, the job will try to use GPU '
                           'acceleration.')

        group.addParam('gpusToUse', params.StringParam, default='',
                      label='Which GPUs to use:', condition='doGpu',
                      help='This argument is not necessary. If left empty, '
                           'the job itself will try to allocate available GPU '
                           'resources. You can override the default '
                           'allocation by providing a list of which GPUs '
                           '(0,1,2,3, etc) to use. MPI-processes are '
                           'separated by ":", threads by ",". '
                           'For example: "0,0:1,1:0,0:1,1"')

        form.addParam('extraParams', params.StringParam, default='',
                      label='Additional arguments:',
                      help='In this box command-line arguments may be provided '
                           'that are not generated by the GUI. This may be '
                           'useful for testing developmental options and/or '
                           'expert use of the program. \n'
                           'The command "relion_autopick" will print a list '
                           'of possible options.')

        form.addSection('Helix')
        form.addParam('fomLabel', params.LabelParam,
                      important=True,
                      label='Helix processing is not implemented still.')

        self._defineStreamingParams(form)

        form.addParallelSection(threads=0, mpi=4)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self.inputStreaming = self.getInputMicrographs().isStreamOpen()

        if self.inputStreaming and not self.isRunOptimize():
            # If the input is in streaming, follow the base class policy
            # about inserting new steps and discovery new input/output
            ProtParticlePickingAuto._insertAllSteps(self)
            self.createOutputStep = self._doNothing
        else:
            # If not in streaming, then we will just insert a single step to
            # pick all micrographs at once since it is much faster
            self._insertInitialSteps()
            self._insertFunctionStep('_pickMicrographsFromStar',
                                     self._getPath('input_micrographs.star'),
                                     *self._getPickArgs())
            self._insertFunctionStep('createOutputStep')

            # Disable streaming functions:
            self._insertFinalSteps = self._doNothing
            self._stepsCheck = self._doNothing

    def _insertInitialSteps(self):
        # Convert the input micrographs and references to
        # the required Relion star files
        inputRefs = self.getInputReferences()
        refId = inputRefs.strId() if self.useInputReferences() else 'Gaussian'
        convertId = self._insertFunctionStep('convertInputStep',
                                             self.getInputMicrographs().strId(),
                                             refId, self.runType.get())
        return [convertId]

    def _doNothing(self, *args):
        pass # used to avoid some streaming functions

    def _loadInputList(self):
        """ This function is re-implemented in this protocol, because it have
         a SetOfCTF as input, so for streaming, we only want to report those
         micrographs for which the CTF is ready.
        """
        micDict, micClose = self._loadMics(self.getInputMicrographs())
        ctfDict, ctfClosed = self._loadCTFs(self.ctfRelations.get())

        # Remove the micrographs that have not CTF
        # and set the CTF property for those who have it
        for micKey, mic in micDict.iteritems():
            if micKey in ctfDict:
                mic.setCTF(ctfDict[micKey])
            else:
                del micDict[micKey]

        # Return the updated micDict and the closed status
        return micDict, micClose and ctfClosed

    # -------------------------- STEPS functions ------------------------------

    def convertInputStep(self, micsId, refsId, runType):
        # runType is passed as parameter to force a re-execute of this step
        # if there is a change in the type

        self._ih = ImageHandler() # used to convert micrographs
        # Match ctf information against the micrographs
        self.ctfDict = {}
        if self.ctfRelations.get() is not None:
            for ctf in self.ctfRelations.get():
                self.ctfDict[ctf.getMicrograph().getMicName()] = ctf.clone()

        micStar = self._getPath('input_micrographs.star')
        writeSetOfMicrographs(self.getMicrographList(), micStar,
                              alignType=ALIGN_NONE,
                              preprocessImageRow=self._preprocessMicrographRow)

        if self.useInputReferences():
            writeReferences(self.getInputReferences(),
                            self._getPath('input_references'), useBasename=True)

        # FIXME: (JMRT-20180523) The following code does not seems to work
        # here it has been worked around by changing the name of the wizard
        # output but this seems to reflect a deeper problem of deleting
        # already existing output objects in a protocol. Maybe when updating
        # from run.db to project.sqlite?

        # Clean up if previously created the outputMicrographs and Coordinates
        # in the wizard - optimization run
        # if self.hasAttribute('outputMicrographs'):
        #     self._deleteChild('outputMicrographs', self.outputMicrographs)
        # if self.hasAttribute('outputCoordinates'):
        #     self._deleteChild('outputCoordinates', self.outputCoordinates)
        # self._store()

    def getAutopickParams(self):
        # Return the autopicking parameters except for the interative ones:
        # - threshold
        # - minDistance
        # - maxStd
        params = ' --pickname autopick'
        params += ' --odir ""'
        params += ' --particle_diameter %d' % self.particleDiameter
        params += ' --angpix %0.3f' % self.getInputMicrographs().getSamplingRate()
        params += ' --shrink %0.3f' % self.shrinkFactor

        if self.doGpu:
            params += ' --gpu "%s"' % self.gpusToUse

        # Now in Relion2.0 autopick can use gassian blobs
        if self.useInputReferences():
            params += ' --ref input_references.star'
            ps = self.getInputReferences().getSamplingRate()
            params += ' --angpix_ref %0.3f' % ps
        else: # Gaussian blobs
            params += ' --ref gauss --gauss_max %0.3f' % self.gaussianPeak

        if self.refsHaveInvertedContrast:
            params += ' --invert'

        if self.refsCtfCorrected:
            params += ' --ctf'

        params += ' --ang %d' % self.angularSampling
        # Negative values for filters means no-filter
        if self.lowpassFilterRefs > 0:
            params += ' --lowpass %d' % self.lowpassFilterRefs
        if self.highpassFilterMics > 0:
            params += ' --highpass %d' % self.highpassFilterMics

        # Add extra params is any
        params += ' %s' % self.extraParams

        return params

    def _getPickArgs(self):
        basicArgs = self.getAutopickParams()
        threshold = self.pickingThreshold.get()
        interDist = self.interParticleDistance.get()
        fomParam = ' --write_fom_maps' if self.isRunOptimize() else ''
        return [basicArgs, threshold, interDist, fomParam]

    def _pickMicrographsFromStar(self, micStarFile, params,
                                 threshold, minDistance, fom):
        """ Launch the 'relion_autopick' for micrographs in the inputStarFile.
         If the input set of complete, the star file will contain all the
         micrographs. If working in streaming, it will be only one micrograph.
        """
        params += ' --i %s' % relpath(micStarFile, self.getWorkingDir())
        params += ' --threshold %0.3f ' % threshold
        params += ' --min_distance %0.3f %s' % (minDistance, fom)

        program = self._getProgram('relion_autopick')

        self.runJob(program, params, cwd=self.getWorkingDir())

    def _pickMicrograph(self, mic, params, threshold, minDistance, fom):
        """ This method should be invoked only when working in streaming mode.
        """
        micRow = md.Row()
        self._preprocessMicrographRow(mic, micRow)
        micrographToRow(mic, micRow)
        self._postprocessMicrographRow(mic, micRow)
        self._pickMicrographsFromStar(self._getMicStarFile(mic), params,
                                      threshold, minDistance, fom)

    def _pickMicrographList(self, micList, params, threshold, minDistance, fom):
        micStar = self._getPath('input_micrographs_%s-%s.star' %
                                (micList[0].strId(), micList[-1].strId()))
        writeSetOfMicrographs(micList, micStar,
                              alignType=ALIGN_NONE,
                              preprocessImageRow=self._preprocessMicrographRow)
        self._pickMicrographsFromStar(micStar, params, threshold, minDistance,
                                      fom)

    def _createSetOfCoordinates(self, micSet, suffix=''):
        """ Override this method to set the box size. """
        coordSet = ProtParticlePickingAuto._createSetOfCoordinates(self, micSet,
                                                                   suffix=suffix)
        coordSet.setBoxSize(self.getBoxSize())

        return coordSet

    def readCoordsFromMics(self, workingDir, micList, coordSet):
        """ Parse back the output star files and populate the SetOfCoordinates.
        """
        template = self._getExtraPath("%s_autopick.star")
        starFiles = [template % pwutils.removeBaseExt(mic.getFileName())
                     for mic in micList]
        readSetOfCoordinates(coordSet, starFiles, micList)

    # -------------------------- STEPS functions -------------------------------
    def autopickStep(self, micStarFile, params, threshold,
                     minDistance, maxStddevNoise, fom):
        """ This method is used from the wizard to optimize the parameters. """
        self._pickMicrographsFromStar(micStarFile, params, threshold,
                                      minDistance, maxStddevNoise, fom)

    def createOutputStep(self):
        micSet = self.getInputMicrographs()
        outputCoordinatesName = 'outputCoordinates'
        outputSuffix = ''

        # If in optimization phase, let's create a subset of the micrographs
        if self.isRunOptimize():
            outputSuffix = '_subset'
            outputCoordinatesName = 'outputCoordinatesSubset'
            micSubSet = self._createSetOfMicrographs(suffix=outputSuffix)
            micSubSet.copyInfo(micSet)
            # Use previously written star file for reading the subset of micrographs,
            for row in md.iterRows(self._getPath('input_micrographs.star')):
                mic = micSet[row.getValue('rlnImageId')]
                micSubSet.append(mic)
            self._defineOutputs(outputMicrographsSubset=micSubSet)
            self._defineTransformRelation(self.getInputMicrographsPointer(),
                                          micSubSet)
            micSet = micSubSet

        coordSet = self._createSetOfCoordinates(micSet)
        template = self._getExtraPath("%s_autopick.star")
        starFiles = [template % pwutils.removeBaseExt(mic.getFileName())
                     for mic in micSet]
        readSetOfCoordinates(coordSet, starFiles, micSet)

        self._defineOutputs(**{outputCoordinatesName: coordSet})
        self._defineSourceRelation(self.getInputMicrographsPointer(),
                                   coordSet)

    # -------------------------- INFO functions --------------------------------
    def _validate(self):
        errors = []
        self.validatePackageVersion('RELION_HOME', errors)

        if self.useInputReferences():
            if self.particleDiameter > self.getInputDimA():
                errors.append('Particle diameter (%d) can not be greater than '
                              'size (%d)' % (self.particleDiameter,
                                             self.getInputDimA()))
            if self.getInputReferences().isOddX():
                errors.append("Relion only works with even values for the "
                              "average dimensions!")
        else:
            if self.particleDiameter <= 0:
                errors.append('When using Gaussian blobs, you need to specify '
                              'the particles diameter manually. ')

        if self.ctfRelations.get() is None and self.refsCtfCorrected:
            errors.append("References CTF corrected parameter must be set to "
                          "False or set ctf relations.")


        errors.extend(self._validateMicSelection())

        return errors

    def _validateMicSelection(self):
        """ Validate the cases when selecting a subset of micrographs
        to optimize.
        """
        inputMics = self.getInputMicrographs()
        inputCTFs = self.ctfRelations.get()

        if self.isRunOptimize():
            if self.micrographsSelection == MICS_AUTO:
                n = self.micrographsNumber.get()
                if n < 3 or n > min(30, inputMics.getSize()):
                    return ['Number of micrographs should be between 3 and '
                            'min(30, input_size)']
            else:
                micSubset = self.micrographsSubset.get()
                if micSubset is None:
                    return ['Select the subset of micrographs']

                def missing(mic):
                    micId = mic.getObjId()
                    return inputMics[micId] is None or inputCTFs[micId] is None

                if any(missing(mic) for mic in micSubset):
                    return ['Some selected micrograph IDs are missing from the '
                            'input micrographs or CTFs.']
        return []

    def _warnings(self):
        if not self.isRunOptimize():
            if not hasattr(self, 'wizardExecuted'):
                return ['It seems that you have not executed the wizard to '
                        'optimize the picking parameters. \n'
                        'Do you want to launch the whole picking anyway?']

        return []

    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methodsMsgs = []
        if self.getOutputsSize() > 0:
            output = self.getCoords()
            methodsMsgs.append("%s: User picked %d particles with a particle "
                               "size of %d px."
                               % (self.getObjectTag(output), output.getSize(),
                                  output.getBoxSize()))
        else:
            methodsMsgs.append(Message.TEXT_NO_OUTPUT_CO)
    
        return methodsMsgs

    # -------------------------- UTILS functions -------------------------------
    def useInputReferences(self):
        return self.referencesType == REF_AVERAGES

    def isRunOptimize(self):
        return self.runType == RUN_OPTIMIZE

    def getInputDimA(self):
        """ Return the dimension of input references in A. """
        inputRefs = self.getInputReferences()
        if inputRefs is None:
            return None
        else:
            return inputRefs.getXDim() * inputRefs.getSamplingRate()

    def getBoxSize(self):
        """ Return a reasonable box-size in pixels. """
        inputRefs = self.getInputReferences()
        inputMics = self.getInputMicrographs()
        micsSampling = inputMics.getSamplingRate()

        if inputRefs is None:
            boxSize = int(self.particleDiameter.get() * 1.25 / micsSampling)
        else:
            # Scale boxsize if the pixel size of the references is not the same
            # of the micrographs
            scale = inputRefs.getSamplingRate() / micsSampling
            boxSize = int(inputRefs.getXDim() * scale)

        if boxSize % 2 == 1:
            boxSize += 1 # Use even box size for relion

        return boxSize
                
    def getInputReferences(self):
        return self.inputReferences.get()

    def getInputMicrographsPointer(self):
        return self.inputMicrographs

    def getInputMicrographs(self):
        return self.getInputMicrographsPointer().get()

    def getMicrographList(self):
        """ Return the list of micrographs (either a subset or the full set)
        that will be used for optimizing the parameters or the picking.
        """
        # Use all micrographs only when going for the full picking
        inputMics = self.getInputMicrographs()

        if not self.isRunOptimize():
            return inputMics

        if self.micrographsSelection == MICS_AUTO:
            mics = getSubsetByDefocus(self.ctfRelations.get(), inputMics,
                                      self.micrographsNumber.get())
        else:  # Subset selection
            mics = [mic.clone() for mic in self.micrographsSubset.get()]

        return mics

    def getCoordsDir(self):
        return self._getTmpPath('xmipp_coordinates')

    def _writeXmippCoords(self, coordSet):
        micSet = self.getInputMicrographs()
        coordPath = self._getTmpPath('xmipp_coordinates')
        pwutils.cleanPath(coordPath)
        pwutils.makePath(coordPath)
        import pyworkflow.em.packages.xmipp3 as xmipp3
        micPath = micSet.getFileName()
        xmipp3.writeSetOfCoordinates(coordPath, coordSet, ismanual=False)
        return micPath, coordPath

    def writeXmippOutputCoords(self):
        return self._writeXmippCoords(self.outputCoordinates)

    def writeXmippCoords(self):
        """ Write the SetOfCoordinates as expected by Xmipp
        to be displayed with its GUI.
        """
        micSet = self.getInputMicrographs()
        coordSet = self._createSetOfCoordinates(micSet)
        coordSet.setBoxSize(self.getBoxSize())
        starFiles = [self._getExtraPath(pwutils.removeBaseExt(mic.getFileName())
                                        + '_autopick.star') for mic in micSet]
        readSetOfCoordinates(coordSet, starFiles)
        return self._writeXmippCoords(coordSet)

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
        # JMRT: The following is not needed since it is set when loading CTF
        # if self.ctfRelations.get() is not None:
        #     img.setCTF(self.ctfDict[img.getMicName()])

    def _postprocessMicrographRow(self, img, imgRow):
        imgRow.writeToFile(self._getMicStarFile(img))

    def _getMicStarFile(self, mic):
        micBase = pwutils.replaceBaseExt(mic.getFileName(), 'star')
        return self._getExtraPath(micBase)