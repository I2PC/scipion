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

"""
This file contains the picking protocols for Relion 2.0 or greater.
"""

import os
from os.path import relpath
import pyworkflow.protocol.params as params
from pyworkflow.protocol.constants import STEPS_PARALLEL
from pyworkflow.em.protocol.protocol_particles import ProtParticlePicking
from pyworkflow.em.constants import RELATION_CTF
from protocol_base import ProtRelionBase
from convert import writeSetOfMicrographs, writeReferences, readSetOfCoordinates
from pyworkflow.em.convert import ImageHandler
import pyworkflow.utils as pwutils

REF_AVERAGES = 0
REF_BLOBS = 1

RUN_OPTIMIZE = 0 # Run only on several micrographs to optimize parameters
RUN_COMPUTE = 1 # Run the picking for all micrographs after optimize


class ProtRelion2Autopick(ProtParticlePicking, ProtRelionBase):
    """
    This Relion protocol uses 2D class averages as templates to run the auto-picking
    job-type. In this first stage, the auto-picking will be run just in few micrographs
    to optimise two of its main parameters ( _Picking threshold_ and _Minimum inter-particle distance_).

    In order to save time, only 2 or 3 micrographs should be used with their CTF
    information. One should use representative micrographs for the entire data set,
    e.g. a high and a low-defocus one, and/or with thin or thick ice.

    The expensive part of this calculation is to calculate a probability-based figure-of-merit
    (related to the cross-correlation coefficient between each rotated reference and all positions
    in the micrographs. This calculation is followed by a much cheaper peak-detection algorithm that
    uses the threshold and minimum distance parameters mentioned above. Because these parameters
    need to be optimised, this first stage of the auto-picking will write out so-called FOM maps.
    These are two large (micrograph-sized) files per reference. To avoid running into hard disc I/O
    problems, the autopicking program can only be run sequentially (hence there is no option to use
    more than one single MPI processor).
    """
    _label = 'auto-picking'

    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        # Convert the input micrographs and references to
        # the required Relion star files
        convertId = self._insertFunctionStep('convertInputStep',
                                             self.getInputMicrographs().strId(),
                                             self.getInputReferences().strId())
        # Insert the picking steps for each micrograph separately 
        # (Relion requires in that way if micrographs have different dimensions)
        allMicSteps = []
        for mic in self.getInputMicrographs():
            autopickId = self._insertAutopickStep(mic, convertId)
            allMicSteps.append(autopickId)
        # Register final coordinates as output
        self._insertFunctionStep('createOutputStep',
                                 prerequisites=allMicSteps)

    #--------------------------- STEPS functions -------------------------------
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

    def _getMicStarFile(self, mic):
        return self._getExtraPath(pwutils.replaceBaseExt(mic.getFileName(),
                                                         'star'))

    def _postprocessMicrographRow(self, img, imgRow):
        imgRow.writeToFile(self._getMicStarFile(img))

    def convertInputStep(self, micsId, refsId):
        self._ih = ImageHandler() # used to convert micrographs
        # Match ctf information against the micrographs
        self.ctfDict = {}
        if self.ctfRelations.get() is not None:
            for ctf in self.ctfRelations.get():
                self.ctfDict[ctf.getMicrograph().getMicName()] = ctf.clone()

        micStar = self._getPath('input_micrographs.star')
        writeSetOfMicrographs(self.getInputMicrographs(), micStar,
                              preprocessImageRow=self._preprocessMicrographRow,
                              postprocessImageRow=self._postprocessMicrographRow)

        writeReferences(self.getInputReferences(),
                        self._getPath('input_references'), useBasename=True)

    def _pickMicrograph(self, micStarFile, params, threshold,
                               minDistance, fom):
        """ Launch the 'relion_autopick' for a micrograph with the given parameters. """
        # Call relion_autopick to allow picking of micrographs with different size
        params += ' --i %s' % relpath(micStarFile, self.getWorkingDir())
        params += ' --threshold %0.3f --min_distance %0.3f %s' % (threshold,
                                                                  minDistance,
                                                                  fom)

        self.runJob(self._getProgram('relion_autopick'), params,
                    cwd=self.getWorkingDir())

    #--------------------------- UTILS functions -------------------------------
    def getInputReferences(self):
        return self.inputReferences.get()

    def getInputMicrographsPointer(self):
        return self.inputMicrographs

    def getInputMicrographs(self):
        return self.getInputMicrographsPointer().get()

    def getCoordsDir(self):
        return self._getTmpPath('xmipp_coordinates')

    def _writeXmippCoords(self, coordSet):
        micSet = self.getInputMicrographs()
        coordPath = self._getTmpPath('xmipp_coordinates')
        pwutils.cleanPath(coordPath)
        pwutils.makePath(coordPath)
        import pyworkflow.em.packages.xmipp3 as xmipp3
        #         micPath = os.path.join(coordPath, 'micrographs.xmd')
        #         xmipp3.writeSetOfMicrographs(micSet, micPath)
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
        coordSet.setBoxSize(self.getInputReferences().getDim()[0])
        starFiles = [self._getExtraPath(pwutils.removeBaseExt(mic.getFileName())
                                        + '_autopick.star') for mic in micSet]
        readSetOfCoordinates(coordSet, starFiles)
        return self._writeXmippCoords(coordSet)

    #--------------------------- DEFINE param functions ------------------------
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
        form.addParam('ctfRelations', params.RelationParam, allowsNull=True,
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

        form.addParam('micrographsList', params.IntParam, default='10',
                      label='Micrographs for optimization:',
                      help='Select the micrographs to be used for parameters '
                           'optimization. If you provide a single number, it '
                           'will be considered as the number of micrographs '
                           'you want to use and they will be picked in a way '
                           'to cover the CTF defocus range. On the other hand, '
                           'you can provide the IDs of each micrograph '
                           'separated by spaces. IDs can be typed or selected '
                           'using the provided wizard. ')

        # From Relion 2.+, it can be picked with gaussian blobs, so we
        # need to add these parameters
        refCondition = 'referencesType==%s' % REF_AVERAGES

        form.addParam('referencesType', params.EnumParam,
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

        form.addParam('gaussianPeak', params.FloatParam, default=0.1,
                      condition='referencesType==%s' % REF_BLOBS,
                      label='Gaussian peak value',
                      help='The peak value of the Gaussian blob. '
                           'Weaker data will need lower values.')

        form.addParam('inputReferences', params.PointerParam,
                      pointerClass='SetOfAverages',
                      condition=refCondition,
                      label='Input references', important=True,
                      help='Input references (SetOfAverages) for auto-pick. \n\n'
                           'Note that the absolute greyscale needs to be correct, \n'
                           'so only use images with proper normalization.')

        form.addParam('particleDiameter', params.IntParam, default=200,
                      label='Mask diameter (A)',
                      help='Diameter of the circular mask that will be applied '
                           'around the templates in Angstroms. When set to a '
                           'negative value, this value is estimated '
                           'automatically from the templates themselves.')

        form.addSection('References')

        form.addParam('lowpassFilterRefs', params.IntParam, default=20,
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
                      condition=refCondition,
                      default=False,
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

        form.addParam('pickingThreshold', params.FloatParam, default=0.25,
                      label='Picking threshold',
                      help='Use lower thresholds to pick more particles '
                           '(and more junk probably)')
        form.addParam('interParticleDistance', params.IntParam, default=100,
                      label='Minimum inter-particle distance (A)',
                      help='Particles closer together than this distance \n'
                           'will be consider to be a single cluster. \n'
                           'From each cluster, only one particle will be '
                           'picked.')

        form.addParam('shrinkFactor', params.FloatParam, default=1,
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

        form.addParam('doGpu', params.BooleanParam, default=True,
                      label='Use GPU acceleration?',
                      help='If set to Yes, the job will try to use GPU '
                           'acceleration.')

        form.addParam('gpusToUse', params.StringParam, default='',
                      label='Which GPUs to use:', condition='doGpu',
                      help='This argument is not necessary. If left empty, '
                           'the job itself will try to allocate available GPU '
                           'resources. You can override the default '
                           'allocation by providing a list of which GPUs '
                           '(0,1,2,3, etc) to use. MPI-processes are '
                           'separated by ":", threads by ",". '
                           'For example: "0,0:1,1:0,0:1,1"')

        form.addSection('Helix')
        form.addParam('fomLabel', params.LabelParam,
                      important=True,
                      label='Helix processing is not implemented still.')

    #--------------------------- INSERT steps functions ------------------------

    def getAutopickParams(self):
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
            params += ' --ref gauss --gaus_max %0.3f' % self.gaussianPeak

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

        return params

    def _insertAutopickStep(self, mic, convertId):
        return self._insertFunctionStep('autopickStep',
                                        self._getMicStarFile(mic),
                                        self.getAutopickParams(), 0.25,
                                        -1, ' --write_fom_maps',
                                        prerequisites=[convertId])

    #--------------------------- STEPS functions -------------------------------
    def autopickStep(self, micStarFile, params, threshold, minDistance, fom):
        """ This method is used from the wizard to optimize the parameters. """
        for mic in self.getInputMicrographs():
            self._pickMicrograph(micStarFile, params,
                                 threshold, minDistance, fom)

    def createOutputStep(self):
        micSet = self.getInputMicrographs()
        coordSet = self._createSetOfCoordinates(micSet)
        coordSet.setBoxSize(self.getInputReferences().getDim()[0])
        starFiles = [
            self._getExtraPath(pwutils.removeBaseExt(mic.getFileName())
                               + '_autopick.star') for mic in micSet]
        readSetOfCoordinates(coordSet, starFiles)

        self._defineOutputs(outputCoordinates=coordSet)
        self._defineSourceRelation(self.getInputMicrographsPointer(),
                                   coordSet)

    #--------------------------- INFO functions --------------------------------
    def _validate(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        errors = []
        self.validatePackageVersion('RELION_HOME', errors)

        if (self.useInputReferences() and
            self.particleDiameter > self.getInputDimA()):
            errors.append('Particle diameter (%d) can not be greater than '
                          'size (%d)' % (self.particleDiameter,
                                         self.getInputDimA()))

        if self.ctfRelations.get() is None and self.refsCtfCorrected:
            errors.append("References CTF corrected parameter must be set to "
                          "False or set ctf relations.")

        return errors

    def _summary(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        summary = []
        return summary

    #--------------------------- UTILS functions -------------------------------
    def useInputReferences(self):
        return self.referencesType == REF_AVERAGES

    def getInputDimA(self):
        """ Return the dimension of input references in A. """
        inputRefs = self.getInputReferences()
        if inputRefs is None:
            return None
        else:
            return inputRefs.getDim()[0] * inputRefs.getSamplingRate()
                
        
        
class ProtRelionAutopick():
    """    
    This Relion protocol uses 2D class averages as templates to run the auto-picking 
    job-type. In this second stage, the protocol reads the FOM maps to optimize
    the 'Threshold' and 'Inter-particle distance'.
    """
    _label = 'auto-picking (step 2)'

    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        
        group = form.addGroup('Autopicking')
           
        group.addParam('fomLabel', params.LabelParam,
                      important=True,
                      label="Modify parameters and click the 'wizard' "
                            "to see output coordinates")

        
        form.addParallelSection(threads=4, mpi=1)
        
    #--------------------------- INSERT steps functions ------------------------
        
    def _insertAutopickStep(self, mic, convertId):
        """ Prepare the command line for calling 'relion_autopick' program """
        params = self.getAutopickParams()
        # Use by default the references size as --min_distance
        return self._insertFunctionStep('autopickStep',
                                        self._getMicStarFile(mic),
                                        params,
                                        self.pickingThreshold.get(),
                                        self.interParticleDistance.get(),
                                        '', # neither read or write fom
                                        prerequisites=[convertId])
                
    #--------------------------- STEPS functions -------------------------------
    

    
