# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (josue.gomez-blanco@mcgill.ca)
# *              Javier Vargas Balbuena (javier.vargasbalbuena@mcgill.ca)
# *
# * Department of Anatomy and Cell Biology, McGill University
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

import re
from glob import glob
from os.path import exists, basename, join

from pyworkflow.protocol.params import (BooleanParam, PointerParam, FloatParam,
                                        IntParam, EnumParam, StringParam,
                                        LabelParam, PathParam)
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.utils.path import cleanPath, replaceBaseExt

import pyworkflow.em as em
import pyworkflow.em.metadata as md
from pyworkflow.em.data import SetOfVolumes
from pyworkflow.em.protocol import EMProtocol

from convert import (writeSetOfParticles, isVersion2, convertMask)

ANGULAR_SAMPLING_LIST = ['30', '15', '7.5', '3.7', '1.8', '0.9', '0.5', '0.2', '0.1']
MASK_FILL_ZERO = 0
MASK_FILL_NOISE = 1


class ProtInitialVolumeSelector(EMProtocol):
    """ This class contains the common functions for all Relion protocols.
    In subclasses there should be little changes about how to create the command
    line and the files produced.

    Most of the Relion protocols, have two modes: NORMAL or CONTINUE. That's why
    some of the function have a template pattern approach to define the behaviour
    depending on the case.
    """
    _label = 'volume selector'

    OUTPUT_TYPE = SetOfVolumes
    FILE_KEYS = ['data', 'optimiser', 'sampling']
    PREFIXES = ['']

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    def _initialize(self):
        """ This function is mean to be called after the
        working dir for the protocol have been set. (maybe after recovery from mapper)
        """
        self._createFilenameTemplates()
        self._createIterTemplates()

        self.ClassFnTemplate = '%(rootDir)s/relion_it%(iter)03d_class%(ref)03d.mrc:mrc'

    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        self.extraIter = self._getExtraPath('relion_it%(iter)03d_')
        myDict = {
            'input_star': self._getPath('input_particles.star'),
            'input_mrcs': self._getPath('input_particles.mrcs'),
            'data_scipion': self.extraIter + 'data_scipion.sqlite',
            'projections': self.extraIter + '%(half)sclass%(ref3d)03d_projections.sqlite',
            'volumes_scipion': self.extraIter + 'volumes.sqlite',
            'data': self.extraIter + 'data.star',
            'model': self.extraIter + 'model.star',
            'shiny': self._getExtraPath('shiny/shiny.star'),
            'optimiser': self.extraIter + 'optimiser.star',
            'angularDist_xmipp': self.extraIter + 'angularDist_xmipp.xmd',
            'all_avgPmax_xmipp': self._getTmpPath(
                'iterations_avgPmax_xmipp.xmd'),
            'all_changes_xmipp': self._getTmpPath(
                'iterations_changes_xmipp.xmd'),
            'selected_volumes': self._getTmpPath('selected_volumes_xmipp.xmd'),
            'movie_particles': self._getPath('movie_particles.star'),
            'dataFinal': self._getExtraPath("relion_data.star"),
            'modelFinal': self._getExtraPath("relion_model.star"),
            'finalvolume': self._getExtraPath(
                "relion_class%(ref3d)03d.mrc:mrc"),
            'final_half1_volume': self._getExtraPath(
                "relion_half1_class%(ref3d)03d_unfil.mrc:mrc"),
            'final_half2_volume': self._getExtraPath(
                "relion_half2_class%(ref3d)03d_unfil.mrc:mrc"),
            'finalSGDvolume': self._getExtraPath(
                "relion_it%(iter)03d_class%(ref3d)03d.mrc:mrc"),
            'preprocess_particles': self._getPath("preprocess_particles.mrcs"),
            'preprocess_particles_star': self._getPath(
                "preprocess_particles.star"),
            'preprocess_particles_preffix': "preprocess_particles"
        }
        # add to keys, data.star, optimiser.star and sampling.star
        for key in self.FILE_KEYS:
            myDict[key] = self.extraIter + '%s.star' % key
            key_xmipp = key + '_xmipp'
            myDict[key_xmipp] = self.extraIter + '%s.xmd' % key
        # add other keys that depends on prefixes
        for p in self.PREFIXES:
            myDict['%smodel' % p] = self.extraIter + '%smodel.star' % p
            myDict[
                '%svolume' % p] = self.extraIter + p + 'class%(ref3d)03d.mrc:mrc'

        self._updateFilenamesDict(myDict)

    def _createIterTemplates(self):
        """ Setup the regex on how to find iterations. """
        self._iterTemplate = self._getFileName('data', iter=0).replace('000',
                                                                       '???')
        # Iterations will be identify by _itXXX_ where XXX is the iteration number
        # and is restricted to only 3 digits.
        self._iterRegex = re.compile('_it(\d{3,3})_')

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      important=True,
                      label="Input particles",
                      help='Select the input images from the project.')
        form.addParam('subsetSize', IntParam, default=1000,
                      label='Subset size',
                      help='Number of individual particles that will be use '
                           'to obtain the best initial volume')
        form.addParam('maskDiameterA', IntParam, default=-1,
                      label='Particle mask diameter (A)',
                      help='The experimental images will be masked with a '
                           'soft circular mask with this <diameter>. '
                           'Make sure this diameter is not set too small '
                           'because that may mask away part of the signal! If '
                           'set to a value larger than the image size no '
                           'masking will be performed.\n\n'
                           'The same diameter will also be used for a '
                           'spherical mask of the reference structures if no '
                           'user-provided mask is specified.')
        form.addParam('targetResol', FloatParam, default=10,
                       label='Target Resolution (A)',
                       help='In order to save time, you could rescale both '
                            'particles and maps to a pisel size = resol/2. '
                            'If set to 0, no rescale will be applied to the '
                            'initial references.')

        group = form.addGroup('Reference 3D map')

        group.addParam('inputVolumes', PointerParam,
                       pointerClass='SetOfVolumes',
                       important=True,
                       label='Input volumes',
                       help='Initial reference 3D map, it should have the same '
                            'dimensions and the same pixel size as your input '
                            'particles.')
        group.addParam('isMapAbsoluteGreyScale', BooleanParam, default=False,
                       label="Is initial 3D map on absolute greyscale?",
                       help='The probabilities are based on squared differences,'
                            ' so that the absolute grey scale is important. \n'
                            'Probabilities are calculated based on a Gaussian '
                            'noise model, which contains a squared difference '
                            'term between the reference and the experimental '
                            'image. This has a consequence that the reference '
                            'needs to be on the same absolute intensity '
                            'grey-scale as the experimental images. RELION and '
                            'XMIPP reconstruct maps at their absolute '
                            'intensity grey-scale. Other packages may perform '
                            'internal normalisations of the reference density, '
                            'which will result in incorrect grey-scales. '
                            'Therefore: if the map was reconstructed in RELION '
                            'or in XMIPP, set this option to Yes, otherwise '
                            'set it to No. If set to No, RELION will use a ('
                            'grey-scale invariant) cross-correlation criterion '
                            'in the first iteration, and prior to the second '
                            'iteration the map will be filtered again using '
                            'the initial low-pass filter. This procedure is '
                            'relatively quick and typically does not '
                            'negatively affect the outcome of the subsequent '
                            'MAP refinement. Therefore, if in doubt it is '
                            'recommended to set this option to No.')

        self.addSymmetry(group)

        group.addParam('initialLowPassFilterA', FloatParam, default=40,
                       label='Initial low-pass filter (A)',
                       help='It is recommended to strongly low-pass filter your '
                            'initial reference map. If it has not yet been '
                            'low-pass filtered, it may be done internally using '
                            'this option. If set to 0, no low-pass filter will '
                            'be applied to the initial reference(s).')

        form.addSection('CTF')
        form.addParam('doCTF', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label='Do CTF-correction?',
                      help='If set to Yes, CTFs will be corrected inside the '
                           'MAP refinement. The resulting algorithm '
                           'intrinsically implements the optimal linear, or '
                           'Wiener filter. Note that input particles should '
                           'contains CTF parameters.')
        form.addParam('hasReferenceCTFCorrected', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Has reference been CTF-corrected?',
                      help='Set this option to Yes if the reference map '
                           'represents CTF-unaffected density, e.g. it was '
                           'created using Wiener filtering inside RELION or '
                           'from a PDB. If set to No, then in the first '
                           'iteration, the Fourier transforms of the reference '
                           'projections are not multiplied by the CTFs.')
        form.addParam('haveDataBeenPhaseFlipped', LabelParam,
                      expertLevel=LEVEL_ADVANCED,
                      label='Have data been phase-flipped?      '
                            '(Don\'t answer, see help)',
                      help='The phase-flip status is recorded and managed by '
                           'Scipion. \n In other words, when you import or '
                           'extract particles, \nScipion will record whether '
                           'or not phase flipping has been done.\n\n'
                           'Note that CTF-phase flipping is NOT a necessary '
                           'pre-processing step \nfor MAP-refinement in '
                           'RELION, as this can be done inside the internal\n'
                           'CTF-correction. However, if the phases have been '
                           'flipped, the program will handle it.')
        form.addParam('ignoreCTFUntilFirstPeak', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Ignore CTFs until first peak?',
                      help='If set to Yes, then CTF-amplitude correction will '
                           'only be performed from the first peak '
                           'of each CTF onward. This can be useful if the CTF '
                           'model is inadequate at the lowest resolution. '
                           'Still, in general using higher amplitude contrast '
                           'on the CTFs (e.g. 10-20%) often yields better '
                           'results. Therefore, this option is not generally '
                           'recommended.')

        form.addSection(label='Optimisation')
        form.addParam('numberOfIterations', IntParam, default=25,
                      expertLevel=LEVEL_ADVANCED,
                      label='Number of iterations',
                      help='Number of iterations to be performed. Note '
                           'that the current implementation does NOT '
                           'comprise a convergence criterium. Therefore, '
                           'the calculations will need to be stopped '
                           'by the user if further iterations do not yield '
                           'improvements in resolution or classes.')
        form.addParam('regularisationParamT', IntParam, default=2,
                      expertLevel=LEVEL_ADVANCED,
                      label='Regularisation parameter T',
                      help='Bayes law strictly determines the relative '
                           'weight between the contribution of the '
                           'experimental data and the prior. '
                           'However, in practice one may need to adjust '
                           'this weight to put slightly more weight on the '
                           'experimental data to allow optimal results. '
                           'Values greater than 1 for this regularisation '
                           'parameter (T in the JMB2011 paper) put more '
                           'weight on the experimental data. Values around '
                           '2-4 have been observed to be useful for 3D '
                           'refinements, values of 1-2 for 2D refinements. '
                           'Too small values yield too-low resolution '
                           'structures; too high values result in '
                           'over-estimated resolutions and overfitting.')
        form.addParam('maskZero', EnumParam, default=0,
                      expertLevel=LEVEL_ADVANCED,
                      choices=['Yes, fill with zeros',
                               'No, fill with random noise'],
                      label='Mask particles with zeros?',
                      help='If set to <Yes>, then in the individual particles, '
                           'the area outside a circle with the radius '
                           'of the particle will be set to zeros prior to '
                           'taking the Fourier transform. '
                           'This will remove noise and therefore increase '
                           'sensitivity in the alignment and classification. '
                           'However, it will also introduce correlations '
                           'between the Fourier components that are not '
                           'modelled. When set to <No>, then the solvent area '
                           'is filled with random noise, which prevents '
                           'introducing correlations.High-resolution '
                           'refinements (e.g. in 3D auto-refine) tend to work '
                           'better when filling the solvent area with random '
                           'noise, some classifications go better when using '
                           'zeros.')
        form.addParam('referenceMask', PointerParam,
                      pointerClass='VolumeMask', expertLevel=LEVEL_ADVANCED,
                      label='Reference mask (optional)', allowsNull=True,
                      help='A volume mask containing a (soft) mask with '
                           'the same dimensions as the reference(s), '
                           'and values between 0 and 1, with 1 being 100% '
                           'protein and 0 being 100% solvent. The '
                           'reconstructed reference map will be multiplied '
                           'by this mask. If no mask is given, a soft '
                           'spherical mask based on the <radius> of the '
                           'mask for the experimental images will be '
                           'applied.\n\n'
                           'In some cases, for example for non-empty '
                           'icosahedral viruses, it is also useful to use '
                           'a second mask. Check _Advaced_ level and '
                           'select another volume mask')
        form.addParam('solventMask', PointerParam,
                      pointerClass='VolumeMask',
                      expertLevel=LEVEL_ADVANCED, allowsNull=True,
                      label='Second reference mask (optional)',
                      help='For all white (value 1) pixels in this second '
                           'mask the corresponding pixels in the '
                           'reconstructed map are set to the average value '
                           'of these pixels. Thereby, for example, the '
                           'higher density inside the virion may be set to '
                           'a constant. Note that this second mask should '
                           'have one-values inside the virion and '
                           'zero-values in the capsid and the solvent '
                           'areas.')
        form.addParam('solventFscMask', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Use solvent-flattened FSCs?',
                      help='If set to Yes, then instead of using '
                           'unmasked maps to calculate the gold-standard '
                           'FSCs during refinement, masked half-maps '
                           'are used and a post-processing-like '
                           'correction of the FSC curves (with '
                           'phase-randomisation) is performed every '
                           'iteration. This only works when a reference '
                           'mask is provided on the I/O tab. This may '
                           'yield higher-resolution maps, especially '
                           'when the mask contains only a relatively '
                           'small volume inside the box.')
        form.addParam('limitResolEStep', FloatParam, default=-1,
                      expertLevel=em.LEVEL_ADVANCED,
                      label='Limit resolution E-step to (A)',
                      help='If set to a positive number, then the '
                           'expectation step (i.e. the alignment) will be '
                           'done only including the Fourier components up '
                           'to this resolution (in Angstroms). This is '
                           'useful to prevent overfitting, as the '
                           'classification runs in RELION are not to be '
                           'guaranteed to be 100% overfitting-free (unlike '
                           'the _3D auto-refine_ with its gold-standard '
                           'FSC). In particular for very difficult data '
                           'sets, e.g. of very small or featureless '
                           'particles, this has been shown to give much '
                           'better class averages. In such cases, values '
                           'in the range of 7-12 Angstroms have proven '
                           'useful.')

        form.addSection('Sampling')
        form.addParam('angularSamplingDeg', EnumParam, default=1,
                      choices=ANGULAR_SAMPLING_LIST,
                      expertLevel=em.LEVEL_ADVANCED,
                      label='Angular sampling interval (deg)',
                      help='There are only a few discrete angular samplings'
                           ' possible because we use the HealPix library to'
                           ' generate the sampling of the first two Euler '
                           'angles on the sphere. The samplings are '
                           'approximate numbers and vary slightly over '
                           'the sphere.')

        form.addSection('Additional')
        form.addParam('useParallelDisk', BooleanParam, default=True,
                      label='Use parallel disc I/O?',
                      help='If set to Yes, all MPI slaves will read '
                           'their own images from disc. Otherwise, only '
                           'the master will read images and send them '
                           'through the network to the slaves. Parallel '
                           'file systems like gluster of fhgfs are good '
                           'at parallel disc I/O. NFS may break with many '
                           'slaves reading in parallel.')
        form.addParam('pooledParticles', IntParam, default=3,
                      label='Number of pooled particles:',
                      help='Particles are processed in individual batches '
                           'by MPI slaves. During each batch, a stack of '
                           'particle images is only opened and closed '
                           'once to improve disk access times. All '
                           'particle images of a single batch are read '
                           'into memory together. The size of these '
                           'batches is at least one particle per thread '
                           'used. The nr_pooled_particles parameter '
                           'controls how many particles are read together '
                           'for each thread. If it is set to 3 and one '
                           'uses 8 threads, batches of 3x8=24 particles '
                           'will be read together. This may improve '
                           'performance on systems where disk access, and '
                           'particularly metadata handling of disk '
                           'access, is a problem. It has a modest cost of '
                           'increased RAM usage.')
        form.addParam('allParticlesRam', BooleanParam, default=False,
                      label='Pre-read all particles into RAM?',
                      help='If set to Yes, all particle images will be '
                           'read into computer memory, which will greatly '
                           'speed up calculations on systems with slow '
                           'disk access. However, one should of course be '
                           'careful with the amount of RAM available. '
                           'Because particles are read in '
                           'float-precision, it will take \n'
                           '( N * (box_size)^2 * 4 / (1024 * 1024 '
                           '* 1024) ) Giga-bytes to read N particles into '
                           'RAM. For 100 thousand 200x200 images, that '
                           'becomes 15Gb, or 60 Gb for the same number of '
                           '400x400 particles. Remember that running a '
                           'single MPI slave on each node that runs as '
                           'many threads as available cores will have '
                           'access to all available RAM.\n\n'
                           'If parallel disc I/O is set to No, then only '
                           'the master reads all particles into RAM and '
                           'sends those particles through the network to '
                           'the MPI slaves during the refinement '
                           'iterations.')
        form.addParam('scratchDir', PathParam,
                      condition='not allParticlesRam',
                      label='Copy particles to scratch directory: ',
                      help='If a directory is provided here, then the job '
                           'will create a sub-directory in it called '
                           'relion_volatile. If that relion_volatile '
                           'directory already exists, it will be wiped. '
                           'Then, the program will copy all input '
                           'particles into a large stack inside the '
                           'relion_volatile subdirectory. Provided this '
                           'directory is on a fast local drive (e.g. an '
                           'SSD drive), processing in all the iterations '
                           'will be faster. If the job finishes '
                           'correctly, the relion_volatile directory will '
                           'be wiped. If the job crashes, you may want to '
                           'remove it yourself.')
        form.addParam('combineItersDisc', BooleanParam, default=False,
                      label='Combine iterations through disc?',
                      help='If set to Yes, at the end of every iteration '
                           'all MPI slaves will write out a large file '
                           'with their accumulated results. The MPI '
                           'master will read in all these files, combine '
                           'them all, and write out a new file with the '
                           'combined results. All MPI salves will then '
                           'read in the combined results. This reduces '
                           'heavy load on the network, but increases load '
                           'on the disc I/O. This will affect the time it '
                           'takes between the progress-bar in the '
                           'expectation step reaching its end (the mouse '
                           'gets to the cheese) and the start of the '
                           'ensuing maximisation step. It will depend on '
                           'your system setup which is most efficient.')
        form.addParam('doGpu', BooleanParam, default=True,
                      label='Use GPU acceleration?',
                      help='If set to Yes, the job will try to use GPU '
                           'acceleration.')
        form.addParam('gpusToUse', StringParam, default='',
                      label='Which GPUs to use:', condition='doGpu',
                      help='This argument is not necessary. If left empty, '
                           'the job itself will try to allocate available '
                           'GPU resources. You can override the default '
                           'allocation by providing a list of which GPUs '
                           '(0,1,2,3, etc) to use. MPI-processes are '
                           'separated by ":", threads by ",". '
                           'For example: "0,0:1,1:0,0:1,1"')

        form.addParam('oversampling', IntParam, default=1,
                      expertLevel=LEVEL_ADVANCED,
                      label="Over-sampling",
                      help="Adaptive oversampling order to speed-up "
                           "calculations (0=no oversampling, 1=2x, 2=4x, etc)")

        form.addParallelSection(threads=1, mpi=3)

    def addSymmetry(self, container):
        container.addParam('symmetryGroup', StringParam, default='c1',
                           label="Symmetry",
                           help='If the molecule is asymmetric, set Symmetry '
                                'group to C1. Note their are multiple '
                                'possibilities for icosahedral symmetry:\n'
                                '* I1: No-Crowther 222 (standard in Heymann,'
                                'Chagoyen  & Belnap, JSB, 151 (2005) 196-207)\n'
                                '* I2: Crowther 222                          \n'
                                '* I3: 52-setting (as used in SPIDER?)       \n'
                                '* I4: A different 52 setting                \n'
                                'The command *relion_refine --sym D2 '
                                '--print_symmetry_ops* prints a list of all '
                                'symmetry operators for symmetry group D2. '
                                'RELION uses MIPP\'s libraries for symmetry '
                                'operations.  Therefore, look at the XMIPP '
                                'Wiki for more details:\n'
                                ' http://xmipp.cnb.csic.es/twiki/bin/view/'
                                'Xmipp/WebHome?topic=Symmetry')

    # -------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self._initialize()
        partsId = self._getInputParticles().getObjId()
        volsId = self.inputVolumes.get().getObjId()
        self._insertFunctionStep('convertInputStep', partsId, volsId,
                                 self.targetResol.get())
        self._insertClassifyStep()
        self._insertFunctionStep('createOutputStep')

    def _insertClassifyStep(self):
        """ Prepare the command line arguments before calling Relion. """
        # Join in a single line all key, value pairs of the args dict
        args = {}

        self._setNormalArgs(args)
        self._setComputeArgs(args)

        params = ' '.join(['%s %s' % (k, str(v)) for k, v in args.iteritems()])

        self._insertFunctionStep('runClassifyStep', params)

    # -------------------------- STEPS functions -------------------------------
    def convertInputStep(self, particlesId, volumesId, tgResol):
        """ Create the input file in STAR format as expected by Relion.
        If the input particles comes from Relion, just link the file.
        Params:
            particlesId, volumesId: use this parameters just to force redo of
            convert if either the input particles and/or input volumes are
            changed.
        """
        self._imgFnList = []
        imgSet = self._getInputParticles()
        imgStar = self._getFileName('input_star')

        subset = em.SetOfParticles(filename=":memory:")

        newIndex = 1
        for img in imgSet.iterItems(orderBy='RANDOM()', direction='ASC'):
            self._scaleImages(newIndex, img)
            newIndex += 1
            subset.append(img)
            subsetSize = self.subsetSize.get()
            minSize = min(subsetSize, imgSet.getSize())
            if subsetSize   > 0 and subset.getSize() == minSize:
                break
        writeSetOfParticles(subset, imgStar, self._getExtraPath(),
                            alignType=em.ALIGN_NONE,
                            postprocessImageRow=self._postprocessParticleRow)
        self._convertInput(subset)
        self._convertRef()

    def runClassifyStep(self, params):
        """ Execute the relion steps with the give params. """
        params += ' --j %d' % self.numberOfThreads.get()
        self.runJob(self._getProgram(), params)

    def createOutputStep(self):
        # create a SetOfVolumes and define its relations
        volumes = self._createSetOfVolumes()
        self._fillVolSetFromIter(volumes, self._lastIter())

        self._defineOutputs(outputVolumes=volumes)
        self._defineSourceRelation(self.inputVolumes, volumes)

    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        self.validatePackageVersion('RELION_HOME', errors)

        if self._getInputParticles().isOddX():
            errors.append("Relion only works with even values for the "
                          "image dimensions!")

        errors += self._validateNormal()

        return errors

    def _validateNormal(self):
        """ Should be overwritten in subclasses to
        return summary message for NORMAL EXECUTION.
        """
        return []

    def _citations(self):
        cites = []
        return cites

    def _summary(self):
        self._initialize()

        lastIter = self._lastIter()

        if lastIter is not None:
            iterMsg = 'Iteration %d' % lastIter
            if self.hasAttribute('numberOfIterations'):
                iterMsg += '/%d' % self._getnumberOfIters()
        else:
            iterMsg = 'No iteration finished yet.'
        summary = [iterMsg]

        if self._getInputParticles().isPhaseFlipped():
            flipMsg = "Your input images are ctf-phase flipped"
            summary.append(flipMsg)

        summary += self._summaryNormal()
        return summary

    def _summaryNormal(self):
        """ Should be overwritten in subclasses to
        return summary message for NORMAL EXECUTION.
        """
        return []



    def _methods(self):
        """ Should be overwritten in each protocol.
        """
        return []

    # -------------------------- UTILS functions ------------------------------
    def _setNormalArgs(self, args):
        maskDiameter = self.maskDiameterA.get()
        if maskDiameter <= 0:
            x, _, _ = self._getInputParticles().getDim()
            maskDiameter = self._getInputParticles().getSamplingRate() * x

        args.update({'--i': self._getFileName('input_star'),
                     '--particle_diameter': maskDiameter,
                     '--angpix': self._getPixeSize(),
                     })
        self._setCTFArgs(args)

        if self.maskZero == MASK_FILL_ZERO:
            args['--zero_mask'] = ''



        args['--K'] = self.inputVolumes.get().getSize()
        if self.limitResolEStep > 0:
            args['--strict_highres_exp'] = self.limitResolEStep.get()

        args['--firstiter_cc'] = ''
        args['--ini_high'] = self.initialLowPassFilterA.get()
        args['--sym'] = self.symmetryGroup.get()

        refArg = self._getRefArg()
        if refArg:
            args['--ref'] = refArg

        self._setBasicArgs(args)

    def _getScratchDir(self):
        """ Returns the scratch dir value without spaces.
         If none, the empty string will be returned.
        """
        scratchDir = self.scratchDir.get() or ''
        return scratchDir.strip()

    def _setComputeArgs(self, args):
        if not self.combineItersDisc:
            args['--dont_combine_weights_via_disc'] = ''

        if not self.useParallelDisk:
            args['--no_parallel_disc_io'] = ''

        if self.allParticlesRam:
            args['--preread_images'] = ''
        else:
            if self._getScratchDir():
                args['--scratch_dir'] = self._getScratchDir()

        args['--pool'] = self.pooledParticles.get()

        if self.doGpu:
            args['--gpu'] = self.gpusToUse.get()

    def _getSamplingFactor(self):
        return 1 if self.oversampling == 0 else 2 * self.oversampling.get()

    def _setBasicArgs(self, args):
        """ Return a dictionary with basic arguments. """
        args.update({'--flatten_solvent': '',
                     '--norm': '',
                     '--scale': '',
                     '--o': self._getExtraPath('relion'),
                     '--oversampling': self.oversampling.get()
                     })

        args['--tau2_fudge'] = self.regularisationParamT.get()
        args['--iter'] = self._getnumberOfIters()

        self._setSamplingArgs(args)
        self._setMaskArgs(args)

    def _setSamplingArgs(self, args):
        """ Set sampling related params. """
        args['--healpix_order'] = self.angularSamplingDeg.get()
        args['--offset_range'] = 5
        args['--offset_step'] = (self._getSamplingFactor())

    def _setCTFArgs(self, args):
        # CTF stuff
        if self.doCTF:
            args['--ctf'] = ''

        if self.hasReferenceCTFCorrected:
            args['--ctf_corrected_ref'] = ''

        if self._getInputParticles().isPhaseFlipped():
            args['--ctf_phase_flipped'] = ''

        if self.ignoreCTFUntilFirstPeak:
            args['--ctf_intact_first_peak'] = ''

    def _setMaskArgs(self, args):
        if self.referenceMask.hasValue():
            mask = convertMask(self.referenceMask.get(), self._getTmpPath())
            args['--solvent_mask'] = mask

        if self.solventMask.hasValue():
            solventMask = convertMask(self.solventMask.get(),
                                      self._getTmpPath())
            args['--solvent_mask2'] = solventMask

        if (isVersion2() and self.referenceMask.hasValue() and
                self.solventFscMask):
            args['--solvent_correct_fsc'] = ''

    def _getProgram(self, program='relion_refine'):
        """ Get the program name depending on the MPI use or not. """
        if self.numberOfMpi > 1:
            program += '_mpi'
        return program

    def _getInputParticles(self):
        return self.inputParticles.get()

    def _getIterNumber(self, index):
        """ Return the list of iteration files, give the iterTemplate. """
        result = None
        files = sorted(glob(self._iterTemplate))
        if files:
            f = files[index]
            s = self._iterRegex.search(f)
            if s:
                result = long(s.group(1))  # group 1 is 3 digits iteration
                # number
        return result

    def _lastIter(self):
        return self._getIterNumber(-1)

    def _firstIter(self):
        return self._getIterNumber(0) or 1

    def _getIterVolumes(self, it, clean=False):
        """ Return a volumes .sqlite file for this iteration.
        If the file doesn't exists, it will be created by
        converting from this iteration data.star file.
        """
        sqlteVols = self._getFileName('volumes_scipion', iter=it)

        if clean:
            cleanPath(sqlteVols)

        if not exists(sqlteVols):
            volSet = self.OUTPUT_TYPE(filename=sqlteVols)
            self._fillVolSetFromIter(volSet, it)
            volSet.write()
            volSet.close()

        return sqlteVols

    def _fillVolSetFromIter(self, volSet, it):
        volSet.setSamplingRate(self._getInputParticles().getSamplingRate())
        modelStar = md.MetaData('model_classes@' +
                                self._getFileName('model', iter=it))
        for row in md.iterRows(modelStar):
            fn = row.getValue('rlnReferenceImage') + ":mrc"
            classDistrib = row.getValue('rlnClassDistribution')
            accurracyRot = row.getValue('rlnAccuracyRotations')
            accurracyTras = row.getValue('rlnAccuracyTranslations')
            resol = row.getValue('rlnEstimatedResolution')

            if classDistrib > 0:
                vol = em.Volume()
                self._invertScaleVol(fn)
                vol.setFileName(self._getOutputVolFn(fn))
                vol._rlnClassDistribution = em.Float(classDistrib)
                vol._rlnAccuracyRotations = em.Float(accurracyRot)
                vol._rlnAccuracyTranslations = em.Float(accurracyTras)
                vol._rlnEstimatedResolution = em.Float(resol)
                volSet.append(vol)

    def _splitInCTFGroups(self, imgStar):
        """ Add a new column in the image star to separate the particles
        into ctf groups """
        from convert import splitInCTFGroups
        splitInCTFGroups(imgStar,
                         self.defocusRange.get(),
                         self.numParticles.get())

    def _getnumberOfIters(self):
        return self.numberOfIterations.get()

    def _getRefArg(self):
        """ Return the filename that will be used for the --ref argument.
        The value will depend if in 2D and 3D or if input references will
        be used.
        It will return None if no --ref should be used. """
        inputObj = self.inputVolumes.get()
        if isinstance(inputObj, em.SetOfVolumes):
            # input SetOfVolumes as references
            return self._getRefStar()
        else:
            return None  # No --ref should be used at this point

    def _convertVolFn(self, inputVol):
        """ Return a new name if the inputFn is not .mrc """
        index, fn = inputVol.getLocation()
        return self._getTmpPath(replaceBaseExt(fn, '%02d.mrc' % index))

    def _convertVol(self, ih, inputVol):
        outputFn = self._convertVolFn(inputVol)

        if outputFn:
            xdim = self._getNewDim()
            img = ih.read(inputVol)
            img.scale(xdim, xdim, xdim)
            img.write(outputFn)

        return outputFn

    def _getRefStar(self):
        return self._getTmpPath("input_references.star")

    def _convertRef(self):
        ih = em.ImageHandler()

        inputObj = self.inputVolumes.get()
        row = md.Row()
        refMd = md.MetaData()
        for vol in inputObj:
            newVolFn = self._convertVol(ih, vol)
            row.setValue(md.RLN_MLMODEL_REF_IMAGE, newVolFn)
            row.addToMd(refMd)
        refMd.write(self._getRefStar())

    def _postprocessImageRow(self, img, imgRow):
        partId = img.getParticleId()
        imgRow.setValue(md.RLN_PARTICLE_ID, long(partId))
        imgRow.setValue(md.RLN_MICROGRAPH_NAME,
                        "%06d@fake_movie_%06d.mrcs"
                        % (img.getFrameId(), img.getMicId()))

    def _postprocessParticleRow(self, part, partRow):
        if part.hasAttribute('_rlnGroupName'):
            partRow.setValue(md.RLN_MLMODEL_GROUP_NAME,
                             '%s' % part.getAttributeValue('_rlnGroupName'))
        else:
            partRow.setValue(md.RLN_MLMODEL_GROUP_NAME,
                             '%s' % part.getMicId())

        ctf = part.getCTF()

        if ctf is not None and ctf.getPhaseShift():
            partRow.setValue(md.RLN_CTF_PHASESHIFT, ctf.getPhaseShift())

    def _getNewDim(self):
        tgResol = self.targetResol.get()
        partSet = self._getInputParticles()
        size = partSet.getXDim()
        nyquist = 2 * partSet.getSamplingRate()

        if tgResol > nyquist:
            newSize = long(round(size * nyquist / tgResol))
            if newSize % 2 == 1:
                newSize += 1
            return newSize
        else:
            return size

    def _getPixeSize(self):
        partSet = self._getInputParticles()
        oldSize = partSet.getXDim()
        newSize  = self._getNewDim()
        pxSize = partSet.getSamplingRate() * oldSize / newSize
        return pxSize

    def _scaleImages(self,indx, img):
        fn = img.getFileName()
        index = img.getIndex()
        newFn = self._getTmpPath('particles_subset.mrcs')
        xdim = self._getNewDim()

        ih = em.ImageHandler()
        image = ih.read((index, fn))
        image.scale(xdim, xdim)

        image.write((indx, newFn))

        img.setFileName(newFn)
        img.setIndex(indx)
        img.setSamplingRate(self._getPixeSize())

    def _convertInput(self, imgSet):
        newDim = self._getNewDim()
        bg = newDim / 2

        args = '--operate_on %s --operate_out %s --norm --bg_radius %d'

        params = args % (self._getFileName('input_star'),
                         self._getFileName('preprocess_particles_star'), bg)
        self.runJob(self._getProgram(program='relion_preprocess'), params)

        from pyworkflow.utils import moveFile

        moveFile(self._getFileName('preprocess_particles'),
                 self._getTmpPath('particles_subset.mrcs'))

    def _invertScaleVol(self, fn):
        xdim = self._getInputParticles().getXDim()
        outputFn = self._getOutputVolFn(fn)
        ih = em.ImageHandler()
        img = ih.read(fn)
        img.scale(xdim, xdim, xdim)
        img.write(outputFn)

    def _getOutputVolFn(self, fn):
        return self._getExtraPath(replaceBaseExt(fn, '_origSize.mrc'))