# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology, MRC-LMB
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

import pyworkflow.em.metadata as md
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        LabelParam, IntParam,
                                        EnumParam, StringParam,
                                        BooleanParam, PathParam,
                                        LEVEL_ADVANCED)
from pyworkflow.em.data import Volume
from pyworkflow.em.protocol import ProtInitialVolume
from pyworkflow.em.packages.relion.protocol_base import ProtRelionBase
from convert import isVersion2
from constants import ANGULAR_SAMPLING_LIST


class ProtRelionInitialModel(ProtInitialVolume, ProtRelionBase):
    """    
    Generate a 3D initial model _de novo_ from 2D particles using
    Relion Stochastic Gradient Descent (SGD) algorithm.
    """
    _label = '3D initial model'
    IS_CLASSIFY = False
    IS_3D_INIT = True
    IS_2D = False
    CHANGE_LABELS = [md.RLN_OPTIMISER_CHANGES_OPTIMAL_ORIENTS,
                     md.RLN_OPTIMISER_CHANGES_OPTIMAL_OFFSETS]

    @classmethod
    def isDisabled(cls):
        return not isVersion2()

    def __init__(self, **args):
        ProtRelionBase.__init__(self, **args)

    def _initialize(self):
        """ This function is mean to be called after the
        working dir for the protocol have been set.
        (maybe after recovery from mapper)
        """
        ProtRelionBase._initialize(self)
        #self.ClassFnTemplate = '%(ref)03d@%(rootDir)s/relion_it%(iter)03d_class001.mrc'
        self.maskZero = False
        self.copyAlignment = False
        self.hasReferenceCTFCorrected = False
        self.doCtfManualGroups = False
        self.realignMovieFrames = False


    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        self.IS_3D = not self.IS_2D
        form.addSection(label='Input')
        form.addParam('doContinue', BooleanParam, default=False,
                      label='Continue from a previous run?',
                      help='If you set to *Yes*, you should select a previous'
                           'run of type *%s* class and most of the input parameters'
                           'will be taken from it.' % self.getClassName())
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      condition='not doContinue',
                      important=True,
                      label="Input particles",
                      help='Select the input images from the project.')
        form.addParam('maskDiameterA', IntParam, default=-1,
                      condition='not doContinue',
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
        form.addParam('continueRun', PointerParam,
                      pointerClass=self.getClassName(),
                      condition='doContinue', allowsNull=True,
                      label='Select previous run',
                      help='Select a previous run to continue from.')
        form.addParam('continueIter', StringParam, default='last',
                      condition='doContinue',
                      label='Continue from iteration',
                      help='Select from which iteration do you want to '
                           'continue. If you use *last*, then the last '
                           'iteration will be used. Otherwise, a valid '
                           'iteration number should be provided.')

        self.addSymmetry(form)

        form.addSection(label='CTF')
        form.addParam('continueMsg', LabelParam, default=True,
                      condition='doContinue',
                      label='CTF parameters are not available in continue mode')
        form.addParam('doCTF', BooleanParam, default=True,
                      label='Do CTF-correction?', condition='not doContinue',
                      help='If set to Yes, CTFs will be corrected inside the '
                           'MAP refinement. The resulting algorithm '
                           'intrinsically implements the optimal linear, or '
                           'Wiener filter. Note that input particles should '
                           'contains CTF parameters.')
        form.addParam('haveDataBeenPhaseFlipped', LabelParam,
                      condition='not doContinue',
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
                      condition='not doContinue',
                      help='If set to Yes, then CTF-amplitude correction will '
                           'only be performed from the first peak '
                           'of each CTF onward. This can be useful if the CTF '
                           'model is inadequate at the lowest resolution. '
                           'Still, in general using higher amplitude contrast '
                           'on the CTFs (e.g. 10-20%) often yields better '
                           'results. Therefore, this option is not generally '
                           'recommended.')

        form.addSection(label='SGD')
        form.addParam('numberOfIterations', IntParam, default=1,
                      label='Number of iterations',
                      help='Number of iterations to be performed. '
                           'Often 1 or 2 iterations with approximately '
                           'ten thousand particles, or 5-10 iterations '
                           'with several thousand particles is enough.')
        form.addParam('sgdSubsetSize', IntParam, default=200,
                      label='SGD subset size',
                      help='How many particles will be processed for each '
                           'SGD step. Often 200 seems to work well.')
        form.addParam('writeSubsets', IntParam, default=10,
                      expertLevel=LEVEL_ADVANCED,
                      label='Write-out frequency subsets',
                      help='Every how many subsets do you want to write the '
                           'model to disk. Negative value means only write '
                           'out model after entire iteration.')
        form.addParam('sgdResLimit', IntParam, default=20,
                      label='Limit resolution SGD to (A)',
                      help='If set to a positive number, then the SGD will '
                           'be done only including the Fourier components '
                           'up to this resolution (in Angstroms). This is '
                           'essential in SGD, as there is very little '
                           'regularisation, i.e. overfitting will start '
                           'to happen very quickly. Values in the range '
                           'of 15-30 Angstroms have proven useful.')
        form.addParam('sgdNoiseVar', IntParam, default=-1,
                      expertLevel=LEVEL_ADVANCED,
                      label='SGD increased noise variance half-life',
                      help='When set to a positive value, the initial '
                           'estimates of the noise variance will internally '
                           'be multiplied by 8, and then be gradually '
                           'reduced, having 50% after this many particles '
                           'have been processed. By default, this option '
                           'is switched off by setting this value to a '
                           'negative number. In some difficult cases, '
                           'switching this option on helps. In such cases, '
                           'values around 1000 have found to be useful. '
                           'Change the factor of eight with the additional '
                           'argument *--sgd_sigma2fudge_ini*')

        form.addSection('Sampling')
        form.addParam('angularSamplingDeg', EnumParam, default=2,
                      choices=ANGULAR_SAMPLING_LIST,
                      label='Angular sampling interval (deg)',
                      help='There are only a few discrete angular samplings'
                           ' possible because we use the HealPix library to'
                           ' generate the sampling of the first two Euler '
                           'angles on the sphere. The samplings are '
                           'approximate numbers and vary slightly over '
                           'the sphere.')
        form.addParam('offsetSearchRangePix', FloatParam, default=5,
                      label='Offset search range (pix)',
                      help='Probabilities will be calculated only for '
                           'translations in a circle with this radius (in '
                           'pixels). The center of this circle changes at '
                           'every iteration and is placed at the optimal '
                           'translation for each image in the previous '
                           'iteration.')
        form.addParam('offsetSearchStepPix', FloatParam, default=1.0,
                      label='Offset search step (pix)',
                      help='Translations will be sampled with this step-size '
                           '(in pixels). Translational sampling is also done '
                           'using the adaptive approach. Therefore, if '
                           'adaptive=1, the translations will first be '
                           'evaluated on a 2x coarser grid.')

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

        form.addParam('extraParams', StringParam, default='',
                      label='Additional parameters',
                      help="In this box command-line arguments may be "
                           "provided that are not generated by the GUI. This "
                           "may be useful for testing developmental options "
                           "and/or expert use of the program, e.g: \n"
                           "--dont_combine_weights_via_disc\n"
                           "--verb 1\n"
                           "--pad 2\n")

        form.addParallelSection(threads=1, mpi=3)

    def addSymmetry(self, container):
        container.addParam('symmetryGroup', StringParam, default='c1',
                           label="Symmetry",
                           help='If the molecule is asymmetric, set Symmetry group '
                                'to C1. Note their are multiple possibilities for '
                                'icosahedral symmetry: \n'
                                '* I1: No-Crowther 222 (standard in Heymann,Chagoyen '
                                '& Belnap, JSB, 151 (2005) 196-207)               \n'
                                '* I2: Crowther 222                                 \n'
                                '* I3: 52-setting (as used in SPIDER?)              \n'
                                '* I4: A different 52 setting                       \n'
                                'The command *relion_refine --sym D2 '
                                '--print_symmetry_ops* prints a list of all symmetry '
                                'operators for symmetry group D2. RELION uses '
                                'XMIPP\'s libraries for symmetry operations. '
                                'Therefore, look at the XMIPP Wiki for more details:'
                                ' http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/'
                                'WebHome?topic=Symmetry')


    #--------------------------- INSERT steps functions -------------------------------------

    #--------------------------- STEPS functions --------------------------------------------
    def createOutputStep(self):
        imgSet = self._getInputParticles()
        vol = Volume()
        fnVol = self._getExtraPath('relion_it%s_class001.mrc') % self._lastIter()
        vol.setFileName(fnVol)
        vol.setSamplingRate(imgSet.getSamplingRate())

        outImgSet = self._createSetOfParticles()
        outImgSet.copyInfo(imgSet)
        self._fillDataFromIter(outImgSet, self._lastIter())

        self._defineOutputs(outputVolume=vol)
        self._defineSourceRelation(self.inputParticles, vol)
        self._defineOutputs(outputParticles=outImgSet)
        self._defineTransformRelation(self.inputParticles, outImgSet)
    
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        """ Should be overwritten in subclasses to
        return summary message for NORMAL EXECUTION. 
        """
        errors = []
        self.validatePackageVersion('RELION_HOME', errors)
        return errors
    
    def _summary(self):
        """ Should be overwritten in subclasses to
        return summary message for NORMAL EXECUTION. 
        """
        return []
    
    #--------------------------- UTILS functions --------------------------------------------
    def _setBasicArgs(self, args):
        """ Return a dictionary with basic arguments. """
        args.update({'--o': self._getExtraPath('relion'),
                     '--oversampling': '1',
                     '--iter': self._getnumberOfIters(),
                     '--sym': self.symmetryGroup.get()
                     })

        self._setSGDArgs(args)
        self._setSamplingArgs(args)

    def _setSGDArgs(self, args):
        args['--sgd'] = ''
        args['--denovo_3dref'] = ''
        args['--subset_size'] = self.sgdSubsetSize.get()
        args['--strict_highres_sgd'] = self.sgdResLimit.get()
        args['--write_subsets'] = self.writeSubsets.get()
        args['--sgd_sigma2fudge_halflife'] = self.sgdNoiseVar.get()

    def _setSamplingArgs(self, args):
        """ Set sampling related params"""
        if not self.doContinue:
            args['--healpix_order'] = self.angularSamplingDeg.get()
            args['--offset_range'] = self.offsetSearchRangePix.get()
            args['--offset_step'] = self.offsetSearchStepPix.get() * 2

    def _fillDataFromIter(self, imgSet, iteration):
        outImgsFn = self._getFileName('data', iter=iteration)
        imgSet.setAlignmentProj()
        imgSet.copyItems(self._getInputParticles(),
                         updateItemCallback=self._createItemMatrix,
                         itemDataIterator=md.iterRows(outImgsFn, sortByLabel=md.RLN_IMAGE_ID))

    def _createItemMatrix(self, item, row):
        from pyworkflow.em.packages.relion.convert import createItemMatrix
        from pyworkflow.em import ALIGN_PROJ

        createItemMatrix(item, row, align=ALIGN_PROJ)
