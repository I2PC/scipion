# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
"""
This module contains the protocol base class for Relion protocols
"""

import re
from glob import glob
from os.path import exists

from pyworkflow.protocol.params import (BooleanParam, PointerParam, FloatParam, 
                                        IntParam, EnumParam, StringParam, 
                                        LabelParam, PathParam)
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.utils.path import cleanPath

import pyworkflow.em as em
import pyworkflow.em.metadata as md
from pyworkflow.em.data import SetOfClasses3D
from pyworkflow.em.protocol import EMProtocol

from constants import ANGULAR_SAMPLING_LIST, MASK_FILL_ZERO, V1_3, V2_0
from convert import (convertBinaryVol, writeSetOfParticles, isVersion2,
                     getVersion, getImageLocation, convertMask)


class ProtRelionBase(EMProtocol):
    """ This class contains the common functions for all Relion protocols.
    In subclasses there should be little changes about how to create the command
    line and the files produced.
    
    Most of the Relion protocols, have two modes: NORMAL or CONTINUE. That's why
    some of the function have a template pattern approach to define the behaviour
    depending on the case.
    """
    IS_CLASSIFY = True
    IS_2D = False
    IS_3D_INIT = False
    OUTPUT_TYPE = SetOfClasses3D
    FILE_KEYS = ['data', 'optimiser', 'sampling'] 
    CLASS_LABEL = md.RLN_PARTICLE_CLASS
    CHANGE_LABELS = [md.RLN_OPTIMISER_CHANGES_OPTIMAL_ORIENTS, 
                     md.RLN_OPTIMISER_CHANGES_OPTIMAL_OFFSETS, 
                     md.RLN_OPTIMISER_CHANGES_OPTIMAL_CLASSES]
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
        if not self.doContinue:
            self.continueRun.set(None)
        else:
            self.referenceVolume.set(None)
    
    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        self.extraIter = self._getExtraPath('relion_it%(iter)03d_')
        myDict = {
                  'input_star': self._getPath('input_particles.star'),
                  'input_mrcs': self._getPath('input_particles.mrcs'),
                  'data_scipion': self.extraIter + 'data_scipion.sqlite',
                  'projections': self.extraIter + '%(half)sclass%(ref3d)03d_projections.sqlite',
                  'classes_scipion': self.extraIter + 'classes_scipion.sqlite',
                  'data': self.extraIter + 'data.star',
                  'model': self.extraIter + 'model.star',
                  'shiny': self._getExtraPath('shiny/shiny.star'),
                  'optimiser': self.extraIter + 'optimiser.star',
                  'angularDist_xmipp': self.extraIter + 'angularDist_xmipp.xmd',
                  'all_avgPmax_xmipp': self._getTmpPath('iterations_avgPmax_xmipp.xmd'),
                  'all_changes_xmipp': self._getTmpPath('iterations_changes_xmipp.xmd'),
                  'selected_volumes': self._getTmpPath('selected_volumes_xmipp.xmd'),
                  'movie_particles': self._getPath('movie_particles.star'),
                  'volume_shiny': self._getExtraPath('shiny/shiny_post.mrc:mrc'),
                  'volume_frame': self._getExtraPath('shiny/frame%(frame)03d_%(halve)sclass%(ref3d)03d_unfil.mrc:mrc'),
                  'guinier_frame': self._getExtraPath('shiny/frame%(frame)03d_guinier.star'),
                  'fsc_shiny': self._getExtraPath('shiny/shiny_post.star'),
                  'bfactors': self._getExtraPath('shiny/bfactors.star'),
                  'dataFinal': self._getExtraPath("relion_data.star"),
                  'modelFinal': self._getExtraPath("relion_model.star"),
                  'finalvolume': self._getExtraPath("relion_class%(ref3d)03d.mrc:mrc"),
                  'final_half1_volume': self._getExtraPath("relion_half1_class%(ref3d)03d_unfil.mrc:mrc"),
                  'final_half2_volume': self._getExtraPath("relion_half2_class%(ref3d)03d_unfil.mrc:mrc"),
                  'finalSGDvolume': self._getExtraPath("relion_it%(iter)03d_class%(ref3d)03d.mrc:mrc"),
                  'preprocess_particles': self._getPath("preprocess_particles.mrcs"),
                  'preprocess_particles_star': self._getPath("preprocess_particles.star"),
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
            myDict['%svolume' % p] = self.extraIter + p + 'class%(ref3d)03d.mrc:mrc'

        self._updateFilenamesDict(myDict)
    
    def _createIterTemplates(self):
        """ Setup the regex on how to find iterations. """
        self._iterTemplate = self._getFileName('data', iter=0).replace('000','???')
        # Iterations will be identify by _itXXX_ where XXX is the iteration number
        # and is restricted to only 3 digits.
        self._iterRegex = re.compile('_it(\d{3,3})_')
        
        
    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        self.IS_3D = not self.IS_2D
        form.addSection(label='Input')
        # Some hidden variables to be used for conditions
        form.addHidden('isClassify', BooleanParam, default=self.IS_CLASSIFY)
        form.addHidden('is2D', BooleanParam, default=self.IS_2D)
        
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
        form.addParam('copyAlignment', BooleanParam, default=True,
                      label='Consider previous alignment?',
                      help='If set to Yes, then alignment information from input'
                           ' particles will be considered.')
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
        form.addParam('numberOfClasses', IntParam, default=3,
                      condition='not doContinue and isClassify',
                      label='Number of classes:',
                      help='The number of classes (K) for a multi-reference '
                           'refinement. These classes will be made in an '
                           'unsupervised manner from a single reference by '
                           'division of the data into random subsets during '
                           'the first iteration.')
        group = form.addGroup('Reference 3D map',
                              condition='not doContinue and not is2D')
        group.addParam('referenceVolume', PointerParam, pointerClass='Volume',
                       important=True,
                       label="Input volume",
                       condition='not doContinue and not is2D',
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

        group.addParam('initialLowPassFilterA', FloatParam, default=60,
                       condition='not is2D',
                       label='Initial low-pass filter (A)',
                       help='It is recommended to strongly low-pass filter your '
                            'initial reference map. If it has not yet been '
                            'low-pass filtered, it may be done internally using '
                            'this option. If set to 0, no low-pass filter will '
                            'be applied to the initial reference(s).')
        
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
        form.addParam('hasReferenceCTFCorrected', BooleanParam, default=False,
                      condition='not is2D and not doContinue',
                      label='Has reference been CTF-corrected?',
                      help='Set this option to Yes if the reference map '
                           'represents CTF-unaffected density, e.g. it was '
                           'created using Wiener filtering inside RELION or '
                           'from a PDB. If set to No, then in the first '
                           'iteration, the Fourier transforms of the reference '
                           'projections are not multiplied by the CTFs.')
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
        form.addParam('doCtfManualGroups', BooleanParam, default=False,
                      label='Do manual grouping ctfs?',
                      condition='not doContinue',
                      help='Set this to Yes the CTFs will grouping manually.')
        form.addParam('defocusRange', FloatParam, default=1000,
                      label='defocus range for group creation (in Angstroms)',
                      condition='doCtfManualGroups and not doContinue',
                      help='Particles will be grouped by defocus.'
                      'This parameter is the bin for an histogram.'
                      'All particles assigned to a bin form a group')
        form.addParam('numParticles', FloatParam, default=10,
                      label='minimum size for defocus group',
                      condition='doCtfManualGroups and not doContinue',
                      help='If defocus group is smaller than this value, '
                           'it will be expanded until number of particles '
                           'per defocus group is reached')
        
        form.addSection(label='Optimisation')
        if self.IS_CLASSIFY:
            form.addParam('numberOfIterations', IntParam, default=25,
                          label='Number of iterations',
                          help='Number of iterations to be performed. Note '
                               'that the current implementation does NOT '
                               'comprise a convergence criterium. Therefore, '
                               'the calculations will need to be stopped '
                               'by the user if further iterations do not yield '
                               'improvements in resolution or classes. '
                               'If continue option is True, you going to do '
                               'this number of new iterations (e.g. if '
                               '*Continue from iteration* is set 3 and this '
                               'param is set 25, the final iteration of the '
                               'protocol will be the 28th.')
            # Default T is 2 for 2D but 4 for 3D in Relion GUI
            form.addParam('regularisationParamT', IntParam,
                          default=2 if self.IS_2D else 4,
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
            if isVersion2() and getVersion() != V2_0:  # version 2.1+ only
                form.addParam('doSubsets', BooleanParam, default=False,
                              condition='not doContinue',
                              expertLevel=LEVEL_ADVANCED,
                              label='Use subsets for initial updates?',
                              help='If set to True, multiple maximization updates '
                                   '(as many as defined by the _Number of subset '
                                   'updates_) will be performed during the first '
                                   'iteration(s): each time after the number of '
                                   'particles in a subset has been processed. By '
                                   'using subsets with much fewer particles than '
                                   'the entire data set, the initial updates '
                                   'will be much faster, while the very low '
                                   'resolution class averages will not be '
                                   'notably worse than with the entire data set. '
                                   '\nThis will greatly speed up 2D '
                                   'classifications with very many (hundreds of '
                                   'thousands or more) particles. A useful '
                                   'subset size is probably in the order of ten '
                                   'thousand particles. If the data set only '
                                   'comprises (tens of) thousands of particles, '
                                   'this option may be less useful.')
                form.addParam('subsetSize', IntParam, default=10000,
                              condition='doSubsets and not doContinue',
                              expertLevel=LEVEL_ADVANCED,
                              label='Initial subset size',
                              help='Number of individual particles after which one '
                                   'will perform a maximization update in the first '
                                   'iteration(s). A useful subset size is probably '
                                   'in the order of ten thousand particles.')
                form.addParam('subsetUpdates', IntParam, default=3,
                              condition='doSubsets and not doContinue',
                              expertLevel=LEVEL_ADVANCED,
                              label='Number of subset updates',
                              help='This option is only used when a positive '
                                   'number is given for the _Initial subset size_. '
                                   'In that case, in the first iteration, '
                                   'maximization updates are performed over '
                                   'a smaller subset of the particles to speed '
                                   'up calculations.Useful values are probably in '
                                   'the range of 2-5 subset updates. Using more '
                                   'might speed up further, but with the risk of '
                                   'affecting the results. If the number of subsets '
                                   'times the subset size is larger than the number '
                                   'of particles in the data set, then more than 1 '
                                   'iteration will be split into subsets.')
        form.addParam('maskZero', EnumParam, default=0,
                      choices=['Yes, fill with zeros',
                               'No, fill with random noise'],
                      label='Mask particles with zeros?',
                      condition='not doContinue',
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
        if self.IS_3D:
            form.addParam('referenceMask', PointerParam,
                          pointerClass='VolumeMask',
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
            form.addParam('solventMask', PointerParam, pointerClass='VolumeMask',
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
            if isVersion2():
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
        else:
            form.addParam('referenceMask', PointerParam, pointerClass='Mask',
                          label='Reference mask (optional)', allowsNull=True,
                          expertLevel=LEVEL_ADVANCED,
                          help='User-provided mask for the references ('
                               'default is to use spherical mask with '
                               'particle_diameter)')
        
        if self.IS_CLASSIFY:
            form.addParam('limitResolEStep', FloatParam, default=-1,
                          label='Limit resolution E-step to (A)',
                          condition="not doContinue",
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
        
        # Change the Sampling section name depending if classify or refine 3D
        if self.IS_CLASSIFY:
            form.addSection('Sampling')
        else:
            form.addSection('Auto-Sampling')
            form.addParam('noteAutoSampling', LabelParam,
                          label='Note that initial sampling rates will be '
                                'auto-incremented!')
        form.addParam('doImageAlignment', BooleanParam, default=True,
                      label='Perform image alignment?', condition="isClassify",
                      help='If set to No, then rather than performing both alignment '
                           'and classification, only classification will be performed. '
                           'This allows the use of very focused masks.This requires '
                           'that the optimal orientations of all particles are already '
                           'calculated.')
        if self.IS_3D:
            form.addParam('angularSamplingDeg', EnumParam, default=2,
                          choices=ANGULAR_SAMPLING_LIST,
                          label='Angular sampling interval (deg)',
                          condition='not isClassify or doImageAlignment',
                          help='There are only a few discrete angular samplings'
                               ' possible because we use the HealPix library to'
                               ' generate the sampling of the first two Euler '
                               'angles on the sphere. The samplings are '
                               'approximate numbers and vary slightly over '
                               'the sphere.')
        else:
            form.addParam('inplaneAngularSamplingDeg', FloatParam, default=5,
                          label='In-plane angular sampling (deg)',
                          condition="doImageAlignment",
                          help='The sampling rate for the in-plane rotation '
                               'angle (psi) in degrees.\n'
                               'Using fine values will slow down the program. '
                               'Recommended value for\n'
                               'most 2D refinements: 5 degrees. \n\n'
                               'If auto-sampling is used, this will be the '
                               'value for the first \niteration(s) only, and '
                               'the sampling rate will be increased \n'
                               'automatically after that.')
        form.addParam('offsetSearchRangePix', FloatParam, default=5,
                      condition='not isClassify or doImageAlignment',
                      label='Offset search range (pix)',
                      help='Probabilities will be calculated only for '
                           'translations in a circle with this radius (in '
                           'pixels). The center of this circle changes at '
                           'every iteration and is placed at the optimal '
                           'translation for each image in the previous '
                           'iteration.')
        form.addParam('offsetSearchStepPix', FloatParam, default=1.0,
                      condition='not isClassify or doImageAlignment',
                      label='Offset search step (pix)',
                      help='Translations will be sampled with this step-size '
                           '(in pixels). Translational sampling is also done '
                           'using the adaptive approach. Therefore, if '
                           'adaptive=1, the translations will first be '
                           'evaluated on a 2x coarser grid.')
        if self.IS_3D: 
            if self.IS_CLASSIFY:
                form.addParam('localAngularSearch', BooleanParam, default=False,
                              condition='not is2D and doImageAlignment',
                              label='Perform local angular search?',
                              help='If set to Yes, then rather than performing '
                                   'exhaustive angular searches, local '
                                   'searches within the range given below will '
                                   'be performed. A prior Gaussian distribution'
                                   ' centered at the optimal orientation in the'
                                   ' previous iteration and with a stddev of '
                                   '1/3 of the range given below will be '
                                   'enforced.')
                form.addParam('localAngularSearchRange', FloatParam,
                              default=5.0,
                              condition='localAngularSearch',
                              label='Local angular search range',
                              help='Local angular searches will be performed '
                                   'within +/- the given amount (in degrees) '
                                   'from the optimal orientation in the '
                                   'previous iteration. A Gaussian prior (also '
                                   'see previous option) will be applied, so '
                                   'that orientations closer to the optimal '
                                   'orientation in the previous iteration will '
                                   'get higher weights than those further away.'
                              )
            else:
                form.addParam('localSearchAutoSamplingDeg', EnumParam,
                              default=4, choices=ANGULAR_SAMPLING_LIST,
                              label='Local search from auto-sampling (deg)',
                              help='In the automated procedure to increase the '
                                   'angular samplings, local angular searches '
                                   'of -6/+6 times the sampling rate will be '
                                   'used from this angular sampling rate '
                                   'onwards.')
                
                form.addSection("Movies")
                form.addParam('realignMovieFrames', BooleanParam, default=False,
                              label='Refine movie particles?',
                              help='If set to Yes, then running averages of '
                                   'the particles from individual frames of '
                                   'recorded movies will be aligned to the '
                                   '3D reference.')
                
                group = form.addGroup('Movie frames alignment',
                                      condition='realignMovieFrames and '
                                                'doContinue')
                group.addParam('inputMovieParticles', PointerParam,
                               pointerClass='SetOfMovieParticles',
                               allowsNull=True, important=True,
                               label='Input movie particles')
                group.addParam('movieAvgWindow', FloatParam, default=5,
                               label='Running average window',
                               help='The individual movie frames will be '
                                    'averaged using a running average window '
                                    'with the specified width. Use an odd '
                                    'number. The optimal value will depend on '
                                    'the SNR in the individual movie frames. '
                                    'For ribosomes, we used a value of 5, '
                                    'where each movie frame integrated '
                                    'approximately 1 electron per squared '
                                    'Angstrom.')
                group.addParam('movieStdTrans', FloatParam, default=1,
                               label='Stddev on the translations (px)',
                               help='A Gaussian prior with the specified '
                                    'standard deviation will be centered at '
                                    'the rotations determined for the '
                                    'corresponding particle where all '
                                    'movie-frames were averaged. For '
                                    'ribosomes, we used a value of 2 pixels. ')
                group.addParam('movieIncludeRotSearch', BooleanParam,
                               default=False,
                               label='Also include rotational searches?',
                               help='If set to Yes, then running averages of '
                                    'the individual frames of recorded movies '
                                    'will also be aligned rotationally. \n'
                                    'If one wants to perform particle '
                                    'polishing, then rotational alignments of '
                                    'the movie frames is NOT necessary and '
                                    'will only take more computing time.')
                group.addParam('movieStdRot', FloatParam, default=1,
                               condition='movieIncludeRotSearch',
                               label='Stddev on the rotations (deg)',
                               help='A Gaussian prior with the specified '
                                    'standard deviation will be centered at '
                                    'the rotations determined for the '
                                    'corresponding particle where all '
                                    'movie-frames were averaged. For '
                                    'ribosomes, we used a value of 1 degree')
        
        form.addSection('Additional')
        if isVersion2():
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
        else:
            form.addParam('memoryPreThreads', IntParam, default=2,
                          label='Memory per Threads',
                          help='Computer memory in Gigabytes that is '
                               'available for each thread. This will only '
                               'affect some of the warnings about required '
                               'computer memory.')
        
        joinHalves = ("--low_resol_join_halves 40 (only not continue mode)"
                      if not self.IS_CLASSIFY else "")

        form.addParam('oversampling', IntParam, default=1,
                      expertLevel=LEVEL_ADVANCED,
                      label="Over-sampling",
                      help="Adaptive oversampling order to speed-up "
                           "calculations (0=no oversampling, 1=2x, 2=4x, etc)")

        form.addParam('extraParams', StringParam,
                      default='',
                      label='Additional parameters',
                      help="In this box command-line arguments may be "
                           "provided that are not generated by the GUI. This "
                           "may be useful for testing developmental options "
                           "and/or expert use of the program, e.g: \n"
                           "--dont_combine_weights_via_disc\n"
                           "--verb 1\n"
                           "--pad 2\n" + joinHalves)
        
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

    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep('convertInputStep',
                                 self._getInputParticles().getObjId(),
                                 self.copyAlignment)
        self._insertRelionStep()
        self._insertFunctionStep('createOutputStep')
    
    def _insertRelionStep(self):
        """ Prepare the command line arguments before calling Relion. """
        # Join in a single line all key, value pairs of the args dict    
        args = {}
        
        if self.doContinue:
            self._setContinueArgs(args)
        else:
            self._setNormalArgs(args)
        if isVersion2():
            self._setComputeArgs(args)
        
        params = ' '.join(['%s %s' % (k, str(v)) for k, v in args.iteritems()])

        if self.extraParams.hasValue():
            params += ' ' + self.extraParams.get()
        
        self._insertFunctionStep('runRelionStep', params)
    
    #--------------------------- STEPS functions -------------------------------
    
    def convertInputStep(self, particlesId, copyAlignment):
        """ Create the input file in STAR format as expected by Relion.
        If the input particles comes from Relion, just link the file.
        Params:
            particlesId: use this parameters just to force redo of convert if 
                the input particles are changed.
        """
        imgSet = self._getInputParticles()
        imgStar = self._getFileName('input_star')

        self.info("Converting set from '%s' into '%s'" %
                  (imgSet.getFileName(), imgStar))

        # Pass stack file as None to avoid write the images files
        # If copyAlignmet is set to False pass alignType to ALIGN_NONE
        if copyAlignment:
            alignType = imgSet.getAlignment()
        else:
            alignType = em.ALIGN_NONE

        writeSetOfParticles(imgSet, imgStar, self._getExtraPath(),
                            alignType=alignType,
                            postprocessImageRow=self._postprocessParticleRow)
        
        if self.doCtfManualGroups:
            self._splitInCTFGroups(imgStar)
        
        if not self.IS_CLASSIFY:
            if self.realignMovieFrames:
                movieParticleSet = self.inputMovieParticles.get()
                movieFn = self._getFileName('movie_particles')
                self.info("Converting set from '%s' into '%s'" %
                          (movieParticleSet.getFileName(), movieFn))
                
                auxMovieParticles = self._createSetOfMovieParticles(suffix='tmp')
                auxMovieParticles.copyInfo(movieParticleSet)
                # Discard the movie particles that are not present in the
                # refinement set
                for movieParticle in movieParticleSet:
                    particle = imgSet[movieParticle.getParticleId()]
                    if particle is not None:
                        auxMovieParticles.append(movieParticle)
                writeSetOfParticles(auxMovieParticles,
                                    movieFn, None, originalSet=imgSet,
                                    postprocessImageRow=self._postprocessImageRow)
                mdMovies = md.MetaData(movieFn)
                continueRun = self.continueRun.get()
                continueIter = self._getContinueIter()
                mdParts = md.MetaData(continueRun._getFileName('data', iter = continueIter))

                if getVersion() == V1_3:
                    mdParts.renameColumn(md.RLN_IMAGE_NAME,
                                         md.RLN_PARTICLE_NAME)
                else:
                    mdParts.renameColumn(md.RLN_IMAGE_NAME,
                                         md.RLN_PARTICLE_ORI_NAME)
                mdParts.removeLabel(md.RLN_MICROGRAPH_NAME)
                
                mag = movieParticleSet.getAcquisition().getMagnification()
                movieSamplingRate = movieParticleSet.getSamplingRate()
                detectorPxSize = mag * movieSamplingRate / 10000
                
                mdAux = md.MetaData()
                mdMovies.fillConstant(md.RLN_CTF_DETECTOR_PIXEL_SIZE,
                                      detectorPxSize)
                mdMovies.fillConstant(md.RLN_CTF_MAGNIFICATION, mag)
                mdAux.join2(mdMovies, mdParts, md.RLN_PARTICLE_ID,
                            md.RLN_IMAGE_ID, md.INNER_JOIN)
                # set priors equal to orig. values
                mdAux.copyColumn(md.RLN_ORIENT_ORIGIN_X_PRIOR, md.RLN_ORIENT_ORIGIN_X)
                mdAux.copyColumn(md.RLN_ORIENT_ORIGIN_Y_PRIOR, md.RLN_ORIENT_ORIGIN_Y)
                mdAux.copyColumn(md.RLN_ORIENT_PSI_PRIOR, md.RLN_ORIENT_PSI)
                mdAux.copyColumn(md.RLN_ORIENT_ROT_PRIOR, md.RLN_ORIENT_ROT)
                mdAux.copyColumn(md.RLN_ORIENT_TILT_PRIOR, md.RLN_ORIENT_TILT)
                mdAux.fillConstant(md.RLN_PARTICLE_NR_FRAMES, self._getNumberOfFrames())
                if isVersion2():
                    # FIXME: set to 1 till frame averaging is implemented in xmipp
                    mdAux.fillConstant(md.RLN_PARTICLE_NR_FRAMES_AVG, 1)

                mdAux.write(movieFn, md.MD_OVERWRITE)
                cleanPath(auxMovieParticles.getFileName())
    
    def runRelionStep(self, params):
        """ Execute the relion steps with the give params. """
        params += ' --j %d' % self.numberOfThreads.get()
        self.runJob(self._getProgram(), params)
    
    def createOutputStep(self):
        pass  # should be implemented in subclasses
    
    #--------------------------- INFO functions --------------------------------

    def _validate(self):
        errors = []
        self.validatePackageVersion('RELION_HOME', errors)
        if getVersion() == V1_3 and not self.doContinue and self.copyAlignment:
            errors.append('In RELION v1.3 a new refinement always starts '
                          'from global search. You cannot consider previous'
                          'alignment unless you run in Continue mode.')

        if self.doContinue:
            continueProtocol = self.continueRun.get()
            if (continueProtocol is not None and
                continueProtocol.getObjId() == self.getObjId()):
                errors.append('In Scipion you must create a new Relion run')
                errors.append('and select the continue option rather than')
                errors.append('select continue from the same run.')
                errors.append('') # add a new line
            errors += self._validateContinue()
        else:
            if self._getInputParticles().isOddX():
                errors.append("Relion only works with even values for the "
                              "image dimensions!")
        
            errors += self._validateNormal()

        if self.IS_CLASSIFY:
            if self._doSubsets():
                total = self._getInputParticles().getSize()
                if total <= self.subsetSize.get():
                    errors.append('Subset size is bigger than the total number '
                                  'of particles!')

        return errors
    
    def _validateNormal(self):
        """ Should be overwritten in subclasses to
        return summary message for NORMAL EXECUTION. 
        """
        return []
    
    def _validateContinue(self):
        """ Should be overwritten in subclasses to
        return summary messages for CONTINUE EXECUTION.
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
        
        if self.doContinue:
            summary += self._summaryContinue()
        summary += self._summaryNormal()
        return summary
    
    def _summaryNormal(self):
        """ Should be overwritten in subclasses to
        return summary message for NORMAL EXECUTION. 
        """
        return []
    
    def _summaryContinue(self):
        """ Should be overwritten in subclasses to
        return summary messages for CONTINUE EXECUTION.
        """
        return []
    
    def _methods(self):
        """ Should be overwritten in each protocol.
        """
        return []
    
    #--------------------------- UTILS functions -------------------------------
    def _setNormalArgs(self, args):
        maskDiameter = self.maskDiameterA.get()
        if maskDiameter <= 0:
            x, _, _ = self._getInputParticles().getDim()
            maskDiameter = self._getInputParticles().getSamplingRate() * x
    
        args.update({'--i': self._getFileName('input_star'),
                     '--particle_diameter': maskDiameter,
                     '--angpix': self._getInputParticles().getSamplingRate(),
                     })
        self._setCTFArgs(args)
    
        if self.maskZero == MASK_FILL_ZERO:
            args['--zero_mask'] = ''
    
        if self.IS_CLASSIFY:
            args['--K'] = self.numberOfClasses.get()
            if self.limitResolEStep > 0:
                args['--strict_highres_exp'] = self.limitResolEStep.get()
    
        if self.IS_3D:
            if not self.IS_3D_INIT:
                args['--ref'] = convertBinaryVol(self.referenceVolume.get(),
                                                 self._getTmpPath())
                if not self.isMapAbsoluteGreyScale:
                    args['--firstiter_cc'] = ''
                args['--ini_high'] = self.initialLowPassFilterA.get()
                args['--sym'] = self.symmetryGroup.get()
        
        if not isVersion2():
            args['--memory_per_thread'] = self.memoryPreThreads.get()
    
        self._setBasicArgs(args)

    def _setContinueArgs(self, args):
        continueRun = self.continueRun.get()
        continueRun._initialize()
    
        if self.IS_CLASSIFY:
            self.copyAttributes(continueRun, 'regularisationParamT')
        self._setBasicArgs(args)
    
        continueIter = self._getContinueIter()
        args['--continue'] = continueRun._getFileName('optimiser',
                                                      iter=continueIter)

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
    
        if self.IS_CLASSIFY:
            args['--tau2_fudge'] = self.regularisationParamT.get()
            args['--iter'] = self._getnumberOfIters()

            if not self.doContinue and isVersion2() and getVersion() != V2_0:
                self._setSubsetArgs(args)
    
        self._setSamplingArgs(args)
        self._setMaskArgs(args)


    def _setCTFArgs(self, args):
        # CTF stuff
        if self.doCTF:
            args['--ctf'] = ''
    
        # this only can be true if is 3D.
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
    
        if self.IS_3D and self.solventMask.hasValue():
            solventMask = convertMask(self.solventMask.get(), self._getTmpPath())
            args['--solvent_mask2'] = solventMask

        if (isVersion2() and self.IS_3D and self.referenceMask.hasValue() and
            self.solventFscMask):
            args['--solvent_correct_fsc'] = ''

    def _setSubsetArgs(self, args):
        if self._doSubsets():
            args['--write_subsets'] = 1
            args['--subset_size'] = self.subsetSize.get()
            args['--max_subsets'] = self.subsetUpdates.get()

    def _getProgram(self, program='relion_refine'):
        """ Get the program name depending on the MPI use or not. """
        if self.numberOfMpi > 1:
            program += '_mpi'
        return program
    
    def _getInputParticles(self):
        if self.doContinue:
            self.inputParticles.set(self.continueRun.get().inputParticles.get())
        return self.inputParticles.get()
    
    def _getIterNumber(self, index):
        """ Return the list of iteration files, give the iterTemplate. """
        result = None
        files = sorted(glob(self._iterTemplate))
        if files:
            f = files[index]
            s = self._iterRegex.search(f)
            if s:
                result = int(s.group(1)) # group 1 is 3 digits iteration number
        return result
        
    def _lastIter(self):
        return self._getIterNumber(-1)

    def _firstIter(self):
        return self._getIterNumber(0) or 1
    
    def _getIterClasses(self, it, clean=False):
        """ Return a classes .sqlite file for this iteration.
        If the file doesn't exists, it will be created by 
        converting from this iteration data.star file.
        """
        data_classes = self._getFileName('classes_scipion', iter=it)
        
        if clean:
            cleanPath(data_classes)
        
        if not exists(data_classes):
            clsSet = self.OUTPUT_TYPE(filename=data_classes)
            clsSet.setImages(self.inputParticles.get())
            self._fillClassesFromIter(clsSet, it)
            clsSet.write()
            clsSet.close()

        return data_classes
    
    def _getIterData(self, it, **kwargs):
        """ Sort the it??.data.star file by the maximum likelihood. """
        data_sqlite = self._getFileName('data_scipion', iter=it)
        
        if not exists(data_sqlite):
            iterImgSet = em.SetOfParticles(filename=data_sqlite)
            iterImgSet.copyInfo(self._getInputParticles())
            self._fillDataFromIter(iterImgSet, it)
            iterImgSet.write()
            iterImgSet.close()
        
        return data_sqlite
    
    def _splitInCTFGroups(self, imgStar):
        """ Add a new column in the image star to separate the particles
        into ctf groups """
        from convert import splitInCTFGroups
        splitInCTFGroups(imgStar,
                         self.defocusRange.get(),
                         self.numParticles.get())
    
    def _getContinueIter(self):
        continueRun = self.continueRun.get()
        
        if continueRun is not None:
            continueRun._initialize()
        
        if self.doContinue:
            if self.continueIter.get() == 'last':
                continueIter = continueRun._lastIter()
            else:
                continueIter = int(self.continueIter.get())
        else:
            continueIter = 0
            
        return continueIter
    
    def _getnumberOfIters(self):
        return self._getContinueIter() + self.numberOfIterations.get()

    def _getNumberOfFrames(self):
        movieProt = self.inputMovieParticles.get().getObjParentId()
        inputMovies = self.getProject().getProtocol(int(movieProt)).inputMovies.get()
        frames = inputMovies.getFirstItem().getNumberOfFrames()

        return frames
    
    def _postprocessImageRow(self, img, imgRow):
        partId = img.getParticleId()
        imgRow.setValue(md.RLN_PARTICLE_ID, long(partId))
        imgRow.setValue(md.RLN_MICROGRAPH_NAME,
                        "%06d@fake_movie_%06d.mrcs"
                        % (img.getFrameId(), img.getMicId()))  # fix relion-2.1

    def _postprocessParticleRow(self, part, partRow):
        if part.hasAttribute('_rlnGroupName'):
            partRow.setValue(md.RLN_MLMODEL_GROUP_NAME,
                             '%s' % part.getAttributeValue('_rlnGroupName'))
        else:
            partRow.setValue(md.RLN_MLMODEL_GROUP_NAME,
                             '%s' % part.getMicId())

        ctf = part.getCTF()

        if ctf is not None:
            partRow.setValue(md.RLN_CTF_PHASESHIFT, ctf.getPhaseShift())

    def _doSubsets(self):
        # Since 'doSubsets' property is only valid for 2.1+ protocols
        # we need provide a default value for backward compatibility
        return self.getAttributeValue('doSubsets', False)
