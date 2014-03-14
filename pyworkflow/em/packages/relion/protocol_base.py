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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
This module contains the protocol base class for Relion protocols
"""

import xmipp
from pyworkflow.protocol.params import BooleanParam, PointerParam, IntParam
from pyworkflow.utils import environAdd
from pyworkflow.em import *
from constants import ANGULAR_SAMPLING_LIST
from convert import createRelionInputParticles, createClassesFromImages


class ProtRelionBase(EMProtocol):
    """ This class cointains the common functionalities for all Relion protocols.
    In subclasses there should be little changes about how to create the command line
    and the files produced.
    
    Most of the Relion protocols, have two modes: NORMAL or CONTINUE. That's why
    some of the function have a template pattern approach to define the behaivour
    depending on the case.
    """
    _label = '2d classify'
    IS_CLASSIFY = True
    IS_2D = False
    FILE_KEYS = ['data', 'optimiser', 'sampling'] 
    CLASS_LABEL = xmipp.MDL_REF
    CHANGE_LABELS = [xmipp.MDL_AVG_CHANGES_ORIENTATIONS, 
                 xmipp.MDL_AVG_CHANGES_OFFSETS, 
                 xmipp.MDL_AVG_CHANGES_CLASSES]
    PREFIXES = ['']
    
    def __init__(self, **args):        
        EMProtocol.__init__(self, **args)
        self.IS_3D = not self.IS_2D
        
    def _initialize(self):
        """ This function is mean to be called after the 
        working dir for the protocol have been set. (maybe after recovery from mapper)
        """
        self._createFilenameTemplates()
        self._createIterTemplates()
        
        self.ClassFnTemplate = '%(rootDir)s/relion_it%(iter)03d_class%(ref)03d.mrc:mrc'
        self.outputClasses = 'classes_ref3D.xmd'
        self.outputVols = 'volumes.xmd'
        
        
    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        self.extraIter = self._getExtraPath('relion_it%(iter)03d_')
        myDict = {
                  'input_star': self._getPath('input_particles.star'),
                  'input_mrcs': self._getPath('input_particles.mrcs'),
                  'data_sorted_xmipp': self.extraIter + 'data_sorted_xmipp.star',
                  'classes_xmipp': self.extraIter + 'classes_xmipp.star',
                  'angularDist_xmipp': self.extraIter + 'angularDist_xmipp.xmd',
                  'all_avgPmax_xmipp': self._getTmpPath('iterations_avgPmax_xmipp.xmd'),
                  'all_changes_xmipp': self._getTmpPath('iterations_changes_xmipp.xmd'),
                  'selected_volumes': self._getTmpPath('selected_volumes_xmipp.xmd')
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

        self._fnDict = myDict
    
    def _createIterTemplates(self):
        """ Setup the regex on how to find iterations. """
        self._iterTemplate = self._getFileName('data', iter=0).replace('000','???')
        # Iterations will be identify by _itXXX_ where XXX is the iteration number
        # and is restricted to only 3 digits.
        self._iterRegex = re.compile('_it(\d{3,3})_')
        
        
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')
        # Some hidden variables to be used for conditions
        form.addHidden('isClassify', BooleanParam, default=self.IS_CLASSIFY)
        form.addHidden('is2D', BooleanParam, default=self.IS_2D)
        
        form.addParam('doContinue', BooleanParam, default=False,
                      label='Continue from a previous run?',
                      help='If you set to *Yes*, you should select a previous'
                      'run of type *%s* class and most of the input parameters'
                      'will be taken from it.' % self.getClassName())
        form.addParam('inputParticles', PointerParam, pointerClass='SetOfParticles',
                      condition='not doContinue',
                      label="Input particles",  
                      help='Select the input images from the project.')   
        form.addParam('continueRun', PointerParam, pointerClass=self.getClassName(),
                      condition='doContinue', allowNull=True,
                      label='Select previous run',
                      help='Select a previous run to continue from.')
        form.addParam('continueIter', StringParam, default='last',
                      condition='doContinue', 
                      label='Continue from iteration',
                      help='Select from which iteration do you want to continue.'
                           'if you use *last*, then the last iteration will be used.'
                           'otherwise, a valid iteration number should be provided.')        
        
        form.addParam('numberOfClasses', IntParam, default=3, 
                      condition='not doContinue and isClassify',
                      label='Number of classes:',
                      help='The number of classes (K) for a multi-reference refinement.'
                           'These classes will be made in an unsupervised manner from a single'
                           'reference by division of the data into random subsets during the'
                           'first iteration.')
        form.addParam('referenceVolume', PointerParam, pointerClass='Volume',
                      condition='not doContinue and not is2D',
                      label="Initial 3D map", 
                      help='Initial reference 3D map, it should have the same '
                           'dimensions and the same pixel size as your input particles.')
        form.addParam('isMapAbsoluteGreyScale', BooleanParam, default=False,
                      condition='not doContinue and not is2D',
                      label="Is initial 3D map on absolute greyscale?", 
                      help='The probabilities are based on squared differences, '
                           'so that the absolute grey scale is important. \n'
                           'Probabilities are calculated based on a Gaussian noise model,'
                           'which contains a squared difference term between the reference and the experimental image.' 
                           'This has a consequence that the reference needs to be on the same absolute intensity '
                           'grey-scale as the experimental images. RELION and XMIPP reconstruct maps at their absolute'
                           'intensity grey-scale. Other packages may perform internal normalisations of the reference' 
                           'density, which will result in incorrect grey-scales. Therefore: if the map was reconstructed'
                           'in RELION or in XMIPP, set this option to Yes, otherwise set it to No. If set to No, RELION '
                           'will use a (grey-scale invariant) cross-correlation criterion in the first iteration, and '
                           'prior to the second iteration the map will be filtered again using the initial low-pass filter.'
                           'This procedure is relatively quick and typically does not negatively affect the outcome of the'
                           'subsequent MAP refinement. Therefore, if in doubt it is recommended to set this option to No.')        
        form.addParam('symmetryGroup', StringParam, default='c1',
                      condition='not doContinue and not is2D',
                      label="Symmetry group", 
                      help='See http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry '
                           'for a description of the symmetry groups format')  
        form.addParam('paddingFactor', FloatParam, default=3,
                      condition='isClassify', expertLevel=LEVEL_ADVANCED,
                      label='',
                      help='The padding factor used for oversampling of the Fourier transform. The default is 3x padding, '
                           'which is combined with nearest-neighbour interpolation. However, for large 3D maps, storing the '
                           'entire 3x oversampled Fourier transform (as doubles) plus the real-space padded map in memory may ' 
                           'be too much. Therefore, for large maps or in cases of many 3D references, in order to fit into memory ' 
                           'one may need to use a smaller padding factor: e.g. 2, or (not recommended) 1. For padding factors smaller '
                           'than 3, (slower) linear interpolation will be used.\n'
                           'The approximate amount of memory (in Gb) required to store K maps of (size x size x size) voxels and a' 
                           'padding factor (pad) may be calculated as: K*2*8*(size*pad)^3/1024/1024/1024\n' 
                           '<Note>: also consider the use of threads if memory is an issue.')
        form.addParam('extraParams', StringParam, default='',
                      expertLevel=LEVEL_EXPERT,
                      label='Additional parameters',
                      help='')
               
        form.addSection(label='CTF')
        form.addParam('doCTF', BooleanParam, default=True,
                      label='Do CTF-amplitud correction?',
                      help='If set to Yes, CTFs will be corrected inside the MAP refinement. '
                           'The resulting algorithm intrinsically implements the optimal linear, ' 
                           'or Wiener filter. Note that input particles should contains CTF parameters.')
        form.addParam('hasReferenceCTFCorrected', BooleanParam, default=False,
                      condition='not is2D',
                      label='Has reference been CTF-corrected?',
                      help='Set this option to Yes if the reference map represents CTF-unaffected density, '
                           'e.g. it was created using Wiener filtering inside RELION or from a PDB. If set to No, ' 
                           'then in the first iteration, the Fourier transforms of the reference projections ' 
                           'are not multiplied by the CTFs.')        
        form.addParam('onlyFlipPhases', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Only flip phases?',
                      help='Set this to Yes to switch CTF-amplitude correction off.'
                           'This option is NOT generally recommended.')          
        form.addParam('haveDataPhaseFlipped', BooleanParam, default=False,
                      label='Have data been phase-flipped?',
                      help='Set this to Yes if the images have been ctf-phase corrected during the pre-processing steps. '
                           'Note that CTF-phase flipping is NOT a necessary pre-processing step for MAP-refinement in RELION, ' 
                           'as this can be done inside the internal CTF-correction. However, if the phases do have been flipped, ' 
                           'one should tell the program about it using this option.'
                           '*TODO*: remove this property if we keep track of phase flipped')   
        form.addParam('ignoreCTFUntilFirstPeak', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Ignore CTFs until first peak?',
                      help='If set to Yes, then CTF-amplitude correction will only be performed from the first peak ' 
                           'of each CTF onward. This can be useful if the CTF model is inadequate at the lowest resolution. ' 
                           'Still, in general using higher amplitude contrast on the CTFs (e.g. 10-20%) often yields better results. '
                           'Therefore, this option is not generally recommended.')    
        form.addParam('doIntensityCorrection', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Do intensity correction?',
                      help='An internal correction for differences in the intensity (grey-scale) of the signal between '
                           'distinct micrographs is applied. This is useful if micrographs have very different ' 
                           'signal-to-noise ratios, e.g. due to different ice thickness or contamination. ' 
                           'Because one typically normalises the noise, this leads to distinct signal intensities in the data, ' 
                           'and this procedure corrects for this. It is quite robust and therefore recommended for the general case.')                                
        
        form.addSection(label='Optimisation')
        form.addParam('initialLowPassFilterA', FloatParam, default=60,
                      condition='not is2D',
                      label='Initial low-pass filter (A)',
                      help='It is recommended to strongly low-pass filter your initial reference map. '
                           'If it has not yet been low-pass filtered, it may be done internally using this option. ' 
                           'If set to 0, no low-pass filter will be applied to the initial reference(s).')
        if self.IS_CLASSIFY:
            form.addParam('numberOfIterations', IntParam, default=25,
                          label='Number of iterations',
                          help='Number of iterations to be performed. Note that the current implementation does NOT '
                               'comprise a convergence criterium. Therefore, the calculations will need to be stopped'
                               'by the user if further iterations do not yield improvements in resolution or classes.') 
            form.addParam('regularisationParamT', IntParam, default=1,
                          label='Regularisation parameter T',
                          help='Bayes law strictly determines the relative weight between the contribution of the '
                               'experimental data and the prior. '
                               'However, in practice one may need to adjust this weight to put slightly more weight on the experimental '
                               'data to allow optimal results. Values greater than 1 for this regularisation parameter '
                               '(T in the JMB2011 paper) put more weight on the experimental data. Values around 2-4 '
                               'have been observed to be useful for 3D refinements, values of 1-2 for 2D refinements. '
                               'Too small values yield too-low resolution structures; too high values result in ' 
                               'over-estimated resolutions and overfitting.') 
        form.addParam('maskRadiusA', IntParam, default=200,
                      label='Particle mask RADIUS (A)',
                      help='The experimental images will be masked with a soft circular mask '
                           'with this <radius> (Note that in Relion GUI the diameter is used). '
                           'Make sure this radius is not set too small because that may mask '
                           'away part of the signal! If set to a value larger than the image '
                           'size no masking will be performed.\n\n'
                           'The same radius will also be used for a spherical mask of the '
                           'reference structures if no user-provided mask is specified.') 
        if self.IS_CLASSIFY:
            form.addParam('maskZero', EnumParam, default=0,
                          choices=['Yes, fill with zeros', 'No, fill with random noise'],
                          label='Mask particles with zeros?',
                          help='If set to <Yes>, then in the individual particles, the area outside a circle with the radius '
                               'of the particle will be set to zeros prior to taking the Fourier transform. '
                               'This will remove noise and therefore increase sensitivity in the alignment and classification. ' 
                               'However, it will also introduce correlations between the Fourier components that are not modelled. ' 
                               'When set to <No>, then the solvent area is filled with random noise, which prevents introducing '
                               'correlations.High-resolution refinements (e.g. in 3D auto-refine) tend to work better when filling ' 
                               'the solvent area with random noise, some classifications go better when using zeros.') 
        form.addParam('referenceMask', PointerParam, pointerClass='Mask',
                      label='Reference mask (optional)', allowNull=True,
                      help='A volume mask containing a (soft) mask with the same dimensions ' 
                           'as the reference(s), and values between 0 and 1, with 1 being 100% protein '
                           'and 0 being 100% solvent. The reconstructed reference map will be multiplied '
                           'by this mask. If no mask is given, a soft spherical mask based on the <radius> '
                           'of the mask for the experimental images will be applied.\n\n'  
                           'In some cases, for example for non-empty icosahedral viruses, it is also useful ' 
                           'to use a second mask. Use <More options> and check <Solvent mask> parameter. ') 
        form.addParam('solventMask', PointerParam, pointerClass='Mask',
                      expertLevel=LEVEL_ADVANCED, allowNull=True,
                      label='Solvent mask (optional)',
                      help='For all white (value 1) pixels in this second mask the '
                           'corresponding pixels in the reconstructed map are set to the average value of '
                           'these pixels. Thereby, for example, the higher density inside the virion may be '
                           'set to a constant. Note that this second mask should have one-values inside the '
                           'virion and zero-values in the capsid and the solvent areas.') 
        
        samplingSection = 'Sampling'
        if not self.IS_CLASSIFY:
            samplingSection += ' (initial value will be autoincremented)'
        form.addSection(samplingSection)
        if self.IS_3D:
            form.addParam('angularSamplingDeg', EnumParam, default=2,
                          choices=ANGULAR_SAMPLING_LIST,
                          label='Angular sampling interval (deg)',
                          help='There are only a few discrete angular samplings possible because '
                           'we use the HealPix library to generate the sampling of the first '
                           'two Euler angles on the sphere. The samplings are approximate numbers ' 
                           'and vary slightly over the sphere.')
        else:
            form.addParam('inplaneAngularSamplingDeg', FloatParam, default=5,
                          label='Angular sampling interval (deg)',
                          help='There are only a few discrete angular samplings possible because '
                           'we use the HealPix library to generate the sampling of the first '
                           'two Euler angles on the sphere. The samplings are approximate numbers ' 
                           'and vary slightly over the sphere.')           
        form.addParam('offsetSearchRangePix', FloatParam, default=5,
                      condition='isClassify or not doContinue',
                      label='Offset search range (pix)',
                      help='Probabilities will be calculated only for translations in a circle '
                           'with this radius (in pixels). The center of this circle changes at '
                           'every iteration and is placed at the optimal translation for each '
                           'image in the previous iteration.')
        form.addParam('offsetSearchStepPix', FloatParam, default=1.0,
                      condition='isClassify or not doContinue',
                      label='Offset search step (pix)',
                      help='Translations will be sampled with this step-size (in pixels). '
                           'Translational sampling is also done using the adaptive approach. '
                           'Therefore, if adaptive=1, the translations will first be evaluated'
                           'on a 2x coarser grid.')
        if self.IS_3D: 
            if self.IS_CLASSIFY:
                form.addParam('localAngularSearch', BooleanParam, default=False,
                              condition='not is2D',
                              label='Perform local angular search?',
                              help='If set to Yes, then rather than performing exhaustive angular searches, '
                                   'local searches within the range given below will be performed. A prior '
                                   'Gaussian distribution centered at the optimal orientation in the previous '
                                   'iteration and with a stddev of 1/3 of the range given below will be enforced.')
                form.addParam('localAngularSearchRange', FloatParam, default=5.0,
                              condition='localAngularSearch',
                              label='Local angular search range',
                              help='Local angular searches will be performed within +/- \n'
                                   'the given amount (in degrees) from the optimal orientation\n'
                                   'in the previous iteration. A Gaussian prior (also see \n'
                                   'previous option) will be applied, so that orientations \n'
                                   'closer to the optimal orientation in the previous iteration \n'
                                   'will get higher weights than those further away.')
            else:
                form.addParam('localSearchAutoSamplingDeg', EnumParam, default=4,
                              choices=ANGULAR_SAMPLING_LIST,
                              label='Local search from auto-sampling (deg)',
                              help='In the automated procedure to increase the angular samplings,\n'
                                   'local angular searches of -6/+6 times the sampling rate will\n'
                                   'be used from this angular sampling rate onwards.')
                
         
    #--------------------------- INSERT steps functions --------------------------------------------  
    def _insertAllSteps(self): 
        self._initialize()
        self._insertFunctionStep('convertInputStep')
        self._insertRelionStep()            
        self._insertFunctionStep('createOutputStep', 2)
 
    def _insertRelionStep(self):
        """ Prepare the command line arguments before calling Relion. """
        # Join in a single line all key, value pairs of the args dict    
        args = {}
        if self.doContinue:
            self._setContinueArgs(args)
        else:
            self._setNormalArgs(args)
        params = ' '.join(['%s %s' % (k, str(v)) for k, v in args.iteritems()])
        
        if self.extraParams.hasValue():
            params += ' ' + self.extraParams.get()

        self._insertFunctionStep('runRelionStep', params)
        
    def _setNormalArgs(self, args):
        args.update({'--i': self._getFileName('input_star'),
                     '--particle_diameter': self.maskRadiusA.get() * 2.0,
                     '--angpix': self.inputParticles.get().getSamplingRate(),
                    })
        self._setBasicArgs(args)
        self._setMaskArgs(args)
        
        if self.IS_CLASSIFY:
            if self.maskZero == MASK_FILL_ZERO:
                args['--zero_mask'] = ''                
            args['--K'] = self.numberOfClasses.get()
        
    def _setContinueArgs(self, args):
        args = {}
        if self.IS_CLASSIFY:
            self.copyAttributes(self.continueRun.get(), 'regularisationParamT')
        self._setBasicArgs(args)
        args['--continue'] = self._getFileName('optimiser', iter=self.continueIter.get())
            
    def _setBasicArgs(self, args):
        """ Return a dictionary with basic arguments. """
        args.update({'--flatten_solvent': '',
                     '--norm': '',
                     '--scale': '',
                     '--o': self._getExtraPath('relion'),
                     '--oversampling': '1'
                    })
        if self.IS_CLASSIFY:
            args['--tau2_fudge'] = self.regularisationParamT.get()
            args['--iter'] = self.numberOfIterations.get()
            
        self._setSamplingArgs(args)
        args['--offset_range'] = self.offsetSearchRangePix.get()
        args['--offset_step']  = self.offsetSearchStepPix.get() * 2
                        
    def _setCTFArgs(self, args):        
        # CTF stuff
        if self.doCTF:
            args['--ctf'] = ''
        
        if self.hasReferenceCTFCorrected:
            args['--ctf_corrected_ref'] = ''
            
        if self.haveDataPhaseFlipped:
            args['--ctf_phase_flipped'] = ''
            
        if self.ignoreCTFUntilFirstPeak:
            args['--ctf_intact_first_peak'] = ''
            
    def _setMaskArgs(self, args):
        if self.referenceMask.hasValue():
            args['--solvent_mask'] = self.referenceMask.get().getFileName() #FIXE: CHANGE BY LOCATION
            
        if self.solventMask.hasValue():
            args['--solvent_mask2'] = self.solventMask.get().getFileName() #FIXME: CHANGE BY LOCATION
 
   
    #--------------------------- STEPS functions --------------------------------------------       
    def convertInputStep(self):
        """ Create the input file in STAR format as expected by Relion.
        If the input particles comes from Relion, just link the file. 
        """
        createRelionInputParticles(self.inputParticles.get(), 
                                   self._getFileName('input_star'),
                                   self._getFileName('input_mrcs'))
        
    def runRelionStep(self, params):
        """ Execute the relion steps with the give params. """
        self._loadEnvironment()
        params += ' --j %d' % self.numberOfThreads.get()
        self.runJob(self._getProgram(), params)
        
    def createOutputStep(self):
        pass # should be implemented in subclasses
        
    
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        errors = []
        if self.doContinue:
            errors += self._validateContinue()
        else:
            errors += self._validateNormal()
        return errors
            
    def _validateNormal(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        return []
    
    def _validateContinue(self):
        """ Should be overriden in subclasses to
        return summary messages for CONTINUE EXECUTION.
        """
        return []
    
    def _citations(self):
        cites = []
        return cites
    
    def _summary(self):
        if self.doContinue:
            return self._summaryContinue()
        return self._summaryNormal()

    def _summaryNormal(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        return []
    
    def _summaryContinue(self):
        """ Should be overriden in subclasses to
        return summary messages for CONTINUE EXECUTION.
        """
        return []


    #--------------------------- UTILS functions --------------------------------------------
    def _getProgram(self):
        """ Get the program name depending on the MPI use or not. """
        program = 'relion_refine'
        if self.numberOfMpi > 1:
            program += '_mpi'
        return program
    
    def _loadEnvironment(self):
        """ Setup the environment variables needed to launch Relion. """
        RELION_BIN = join(os.environ['RELION_HOME'], 'bin')
        RELION_LD = join(os.environ['RELION_HOME'], 'lib64')
        environAdd('PATH', RELION_BIN)
        environAdd('LD_LIBRARY_PATH', RELION_LD)
    
    def _getFileName(self, key, **args):
        """ Retrieve a filename from the templates. """
        return self._fnDict[key] % args
    
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
    
    def _getIterClasses(self, it):
        """ Return the .star file with the classes for this iteration.
        If the file doesn't exists, it will be created. 
        """
        data_star = self._getFileName('data', iter=it)
        data_classes = self._getFileName('classes_xmipp', iter=it)
        
        if not exists(data_classes):
            createClassesFromImages(data_star, data_classes, it, 
                                    self.CLASS_LABEL, self.ClassFnTemplate)
        return data_classes

    def _getIterSortedData(self, it):
        """ Sort the it??.data.star file by the maximum likelihood. """
        data_sorted = self._getFileName('data_sorted_xmipp', iter=it)
        
        if not exists(data_sorted):
            data = self._getFileName('data', iter=it)
            from convert import sortImagesByLL
            sortImagesByLL(data, data_sorted)
        
        return data_sorted