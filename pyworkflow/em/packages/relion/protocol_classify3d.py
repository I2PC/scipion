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
This module contains the protocol for 3d classification with relion.
"""


from pyworkflow.em import *  
from pyworkflow.utils.which import which
from pyworkflow.utils.path import makePath, replaceBaseExt, join, basename
from convert import createRelionInputImages, createRelionInputVolume

class ProtRelionBase(EMProtocol):
    """ Base class for all relion protocols"""
    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        self.ParamsStr = ''
        if self.numberOfMpi > 1:
            self.program = 'relion_refine_mpi'
        else:
            self.program = 'relion_refine'

    def loadEnvironment(self):
        """ Load the environment variables needed for use RELION tools. """
        RELION_DIR = os.environ['RELION_HOME']
        os.environ['PATH'] = os.environ['PATH'] + os.pathsep + join(RELION_DIR, 'bin', '')
        os.environ['LD_LIBRARY_PATH'] = os.environ['LD_LIBRARY_PATH'] + os.pathsep + join(RELION_DIR, 'lib64', '')
        os.environ['RELION_QSUB_TEMPLATE'] = join(RELION_DIR, 'bin', 'qsub.csh')

            
    def defineSteps(self): 
        
        self.loadEnvironment()
        
        self.imgStar    = createRelionInputImages(self, self.inputParticles.get())        
        self.volumeName = createRelionInputVolume(self, self.input3DReferences.get())


    def defineSteps2(self, firstIteration
                         , lastIteration
                         , NumberOfClasses):
        
        inputs=[]

        #data to store in Working dir so it can be easily accessed by users
        lastIterationVolumeFns = []
        standardOutputClassFns = []
        extraIter = self._getExtraPath('relion_it%03d_' % lastIteration)
        extraIter2 = self._getExtraPath('relion_')
        for ref3d in range (1,self.numberOfClasses.get()+1):
            try:
                lastIterationVolumeFns += [extraIter + "class%03d.mrc:mrc" % ref3d]
                #lastIterationVolumeFns += [self.getFilename('volume', iter=lastIteration, ref3d=ref3d, workingDir=self.WorkingDir )]
            except:
                lastIterationVolumeFns += [extraIter2 + "class%03d.mrc:mrc" % ref3d]
                #lastIterationVolumeFns += [self.getFilename('volumeFinal', iter=lastIteration, ref3d=ref3d, workingDir=self.WorkingDir )]
            
            standardOutputClassFns += ["images_ref3d%06d@"%ref3d + self._getPath("classes_ref3D.xmd")]
        #lastIterationMetadata = "images@"+self.getFilename('data'+'Xm', iter=lastIteration, workingDir=self.WorkingDir )
        lastIterationMetadata = "images@" + extraIter +'data.xmd'
        
        inputMetadatas = []
        for it in range (firstIteration,self.numberOfIterations.get()+1): #always list all iterations
            inputMetadatas += ["images@" + self._getExtraPath('relion_it%03d_' % it) +'data.xmd']

#        self.insertStep('convertRelionMetadata',
#                        inputs  = inputs,
#                        lastIterationVolumeFns = lastIterationVolumeFns,
#                        lastIterationMetadata  = lastIterationMetadata,
#                        NumberOfClasses        = self.numberOfClasses.get(),
#                        standardOutputClassFns = standardOutputClassFns,
#                        standardOutputImageFn  = "images@" + self._getPath("images.xmd"),
#                        standardOutputVolumeFn = "volumes@" + self._getPath("volumes.xmd"),
#                        inputMetadatas=inputMetadatas,
#                        imagesAssignedToClassFn='imagesAssignedToClass@'+self._getExtraPath('/dataForVisualize.xmd')
#                        )

                        
class Relion3DClassification(ProtClassify3D, ProtRelionBase):
    """Protocol to perform CTF estimation on a set of micrographs
    using the ctffind3 program"""
    _label = '3D classify'
    
    def __init__(self, **args):
        ProtRelionBase.__init__(self, **args)
        self.Import = 'from protocol_relion_classify import *'
        self.relionType='classify'
            
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam, label="Input images", important=True, 
                      pointerClass='SetOfParticles',
                      help='Provide a set of images from a stack (Spider/MRC) or metadata file that make up your data set.'
                      'Note that in the case of a stack, no metadata can be included and thus no CTF correction '
                      'can be performed.Nor will it be possible to perform noise spectra estimation or intensity '
                      'scale corrections in image groups.')  
        form.addParam('numberOfClasses', IntParam, default=1,
                      label='Number of classes',
                      help='The number of classes (K) for a multi-reference refinement. '
                      'These classes will be made in an unsupervised manner from a single reference '
                      'by division of the data into random subsets during the first iteration.')
        form.addParam('input3DReferences', PointerParam,label='Initial 3D reference volume',important=True, 
                      pointerClass='SetOfVolumes',
                      help='Initial 3D density map with the same dimensions as your particles.')
        form.addParam('isRefMapGreyScale', BooleanParam, default=True,
                      label="Ref. Map is on absolute greyscale?", 
                      help=' The probabilities are based on squared differences, so that the absolute grey scale is important.'
                      'Probabilities are calculated based on a Gaussian noise model,'
                      'which contains a squared difference term between the reference and the experimental image.'
                      'This has a consequence that the reference needs to be on the same absolute intensity '
                      'grey-scale as the experimental images. RELION and XMIPP reconstruct maps at their absolute '
                      'intensity grey-scale. Other packages may perform internal normalisations of the reference '
                      'density, which will result in incorrect grey-scales. Therefore: if the map was reconstructed '
                      'in RELION or in XMIPP, set this option to Yes, otherwise set it to No. If set to No, RELION '
                      'will use a (grey-scale invariant) cross-correlation criterion in the first iteration, and '
                      'prior to the second iteration the map will be filtered again using the initial low-pass filter.'
                      'This procedure is relatively quick and typically does not negatively affect the outcome of the '
                      'subsequent MAP refinement. Therefore, if in doubt it is recommended to set this option to No.')
#        form.addParam('doNormalize', BooleanParam, default=False,
#                      label="Normalize input images?", 
#                      help='Average background value must be 0 and a stddev value must be 1.'
#                      'Note that the average and stddev values for the background are '
#                      'calculated outside a circle with the particle diameter.')
        form.addParam('symmetryGroup', TextParam, default='c1',
                      label='Symmetry group',
                      help='See http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry '
                      'for a description of the symmetry groups format.'
                      'If no symmetry is present, give c1.')   
        form.addSection(label='CTF', questionParam='doCtf')        
        form.addParam('doCtf', BooleanParam, default=False,
                      label='Use CTF-amplitude correction?')
        form.addParam('isReferenceCtfCorrected', BooleanParam, default=False,
                      label='Has reference been CTF-corrected?',
                      help='Set this option to Yes if the reference map represents CTF-unaffected density,'
                      'e.g. it was created using Wiener filtering inside RELION or from a PDB. If set to No,'
                      'then in the first iteration, the Fourier transforms of the reference projections '
                      'are not multiplied by the CTFs.')   
        form.addParam('onlyFlipPhases', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Only flip phases?',
                      help='Set this to Yes to switch CTF-amplitude correction off.'
                      'This option is NOT generally recommended.')   
        form.addParam('areImagesPhaseFlipped', BooleanParam, default=False,
                      label='Are the images CTF phase flipped?',
                      help='Set this to Yes if the images have been ctf-phase corrected during the pre-processing steps.'
                      'Note that CTF-phase flipping is NOT a necessary pre-processing step for MAP-refinement in RELION,'
                      'as this can be done inside the internal CTF-correction. However, if the phases do have been flipped,'
                      'one should tell the program about it using this option.')     
        form.addParam('ignoreFirstPeakCTF', BooleanParam, default=False,
                      label='Ignore CTFs until first peak?',
                      help='If set to Yes, then CTF-amplitude correction will only be performed from the first peak '
                      'of each CTF onward. This can be useful if the CTF model is inadequate at the lowest resolution.'
                      'Still, in general using higher amplitude contrast on the CTFs (e.g. 10-20%) often yields better results.'
                      'Therefore, this option is not generally recommended.')   
        form.addParam('doIntensityCorrection', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Do intensity correction?',
                      help='An internal correction for differences in the intensity (grey-scale) of the signal between '
                      'distinct micrographs is applied. This is useful if micrographs have very different '
                      'signal-to-noise ratios, e.g. due to different ice thickness or contamination.'
                      'Because one typically normalises the noise, this leads to distinct signal intensities in the data,'
                      'and this procedure corrects for this. It is quite robust and therefore recommended for the general case.')   
        form.addSection(label='Optimisation')
        form.addParam('iniLowPassFilter', FloatParam, default=60.0,
                      label='Initial low-pass filter (in Angstroms)',
                      help='It is recommended to strongly low-pass filter your initial reference map.'
                      'If it has not yet been low-pass filtered, it may be done internally using this option.'
                      'If set to 0, no low-pass filter will be applied to the initial reference(s).')    
        form.addParam('numberOfIterations', IntParam, default=25,
                      label='Number of iterations',
                      help='Number of iterations to be performed. Note that the current implementation does NOT comprise a convergence criterium. '
                      'Therefore, the calculations will need to be stopped by the user if further iterations do not yield improvements in resolution or classes.')    
        form.addParam('regularParT', IntParam, default=1,
                      label='Regularisation parameter T',
                      help='Bayes law strictly determines the relative weight between the contribution of the experimental data and the prior.'
                      'However, in practice one may need to adjust this weight to put slightly more weight on the experimental '
                      'data to allow optimal results. Values greater than 1 for this regularisation parameter '
                      '(T in the JMB2011 paper) put more weight on the experimental data. Values around 2-4 '
                      'have been observed to be useful for 3D refinements, values of 1-2 for 2D refinements. '
                      'Too small values yield too-low resolution structures; too high values result in over-estimated resolutions and overfitting.')
        form.addParam('maskRadius', IntParam, default=200,
                      label='Particles mask RADIUS (in Angstroms)',
                      help='The experimental images will be masked with a soft circular mask with this diameter. '
                      'Make sure this radius is not set too small because that may mask away part of the signal! '
                      'If set to a value larger than the image size no masking will be performed.')    
        form.addParam('doMaskReferences', BooleanParam, default=True, 
                      label='Mask references structures?',
                      help='If set to yes, a mask will also be applied to the reconstructed references. '
                      'This is useful to set the solvent region of your reconstruction to 0. Either a soft spherical mask '
                      '(based on the diameter of the experimental image mask given above) or a user-provided mask (next option) may be used. '
                      'The user-provided mask should have values between 0 and 1 only. Solvent flattening is recommended, but make sure not to mask any signal away.')
        #TODO: the following parameter should be a Pointer to an Image and not a string containing the path  
        form.addParam('maskFile', StringParam, default='', condition='doMaskReferences',
                      label='Reference mask file', 
                      help='A Spider/mrc map containing a (soft) mask with the same dimensions as the reference(s), and values between 0 and 1, '
                      'with 1 being 100% protein and 0 being 100% solvent. The reconstructed reference map will be multiplied by this mask.'
                      'If no mask is given, a soft spherical mask based on the diameter of the mask for the experimental images will be applied.'
                      'In some cases, for example for non-empty icosahedral viruses, it is also useful to use a second mask. '
                      'For all white (value 1) pixels in this second mask the corresponding pixels in the reconstructed map are set to the average value of these pixels. '
                      'Thereby, for example, the higher density inside the virion may be set to a constant. '
                      'Note that this second mask should have one-values inside the virion and zero-values in the capsid and the solvent areas. '
                      'To use a second mask, use the additional option --solvent_mask2, which may given in the Additional arguments line (in the Sampling tab).')     
        form.addSection(label='Sampling')
        form.addParam('angSampling', EnumParam, choices=['30', '15', '7.5', '3.7', '1.8', '0.9', '0.5', '0.2', '0.1'], 
                      default=2, label='Angular sampling interval (deg)', display=EnumParam.DISPLAY_COMBO,
                      help='There are only a few discrete angular samplings possible because we use the HealPix library ' 
                      'to generate the sampling of the first two Euler angles on the sphere. '
                      'The samplings are approximate numbers and vary slightly over the sphere.')
        form.addParam('offsetSearchRange', IntParam, default=5,
                      label='Offset search range (px)',
                      help='Probabilities will be calculated only for translations in a circle with this radius (in pixels). '
                      'The center of this circle changes at every iteration and is placed at the optimal translation for each image in the previous iteration.')
        form.addParam('offsetSearchStep', IntParam, default=1,
                      label='Offset search step (px)',
                      help='Translations will be sampled with this step-size (in pixels).Translational sampling is also done using the adaptive approach.' 
                      'Therefore, if adaptive=1, the translations will first be evaluated on a 2x coarser grid.')
        form.addParam('doLocalAngSearch', BooleanParam, default=False,
                      label='Perform local angular search?',
                      help='If set to Yes, then rather than performing exhaustive angular searches, local searches within the range given below will be performed. '
                      'A prior Gaussian distribution centered at the optimal orientation in the previous iteration '
                      'and with a stddev of 1/3 of the range given below will be enforced.')
        form.addParam('additionalArgs', TextParam, expertLevel=LEVEL_ADVANCED,
                      label='Additional arguments', default="",
                      help='In this box command-line arguments may be provided that are not generated by the GUI. '
                      'This may be useful for testing developmental options and/or expert use of the program. '
                      'The command \'relion_refine\' will print a list of possible options.')  
        
    def _defineSteps(self):
        
        ProtRelionBase.defineSteps(self)
        # launch relion program
        self.insertRelionClassify()
        
        ProtRelionBase.defineSteps2(self, 1,self.numberOfIterations.get(), self.numberOfClasses.get())       
        
        self._insertFunctionStep('createOutput') 
       
    def insertRelionClassify(self):
        args = {'--iter': self.numberOfIterations.get(),
                '--tau2_fudge': self.regularParT.get(),
                '--flatten_solvent': '',
                #'--zero_mask': '',# this is an option but is almost always true
                '--norm': '',
                '--scale': '',
                '--o': self._getExtraPath('relion')
                }
        if len(self.maskFile.get()):
            args['--solvent_mask'] = self.maskFile.get()

        args.update({'--i': self.imgStar,
                     '--particle_diameter': self.maskRadius.get() * 2.0 ,
                     '--angpix': self.inputParticles.get().getSamplingRate(),
                     '--ref': self.volumeName,
                     '--oversampling': '1'
                     })
        
        if not self.isRefMapGreyScale:
            args[' --firstiter_cc'] = '' 
            
        if self.iniLowPassFilter.get() > 0:
            args['--ini_high'] = self.iniLowPassFilter.get()
            
        # CTF stuff
        if self.doCtf:
            args['--ctf'] = ''
        
        if self.isReferenceCtfCorrected:
            args['--ctf_corrected_ref'] = ''
            
        if self.areImagesPhaseFlipped:
            args['--ctf_phase_flipped'] = ''
            
        if self.ignoreFirstPeakCTF:
            args['--ctf_intact_first_peak'] = ''
            
        args['--sym'] = self.symmetryGroup.get().upper()
        
        args['--K'] = self.numberOfClasses.get()
            
        # Sampling stuff
        # Find the index(starting at 0) of the selected
        # sampling rate, as used in relion program
        iover = 1 #TODO: check this DROP THIS
#        index = ['30','15','7.5','3.7','1.8',
#                 '0.9','0.5','0.2','0.1'].index(self.AngularSamplingDeg)
        index = self.angSampling.get()
        #Chequear con Roberto porque si iover es 1 es tonteria lo que esta haciendo
        args['--healpix_order'] = float(index + 1 - iover)
        
        if self.doLocalAngSearch:
            args['--sigma_ang'] = self.offsetSearchRange.get() / 3.
            
        args['--offset_range'] = self.offsetSearchRange.get()
        args['--offset_step']  = self.offsetSearchStep.get() * pow(2, iover)

        args['--j'] = self.numberOfThreads.get()
        
        # Join in a single line all key, value pairs of the args dict    
        params = ' '.join(['%s %s' % (k, str(v)) for k, v in args.iteritems()])
        if self.additionalArgs.get() is not None:
            params += ' ' + self.additionalArgs.get()

        self._insertRunJobStep(self.program, params)

    #TODO
    def createOutput(self):
        pass

    #TODO
    def _summary(self):
        summary = []
        summary.append("Input images:  %s" % self.inputParticles.get().getNameId())
        return summary

    #TODO
    def _validate(self):
        errors = []
        return errors
            
