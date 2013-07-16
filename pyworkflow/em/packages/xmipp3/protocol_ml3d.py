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
This sub-package contains wrapper around ML3D Xmipp program
"""


from pyworkflow.em import *  
from pyworkflow.utils import *  
import xmipp
from data import *
from xmipp3 import XmippProtocol


class XmippDefML3D(Form):
    """Create the definition of parameters for
    the XmippProtML2D protocol"""
    def __init__(self):
        Form.__init__(self)
    
        self.addSection(label='Input')
        self.addParam('inputImages', PointerParam, label="Input images", important=True, 
                      pointerClass='SetOfParticles',
                      help='Select the input images from the project.'
                           'It should be a SetOfImages class')  
        #TODO: the following parameter should be a Pointer to a Volume and not a string containing the path      
        self.addParam('ini3DrefVolumes', TextParam,
                      label='Initial 3D reference volumes', 
                      help='Initial 3D density maps with the same dimensions as your particles.')
        self.addParam('numberOfSeedsPerRef', IntParam, default=1,
                      label='Number of seeds per reference',
                      help='The total number of seeds generated will be the number of provided '
                      'reference volumes times the number of seeds per reference. '
                      'If you provide 2 initial volumes and 3 seeds per referece you will '
                      'produce 6 3D maps.')
        self.addParam('doCorrectGreyScale', BooleanParam, default=False,
                      label="Correct the absolute grey-scale of initial references?", 
                      help='The probabilities are based on squared differences, so that the '
                      'absolute grey scale is important.')
        self.addParam('projMatchSampling', FloatParam, default=15.0, condition='doCorrectGreyScale',
                      label='Sampling for projection matching',
                      help='Angular sampling for a quick projection matching '
                      'to obtain right grey scale. As the resolution of the intial reference '
                      'should be low, this sampling can be relatively crude, e.g. 15')          
        self.addParam('doLowPassFilter', BooleanParam, default=True,
                      label="Low-pass filter initial references?", 
                      help='It is highly recommended to low-pass filter your initial reference '
                      'volume as much as you can.')   
        self.addParam('lowPassFilter', IntParam, default=50, condition='doLowPassFilter',
                      label='Resolution of the low-pass filter (Ang)',
                      help='Resolution of the low-pass filter in Angstroms.')       
        self.addSection(label='ML3D classification')
        self.addParam('angularSampling', IntParam, default=10,
                      label='Angular sampling for classification',
                      help='Fine samplings take huge amounts of CPU and memory. '
                      'Therefore, in general, dont use samplings finer than 10 degrees.')    
        self.addParam('numberOfIterations', IntParam, default=25,
                      label='Number of ML(F)3D iterations to perform',
                      help='Number of ML(F)3D iterations to perform.')
        self.addParam('symmetry', TextParam, default='c1',
                      label='Point group symmetry',
                      help='Number of ML(F)3D iterations to perform.')        
        self.addParam('doNorm', BooleanParam, default=False,
                      label="Refine the normalization for each image?", 
                      help='This variant of the algorithm deals with normalization errors.')    
        self.addParam('restartIter', IntParam, default=0,
                      expertLevel=LEVEL_ADVANCED,
                      label='Restart after iteration',
                      help='For previous runs that stopped before convergence, '
                      'resume the calculations after the completely finished iteration, '
                      'i.e. including all 3D reconstructions. '
                      'Note that all flags about grey-scale correction, filtering and '
                      'seed generation will be ignored if a value larger than 0 is given, '
                      'since this option only concerns the ML3D classification part.')   
        self.addParam('extraParams', TextParam,
                      expertLevel=LEVEL_ADVANCED,
                      label='Additional parameters',
                      help='Additional xmipp_ml(f)_refine3d parameters.')                  
        self.addSection(label='MLF parameters', questionParam='doMlf')        
        self.addParam('doMlf', BooleanParam, default=False,
                      label='Use MLF2D instead of ML2D')
        self.addParam('doCorrectAmplitudes', BooleanParam, default=True,
                      label='Use CTF-amplitude correction inside MLF?',
                      help='If set to <Yes>, the input images file should contain '
                           'the CTF information for each image.')
        self.addParam('highResLimit', IntParam, default=20,
                      label='High-resolution limit (in Angstroms)',
                      help='No frequencies higher than this limit will be taken into account. '
                      'If zero is given, no limit is imposed.')     
        self.addParam('areImagesPhaseFlipped', BooleanParam, default=True,
                      label='Are the images CTF phase flipped?',
                      help='You can run MLF with or without having phase flipped the images.')     
        self.addParam('initialMapIsAmplitudeCorrected', BooleanParam, default=False,
                      label='Are initial references CTF-amplitude corrected?',
                      help='If coming from programs other than xmipp_mlf_refine3d this is '
                      'usually not the case. If you will perform a grey-scale correction, '
                      'this parameter becomes irrelevant as the output maps never have the '
                      'CTF-amplitudes corrected.')     
        self.addParam('seedsAreAmplitudeCorrected', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Are the seeds CTF-amplitude corrected?',
                      help='This option is only relevant if you provide your own seeds! '
                      'If the seeds are generated automatically, this parameter becomes '
                      'irrelevant as they will always be amplitude-corrected.') 
        self.addSection(label='3D Reconstruction', expertLevel=LEVEL_ADVANCED)    
        self.addParam('reconstructionMethod', EnumParam, choices=['fourier', 'wlsART'], 
                      default=0, label='Reconstruction method', display=EnumParam.DISPLAY_LIST,
                      expertLevel=LEVEL_ADVANCED, 
                      help='Choose between wslART or fourier.')
        self.addParam('aRTExtraParams', TextParam,
                      condition='reconstructionMethod==1', expertLevel=LEVEL_ADVANCED,
                      label='Extra parameters',
                      help='Additional reconstruction parameters for ART.')  
        self.addParam('fourierExtraParams', TextParam, 
                      condition='reconstructionMethod==0', expertLevel=LEVEL_ADVANCED,
                      label='Extra parameters',
                      help='The Fourier-interpolation reconstruction method is much faster than wlsART '
                      'and may give similar results. It however is not guaranteed to optimize the '
                      'likelihood function. This is an experimental feature. One may limit the '
                      'maximum resolution of the fourier-interpolation using -max_resolution 0.3 '
                      '(to 0.3 digital frequency). Use the extra parameter entry below for that.')  
        
        self.addParallelSection(threads=1, mpi=8)
        
        
class XmippProtML3D(XmippProtocol, ProtRefine3D, ProtClassify3D):
    """ Protocol for Xmipp-based ML3D/MLF3D classification. """
    _definition = XmippDefML3D()
    _label = 'Xmipp ML3D'
    
    # Reconstruction method constants
    FOURIER = 0
    WLSART = 1

    def getProgramId(self):
        progId = "ml"
        if self.doMlf:
            progId += "f" 
        return progId   
         
    def _defineSteps(self):

        self.ParamsDict = {}
        self.ParamsDict['ProgId'] = self.getProgramId()
        
        #TODO: Retrieve initial volume path from object (now is a text) and convert if needed
        refMd = self.ini3DrefVolumes.get()
        initVols = self.ParamsDict['InitialVols'] = self._getExtraPath('initial_volumes.stk')
        self.mdVols = xmipp.MetaData(refMd)
        
        # Convert input images if necessary
        self.inputImgs = self.inputImages.get()        
        imgsFn = self._insertConvertStep('inputImgs', XmippSetOfImages,
                                         self._getPath('input_images.xmd'))
        self.imgMd = self.ParamsDict['ImgMd'] = imgsFn
        
        #Get sampling rate from input images
        self.samplingRate = self.inputImages.get().samplingRate.get()
        
        self._insertFunctionStep('copyVolumes', refMd, initVols)
        
        if self.doCorrectGreyScale.get():
            self.insertCorrectGreyScaleSteps()
                    
        if self.doLowPassFilter.get():
            self.insertFilterStep()
            
        if self.numberOfSeedsPerRef.get() > 1:
            self.insertGenerateRefSteps()
            
        self.insertML3DStep(self.imgMd, self.workingDir.get() + '/', self.ParamsDict['InitialVols'], 
                            self.numberOfIterations.get(), self.seedsAreAmplitudeCorrected.get())
        
        self._insertFunctionStep('renameOutput', self.workingDir.get(), self.getProgramId())
                
        self._insertFunctionStep('createOutput')

    def copyVolumes(self, inputMd, outputStack):
        ''' This function will copy input references into a stack in working directory'''
        if exists(outputStack):
            os.remove(outputStack)
        md = xmipp.MetaData(inputMd)
        img = xmipp.Image()
        for i, idx in enumerate(md):
            img.read(md.getValue(xmipp.MDL_IMAGE, idx))
            img.write('%d@%s' % (i + 1, outputStack))

        
    def insertCorrectGreyScaleSteps(self):
        ''' Correct the initial reference greyscale '''
        cgsDir = self._getPath('CorrectGreyscale')
        makePath(cgsDir)
        volStack = self.ParamsDict['InitialVols'] = self._getExtraPath('corrected_volumes.stk')
        # Grey-scale correction always leads to an amplitude uncorrected map
        self.initialMapIsAmplitudeCorrected.set(False)
        index = 1
        outputVol = ''
        for idx in self.mdVols:
            volDir = join(cgsDir, 'vol%03d' % index)
            projs = join(volDir, 'projections')
            makePath(volDir)
            outputVol = "%(index)d@%(volStack)s" % locals()
            corrRefsRoot = join(volDir, 'corrected_refs')
            self.ParamsDict.update({
                'inputVol': self.mdVols.getValue(xmipp.MDL_IMAGE, idx),
                'outputVol': outputVol,
                'projRefs': projs + ".stk",
                'docRefs': projs + ".doc",
                'corrRefsRoot':corrRefsRoot ,
                'corrRefs': corrRefsRoot + '_Ref3D_001.stk',
                'projMatch': join(volDir, "proj_match.doc"),
                'projMatchSampling': self.projMatchSampling.get(),
                'symmetry': self.symmetry.get(),
                'numberOfThreads': self.numberOfThreads.get()
                })
            self.mdVols.setValue(xmipp.MDL_IMAGE, outputVol, idx)
            self.ParamsStr = ' -i %(inputVol)s --experimental_images %(ImgMd)s -o %(projRefs)s' + \
                    ' --sampling_rate %(projMatchSampling)f --sym %(symmetry)s' + \
                    'h --compute_neighbors --angular_distance -1' 
                       
            self._insertRunJobStep('xmipp_angular_project_library', self.ParamsStr % self.ParamsDict)

            self.ParamsStr = '-i %(ImgMd)s -o %(projMatch)s --ref %(projRefs)s' 
            self._insertRunJobStep('xmipp_angular_projection_matching', self.ParamsStr % self.ParamsDict)
 
            self.ParamsStr = '-i %(projMatch)s --lib %(docRefs)s -o %(corrRefsRoot)s'
            self._insertRunJobStep('xmipp_angular_class_average', self.ParamsStr % self.ParamsDict)

            self.ParamsStr = '-i %(projMatch)s -o %(outputVol)s --sym %(symmetry)s --weight --thr %(numberOfThreads)d'
            self._insertRunJobStep('xmipp_reconstruct_fourier', self.ParamsStr % self.ParamsDict)
            index += 1

    def insertFilterStep(self):
        volStack = self.ParamsDict['FilteredVols'] = self._getExtraPath('filtered_volumes.stk')
        index = 1
        outputVol = ''
        for idx in self.mdVols:
            outputVol = "%(index)d@%(volStack)s" % locals()
            self.mdVols.setValue(xmipp.MDL_IMAGE, outputVol, idx)
            index += 1
        self.ParamsStr = '-i %(InitialVols)s -o %(FilteredVols)s --fourier low_pass %(lowPassFilter)f --sampling %(samplingRate)f'
        self.ParamsDict.update({
                                'lowPassFilter':self.lowPassFilter.get(),
                                'samplingRate':self.samplingRate
                                })
        self._insertRunJobStep('xmipp_transform_filter', self.ParamsStr % self.ParamsDict)
        self.ParamsDict['InitialVols'] = self.ParamsDict['FilteredVols']

    def insertGenerateRefSteps(self):
        ''' Generate more reference volumes than provided in input reference '''
        grDir = self._getPath('GeneratedReferences')
        # Create dir for seeds generation
        makePath(grDir)
        # Split images metadata
        nvols = self.ParamsDict['NumberOfVols'] = self.mdVols.size() * self.numberOfSeedsPerRef.get()
        sroot = self.ParamsDict['SplitRoot'] = join(grDir, 'images')
        self.ParamsStr = '-i %(ImgMd)s -n %(NumberOfVols)d --oroot %(SplitRoot)s'
        files = ['%s%06d.xmd' % (sroot, i) for i in range(1, nvols+1)]        
        self._insertRunJobStep('xmipp_metadata_split', self.ParamsStr % self.ParamsDict, numberOfMpi=1)
        
        volStack = self.ParamsDict['InitialVols'] = self._getExtraPath('generated_volumes.stk') 
        index = 1
        copyVols = []
        for idx in self.mdVols:
            for i in range(self.numberOfSeedsPerRef.get()):
                outputVol = "%d@%s" % (index, volStack)
                generatedVol = join(grDir, "vol%03dextra/iter%03d/vol%06d.vol" % (index, 1, 1))
                copyVols.append((outputVol, generatedVol))
                self.insertML3DStep(files[index-1], join(grDir, 'vol%03d' % index), self.mdVols.getValue(xmipp.MDL_IMAGE, idx), 1, 
                                    self.initialMapIsAmplitudeCorrected)
                #self.mdVols.setValue(MDL_IMAGE, outputVol, idx)
                index += 1
                
        for outVol, genVol in copyVols:
            self.ParamsDict.update({'outVol': outVol, 'genVol':genVol})
            self.ParamsStr = '-i %(genVol)s -o %(outVol)s'
            self._insertRunJobStep('xmipp_image_convert', self.ParamsStr % self.ParamsDict, numberOfMpi=1)
            
        # Seed generation with MLF always does amplitude correction
        self.seedsAreAmplitudeCorrected.set(True)

    def insertML3DStep(self, inputImg, oRoot, initialVols, numberOfIters, amplitudCorrected):
        self.ParamsDict.update({
                         '_ImgMd': inputImg,
                         '_ORoot': oRoot,
                         '_InitialVols': initialVols,
                         '_NumberOfIterations': numberOfIters,
                         'symmetry':self.symmetry.get(),
                         'angularSampling': self.angularSampling.get(),
                         'extraParams': self.extraParams.get(),
                         'numberOfThreads': self.numberOfThreads.get(),
                         'highResLimit': self.highResLimit.get(),
                         'samplingRate': self.samplingRate,
                         'reconstructionMethod': self.getEnumText('reconstructionMethod'),
                         'aRTExtraParams': self.aRTExtraParams.get(),
                         'fourierExtraParams': self.fourierExtraParams.get()
                        })
        self.ParamsStr = "-i %(_ImgMd)s --oroot %(_ORoot)s --ref %(_InitialVols)s --iter %(_NumberOfIterations)d " + \
                         "--sym %(symmetry)s --ang %(angularSampling)s"
        if self.aRTExtraParams.hasValue(): self.ParamsStr += " %(extraParams)s"
#        if self.NumberOfReferences > 1:
#            self.ParamsStr += " --nref %(NumberOfReferences)s"
        if self.numberOfThreads.get() > 1:
            self.ParamsStr += " --thr %(numberOfThreads)d"
        if self.doNorm.get():
            self.ParamsStr += " --norm"
        
        if self.doMlf.get():
            if not self.doCorrectAmplitudes.get():
                self.ParamsStr += " --no_ctf"
            if not self.areImagesPhaseFlipped.get():
                self.ParamsStr += " --not_phase_flipped"
            if not amplitudCorrected:
                self.ParamsStr += " --ctf_affected_refs"
            if self.highResLimit > 0:
                self.ParamsStr += " --limit_resolution 0 %(highResLimit)f"
            self.ParamsStr += ' --sampling_rate %(samplingRate)f'

        self.ParamsStr += " --recons %(reconstructionMethod)s "
        
        if self.reconstructionMethod.get() == self.WLSART:
            if self.aRTExtraParams.hasValue(): self.ParamsStr += " %(aRTExtraParams)s"
        else:
            if self.fourierExtraParams.hasValue(): self.ParamsStr += " %(fourierExtraParams)s" 
            
        self._insertRunJobStep('xmipp_%s_refine3d' % self.getProgramId(), self.ParamsStr % self.ParamsDict)
        
    def renameOutput(self, WorkingDir, ProgId):
        ''' Remove ml2d prefix from:
            ml2dclasses.stk, ml2dclasses.xmd and ml2dimages.xmd'''
        prefix = '%s2d' % ProgId
        for f in ['%sclasses.stk', '%sclasses.xmd', '%simages.xmd']:
            f = join(WorkingDir, f % prefix)
            nf = f.replace(prefix, '')
            shutil.move(f, nf)
                                                    
    def createOutput(self):
        lastIter = 'iter%03d' % self.numberOfIterations.get()
        volumes = XmippSetOfVolumes(self._getExtraPath(lastIter + '/' + 'iter_volumes.xmd'))
        self._defineOutputs(outputVolumes=volumes)

    def _summary(self):
        summary = []
        summary.append("Input images:  %s" % self.inputImages.get().getNameId())
        if self.doMlf.get():
            if self.doCorrectAmplitudes.get():
                suffix = "with CTF correction "
            else:
                suffix = "ignoring CTF effects "
            summary.append("Using a ML in <Fourier-space> " + suffix)
         
        summary.append("Reference volumes(s): [%s]" % self.ini3DrefVolumes.get())

        if self.numberOfSeedsPerRef.get() > 1:
            summary.append("Number of references per volume: <%d>" % self.numberOfSeedsPerRef.get())
           
        # TODO: Add information at info@iter_classes.xmd from last iteration
        
        if not hasattr(self, 'outputVolumes'):
            summary.append("Output volumes not ready yet.")
        else:
            summary.append("Output volumes: %s" % self.outputVolumes.get())
        return summary
    
    def _validate(self):
        validateMsgs = []
        # Volume references cannot be empty (REMOVE this when we use a pointer to a set instead of a path)
        if not self.ini3DrefVolumes.hasValue():
            validateMsgs.append('Please provide an initial reference volume(s).')
        # Correct grey scale needs mpi to b e run
        if self.doCorrectGreyScale.get() and self.numberOfMpi < 2:
            validateMsgs.append('Correct grey scale needs mpi to be run.')
            
        #TODO: Check images dimension when it is implemented on SetOFImages class
        return validateMsgs