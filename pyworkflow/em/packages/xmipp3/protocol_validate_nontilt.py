# **************************************************************************
# *
# * Authors:     Javier Vargas (jvargas@cnb.csic.es)
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

from os.path import join
from pyworkflow.utils.process import runJob
from pyworkflow.object import Float, String
from pyworkflow.protocol.params import (PointerParam, FloatParam, STEPS_PARALLEL,
                                        StringParam, BooleanParam, EnumParam, LEVEL_ADVANCED)
from pyworkflow.em.data import Volume
from pyworkflow.em import Viewer
from pyworkflow.em.protocol import ProtAnalysis3D
from pyworkflow.utils.path import moveFile, makePath
from pyworkflow.em.packages.xmipp3.convert import (writeSetOfParticles,
                                                   writeSetOfVolumes,
                                                   getImageLocation)
import xmipp


class XmippProtValidateNonTilt(ProtAnalysis3D):
    """    
    Ranks a set of volumes according to their alignment reliability obtained from a clusterability test.
    """
    _label = 'validate_nontilt'
    PROJECTION_MATCHING = 0
    SIGNIFICANT = 1
    WEB = 0
    
    def __init__(self, *args, **kwargs):
        ProtAnalysis3D.__init__(self, *args, **kwargs)
        
        if (self.WEB == 1):
            self.stepsExecutionMode = STEPS_PARALLEL
        
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputVolumes', PointerParam, pointerClass='SetOfVolumes, Volume',
                      label="Input volumes",  
                      help='Select the input volumes.')     

        form.addParam('inputParticles', PointerParam, pointerClass='SetOfParticles, SetOfClasses2D', 
                      label="Input particles",  
                      help='Select the input projection images .') 
            
        form.addParam('symmetryGroup', StringParam, default='c1',
                      label="Symmetry group", 
                      help='See [[Xmipp Symmetry][http://www2.mrc-lmb.cam.ac.uk/Xmipp/index.php/Conventions_%26_File_formats#Symmetry]] page '
                           'for a description of the symmetry format accepted by Xmipp') 
        
        form.addParam('alignmentMethod', EnumParam, label='Image alignment', choices=['Projection_Matching','Significant'], default=self.PROJECTION_MATCHING)

        form.addParam('highPassFilter', FloatParam, label='Volume high-pass the (A)', default=15)
        form.addParam('lowPassFilter' , FloatParam, label='Volume low-pass the (A)', default=50)        
        
        form.addParam('angularSampling', FloatParam, default=5, expertLevel=LEVEL_ADVANCED,
                      label="Angular Sampling (degrees)",  
                      help='Angular distance (in degrees) between neighboring projection points ')

        form.addParam('numOrientations', FloatParam, default=6, expertLevel=LEVEL_ADVANCED,
                      label="Number of orientations per particle",  
                      help='Number of possible orientations in which a particle can be \n')
        
        #form.addParallelSection(threads=1, mpi=1)
        form.addParallelSection(threads=0, mpi=4)


    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):        
        convertId = self._insertFunctionStep('convertInputStep', 
                                             self.inputParticles.get().getObjId())
        deps = [] # store volumes steps id to use as dependencies for last step
        
        if ( self.alignmentMethod == self.SIGNIFICANT):
            commonParams = self._getCommonParamsSignificant()
        else:
            commonParams = self._getCommonParamsProjection()
            
        sym = self.symmetryGroup.get()       
        
        for i, vol in enumerate(self._iterInputVols()):
            
            volName = getImageLocation(vol)            
            volDir = self._getVolDir(i+1)  
            volNameFilt = (volDir + '_filt.vol')
          
            filterId=self._insertFunctionStep('filterVolumeStep',volName,volNameFilt,volDir,prerequisites=[convertId])
                                    
            if (self.alignmentMethod == self.SIGNIFICANT):            
                sigStepId = self._insertFunctionStep('significantStep', 
                                                     volNameFilt, volDir,
                                                     commonParams, 
                                                     prerequisites=[filterId])
                
            else:            

                pmStepId = self._insertFunctionStep('projectionLibraryStep', 
                                                     volNameFilt, volDir,
                                                     prerequisites=[filterId])
                
                
                sigStepId = self._insertFunctionStep('projectionMatchingStep', 
                                                     volNameFilt, volDir,
                                                     commonParams, 
                                                     prerequisites=[pmStepId])
                            

            
            volStepId = self._insertFunctionStep('validationStep', 
                                                 volNameFilt, volDir,
                                                 sym, 
                                                 prerequisites=[sigStepId])

            
            deps.append(volStepId)
        self._insertFunctionStep('createOutputStep', 
                                 prerequisites=deps)
        
    def convertInputStep(self, particlesId):
        """ Write the input images as a Xmipp metadata file. 
        particlesId: is only need to detect changes in
        input particles and cause restart from here.
        """
        writeSetOfParticles(self.inputParticles.get(), 
                            self._getPath('input_particles.xmd'))
        
        
        
    def filterVolumeStep(self, volName, volNameFilt, volDir):
        
        params =  '  -i %s' % volName
        params +=  ' -o %s' % volNameFilt
        params += '  --fourier '
        params += '  band_pass %f %f'  % (self.lowPassFilter.get(), self.highPassFilter.get())
        params += '  --sampling %f ' % self.inputParticles.get().getSamplingRate()
        
        self.runJob('xmipp_transform_filter', 
                    params, numberOfMpi=1,numberOfThreads=1)
        
    def _getCommonParamsSignificant(self):
        
        params =  '  -i %s' % self._getPath('input_particles.xmd')        
        params += ' --sym %s' % self.symmetryGroup.get()
        params += ' --angularSampling %0.3f' % self.angularSampling.get()
        params += ' --dontReconstruct'
        params += ' --useForValidation %0.3f' % (self.numOrientations.get())        
        return params


    def _getCommonParamsProjection(self):
        params =  '  -i %s' % self._getPath('input_particles.xmd')        
        params += ' --Ri 0.0'
        params += ' --Ro %0.3f' % ((self.inputParticles.get().getDimensions()[0])/2)
        params += ' --max_shift %0.3f' % ((self.inputParticles.get().getDimensions()[0])/10)
        params += ' --append' 
        params += ' --search5d_shift 5'
        params += ' --number_orientations %0.3f' % self.numOrientations.get()
                     
        return params
    
    def projectionLibraryStep(self, volName, volDir):
        
        # Generate projections from this reconstruction        
        nproc = self.numberOfMpi.get()
        nT=self.numberOfThreads.get() 
        
        makePath(volDir)
        fnGallery= (volDir+'/gallery.stk')
        params = '-i %s -o %s --sampling_rate %f --sym %s --method fourier 1 0.25 bspline --compute_neighbors --angular_distance %f --experimental_images %s --max_tilt_angle 90'\
                    %(volName,fnGallery,self.angularSampling.get(),self.symmetryGroup.get(), -1, self._getPath('input_particles.xmd'))
        
        print params
        self.runJob("xmipp_angular_project_library", params, numberOfMpi=nproc, numberOfThreads=nT)                    
        #self.runJob("ldd","`which xmipp_mpi_angular_project_library`",numberOfMpi=1)                    

    def projectionMatchingStep(self, volName, volDir, params):

        nproc = self.numberOfMpi.get()
        nT=self.numberOfThreads.get() 
        params += '  -o %s' % (volDir+'/angles_iter001_00.xmd')
        params += ' --ref %s' % (volDir+'/gallery.stk')
        self.runJob('xmipp_angular_projection_matching', 
                    params, numberOfMpi=nproc,numberOfThreads=nT)
        
    
    def significantStep(self, volName, volDir, params):

        nproc = self.numberOfMpi.get()
        nT=self.numberOfThreads.get() 
        makePath(volDir)  
        params += '  --initvolumes %s' % volName  
        params += ' --odir %s' % volDir
        params += ' --iter %d' % 1
        self.runJob('xmipp_reconstruct_significant', 
                    params, numberOfMpi=nproc,numberOfThreads=nT)
        
    def validationStep(self, volName, volDir,sym):
        makePath(volDir)                  
        nproc = self.numberOfMpi.get()
        nT=self.numberOfThreads.get() 
        
        params  = '  --i %s' % (volDir+'/angles_iter001_00.xmd')
        params += '  --volume %s' % volName  
        params += ' --odir %s' % volDir
        params += ' --sym %s' % sym
        
        if (self.alignmentMethod == self.SIGNIFICANT):
            params += ' --useSignificant '
        
        self.runJob('xmipp_validation_nontilt', params,numberOfMpi=nproc,numberOfThreads=nT)
        
    def createOutputStep(self):
        outputVols = self._createSetOfVolumes()
        
        for i, vol in enumerate(self._iterInputVols()):
            volume = vol.clone()
            volDir = self._getVolDir(i+1)
            volPrefix = 'vol%03d_' % (i+1)
            validationMd = self._getExtraPath(volPrefix + 'validation.xmd')
            moveFile(join(volDir, 'validation.xmd'), 
                     validationMd)
            clusterMd = self._getExtraPath(volPrefix + 'clusteringTendency.xmd')
            moveFile(join(volDir, 'clusteringTendency.xmd'), clusterMd)
            
            md = xmipp.MetaData(validationMd)
            weight = md.getValue(xmipp.MDL_WEIGHT, md.firstObject())
            volume.weight = Float(weight)
            volume.clusterMd = String(clusterMd)
            volume.cleanObjId() # clean objects id to assign new ones inside the set
            outputVols.append(volume)                
        
        outputVols.setSamplingRate(volume.getSamplingRate())
        self._defineOutputs(outputVolumes=outputVols)
        self._defineTransformRelation(self.inputVolumes, outputVols)
        
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        validateMsgs = []
        # if there are Volume references, it cannot be empty.
        if self.inputVolumes.get() and not self.inputVolumes.hasValue():
            validateMsgs.append('Please provide an input reference volume.')
        if self.inputParticles.get() and not self.inputParticles.hasValue():
            validateMsgs.append('Please provide input particles.')            
        return validateMsgs
        
    
    def _summary(self):
        summary = []

        if  (not hasattr(self,'outputVolumes')):
            summary.append("Output volumes not ready yet.")
        else:
            size = 0
            for i, vol in enumerate(self._iterInputVols()):
                size +=1
            summary.append("Volumes to validate: *%d* " % size)
            summary.append("Angular sampling: %s" % self.angularSampling.get())
            summary.append("Significance value: %s" % self.alpha.get())

        return summary
    
    def _methods(self):
        messages = []
        if (hasattr(self,'outputVolumes')):
            messages.append('The quality parameter(s) has been obtained using the approach [Vargas2014a] with angular sampling of %f and significant value of %f' % (self.angularSampling.get(), self.alpha.get()))
        return messages
    
    def _citations(self):
        return ['Vargas2014a']
    
    #--------------------------- UTILS functions --------------------------------------------
    def _getVolDir(self, volIndex):
        return self._getExtraPath('vol%03d' % volIndex)
    
    def _iterInputVols(self):
        """ In this function we will encapsulate the logic
        to iterate throught the input volumes.
        This give the flexibility of having Volumes, SetOfVolumes or 
        a combination of them as input and the protocol code
        remain the same.
        """
        inputVols = self.inputVolumes.get()
        
        if isinstance(inputVols, Volume):
            yield inputVols
        else:
            for vol in inputVols:
                yield vol
        
        
    def _defineMetadataRootName(self, mdrootname,volId):
        
        if mdrootname=='P':
            VolPrefix = 'vol%03d_' % (volId)
            return self._getExtraPath(VolPrefix+'clusteringTendency.xmd')
        if mdrootname=='Volume':

            VolPrefix = 'vol%03d_' % (volId)
            return self._getExtraPath(VolPrefix+'validation.xmd')
            
    def _definePName(self):
        fscFn = self._defineMetadataRootName('P')
        return fscFn
    
    def _defineVolumeName(self,volId):
        fscFn = self._defineMetadataRootName('Volume',volId)
        return fscFn

