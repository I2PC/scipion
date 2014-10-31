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

from pyworkflow.object import Float, String
from pyworkflow.protocol.params import (PointerParam, FloatParam, STEPS_PARALLEL,
                                        StringParam, BooleanParam, LEVEL_ADVANCED)
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
    Reconstruct a volume using Xmipp_reconstruct_fourier from a given set of particles.
    The alignment parameters will be converted to a Xmipp xmd file
    and used as direction projections to reconstruct.
    """
    _label = 'validate_nontilt'
    
    def __init__(self, *args, **kwargs):
        ProtAnalysis3D.__init__(self, *args, **kwargs)
        
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputVolumes', PointerParam, pointerClass='SetOfVolumes,Volume',
                      label="Input volumes",  
                      help='Select the input volumes.')     

        form.addParam('inputParticles', PointerParam, pointerClass='SetOfParticles, SetOfClasses2D', 
                      label="Input particles",  
                      help='Select the input projection images .') 
            
        form.addParam('symmetryGroup', StringParam, default='c1',
                      label="Symmetry group", 
                      help='See [[Xmipp Symmetry][http://www2.mrc-lmb.cam.ac.uk/Xmipp/index.php/Conventions_%26_File_formats#Symmetry]] page '
                           'for a description of the symmetry format accepted by Xmipp') 
        
        form.addParam('angularSampling', FloatParam, default=10,
                      label="Angular Sampling (degrees)",  
                      help='Angular distance (in degrees) between neighboring projection points ')

        form.addParam('alpha', FloatParam, default=0.05,
                      label="Significance value",  
                      help='Parameter to define the corresponding most similar volume \n' 
                      '    projected images for each projection image')
        
        form.addParallelSection(threads=1, mpi=1)

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):        
        convertId = self._insertFunctionStep('convertInputStep', 
                                             self.inputParticles.get().getObjId())
        deps = [] # store volumes steps id to use as dependencies for last step
        commonParams = self._getCommonParams()
        sym = self.symmetryGroup.get()
        for i, vol in enumerate(self._iterInputVols()):
            volName = getImageLocation(vol)
            volDir = self._getVolDir(i+1)
            
            sigStepId = self._insertFunctionStep('significantStep', 
                                                 volName, volDir,
                                                 commonParams, 
                                                 prerequisites=[convertId])            

            volStepId = self._insertFunctionStep('validationStep', 
                                                 volName, volDir,
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
                            self._getTmpPath('input_particles.xmd'))
                    
    def _getCommonParams(self):
        params =  '  -i %s' % self._getTmpPath('input_particles.xmd')        
        params += ' --sym %s' % self.symmetryGroup.get()
        params += ' --alpha0 %0.3f' % self.alpha.get()
        params += ' --angularSampling %0.3f' % self.angularSampling.get()
        return params
    
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
        params = '  --volume %s' % volName  
        params += ' --odir %s' % volDir
        params += ' --sym %s' % sym
        self.runJob('xmipp_validation_nontilt', params,numberOfMpi=1,numberOfThreads=1)
        
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
        #self._defineTransformRelation(self.inputVolumes.get(), volume)
        
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

        summary.append("Input particles:  %s" % self.inputParticles.get().getNameId())

        summary.append("-----------------")
        if self.inputVolumes.get():
            for i, vol in enumerate(self._iterInputVols()):
                summary.append("Input volume(s)_%d: [%s]" % (i+1,vol))

        summary.append("-----------------")
        if  (not hasattr(self,'outputVolumes')):
            summary.append("Output volumes not ready yet.")
        else:
            for i, vol in enumerate(self._iterInputVols()):
                            
                VolPrefix = 'vol%03d_' % (i+1)
                md = xmipp.MetaData(self._getExtraPath(VolPrefix+'validation.xmd'))                
                weight = md.getValue(xmipp.MDL_WEIGHT, md.firstObject())
                summary.append("Output volume(s)_%d : %s" % (i+1,self.outputVolumes.getNameId()))
                summary.append("Quality parameter_%d : %f" % (i+1,weight))
                summary.append("-----------------")        
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
        return self._getTmpPath('vol%03d' % volIndex)
    
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

