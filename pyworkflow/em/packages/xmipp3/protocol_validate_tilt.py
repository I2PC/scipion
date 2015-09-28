# **************************************************************************
# *
# * Authors:     Jose Luis Vilas (jlvilas@cnb.csic.es)
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

import xmipp
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



class XmippProtValidateTilt(ProtAnalysis3D):
    """    
    Validates initial volumes by using a set of projections/classes from a tilted-pair picking process and using RCT algorithm
    """
    _label = 'validate tilt'
    
    def __init__(self, *args, **kwargs):
        ProtAnalysis3D.__init__(self, *args, **kwargs)
        
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputVolume', PointerParam, pointerClass='Volume',
                      label="Input volume",  
                      help='Select the input volume.')     

        form.addParam('inputTiltPair', PointerParam, label="Input angles tilt pair", 
                      pointerClass='ParticlesTiltPair',
                      help='Select the input particles tilt pair file that will be used.  file. This file is used to associate each micrograph with its tilted equivalent.')
           
        form.addParam('symmetryGroup', StringParam, default='c1',
                      label="Symmetry group", 
                      help='See [[Xmipp Symmetry][http://www2.mrc-lmb.cam.ac.uk/Xmipp/index.php/Conventions_%26_File_formats#Symmetry]] page '
                           'for a description of the symmetry format accepted by Xmipp') 
        
        form.addParam('angularSampling', FloatParam, default=5,
                      label="Angular Sampling (degrees)",  
                      help='Angular distance (in degrees) between neighboring projection points ')

        form.addParam('alpha', FloatParam, default=0.01,
                      label="Significance value",  
                      help='Parameter to define the corresponding most similar volume \n' 
                      '    projected images for each projection image')
        
        form.addParallelSection(threads=1, mpi=8)
        
       

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep', self.inputTiltPair.get().getUntilted().getObjId())
        self._insertFunctionStep('significantStep', 1)      
        self._insertFunctionStep('significantStep', 0)
            
        #self._insertFunctionStep('createOutputStep', s
        #                         prerequisites=deps)
        
    def convertInputStep(self, particlesId):
        """ Write the input images as a Xmipp metadata file. 
        particlesId: is only need to detect changes in
        input particles and cause restart from here.<
        """
        writeSetOfParticles(self.inputTiltPair.get().getUntilted(), self._getPath('input_untilted_particles.xmd'))
        writeSetOfParticles(self.inputTiltPair.get().getTilted(), self._getPath('input_tilted_particles.xmd'))
                    
    def significantStep(self, isTilt):
        if isTilt == 1:
            params =  '  -i %s' % self._getPath('input_tilted_particles.xmd')
            outputDir = self._getExtraPath("tilted")
        elif isTilt == 0:
            params =  '  -i %s' % self._getPath('input_untilted_particles.xmd')
            outputDir = self._getExtraPath("untilted")
                    
        firstImage = self.inputTiltPair.get().getUntilted().getFirstItem()
        maxShift = 0.35 * firstImage.getDim()[0]
        
        params += ' --sym %s' % self.symmetryGroup.get()
        params += ' --alpha0 %f --alphaF %f' % ((1-self.alpha.get())/100,(1-self.alpha.get())/100)
        params += ' --angularSampling %f' % self.angularSampling.get()
        params += ' --dontReconstruct'
        params += ' --iter 1'
        params += ' --initvolumes %s' % getImageLocation(self.inputVolume.get())  
        params += ' --useForValidation 3'  
        params += ' --maxShift %f' % maxShift
        params += ' --odir %s' % outputDir

        makePath(outputDir)  
        self.runJob('xmipp_reconstruct_significant', params)
        
    def validationStep(self, volName, volDir,sym):
        makePath(volDir)                  
        nproc = self.numberOfMpi.get()
        nT=self.numberOfThreads.get() 
        params = '  --volume %s' % volName  
        params += ' --odir %s' % volDir
        params += ' --sym %s' % sym
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
        #self._defineTransformRelation(self.inputVolumes.get(), volume)
       
                    
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        validateMsgs = []
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
