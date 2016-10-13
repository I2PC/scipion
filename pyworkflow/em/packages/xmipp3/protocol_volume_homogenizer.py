# **************************************************************************
# *
# * Authors:     Mohsen Kazemi  (mkazemi@cnb.csic.es)
# *              
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

import pyworkflow.protocol.params as params
from pyworkflow.em.protocol import ProtProcessParticles
from pyworkflow.em.packages.xmipp3.convert import (writeSetOfParticles, 
                                                   readSetOfParticles)
import numpy as np


class XmippProtVolumeHomogenizer(ProtProcessParticles):
    
    """    
    Method to get two volume from different classes (with different conformation)
    and correcting (deforming) all images of one of the volumes (input volume) 
    with respect to the another one as a reference, using optical flow algorithm.
    The output is a setOfParticles contaied deformed particles merged with 
    reference particles.
    """
    
    _label = 'volume homogenizer'    
    #--------------------------- DEFINE param functions --------------------------------------------   
   
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('referenceVolume', params.PointerParam,
                 pointerClass='Volume',
                 label='Reference volume', 
                 help="This is the volume that will be used as the reference "
                      "in OF algorithm.")
        form.addParam('referenceParticles', params.PointerParam, 
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      label="Reference particles",  
                      help="Aligned particles related to the reference "
                           "volume.")
        form.addParam('inputVolume', params.PointerParam,
                 pointerClass='Volume',
                 label='Input volume', 
                 help="Volume that we want to deform its related particles.")
        form.addParam('inputParticles', params.PointerParam, 
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      label="Input particles",  
                      help="Aligned particles related to the input volume. "
                           "These particles will be deformed (corrected) "
                           "based on the reference volume using OF algorithm. "
                           "Deformed particles will be merged with "
                           "reference particles inside the protocol to "
                           "create outputParticles.")        
        form.addParam('doAlignment', params.BooleanParam, default=False,
                      label='Reference and input volumes need to be aligned?',
                      help="Input and reference volumes must be aligned. "
                           "If you have not aligned them before choose this "
                           "option, so protocol will handle it internally.")
        form.addParam('cutOffFrequency', params.FloatParam, default = -1,
                      label="Cut-off frequency",
                      help="This digital frequency is used to filter both "
                           "input and reference voluem."
                           "IF it is (-1), cut-off frequency will be based on "
                           "Nyquist theorem.\n"
                           "Note:\n"
                           "Based on the experimental results, the best value "
                           "for cut-off frequency is 20A "
                           "(digitalFrequency = samplingRate/20)")
        form.addParam('winSize', params.IntParam, default=50,
                       label="Window size",
                       expertLevel=params.LEVEL_ADVANCED,
                       help="Size of the search window at each pyramid level "
                            "(shifts are assumed to be constant "
                            "within this window).")          
                      
        form.addParallelSection(threads=1, mpi=2)
    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):            
        inputParticlesMd = self._getExtraPath('input_particles.xmd')
        inputParticles = self.inputParticles.get()
        writeSetOfParticles(inputParticles, inputParticlesMd)        
        inputVol = self.inputVolume.get().getFileName()
        referenceVol = self.referenceVolume.get().getFileName()
        
        if not self.doAlignment.get():
            self._insertFunctionStep('opticalFlowAlignmentStep', 
                                     inputVol, referenceVol, inputParticlesMd)
        else:
            fnAlignedVolFf = self._getExtraPath('aligned_inputVol_to_refVol_FF.vol')
            fnAlnVolFfLcl = self._getExtraPath('aligned_FfAlnVol_to_refVol_lcl.vol')
            fnInPartsNewAng = self._getExtraPath("inputParts_anglesModified.xmd")
            self._insertFunctionStep('volumeAlignmentStep', 
                                     referenceVol, inputVol, fnAlignedVolFf, 
                                     fnAlnVolFfLcl, fnInPartsNewAng)
            self._insertFunctionStep('opticalFlowAlignmentStep', 
                                     fnAlnVolFfLcl, referenceVol, fnInPartsNewAng)
                        
        self._insertFunctionStep('createOutputStep')        
    #--------------------------- STEPS functions --------------------------------------------
    
    def volumeAlignmentStep (self, referenceVol, inputVol, fnAlignedVolFf, 
                                   fnAlnVolFfLcl, fnInPartsNewAng):
        fnTransformMatFF = self._getExtraPath('transformation-matrix-FF.txt')
        alignArgsFf = "--i1 %s --i2 %s --apply %s" % (referenceVol, 
                                                      inputVol, 
                                                      fnAlignedVolFf)
        alignArgsFf += " --frm --copyGeo %s" % fnTransformMatFF
        self.runJob('xmipp_volume_align', alignArgsFf, numberOfMpi = 1)        
        
        fnTransformMatLocal = self._getExtraPath('transformation-matrix-local.txt')
        alignArgsLocal = "--i1 %s --i2 %s --apply %s" % (referenceVol, 
                                                         fnAlignedVolFf, 
                                                         fnAlnVolFfLcl)
        alignArgsLocal += " --local --rot 0 0 1 --tilt 0 0 1"
        alignArgsLocal += " --psi 0 0 1 -x 0 0 1 -y 0 0 1"
        alignArgsLocal += " -z 0 0 1 --scale 1 1 0.005 --copyGeo %s" % fnTransformMatLocal        
        self.runJob('xmipp_volume_align', alignArgsLocal, numberOfMpi = 1)
                 
        #apply transformation matrix to input ratricles
        transMatFromFileFF = np.loadtxt(fnTransformMatFF)
        transformationArrayFF = np.reshape(transMatFromFileFF,(4,4))
        transformMatrixFF = np.matrix(transformationArrayFF)
        
        transMatFromFileLocal = np.loadtxt(fnTransformMatLocal)
        transformationArrayLocal = np.reshape(transMatFromFileLocal,(4,4))
        transformMatrixLocal = np.matrix(transformationArrayLocal)
        
        transformMatrix = transformMatrixFF * transformMatrixLocal
        print "Transformation matrix used to apply to the input particles =\n"
        print transformMatrix
        inputParts = self.inputParticles.get()
        outputSet = self._createSetOfParticles()
        for part in inputParts.iterItems():        
            partTransformMat = part.getTransform().getMatrix()            
            partTransformMatrix = np.matrix(partTransformMat)            
            newTransformMatrix = transformMatrix * partTransformMatrix                      
            part.getTransform().setMatrix(newTransformMatrix)
            outputSet.append(part)       
        outputSet.copyInfo(inputParts)                
        writeSetOfParticles(outputSet, fnInPartsNewAng)             
        
    def opticalFlowAlignmentStep(self, inputVol, referenceVol, inputParticlesMd):
        winSize = self.winSize.get()
        if self.cutOffFrequency.get() == -1:
            cutFreq = 0.5
        else:
            cutFreq = self.cutOffFrequency.get()
        fnOutput = self._getExtraPath('deformed-particles')        
          
        self.runJob("xmipp_volume_homogenizer", 
                    "-i %s -ref %s -img %s -o %s --winSize %d --cutFreq %f" % (
                    inputVol, referenceVol, inputParticlesMd, 
                    fnOutput, winSize, cutFreq), 
                    numberOfMpi=self.numberOfMpi.get()*self.numberOfThreads.get())        
    
    def createOutputStep(self):        
        refernceParticlesMd = self._getExtraPath('reference_particles.xmd')
        referenceParticles = self.referenceParticles.get()
        writeSetOfParticles(referenceParticles, refernceParticlesMd)
        inputParticles = self.inputParticles.get()
        
        fnDeformedParticles = self._getExtraPath('deformed-particles.xmd')
        fnOutputParticles = self._getExtraPath('OutputParticles_merged.xmd')
        self.runJob("xmipp_metadata_utilities", 
                    '-i %s -o %s -s union_all %s' % (
                    refernceParticlesMd, fnOutputParticles, fnDeformedParticles),
                    numberOfMpi = 1)
        self.runJob("xmipp_metadata_utilities", 
                    '-i %s -o %s -l itemId lineal 1 1' % (
                    fnOutputParticles, fnOutputParticles), numberOfMpi = 1)

        outputSetOfParticles = self._createSetOfParticles()
        readSetOfParticles(fnOutputParticles, outputSetOfParticles)        
        outputSetOfParticles.copyInfo(inputParticles)        
        self._defineOutputs(outputParticles=outputSetOfParticles)              
    #--------------------------- INFO functions -------------------------------------------- 
    
    def _summary(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        summary = []
        inputSize = self.inputParticles.get().getSize()
        referenceSize = self.referenceParticles.get().getSize()        
        if not hasattr(self, 'outputParticles'):
            summary.append("Output particles not ready yet.")
        else:
            summary.append("Applied OF to deform %d particles.\n" % inputSize)
            summary.append("Output has %d particles as it is "
                           "the result of merging deformed input "
                           "particles and reference ones." %
                           (inputSize+referenceSize))
        return summary
    
    def _methods(self):
        messages = []
        if not hasattr(self, 'outputParticles'):
            messages.append("Output particles not ready yet.")
        else:
            messages.append("We deformed %s particles from %s and produced %s."
                    %(self.inputParticles.get().getSize(), 
                      self.getObjectTag('inputParticles'), 
                      self.getObjectTag('outputParticles')))
        return messages
    
    def _citations(self):
        return ['**********????????????????????************']
    
    def _validate(self):
        errors=[]
        inputVolDim = self.inputVolume.get().getDim()[0]
        inputParticlesDim = self.inputParticles.get().getDim()[0]
        referenceVolDim = self.referenceVolume.get().getDim()[0]
        if inputVolDim != referenceVolDim:
            errors.append("Input and reference maps must have the "
                          "same dimensions!!!") 
        if inputParticlesDim != inputVolDim:
            errors.append("Input particles and input map do not have "
                          "the same dimensions!!!")
        return errors              
    