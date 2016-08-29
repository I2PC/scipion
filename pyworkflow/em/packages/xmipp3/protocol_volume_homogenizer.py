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
from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles, readSetOfParticles
from pyworkflow.em.data import SetOfParticles



class XmippProtVolumeHomogenizer(ProtProcessParticles):
    """    
    Method to get two volume from different classes (with different conformation)
    and correcting all images of one of the volume (input volume) with respect to
    the another one as a reference, using optical flow algorithm.
    This is to later merging the corrected images (outputParticles) to the images
    of the reference map to reconstruct a volume with better resolution
    """
    _label = 'volume homogenizer'    
    #--------------------------- DEFINE param functions --------------------------------------------   
   
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('referenceVolume', params.PointerParam,
                 pointerClass='Volume',
                 label='Reference volume', 
                 help="This is the volume that its resolution "
                      "is going to be improved.")
        form.addParam('inputVolume', params.PointerParam,
                 pointerClass='Volume',
                 label='Input volume', 
                 help="volume with the better resolution. "
                      "Input and reference volume are from "
                      "two different class.")
        form.addParam('inputParticles', params.PointerParam, 
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      label="Input particles",  
                      help="Aligned particles related to the input volume "
                           "(not the reference volume!).")
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
        self._insertFunctionStep('opticalFlowAlignmentStep')                   
        self._insertFunctionStep('createOutputStep')        
    #--------------------------- STEPS functions --------------------------------------------
    
    def opticalFlowAlignmentStep(self):
        particlesMd = self._getExtraPath('input_particles.xmd')
        inputParticles = self.inputParticles.get()
        writeSetOfParticles(inputParticles, particlesMd)
        
        inputVol = self.inputVolume.get().getFileName()
        referenceVol = self.referenceVolume.get().getFileName()
        winSize = self.winSize.get()
        if self.cutOffFrequency.get() == -1:
            cutFreq = 0.5
        else:
            cutFreq = self.cutOffFrequency.get()
        fnOutput = self._getExtraPath('OutputParticles_modified')
            
        self.runJob("xmipp_volume_homogenizer", 
                    "-i %s -ref %s -img %s -o %s --winSize %d --cutFreq %f" % (
                    inputVol, referenceVol, particlesMd, 
                    fnOutput, winSize, cutFreq), 
                    numberOfMpi=self.numberOfMpi.get()*self.numberOfThreads.get())
           
    def createOutputStep(self):        
        inputParticles = self.inputParticles.get()
        fnOutputParticles = self._getExtraPath('OutputParticles_modified.xmd')
        
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
        if not hasattr(self, 'outputParticles'):
            summary.append("Output particles not ready yet.")
        return summary
    
    def _methods(self):
        messages = []
        messages.append('********************')
        return messages
    
    def _citations(self):
        return ['********']
    
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
    