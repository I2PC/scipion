# **************************************************************************
# *
# * Authors:     Mohsen Kazemi  (mkazemi@cnb.csic.es)
# *              Javier Vargas  (javier.vargasbalbuena@mcgill.ca)  
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import pyworkflow.protocol.params as params
from pyworkflow import VERSION_1_1
from pyworkflow.em.protocol import ProtProcessParticles
from pyworkflow.em.packages.xmipp3.convert import (writeSetOfParticles, 
                                                   readSetOfParticles,
                                                   geometryFromMatrix,
                                                   SetOfParticles)
from pyworkflow.utils import getExt
import numpy as np

class XmippProtVolumeHomogenizer(ProtProcessParticles):
    
    """    
    Method to get two volume from different classes (with different conformation)
    and correcting (deforming) all images of one of the volumes (input volume) 
    with respect to the another one as a reference, using optical flow algorithm.
    The output is a setOfParticles contaied deformed reference particles.
    """    
    _label = 'volume homogenizer'
    _lastUpdateVersion = VERSION_1_1
    #--------------------------- DEFINE param functions --------------------------------------------   
   
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('doGoldStandard', params.BooleanParam, default=False,
                      label='do Gold Standard?',
                      help="If YES provide half1 and half2 maps for reference"
                           "and for input volumes.")
        form.addParam('referenceVolume', params.PointerParam,
                      condition='not doGoldStandard',
                      pointerClass='Volume',
                      label='Reference volume',
                      help="This is the volume that will be used as the "
                           "reference in OF algorithm. If you want to "
                           "use Gold-Standard provide here half1 map")
        form.addParam('referenceVolume1', params.PointerParam,
                      condition='doGoldStandard',
                      pointerClass='Volume',
                      label='Reference volume half1',
                      help="This is the half1 volume that will be used as the "
                           "reference for half1 in OF algorithm.")
        form.addParam('referenceVolume2', params.PointerParam,
                      condition='doGoldStandard',
                      pointerClass='Volume',
                      label='Reference volume half2',
                      help="This is half2 volume that will be used as the "
                           "reference for half2 in OF algorithm.")
        form.addParam('inputVolume', params.PointerParam,
                      condition='not doGoldStandard',
                      pointerClass='Volume',
                      label='Input volume',
                      help="Volume that we want to process its related "
                           "particles.")
        form.addParam('inputVolume1', params.PointerParam,
                      condition='doGoldStandard',
                      pointerClass='Volume',
                      label='Input volume half1', 
                      help="Volume that we want to process its related"
                           " particles. It should represent half1 map.")
        form.addParam('inputVolume2', params.PointerParam,
                      condition='doGoldStandard',
                      pointerClass='Volume',
                      label='Input volume half2', 
                      help="Volume that we want to process its related "
                           "particles. It should represent half2 map.")                
        form.addParam('inputParticles', params.PointerParam,
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      label="Input particles",  
                      help="Aligned particles related to the input volume. "
                           "These particles will be processed (deformed) "
                           "based on the reference volume using OF algorithm."
                           "If selected doGoldStandard True the particles have"
                           " to have  information about the halfId they "
                           "belong.")        
        form.addParam('doAlignment', params.BooleanParam, default=False,
                      label='Reference and input volumes need to be aligned?',
                      help="Input and reference volumes must be aligned. "
                           "If you have not aligned them before choose this "
                           "option, so protocol will handle it internally.")        
        form.addParam('dofrm', params.BooleanParam, default=True,
                      condition='doAlignment',
                      label='Use Fast Rotational Matching.',
                      help="Use Fast Rotational Matching. Before use it you "
                           "have to install it by scipion install frm"
                           "This method for volume alignment is much more "
                           "fast than exhaustive search")        
        form.addParam('resLimit', params.FloatParam, default=20,
                      label="Resolution Limit (A)",
                      help="Resolution limit used to low pass filter both "
                           "input and reference map(s)."
                           "Based on previous experimental results, a good "
                           "value for  seems to be 20A ")                
        form.addSection(label='OF')
        form.addParam('winSize', params.IntParam, default=50,
                      label="Window size",
                      help="Size of the search window at each pyramid level "
                            "(shifts are assumed to be constant within this "
                            "window).")                   
        form.addParam('pyrScale', params.FloatParam, default=0.5,
                      label="pyramid Scale",
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Parameter, specifying the image scale (<1) to "
                           "build pyramids for each image. pyrScale=0.5 " 
                            "means a classical pyramid, where each next layer" 
                            " is twice smaller than the previous one.")        
        form.addParam('levels', params.IntParam, default=2,
                      label="Number of Levels",
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Number of pyramid layers including the initial "
                           "image; levels=1 means that no extra layers are "
                           "created and only the original images are used.")
        form.addParam('iterations', params.IntParam, default=10,
                      label="Iterations",
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Number of iterations the algorithm does at "                             
                           "each pyramid level.")                      
        form.addParallelSection(threads=1, mpi=2)
        
    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):              
        #If doGoldStandard then we have two halves
        inputPart = self.inputParticles.get()        
        if self.doGoldStandard.get():            
            inputParticlesMd1 = self._getExtraPath('input_particles_half1.xmd')
            
            inputParticlesHalf1 = SetOfParticles(filename=':memory:')            
            inputParticlesHalf1.copyInfo(inputPart)
            inputParticlesHalf1.copyItems(inputPart,
                                 updateItemCallback=self._setHalf1)                                
            writeSetOfParticles(inputParticlesHalf1, inputParticlesMd1)
    
            inputParticlesMd2 = self._getExtraPath('input_particles_half2.xmd')
            inputParticlesHalf2 = SetOfParticles(filename=':memory:')          
            inputParticlesHalf2.copyInfo(inputPart)
            inputParticlesHalf2.copyItems(inputPart,
                                 updateItemCallback=self._setHalf2)                                
            writeSetOfParticles(inputParticlesHalf2, inputParticlesMd2)
            
            inputVol1 = self.changeExtension(self.inputVolume1.get().getFileName())
            inputVol2 = self.changeExtension(self.inputVolume2.get().getFileName())
            referenceVol1 = self.changeExtension(self.referenceVolume1.get().getFileName())
            referenceVol2 = self.changeExtension(self.referenceVolume2.get().getFileName())

            if not self.doAlignment.get():         #No alignment
                fnOutputHalf1 = self._getExtraPath('deformed_particles_half1')
                fnOutputHalf2 = self._getExtraPath('deformed_particles_half2')
                self._insertFunctionStep('opticalFlowAlignmentStep', 
                                         inputVol1, referenceVol1, inputParticlesMd1, fnOutputHalf1)                
                self._insertFunctionStep('opticalFlowAlignmentStep', 
                                         inputVol2, referenceVol2, inputParticlesMd2, fnOutputHalf2)

            else:
                               
                fnInputToRefLocal1 = self._getExtraPath('aligned_inputVol_to_refVol_Local1.vol')
                fnInputToRefLocal2 = self._getExtraPath('aligned_inputVol_to_refVol_Local2.vol')
                
                fnInPartsNewAng1 = self._getExtraPath("inputparts_anglesModified1.xmd")
                fnInPartsNewAng2 = self._getExtraPath("inputparts_anglesModified2.xmd")
                
                self._insertFunctionStep('volumeAlignmentStep', 
                                         inputVol1, referenceVol1, inputParticlesMd1, 
                                         fnInputToRefLocal1, fnInPartsNewAng1)      

                self._insertFunctionStep('volumeAlignmentStep', 
                                         inputVol2, referenceVol2, inputParticlesMd2,
                                         fnInputToRefLocal2, fnInPartsNewAng2)                
                
                fnOutputHalf1 = self._getExtraPath('deformed_particles_half1')
                fnOutputHalf2 = self._getExtraPath('deformed_particles_half2')
                
                self._insertFunctionStep('opticalFlowAlignmentStep', 
                                         fnInputToRefLocal1, referenceVol1, fnInPartsNewAng1,fnOutputHalf1)                
                self._insertFunctionStep('opticalFlowAlignmentStep', 
                                         fnInputToRefLocal2, referenceVol2, fnInPartsNewAng2,fnOutputHalf2)
      
        else:
            
            inputParticlesMd = self._getExtraPath('input_particles.xmd')
            writeSetOfParticles(inputPart, inputParticlesMd)
            
            inputVol = self.changeExtension(self.inputVolume.get().getFileName())
            referenceVol = self.changeExtension(self.referenceVolume.get().getFileName())                
            fnOutput = self._getExtraPath('deformed_particles')
            
            if not self.doAlignment.get():
                self._insertFunctionStep('opticalFlowAlignmentStep', 
                                         inputVol, referenceVol, inputParticlesMd, fnOutput)
            else:
                
                fnInputToRefLocal = self._getExtraPath('aligned_inputVol_to_refVol_Local.vol')                
                fnInPartsNewAng = self._getExtraPath("inputparts_anglesModified.xmd")
                self._insertFunctionStep('volumeAlignmentStep', 
                                         inputVol, referenceVol, inputParticlesMd,  
                                         fnInputToRefLocal, fnInPartsNewAng)                
                self._insertFunctionStep('opticalFlowAlignmentStep', 
                                         fnInputToRefLocal, referenceVol, fnInPartsNewAng, fnOutput)                

                        
        self._insertFunctionStep('createOutputStep')
                
    #--------------------------- STEPS functions --------------------------------------------
        
    def changeExtension(self, vol):
            extVol = getExt(vol)
            if (extVol == '.mrc') or (extVol == '.map'):
                vol = vol + ':mrc'
            return vol
            
    def volumeAlignmentStep (self, inputVol, referenceVol, inputPartsMD,
                                   fnInputToRefLocal, fnInPartsNewAng):
        '''The input vol is modified towards the reference vol first by a global alignment and then by a local one'''
        '''The particles orientations are changed accordingly'''               
        
        fnTransformMatGlobal = self._getExtraPath('transformation-matrix-Global.txt')
        fnTransformMatLocal = self._getExtraPath('transformation-matrix-Local.txt')
        alignArgsGlobal = "--i1 %s --i2 %s --dontScale  --copyGeo %s" % (referenceVol,       
                                                                                            inputVol,                                                                                                                                                                                                                               
                                                                                            fnTransformMatGlobal)
        if (self.dofrm):
            alignArgsGlobal += " --frm --show_fit "                        
        else:
            alignArgsGlobal += " --rot 0.000000 360.000000 5.000000 --tilt 0.000000 180.000000 5.000000 --psi 0.000000 360.000000 5.000000 -x 0.000000 0.000000 1.000000 -y 0.000000 0.000000 1.000000 -z 0.000000 0.000000 1.000000"
                                       
        self.runJob('xmipp_volume_align', alignArgsGlobal, numberOfMpi=1, numberOfThreads=1)
        transMatFromFileFF = np.loadtxt(fnTransformMatGlobal)
        transformationArrayFF = np.reshape(transMatFromFileFF, (4, 4))
        transformMatrixFF = np.matrix(transformationArrayFF)
        shifts, angles = geometryFromMatrix(transformMatrixFF, False) 
        print("Global transformation to be applied: ")
        print(transformMatrixFF)
        print("Shifts and angles to be applied: ")
        print(shifts, angles)
        print("\n ")                
        
        #We calculate one local alignment
        alignArgsLocal = "--i1 %s --i2 %s --apply %s --show_fit  --local --dontScale --copyGeo %s " % (referenceVol, 
                                                                                                       inputVol, 
                                                                                                       fnInputToRefLocal,
                                                                                                       fnTransformMatLocal)
        alignArgsLocal += "--rot %s --tilt %s --psi %s -x %s -y %s -z %s" % (angles[0], angles[1], angles[2],
                                                                             shifts[0],shifts[1],shifts[2])                                                        

        self.runJob('xmipp_volume_align', alignArgsLocal, numberOfMpi=1)
        transMatFromFileLocal = np.loadtxt(fnTransformMatLocal)
        transformationArrayLocal = np.reshape(transMatFromFileLocal, (4, 4))
        transformMatrixLocal = np.matrix(transformationArrayLocal)
        shifts, angles = geometryFromMatrix(transformMatrixLocal, False)        
        print(shifts, angles)        
        print("\n ")
        
        #We calculate another local alignment
        alignArgsLocal = "--i1 %s --i2 %s --apply %s --show_fit  --local --dontScale --copyGeo %s " % (referenceVol, 
                                                                                                       inputVol, 
                                                                                                       fnInputToRefLocal,
                                                                                                       fnTransformMatLocal)        
        alignArgsLocal += "--rot %s --tilt %s --psi %s -x %s -y %s -z %s" % (angles[0], angles[1], angles[2],
                                                                             shifts[0],shifts[1],shifts[2])                                                        
        self.runJob('xmipp_volume_align', alignArgsLocal, numberOfMpi=1)
        transMatFromFileLocal = np.loadtxt(fnTransformMatLocal)
        transformationArrayLocal = np.reshape(transMatFromFileLocal, (4, 4))
        transformMatrixLocal = np.matrix(transformationArrayLocal)
        shifts, angles = geometryFromMatrix(transformMatrixLocal, False)
        print(shifts, angles)        
        print("\n ")

        #And one more.        
        alignArgsLocal = "--i1 %s --i2 %s --apply %s --show_fit  --local --dontScale --copyGeo %s " % (referenceVol, 
                                                                                                       inputVol, 
                                                                                                       fnInputToRefLocal,
                                                                                                       fnTransformMatLocal)        
        alignArgsLocal += "--rot %s --tilt %s --psi %s -x %s -y %s -z %s" % (angles[0], angles[1], angles[2],
                                                                             shifts[0],shifts[1],shifts[2])                                                        
        self.runJob('xmipp_volume_align', alignArgsLocal, numberOfMpi=1)
        transMatFromFileLocal = np.loadtxt(fnTransformMatLocal)
        transformationArrayLocal = np.reshape(transMatFromFileLocal, (4, 4))
        transformMatrixLocal = np.matrix(transformationArrayLocal)
        shifts, angles = geometryFromMatrix(transformMatrixLocal, False)
        print(shifts, angles)        
            
        print("Local transformation to be applied: ")
        print(transformMatrixLocal)
        print("Shifts and angles to be applied: ")
        print(shifts, angles)
        print("\n ")
                
        inputParts = SetOfParticles(filename=':memory:')
        inputParts.copyInfo(self.inputParticles.get())
        readSetOfParticles(inputPartsMD, inputParts)

        outputSet = SetOfParticles(filename=':memory:')
                                          
        for part in inputParts.iterItems():
            part.getTransform().composeTransform(np.matrix(transformMatrixLocal))
            outputSet.append(part)   
                                                   
        outputSet.copyInfo(inputParts)                
        writeSetOfParticles(outputSet, fnInPartsNewAng)             
        
    def opticalFlowAlignmentStep(self, inputVol, referenceVol, inputParticlesMd, fnOutput):
        winSize = self.winSize.get()
        pyrScale = self.pyrScale.get()
        levels = self.levels.get()
        iterations = self.iterations.get()
        
        sampling_rate = self.inputParticles.get().getSamplingRate()
        resLimitDig = sampling_rate/self.resLimit.get()
        
        nproc = self.numberOfMpi.get()
        nT=self.numberOfThreads.get()
        
        self.runJob("xmipp_volume_homogenizer", 
                    "-i %s -ref %s -img %s -o %s --winSize %d --cutFreq %f --pyr_scale %f --levels %d --iterations %d" % (
                    inputVol, referenceVol, inputParticlesMd, 
                    fnOutput, winSize, resLimitDig,
                    pyrScale, levels, iterations), 
                    numberOfMpi=nproc,numberOfThreads=nT)
            
    def createOutputStep(self):
        inputParticles = self.inputParticles.get()        
        inputClassName = str(inputParticles.getClassName())
        key = 'output' + inputClassName.replace('SetOf', '') + '%02d'  

        
        if not self.doGoldStandard.get():            
            fnDeformedParticles = self._getExtraPath('deformed_particles.xmd')
            outputSetOfParticles = self._createSetOfParticles()
            outputSetOfParticles.copyInfo(inputParticles)
            readSetOfParticles(fnDeformedParticles, outputSetOfParticles)
            self._defineOutputs(outputParticles=outputSetOfParticles)              
        else:
            fnDeformedParticlesHalf1 = self._getExtraPath('deformed_particles_half1.xmd')
            outputSetOfParticlesHalf1 = self._createSetOfParticles(suffix="1")            
            outputSetOfParticlesHalf1.copyInfo(inputParticles)
            readSetOfParticles(fnDeformedParticlesHalf1, outputSetOfParticlesHalf1)

            self._defineOutputs(**{key % 1: outputSetOfParticlesHalf1})            
            self._defineTransformRelation(inputParticles, outputSetOfParticlesHalf1)
            
            fnDeformedParticlesHalf2 = self._getExtraPath('deformed_particles_half2.xmd')
            outputSetOfParticlesHalf2 = self._createSetOfParticles(suffix="2")                                  
            outputSetOfParticlesHalf2.copyInfo(inputParticles)
            readSetOfParticles(fnDeformedParticlesHalf2, outputSetOfParticlesHalf2)

            self._defineOutputs(**{key % 2: outputSetOfParticlesHalf2})  
            self._defineTransformRelation(inputParticles, outputSetOfParticlesHalf2)
     
            
    #--------------------------- INFO functions -------------------------------------------- 
    
    def _summary(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        summary = []
        inputSize = self.inputParticles.get().getSize()        
        if not hasattr(self, 'outputParticles'):
            summary.append("Output particles are not ready yet.")
        else:
            summary.append("Applied OF to deform %d particles.\n" % inputSize)            
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
    
    def _setHalf1(self, item, row):
        if item._rlnRandomSubset == 1:
            item._appendItem = False

    def _setHalf2(self, item, row):
        if item._rlnRandomSubset == 2:
            item._appendItem = False
            
    def _validate(self):
        errors=[]
                 
        if self.doGoldStandard.get():
            
            inputVolDim1 = self.inputVolume1.get().getDim()[0]
            inputVolDim2 = self.inputVolume2.get().getDim()[0]
            
            inputParticlesDim = self.inputParticles.get().getDim()[0]
            referenceVolDim1 = self.referenceVolume1.get().getDim()[0]
            referenceVolDim2 = self.referenceVolume2.get().getDim()[0]
            
            if ( (inputVolDim1 != inputVolDim2) or  (referenceVolDim1 != referenceVolDim2) or  (inputVolDim1 != referenceVolDim1) or (inputParticlesDim != inputVolDim1)):
                errors.append("Incorrect dimensions of the input data") 
        else:
            
            inputVolDim = self.inputVolume.get().getDim()[0]
            inputParticlesDim = self.inputParticles.get().getDim()[0]
            referenceVolDim = self.referenceVolume.get().getDim()[0]
            
            if ( (inputVolDim != referenceVolDim) or (inputParticlesDim != inputVolDim)):
                errors.append("Incorrect dimensions of the input data") 

        return errors
    
