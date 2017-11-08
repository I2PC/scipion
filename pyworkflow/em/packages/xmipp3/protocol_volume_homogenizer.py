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
                                                   matrixFromGeometry)
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
                 help="This is the volume that will be used as the reference "
                      "in OF algorithm. If you want to use Gold-Standard "
                      "provide here half1 map")
        
        form.addParam('referenceVolume1', params.PointerParam, 
                 condition='doGoldStandard',
                 pointerClass='Volume',
                 label='Reference volume half1', 
                 help="This is the half1 volume that will be used as the reference"
                      "for half1 in OF algorithm.")
        
        form.addParam('referenceVolume2', params.PointerParam,
                 condition='doGoldStandard',
                 pointerClass='Volume',
                 label='Reference volume half2', 
                 help="This is half2 volume that will be used as the reference"
                      "for half2 in OF algorithm.")
        
        form.addParam('inputVolume', params.PointerParam,
                 condition='not doGoldStandard',
                 pointerClass='Volume',
                 label='Input volume', 
                 help="Volume that we want to process its related particles.")

        form.addParam('inputVolume1', params.PointerParam,
                 condition='doGoldStandard',
                 pointerClass='Volume',
                 label='Input volume half1', 
                 help="Volume that we want to process its related particles."
                       "It should represent half1 map.")

        form.addParam('inputVolume2', params.PointerParam,
                 condition='doGoldStandard',
                 pointerClass='Volume',
                 label='Input volume half2', 
                 help="Volume that we want to process its related particles."
                       "It should represent half2 map.")
                
        form.addParam('inputParticles', params.PointerParam,
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      label="Input particles",  
                      help="Aligned particles related to the input volume. "
                           "These particles will be processed (deformed) "
                           "based on the reference volume using OF algorithm."
                           "If selected doGoldStandard True the particles have to have"
                           "information about the halfId they belong.")
        
        form.addParam('doAlignment', params.BooleanParam, default=False,
                      label='Reference and input volumes need to be aligned?',
                      help="Input and reference volumes must be aligned. "
                           "If you have not aligned them before choose this "
                           "option, so protocol will handle it internally.")

        form.addParam('resLimit', params.FloatParam, default = 20,
                      label="Resolution Limit (A)",
                      help="Resolution limit used to low pass filter both "
                           "input and reference map(s)."
                           "Based on previous experimental results, a good value "
                           "for  seems to be 20A ")
        
        form.addParam('winSize', params.IntParam, default=50,
                       label="Window size",
                       expertLevel=params.LEVEL_ADVANCED,
                       help="Size of the search window at each pyramid level "
                            "(shifts are assumed to be constant "
                            "within this window).")
                      
        form.addParallelSection(threads=1, mpi=2)
    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
              
        #If doGoldStandard then we have two halves
        if self.doGoldStandard.get():            
            inputParticlesMd1 = self._getExtraPath('input_particles_half1.xmd')
            inputParticles = self.inputParticles.get()        
            inputParticlesHalf1 = self._createSetOfParticles()            
            inputParticlesHalf1.copyInfo(inputParticles)
            inputParticlesHalf1.copyItems(inputParticles,
                                 updateItemCallback=self._setHalf1)                                
            writeSetOfParticles(inputParticlesHalf1, inputParticlesMd1)
    
            inputParticlesMd2 = self._getExtraPath('input_particles_half2.xmd')
            inputParticles = self.inputParticles.get()        
            inputParticlesHalf2 = self._createSetOfParticles()            
            inputParticlesHalf2.copyInfo(inputParticles)
            inputParticlesHalf2.copyItems(inputParticles,
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
                                         inputVol1, referenceVol1, inputParticlesMd1,fnOutputHalf1)                
                self._insertFunctionStep('opticalFlowAlignmentStep', 
                                         inputVol2, referenceVol2, inputParticlesMd2,fnOutputHalf2)

            else:
                
                print("Nothing")
                '''
                fnAlignedVolFf = self._getExtraPath('aligned_inputVol_to_refVol_FF.vol')
                fnAlnVolFfLcl = self._getExtraPath('aligned_FfAlnVol_to_refVol_lcl.vol')                
                fnInPartsNewAng = self._getExtraPath("inputParts_anglesModified.xmd")
                self._insertFunctionStep('volumeAlignmentStep', 
                                         inputVol1, referenceVol1, fnAlignedVolFf, 
                                         fnAlnVolFfLcl, fnInPartsNewAng)                
                self._insertFunctionStep('opticalFlowAlignmentStep', 
                                         inputVol1, referenceVol1, inputParticlesMd1)                
                self._insertFunctionStep('opticalFlowAlignmentStep', 
                                         inputVol2, referenceVol2, inputParticlesMd2)

                '''
        else:
            
            inputParticlesMd = self._getExtraPath('input_particles.xmd')
            inputParticles = self.inputParticles.get()                
            writeSetOfParticles(inputParticles, inputParticlesMd)
            
            inputVol = self.changeExtension(self.inputVolume.get().getFileName())
            referenceVol = self.changeExtension(self.referenceVolume.get().getFileName())                
            fnOutput = self._getExtraPath('deformed_particles')
            
            if not self.doAlignment.get():
                self._insertFunctionStep('opticalFlowAlignmentStep', 
                                         inputVol, referenceVol, inputParticlesMd, fnOutput)
            else:
                
                fnInputToRefGlobal = self._getExtraPath('aligned_inputVol_to_refVol_Gobal.vol')
                fnInputToRefLocal = self._getExtraPath('aligned_inputVol_to_refVol_Local.vol')                
                fnInPartsNewAng = self._getExtraPath("inputparts_anglesModified.xmd")
                self._insertFunctionStep('volumeAlignmentStep', 
                                         inputVol, referenceVol, fnInputToRefGlobal, 
                                         fnInputToRefLocal, fnInPartsNewAng)                
                self._insertFunctionStep('opticalFlowAlignmentStep', 
                                         fnInputToRefLocal, referenceVol, fnInPartsNewAng, fnOutput)                
#TODO Ubset of particles according to their HalfID!!!! 
#Now this code does not work!!!!!

                        
        self._insertFunctionStep('createOutputStep')        
    #--------------------------- STEPS functions --------------------------------------------
        
    def changeExtension(self, vol):
            extVol = getExt(vol)
            if (extVol == '.mrc') or (extVol == '.map'):
                vol = vol + ':mrc'
            return vol
            
    def volumeAlignmentStep (self, inputVol, referenceVol, fnInputToRefGlobal, 
                                   fnInputToRefLocal, fnInPartsNewAng):
        '''The input vol is modified towards the reference vol first by a global alignment and then by a local one'''
        '''The particles orientations are changed accordingly'''               
        
        fnTransformMatGlobal = self._getExtraPath('transformation-matrix-Global.txt')
        fnTransformMatLocal = self._getExtraPath('transformation-matrix-Local.txt')

        alignArgsGlobal = "--i1 %s --i2 %s --apply %s --dontScale  --frm --show_fit --copyGeo %s" % (referenceVol,
                                                                                                     inputVol,                                                              
                                                                                                     fnInputToRefGlobal,                                                                         
                                                                                                     fnTransformMatGlobal)
 
        alignArgsLocal = "--i1 %s --i2 %s --apply %s --show_fit  --local --dontScale --copyGeo %s " % (referenceVol, 
                                                                                                       fnInputToRefGlobal, 
                                                                                                       fnInputToRefLocal,
                                                                                                       fnTransformMatLocal)

        inputParts = self.inputParticles.get()
        outputSet = self._createSetOfParticles()
                                
        # Obtain information about the Global transformation                                
        self.runJob('xmipp_volume_align', alignArgsGlobal, numberOfMpi=1, numberOfThreads=1)
        transMatFromFileFF = np.loadtxt(fnTransformMatGlobal)
        transformationArrayFF = np.reshape(transMatFromFileFF, (4, 4))
        transformMatrixFF = np.matrix(transformationArrayFF)
        shifts, angles = geometryFromMatrix(transformMatrixFF, False) 
        print("Global transformation to be applied: ")
        print(transformMatrixFF)
        print("Shifts and angles to be applied: ")
        print(shifts, angles)
        matG = matrixFromGeometry(shifts, angles,False)
        print matG
        print("\n ")        
        #alignArgsLocal += "--rot %s --tilt %s --psi %s -x %s -y %s -z %s" % (angles[0], angles[1], angles[2],
        #                                                                     shifts[0],shifts[1],shifts[2])
                                                        
        # Obtain information about the Local transformation       
        self.runJob('xmipp_volume_align', alignArgsLocal, numberOfMpi=1)
        transMatFromFileLocal = np.loadtxt(fnTransformMatLocal)
        transformationArrayLocal = np.reshape(transMatFromFileLocal, (4, 4))
        transformMatrixLocal = np.matrix(transformationArrayLocal)
        shifts, angles = geometryFromMatrix(transformMatrixLocal, False)
        print("Local transformation to be applied: ")
        print(transformMatrixLocal)
        print("Shifts and angles to be applied: ")
        print(shifts, angles)
        matL = matrixFromGeometry(shifts, angles,False)
        print matL 
        print("\n ")
            
        for part in inputParts.iterItems():  
            part.getTransform().composeTransform(transformMatrixFF)
            part.getTransform().composeTransform(transformMatrixLocal)
            outputSet.append(part)   
                                                   
        outputSet.copyInfo(inputParts)                
        writeSetOfParticles(outputSet, fnInPartsNewAng)             
        
    def opticalFlowAlignmentStep(self, inputVol, referenceVol, inputParticlesMd, fnOutput):
        winSize = self.winSize.get()
        
        sampling_rate = self.inputParticles.get().getSamplingRate()
        resLimitDig = sampling_rate/self.resLimit.get()
        
        nproc = self.numberOfMpi.get()
        nT=self.numberOfThreads.get()
        
        self.runJob("xmipp_volume_homogenizer", 
                    "-i %s -ref %s -img %s -o %s --winSize %d --cutFreq %f" % (
                    inputVol, referenceVol, inputParticlesMd, 
                    fnOutput, winSize, resLimitDig), 
                    numberOfMpi=nproc,numberOfThreads=nT)
            
    def createOutputStep(self):
        inputParticles = self.inputParticles.get()        
        inputClassName = str(inputParticles.getClassName())
        
        if not self.doGoldStandard.get():            
            fnDeformedParticles = self._getExtraPath('deformed_particles.xmd')
            outputSetOfParticles = self._createSetOfParticles()
            readSetOfParticles(fnDeformedParticles, outputSetOfParticles)        
            outputSetOfParticles.copyInfo(inputParticles)        
            self._defineOutputs(outputParticles=outputSetOfParticles)              
        else:
            fnDeformedParticlesHalf1 = self._getExtraPath('deformed_particles_half1.xmd')
            outputSetOfParticlesHalf1 = self._createSetOfParticles()            
            readSetOfParticles(fnDeformedParticlesHalf1, outputSetOfParticlesHalf1)                                
            outputSetOfParticlesHalf1.copyInfo(inputParticles)                      
            key = 'output' + inputClassName.replace('SetOf', '') + '%02d'  
            self._defineOutputs(**{key % 1: outputSetOfParticlesHalf1})
            
            fnDeformedParticlesHalf2 = self._getExtraPath('deformed_particles_half2.xmd')
            outputSetOfParticlesHalf2 = self._createSetOfParticles()                                  
            readSetOfParticles(fnDeformedParticlesHalf2, outputSetOfParticlesHalf2)
            outputSetOfParticlesHalf2.copyInfo(inputParticles)            
            self._defineOutputs(**{key % 2: outputSetOfParticlesHalf2})        
            
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
        if (item._rln_halfId == 1):
            item._appendItem=False

    def _setHalf2(self, item, row):
        if (item._rln_halfId == 2):
            item._appendItem=False
            
    def _validate(self):
        errors=[]
        '''inputVolDim = self.inputVolume.get().getDim()[0]
        inputParticlesDim = self.inputParticles.get().getDim()[0]
        referenceVolDim = self.referenceVolume.get().getDim()[0]
        if inputVolDim != referenceVolDim:
            errors.append("Input and reference maps must have the "
                          "same dimensions!!!") 
        if inputParticlesDim != inputVolDim:
            errors.append("Input particles and input map do not have "
                          "the same dimensions!!!")
        '''
        return errors              
#TODO: Validate that if use gold-standard the particles have halfID metadata info
#TODO: Validate that if use gold-standard the particles have halfID metadata.

    
