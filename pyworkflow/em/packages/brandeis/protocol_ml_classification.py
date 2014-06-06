# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (jgomez@cnb.csic.es)
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
This module contains the protocol to obtain a refined 3D reconstruction from a set of particles using Frealign
"""
import os
from pyworkflow.utils import *
from pyworkflow.em import *
from data import *
from brandeis import *
from constants import *
from protocol_frealign_base import ProtFrealignBase


class ProtFrealignClassify(ProtFrealignBase, ProtClassify3D):
    """ Protocol to classify 3D using Frealign. Frealign employs
maximum likelihood classification for single particle electron
cryo-microscopy.
Particle alignment parameters are determined by maximizing a
joint likelihood that can include hierarchical priors, while
classification is performed by expectation maximization of a
marginal likelihood.
    """
    _label = 'frealign classify'
    IS_REFINE = False
    
    def __init__(self, **args):
        ProtFrealignBase.__init__(self, **args)
    
  #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        """Insert the steps to refine orientations and shifts of the SetOfParticles
        """
        numberOfBlocks = self.numberOfThreads.get() - 1
        depsOcc = []
        
        if self.doContinue:
            continueRun = self.continueRun.get()
            self.inputParticles.set(continueRun.inputParticles.get())
            self.symmetry.set(continueRun.symmetry.get())
            self.input3DReference.set(None)
            self.numberOfRef = continueRun.numberOfClasses.get()
            if self.continueIter.get() == 'last':
                initIter = continueRun._getLastIter() + 1
            else:
                initIter = int(self.continueIter.get()) + 1
            self._setLastIter(initIter-1)
            self._insertFunctionStep('continueStep', initIter)
        else:
            initIter = 1
            self.numberOfRef = self.numberOfClasses.get()
        
        lastIter = initIter + self.numberOfIterations.get()
        cpuList = self._cpusPerClass(numberOfBlocks, self.numberOfRef)
        
        for iter in range(initIter, lastIter):
            depsRecons = []
            initId = self._insertFunctionStep('initIterStep', iter, numberOfBlocks, prerequisites=depsOcc)
            depsRefine = self._insertRefineIterStep(iter, numberOfBlocks, [initId])
            firstOccId = self._insertFunctionStep("calculateOCCStep", iter, False, numberOfBlocks, prerequisites=depsRefine)
            for ref in range(1, self.numberOfRef + 1):
                reconsId = self._insertFunctionStep("reconstructVolumeStep", iter, ref, cpuList, prerequisites=[firstOccId])
                depsRecons.append(reconsId)
            secondOccId = self._insertFunctionStep("calculateOCCStep", iter, True, numberOfBlocks, prerequisites=depsRecons)
            depsOcc = [secondOccId]
        self._insertFunctionStep("createOutputStep", lastIter-1, prerequisites=depsOcc)
    
    def _insertRefineIterStep(self, iter, numberOfBlocks, depsInitId):
        """ execute the refinement for the current iteration """
        
        depsRefine = []
        iterDir = self._iterWorkingDir(iter)
        
        if iter == 1:
            if not self.useInitialAngles.get():
                stepConstructId = self._insertFunctionStep("constructParamFilesStep", numberOfBlocks, iterDir, iter, prerequisites=depsInitId)
                depsConstruct = [stepConstructId]
                for block in range(1, numberOfBlocks + 1):
                    refineId = self._insertFunctionStep("refineBlockStep", iterDir, block, prerequisites=depsConstruct)
                    depsRefine.append(refineId)
            else:
                # ToDo: Construct the function to extract the coordinates and euler angles for each particle from SetOfParticles.
                pass
        else:
            blocks = self._cpusPerClass(numberOfBlocks, self.numberOfRef)
            cpuList = self._cpusPerClass(numberOfBlocks, self.numberOfRef)
            for ref in range(1, self.numberOfRef + 1):
                for block in range(1, blocks[ref-1] + 1):
                    restIterAngle = iter % self.itRefineAngles.get()
                    restIterShifts = iter % self.itRefineShifts.get()
                    if restIterAngle == 0 and restIterShifts == 0:
                        parRefine = 0
                    elif restIterAngle == 0:
                        parRefine = 1
                    elif restIterShifts == 0:
                        parRefine = 2
                    else:
                        parRefine = 3
                    refineId = self._insertFunctionStep("refineParticlesStep", iter, ref, block, parRefine, cpuList[ref-1], prerequisites=depsInitId)
                    depsRefine.append(refineId)
        return depsRefine
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def initIterStep(self, iter, numberOfBlocks):
        """ Prepare files and directories for the current iteration """
        
        self._createIterWorkingDir(iter) # create the working directory for the current iteration.
        prevIter = iter - 1
        self._createFilenameTemplates(iter)
        
        cpusRef = self._cpusPerClass(numberOfBlocks, self.numberOfRef)
            
        if iter==1:
            imgSet = self.inputParticles.get()
            vol = self.input3DReference.get()
            
            imgFn = self._getFileName('particles')
            volFn = self._getFileName('init_vol')
            refVol = self._getFileName('ref_vol', iter=iter) # reference volume of the step.
            imgSet.writeStack(imgFn) # convert the SetOfParticles into a mrc stack.
            ImageHandler().convert(vol.getLocation(), volFn) # convert the reference volume into a mrc volume
            copyFile(volFn, refVol)  #Copy the initial volume in the current directory.
            
        for ref in range(1, self.numberOfRef + 1):
            refVol = self._getFileName('ref_vol_class', iter=iter, ref=ref) # reference volume of the step.
            iterVol =  self._getFileName('iter_vol_class', iter=iter, ref=ref) # refined volumes of the step
            if iter == 1:
                copyFile(volFn, iterVol)  #Copy the initial volume in the current directory.
            else:
                self._splitParFile(iter, ref, cpusRef[ref-1])
                prevIterVol = self._getFileName('prev_vol_class', iter=prevIter, ref=ref) # volumes of the previous iteration
                copyFile(prevIterVol, refVol)   #Copy the reference volume as refined volume.
                copyFile(refVol, iterVol)   #Copy the reference volume as refined volume.
    
    def refineParticlesStep(self, iter, ref, block, parRefine, numberOfBlocks):
        """Only refine the parameters of the SetOfParticles
        """
        param = {}
        imgSet = self.inputParticles.get()
        
        iterDir = self._iterWorkingDir(iter)
        if ref ==1 and block==1:
            self._enterDir(iterDir) # enter to the working directory for the current iteration.
        
        iniPart, lastPart = self._particlesInBlock(block, numberOfBlocks)
        prevIter = iter - 1
        param['inputParFn'] = self._getFile('input_par_block_class', iter=prevIter, ref=ref, block=block)
        param['initParticle'] = iniPart
        param['finalParticle'] = lastPart
        
        paramDic = self._setParamsClassRefineParticles(iter, ref, block)
        initParamsDict = self._getParamsIteration(imgSet, iter)
        
        paramsRefine = dict(initParamsDict.items() + paramDic.items() + param.items())
        
        if parRefine == 1:
            paramsRefine['mode'] = 1
            paramsRefine['paramRefine'] = '1, 1, 1, 0, 0'
        elif parRefine == 2:
            paramsRefine['mode'] = 1
            paramsRefine['paramRefine'] = '0, 0, 0, 1, 1'
        elif parRefine == 3:
            paramsRefine['mode'] = 1
            paramsRefine['paramRefine'] = '0, 0, 0, 0, 0'
        
        args = self._prepareCommand()
        
        # frealign program is already in the args script, that's why runJob('')
        self.runJob('', args % paramsRefine)
    
    def reconstructVolumeStep(self, iter, ref, cpuList):
        """Reconstruct a volume from a SetOfParticles with its current parameters refined
        """
        imgSet = self.inputParticles.get()
        initParticle = 1
        finalParticle = imgSet.getSize()
        params = self._getParamsIteration(imgSet, iter)
        
        os.environ['NCPUS'] = str(cpuList[ref-1])
        params['frealign'] = FREALIGNMP_PATH
        params['outputParFn'] = self._getFileName('output_vol_par_class', iter=iter, ref=ref)
        params['initParticle'] = initParticle
        params['finalParticle'] = finalParticle

        params2 = self._setParams3DR(iter, ref)
        
        params3DR = dict(params.items() + params2.items())
        
        args = self._prepareCommand()
        # frealign program is already in the args script, that's why runJob('')
        self.runJob('', args % params3DR)
    
    def calculateOCCStep(self, iter, leaveDir, numberOfBlocks):
        imgSet = self.inputParticles.get()
        numberOfClasses = self.numberOfRef
        cpusRef = self._cpusPerClass(numberOfBlocks, numberOfClasses)
        
        if iter == 1 and not leaveDir:
            ProtFrealignBase._mergeAllParFiles(self, iter, numberOfBlocks)
            parFile = self._getFile('output_par', iter=iter)
            samplingRate = imgSet.getSamplingRate()
            rootFn = self._getFile('output_par_class_tmp', iter=iter)
            args  = self._rsampleCommand()
            program = RSAMPLE_PATH
        else:
            args = self._occCommand()
            tmp = ''
            for ref in range(1, numberOfClasses + 1):
                if not leaveDir:
                    self._mergeAllParFiles(iter, ref, cpusRef[ref-1])
                args += '%s\n' % self._getFile('output_par_class', iter=iter, ref=ref)
                tmp += '%s\n' % self._getFile('output_par_class', iter=iter, ref=ref)
            args = args + tmp + 'eot'
            program = CALC_OCC_PATH
        
        self.runJob(program, args % locals())
        
        if leaveDir:
            self._setLastIter(iter)
            self._leaveDir()
    
    def createOutputStep(self, lastIter):
        from convert import readSetOfClasses3D
        self._createFilenameTemplates(lastIter)
        numberOfClasses = self.numberOfRef
        imgSet = self.inputParticles.get()
        fileparList = []
        volumeList = []
        
        for ref in range(1, numberOfClasses + 1):
            filepar = self._getFileName('output_par_class', iter=lastIter, ref=ref)
            volFn = self._getFileName('iter_vol_class', iter=lastIter, ref=ref)
            fileparList.append(filepar)
            volumeList.append(volFn)
        
        classes = self._createSetOfClasses3D(imgSet)
        readSetOfClasses3D(classes, fileparList, volumeList)
        self._defineOutputs(outputClasses=classes)
        self._defineSourceRelation(imgSet, outputClasses)
        self._defineSourceRelation(self.input3DReference.get(), outputClasses)
    
    #--------------------------- INFO functions ----------------------------------------------------
    def _validate(self):
        errors = ProtFrealignBase._validate(self)
        
        if self.numberOfClasses.get() <= 1:
            errors.append('The number of classes must be at least 2')
        
        return errors
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _setParamsClassRefineParticles(self, iter, ref, block):
        paramDics = {}
        paramDics['stopParam'] = -100
        paramDics['volume'] = self._getFile('ref_vol_class', iter=iter, ref=ref)
        paramDics['outputParFn'] = self._getFile('output_par_block_class', iter=iter, ref=ref, block=block)
        paramDics['inputParFn'] = paramDics['outputParFn']
        paramDics['imgFnMatch'] = self._getFileName('match_block_class', block=block, iter=iter, ref=ref)
        paramDics['outputShiftFn'] = self._getFileName('shift_block_class', block=block, iter=iter, ref=ref)
        paramDics['3Dweigh'] = self._getFileName('weight_block_class', block=block, iter=iter, ref=ref)
        paramDics['FSC3DR1'] = self._getFileName('vol1_block_class', block=block, iter=iter, ref=ref)
        paramDics['FSC3DR2'] = self._getFileName('vol2_block_class', block=block, iter=iter, ref=ref)
        paramDics['VolPhResidual'] = self._getFileName('phase_block_class', block=block, iter=iter, ref=ref)
        paramDics['VolpointSpread'] = self._getFileName('spread_block_class', block=block, iter=iter, ref=ref)
        return paramDics
    
    def _setParams3DR(self, iter, ref):
        """ Setting the parameters to reconstruct a new 3DR"""
        paramDics = {}
        paramDics['mode'] = 0
        paramDics['stopParam'] = 0   #The stopParam must be 0 if you want obtain a 3D reconstruction.
        paramDics['volume'] = self._getFile('iter_vol_class', iter=iter, ref=ref)
        paramDics['inputParFn'] = self._getFile('output_par_class', iter=iter, ref=ref)
        paramDics['imgFnMatch'] = self._getFileName('match_class', iter=iter, ref=ref)
        paramDics['outputShiftFn'] = self._getFileName('shift_class', iter=iter, ref=ref)
        paramDics['3Dweigh'] = self._getFileName('weight_class', iter=iter, ref=ref)
        paramDics['FSC3DR1'] = self._getFileName('vol1_class', iter=iter, ref=ref)
        paramDics['FSC3DR2'] = self._getFileName('vol2_class', iter=iter, ref=ref)
        paramDics['VolPhResidual'] = self._getFileName('phase_class', iter=iter, ref=ref)
        paramDics['VolpointSpread'] = self._getFileName('spread_class', iter=iter, ref=ref)
        return paramDics
    
    def _cpusPerClass(self, numberOfCpus, numberOfClasses):
        """ Return a list with cpusPerClass values, each value will be
        the number of CPUs assigned to each class. The number of CPUs
        will be distributed equally between each class as possible.
        """
        restBlock = numberOfCpus % numberOfClasses
        colBlock = numberOfCpus / numberOfClasses
        blockCpus = [colBlock] * numberOfClasses
        
        for i, v in enumerate(blockCpus):
            if i < restBlock:
                blockCpus[i] += 1
            if blockCpus[i] == 0:
                blockCpus[i] = 1
        return blockCpus
    
    def _mergeAllParFiles(self, iter, ref, numberOfBlocks):
        """ This method merge all parameters files that has been created in a refineIterStep """
        
        file2 = self._getFile('output_par_class', iter=iter, ref=ref)
        if numberOfBlocks != 1:
            f2 = open(file2, 'w+')
            
            for block in range(1, numberOfBlocks + 1):
                file1 = self._getFile('output_par_block_class', block=block, iter=iter, ref=ref)
                f1 = open(file1)
                
                for l in f1:
                    if not l.startswith('C'):
                        f2.write(l)
                f1.close()
            f2.close()
        else:
            file1 = self._getFile('output_par_block_class', block=1, iter=iter, ref=ref)
            copyFile(file1, file2)
    
    def _splitParFile(self, iter, ref, numberOfBlocks):
        """ This method split the parameter files that has been previosuly merged """
        
        prevIter = iter -1
        file1 = self._getFileName('input_par_class', iter=prevIter, ref=ref)
        if numberOfBlocks != 1:
            for block in range(1, numberOfBlocks + 1):
                f1 = open(file1)
                file2 = self._getFileName('output_par_block_class', block= block, iter=prevIter, ref=ref)
                f2 = open(file2, 'w+')
                initpart, finalPart = self._particlesInBlock(block, numberOfBlocks)
                
                for l in f1:
                    
                    if l.startswith('C'):
                        f2.write(l)
                    else:
                        split = l.split()
                        numPart = int(''.join(split[:1]))
                        
                        if numPart <= finalPart:
                            f2.write(l)
                        else:
                            break
                f2.close()
                f1.close()
        else:
            file2 = self._getFileName('output_par_block_class', block=1, iter=prevIter, ref=ref)
            copyFile(file1, file2)
    
    def _rsampleCommand(self):
        args = """ << eot
%(parFile)s
%(samplingRate)f
%(numberOfClasses)d
%(rootFn)s
eot
"""
        return args
    
    def _occCommand(self):
        args = """ << eot
%(numberOfClasses)d
1.0
"""
        return args