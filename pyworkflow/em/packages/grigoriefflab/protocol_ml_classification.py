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
from pyworkflow.utils import copyFile
from pyworkflow.em import ProtClassify3D, ImageHandler
from data import Volume
from pyworkflow.em.packages.grigoriefflab.grigoriefflab import FREALIGNMP_PATH, RSAMPLE_PATH, CALC_OCC_PATH
# from constants import *
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
    def _insertContinueStep(self):
        if self.doContinue:
            continueRun = self.continueRun.get()
            self.inputParticles.set(continueRun.inputParticles.get())
            self.symmetry.set(continueRun.symmetry.get())
            self.input3DReference.set(None)
            self.numberOfRef = continueRun.numberOfClasses.get()
            if self.continueIter.get() == 'last':
                self.initIter = continueRun._getCurrIter()
            else:
                self.initIter = int(self.continueIter.get()) + 1
            self._insertFunctionStep('continueStep', self.initIter)
        else:
            self.initIter = 1
            self.numberOfRef = self.numberOfClasses.get()
        self.finalIter = self.initIter + self.numberOfIterations.get()
        self.cpuList = self._cpusPerClass(self.numberOfBlocks, self.numberOfRef)
    
    def _insertItersSteps(self):
        """ Insert the steps for all iters """
        
        for iterN in self._allItersN():
            depsRecons = []
            initId = self._insertFunctionStep('initIterStep', iterN)
            paramsDic = self._getParamsIteration(iterN)
            depsRefine = self._insertRefineIterStep(iterN, paramsDic, [initId])
            firstOccId = self._insertFunctionStep("calculateOCCStep", iterN, False, prerequisites=depsRefine)
            for ref in self._allRefs():
                reconsId = self._insertFunctionStep("reconstructVolumeStep", iterN, ref, paramsDic, prerequisites=[firstOccId])
                depsRecons.append(reconsId)
            self._insertFunctionStep("calculateOCCStep", iterN, True, prerequisites=depsRecons)
#             depsOcc = [secondOccId]
    
    def _insertRefineIterStep(self, iterN, paramsDic, depsInitId):
        """ execute the refinement for the current iteration """
        
        depsRefine = []
        
        if iterN == 1:
            if not self.useInitialAngles.get():
                stepConstructId = self._insertFunctionStep("constructParamFilesStep", paramsDic, prerequisites=depsInitId)
                depsConstruct = [stepConstructId]
                for block in self._allBlocks():
                    refineId = self._insertFunctionStep("refineBlockStep", block, prerequisites=depsConstruct)
                    depsRefine.append(refineId)
            else:
                initAngStepId = self._insertFunctionStep("writeInitialAnglesStep", prerequisites=depsInitId)
                paramsDic['paramRefine'] = '0, 0, 0, 1, 1'
                for block in self._allBlocks():
                    refineId = self._insertFunctionStep("refineParticlesStep", iterN, block, paramsDic, prerequisites=[initAngStepId])
                    depsRefine.append(refineId)
        else:
            
            for ref in self._allRefs():
                for block in self._selBlocks(self.cpuList[ref-1]):
                    restIterAngle = iterN % self.itRefineAngles.get()
                    restIterShifts = iterN % self.itRefineShifts.get()
 
                    if restIterAngle == 0 and restIterShifts == 0:
                        pass # no change anything
                    elif restIterAngle == 0:
                        paramsDic['mode'] = 1
                        paramsDic['paramRefine'] = '1, 1, 1, 0, 0'
                    elif restIterShifts == 0:
                        paramsDic['mode'] = 1
                        paramsDic['paramRefine'] = '0, 0, 0, 1, 1'
                    else:
                        paramsDic['mode'] = 1
                        paramsDic['paramRefine'] = '0, 0, 0, 0, 0'
                    
                    refineId = self._insertFunctionStep("refineClassParticlesStep", iterN, ref, block, paramsDic, prerequisites=depsInitId)
                    depsRefine.append(refineId)
        return depsRefine
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def initIterStep(self, iterN):
        """ Prepare files and directories for the current iteration """
        
        self._createIterWorkingDir(iterN) # create the working directory for the current iteration.
        prevIter = iterN - 1
        
        if iterN==1:
            imgSet = self.inputParticles.get()
            vol = self.input3DReference.get()
            
            imgFn = self._getFileName('particles')
            volFn = self._getFileName('init_vol')
            refVol = self._getFileName('ref_vol', iter=iterN) # reference volume of the step.
            imgSet.writeStack(imgFn) # convert the SetOfParticles into a mrc stack.
            ImageHandler().convert(vol.getLocation(), volFn) # convert the reference volume into a mrc volume
            copyFile(volFn, refVol)  #Copy the initial volume in the current directory.
            
        for ref in self._allRefs():
            refVol = self._getFileName('ref_vol_class', iter=iterN, ref=ref) # reference volume of the step.
            iterVol =  self._getFileName('iter_vol_class', iter=iterN, ref=ref) # refined volumes of the step
            if iterN == 1:
                copyFile(volFn, iterVol)  #Copy the initial volume in current directory.
            else:
                self._splitParFile(iterN, ref, self.cpuList[ref-1])
                prevIterVol = self._getFileName('iter_vol_class', iter=prevIter, ref=ref) # volumes of the previous iteration
                copyFile(prevIterVol, refVol)   #Copy the reference volume as refined volume.
                copyFile(refVol, iterVol)   #Copy the reference volume as refined volume.
    
    def refineClassParticlesStep(self, iterN, ref, block, paramsDic):
        """Only refine the parameters of the SetOfParticles
        """
        iterDir = self._iterWorkingDir(iterN)
        
        iniPart, lastPart = self._particlesInBlock(block, self.cpuList[ref-1])
        prevIter = iterN - 1
        
        inOutParam = {'inputParFn' : self._getBaseName('input_par_block_class',prevIter=prevIter, iter=iterN, ref=ref, block=block),
                      'initParticle' : iniPart,
                      'finalParticle' : lastPart
                      }
        
        paramClassRefDic = self._setParamsClassRefineParticles(iterN, ref, block)
        
        paramsRefine = dict(paramsDic.items() + paramClassRefDic.items() + inOutParam.items())
        
        args = self._prepareCommand()
        
        # frealign program is already in the args script, that's why runJob('')
        self.runJob('', args % paramsRefine, cwd=iterDir)
    
    def reconstructVolumeStep(self, iterN, ref, paramsDic):
        """Reconstruct a volume from a SetOfParticles with its current parameters refined
        """
        imgSet = self.inputParticles.get()
        initParticle = 1
        finalParticle = imgSet.getSize()
        iterDir = self._iterWorkingDir(iterN)
        
        os.environ['NCPUS'] = str(self.cpuList[ref-1])
        paramsDic['frealign'] = FREALIGNMP_PATH
        paramsDic['outputParFn'] = self._getFileName('output_vol_par_class', iter=iterN, ref=ref)
        paramsDic['initParticle'] = initParticle
        paramsDic['finalParticle'] = finalParticle
#         paramsDic['paramRefine'] = '0, 0, 0, 0, 0'

        params2 = self._setParams3DR(iterN, ref)
        
        params3DR = dict(paramsDic.items() + params2.items())
        
        args = self._prepareCommand()
        # frealign program is already in the args script, that's why runJob('')
        self.runJob('', args % params3DR, cwd=iterDir)
    
    def calculateOCCStep(self, iterN, isLastIterStep):
        imgSet = self.inputParticles.get()
        numberOfClasses = self.numberOfRef
        cpusRef = self._cpusPerClass(self.numberOfBlocks, numberOfClasses)
        iterDir = self._iterWorkingDir(iterN)
        
        if iterN == 1 and not isLastIterStep:
            ProtFrealignBase._mergeAllParFiles(self, iterN, self.numberOfBlocks)
            parFile = self._getBaseName('output_par', iter=iterN)
            samplingRate = imgSet.getSamplingRate()
            rootFn = self._getBaseName('output_par_class_tmp', iter=iterN)
            args  = self._rsampleCommand()
            program = RSAMPLE_PATH
        else:
            args = self._occCommand()
            tmp = ''
            for ref in self._allRefs():
                if not isLastIterStep:
                    self._mergeAllParFiles(iterN, ref, self.cpuList[ref-1])
                args += '%s\n' % self._getBaseName('output_par_class', iter=iterN, ref=ref)
                tmp += '%s\n' % self._getBaseName('output_par_class', iter=iterN, ref=ref)
            args = args + tmp + 'eot'
            program = CALC_OCC_PATH
        
        self.runJob(program, args % locals(), cwd=iterDir)
        
        if isLastIterStep:
            self._setLastIter(iterN)
    
    def createOutputStep(self):
        from convert import readSetOfClasses3D
        numberOfClasses = self.numberOfRef
        imgSet = self.inputParticles.get()
        volumes = self._createSetOfVolumes()
        volumes.setSamplingRate(imgSet.getSamplingRate())
        
        fileparList = []
        volumeList = []
        
        for ref in range(1, numberOfClasses + 1):
            vol = Volume()
            filepar = self._getFileName('output_par_class', iter=self._getLastIter(), ref=ref)
            volFn = self._getFileName('iter_vol_class', iter=self._getLastIter(), ref=ref)
            vol.setFileName(volFn)
            volumes.append(vol)
            fileparList.append(filepar)
            volumeList.append(volFn)
        
        classes = self._createSetOfClasses3D(imgSet)
        readSetOfClasses3D(classes, fileparList, volumeList)
        
        # Define the outputs and relations
        self._defineOutputs(outputClasses=classes)
        self._defineOutputs(outputVolumes=volumes)
        #TODO: save alignment

        self._defineSourceRelation(imgSet, classes)
        self._defineSourceRelation(self.input3DReference.get(), classes)
        self._defineSourceRelation(imgSet, volumes)
        self._defineSourceRelation(self.input3DReference.get(), volumes)
    
    #--------------------------- INFO functions ----------------------------------------------------
    def _validate(self):
        errors = ProtFrealignBase._validate(self)
        
        if self.numberOfClasses.get() <= 1:
            errors.append('The number of classes must be at least 2')
        
        return errors
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _setParamsClassRefineParticles(self, iterN, ref, block):
        paramDics = {}
        paramDics['stopParam'] = -100
        paramDics['volume'] = self._getBaseName('ref_vol_class', iter=iterN, ref=ref)
        paramDics['outputParFn'] = self._getBaseName('output_par_block_class', iter=iterN, ref=ref, block=block)
        paramDics['inputParFn'] = paramDics['outputParFn']
        paramDics['imgFnMatch'] = self._getFileName('match_block_class', block=block, iter=iterN, ref=ref)
        paramDics['outputShiftFn'] = self._getFileName('shift_block_class', block=block, iter=iterN, ref=ref)
        paramDics['3Dweigh'] = self._getFileName('weight_block_class', block=block, iter=iterN, ref=ref)
        paramDics['FSC3DR1'] = self._getFileName('vol1_block_class', block=block, iter=iterN, ref=ref)
        paramDics['FSC3DR2'] = self._getFileName('vol2_block_class', block=block, iter=iterN, ref=ref)
        paramDics['VolPhResidual'] = self._getFileName('phase_block_class', block=block, iter=iterN, ref=ref)
        paramDics['VolpointSpread'] = self._getFileName('spread_block_class', block=block, iter=iterN, ref=ref)
        paramDics['logFile'] = self._getFileName('logFileRefine', block=block, iter=iterN, ref=ref)
        return paramDics
    
    def _setParams3DR(self, iterN, ref):
        """ Setting the parameters to reconstruct a new 3DR"""
        paramDics = {}
        paramDics['mode'] = 0
        paramDics['stopParam'] = 0   #The stopParam must be 0 if you want obtain a 3D reconstruction.
        paramDics['volume'] = self._getBaseName('iter_vol_class', iter=iterN, ref=ref)
        paramDics['inputParFn'] = self._getBaseName('output_par_class', iter=iterN, ref=ref)
        paramDics['imgFnMatch'] = self._getFileName('match_class', iter=iterN, ref=ref)
        paramDics['outputShiftFn'] = self._getFileName('shift_class', iter=iterN, ref=ref)
        paramDics['3Dweigh'] = self._getFileName('weight_class', iter=iterN, ref=ref)
        paramDics['FSC3DR1'] = self._getFileName('vol1_class', iter=iterN, ref=ref)
        paramDics['FSC3DR2'] = self._getFileName('vol2_class', iter=iterN, ref=ref)
        paramDics['VolPhResidual'] = self._getFileName('phase_class', iter=iterN, ref=ref)
        paramDics['VolpointSpread'] = self._getFileName('spread_class', iter=iterN, ref=ref)
        paramDics['logFile'] = self._getFileName('logFileRecons', iter=iterN, ref=ref)
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
    
    def _mergeAllParFiles(self, iterN, ref, numberOfBlocks):
        """ This method merge all parameters files that has been created in a refineIterStep """
        
        file2 = self._getFileName('output_par_class', iter=iterN, ref=ref)
        if numberOfBlocks != 1:
            f2 = open(file2, 'w+')
            
            for block in range(1, numberOfBlocks + 1):
                file1 = self._getFileName('output_par_block_class', block=block, iter=iterN, ref=ref)
                f1 = open(file1)
                
                for l in f1:
                    if not l.startswith('C'):
                        f2.write(l)
                f1.close()
            f2.close()
        else:
            file1 = self._getFileName('output_par_block_class', block=1, iter=iterN, ref=ref)
            copyFile(file1, file2)
    
    def _splitParFile(self, iterN, ref, numberOfBlocks):
        """ This method split the parameter files that has been previosuly merged """
        
        prevIter = iterN -1
        file1 = self._getFileName('output_par_class', iter=prevIter, ref=ref)
        if numberOfBlocks != 1:
            for block in range(1, numberOfBlocks + 1):
                f1 = open(file1)
                file2 = self._getFileName('input_par_block_class',prevIter=prevIter, iter=iterN, ref=ref, block=block)
                f2 = open(file2, 'w+')
                _, finalPart = self._particlesInBlock(block, numberOfBlocks)
                
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
            file2 = self._getFileName('input_par_block_class',prevIter=prevIter, iter=iterN, ref=ref, block=block)
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
    
    
    def _allRefs(self):
        """ Iterate over all references. """
        for i in range(1, self.numberOfRef+1):
            yield i

    def _selBlocks(self, blocks):
        """ Iterate over all numberOfCPUs. """
        for i in range(1, blocks+1):
            yield i
    
    def _particlesInBlock(self, block, numberOfBlocks):
        """calculate the initial and final particles that belongs to this block"""
        
        imgSet = self.inputParticles.get()
        
        blockParticles = self._particlesPerBlock(numberOfBlocks, imgSet.getSize())
        initPart = 0
        lastPart = 0
        for i in range(block):
            initPart = lastPart + 1
            lastPart = lastPart + blockParticles[i]
        
        return initPart, lastPart

