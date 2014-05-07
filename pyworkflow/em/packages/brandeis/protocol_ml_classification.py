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
# *  e-mail address 'jgomez@cnb.csic.es'
# *
# **************************************************************************
"""
This module contains the protocol to obtain a refined 3D recontruction from a set of particles using Frealign
"""
import os
from pyworkflow.utils import *
from pyworkflow.em import *
from data import *
from brandeis import *
from constants import *
from protocol_frealign_base import ProtFrealignBase


class ProtFrealignClassify(ProtFrealignBase, ProtClassify3D):
    """ This class implements the wrapper to single particle refinement protocol with frealign."""
    _label = 'frealign classify'
    IS_REFINE = False
    
    def __init__(self, **args):
        ProtFrealignBase.__init__(self, **args)
    
    def _createFilenameTemplates(self, iter):
        """ Centralize how files are called for iterations and references. """
        prevIter = iter - 1
        myDict = {
                  'particles': self._getTmpPath('particles.mrc'),
                  'init_vol': self._getTmpPath('volume.mrc'),
                  # Volumes for the iteration
                  'ref_vol': self._iterWorkingDir(iter, 'reference_volume_iter_%(iter)03d.mrc'),
                  'iter_vol': self._iterWorkingDir(iter, 'volume_iter_%(iter)03d.mrc'),
                  'prev_vol': self._iterWorkingDir(prevIter, 'volume_iter_%(iter)03d.mrc'),
                  'output_vol_par': 'output_vol_iter_%(iter)03d.par',
                  # dictionary for all set
                  'input_par': self._iterWorkingDir(prevIter, 'particles_iter_%(iter)03d.par'),
                  'output_par': self._iterWorkingDir(iter, 'particles_iter_%(iter)03d.par'),
                  'shift' : 'particles_shifts_iter_%(iter)03d.shft',
                  'match' : 'particles_match_iter_%(iter)03d.mrc',
                  'weight' : 'volume_weights_iter_%(iter)03d.mrc',
                  'vol1' : 'volume_1_iter_%(iter)03d.mrc',
                  'vol2' : 'volume_2_iter_%(iter)03d.mrc',
                  'phase' : 'volume_phasediffs_iter_%(iter)03d.mrc',
                  'spread' : 'volume_pointspread_iter_%(iter)03d.mrc',
                  # each class volumes for the iteration
                  'ref_vol_class': self._iterWorkingDir(iter, 'reference_volume_iter_%(iter)03d_class_%(ref)02d.mrc'),
                  'iter_vol_class': self._iterWorkingDir(iter, 'volume_iter_%(iter)03d_class_%(ref)02d.mrc'),
                  'prev_vol_class': self._iterWorkingDir(prevIter, 'volume_iter_%(iter)03d_class_%(ref)02d.mrc'),
                  'output_vol_class_par': 'output_vol_iter_%(iter)03d_class_%(ref)02d.par',
                  # dictionary for each class
                  'input_par_class': self._iterWorkingDir(prevIter, 'particles_iter_%(iter)03d_class_%(ref)02d.par'),
                  'output_par_class': self._iterWorkingDir(iter, 'particles_iter_%(iter)03d_class_%(ref)02d.par'),
                  'output_par_class_tmp': self._iterWorkingDir(iter, 'particles_iter_%(iter)03d_class_0.par'),
                  'shift_class' : 'particles_shifts_iter_%(iter)03d_class_%(ref)02d.shft',
                  'match_class' : 'particles_match_iter_%(iter)03d_class_%(ref)02d.mrc',
                  'weight_class' : 'volume_weights_iter_%(iter)03d_class_%(ref)02d.mrc',
                  'vol1_class' : 'volume_1_iter_%(iter)03d_class_%(ref)02d.mrc',
                  'vol2_class' : 'volume_2_iter_%(iter)03d_class_%(ref)02d.mrc',
                  'phase_class' : 'volume_phasediffs_iter_%(iter)03d_class_%(ref)02d.mrc',
                  'spread_class' : 'volume_pointspread_iter_%(iter)03d_class_%(ref)02d.mrc',
                  # dictionary for each processing block and class
                  'input_par_block': self._iterWorkingDir(prevIter, 'particles_iter_%(iter)03d_class_%(ref)02d_%(block)02d.par'),
                  'output_par_block': self._iterWorkingDir(iter, 'particles_iter_%(iter)03d_class_%(ref)02d_%(block)02d.par'),
                  'shift_block' : 'particles_shifts_iter_%(iter)03d_class_%(ref)02d_%(block)02d.shft',
                  'match_block' : 'particles_match_iter_%(iter)03d_class_%(ref)02d_%(block)02d.mrc', 
                  'weight_block' : 'volume_weights_iter_%(iter)03d_class_%(ref)02d_%(block)02d',
                  'vol1_block' : 'volume_1_%(iter)03d_class_%(ref)02d_%(block)02d_iter',
                  'vol2_block' : 'volume_2_iter_%(iter)03d_class_%(ref)02d_%(block)02d',
                  'phase_block' : 'volume_phasediffs_iter_%(iter)03d_class_%(ref)02d_%(block)02d',
                  'spread_block' : 'volume_pointspread_iter_%(iter)03d_class_%(ref)02d_%(block)02d',
                  'output_vol_par': 'output_vol_iter_%(iter)03d.par',
                  }
        
        self._fnDict = myDict
    
  #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        """Insert the steps to refine orientations and shifts of the SetOfParticles
        """
        numberOfBlocks = self.numberOfThreads.get() - 1
        self.numberOfRef = self.numberOfClasses.get()
        cpuList = self._cpusPerClass(numberOfBlocks, self.numberOfRef)
        depsOcc = []
        
        for iter in range(1, self.numberOfIterations.get() + 1):
            depsRecons = []
            initId = self._insertFunctionStep('initIterStep', iter, numberOfBlocks, prerequisites=depsOcc)
            depsRefine = self._insertRefineIterStep(iter, numberOfBlocks, [initId])
            firstOccId = self._insertFunctionStep("calculateOCCStep", iter, False, numberOfBlocks, prerequisites=depsRefine)
            for ref in range(1, self.numberOfRef + 1):
                reconsId = self._insertFunctionStep("reconstructVolumeStep", iter, ref, cpuList, prerequisites=[firstOccId])
                depsRecons = depsRecons + [reconsId]
            secondOccId = self._insertFunctionStep("calculateOCCStep", iter, True, numberOfBlocks, prerequisites=depsRecons)
            depsOcc = [secondOccId]
        self._insertFunctionStep("createOutputStep", prerequisites=depsOcc)
    
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
                    refineId = self._insertFunctionStep("refineParticlesStep", iter, ref, block, parRefine, prerequisites=depsInitId)
                    depsRefine.append(refineId)
        return depsRefine
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def initIterStep(self, iter, numberOfBlocks):
        """ Prepare files and directories for the current iteration """
        
        self._createIterWorkingDir(iter) # create the working directory for the current iteration.
        prevIter = iter - 1
        self._createFilenameTemplates(iter)
        
        cpusRefs = self._cpusPerClass(numberOfBlocks, self.numberOfRef)
            
            
        if iter==1:
            imgSet = self.inputParticles.get()
            vol = self.input3DReference.get()
            
            imgFn = self._getFileName('particles')
            volFn = self._getFileName('init_vol')
            imgSet.writeStack(imgFn) # convert the SetOfParticles into a mrc stack.
            ImageHandler().convert(vol.getLocation(), volFn) # convert the reference volume into a mrc volume
            
        for ref in range(1, self.numberOfRef + 1):
            refVol = self._getFileName('ref_vol_class', iter=iter, ref=ref) # reference volume of the step.
            iterVol =  self._getFileName('iter_vol_class', iter=iter, ref=ref) # refined volumes of the step
            if iter==1:
                copyFile(volFn, refVol)  #Copy the initial volume in the current directory.
            else:
                self._splitParFile(iter, ref, cpusRef[ref-1])
                prevIterVol = self._getFileName('prev_vol_class', iter=prevIter, ref=ref) # volumes of the previous iteration
                copyFile(prevIterVol, refVol)   #Copy the reference volume as refined volume.
            copyFile(refVol, iterVol)   #Copy the reference volume as refined volume.
    
    def refineParticlesStep(self, iter, ref, block, parRefine):
        """Only refine the parameters of the SetOfParticles
        """
        param = {}
        imgSet = self.inputParticles.get()
        
        iterDir = self._iterWorkingDir(iter)
        if block==1:
            self._enterDir(iterDir) # enter to the working directory for the current iteration.
        
        iniPart, lastPart = self._particlesInBlock(block)
        prevIter = iter - 1
        param['inputParFn'] = self._getFile('input_par_block', iter=prevIter, ref=ref, block=block)
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
        params['outputParFn'] = self._getFileName('output_vol_par', iter=iter)
        params['initParticle'] = initParticle
        params['finalParticle'] = finalParticle

        params2 = self._setParams3DR(iter, ref)
        
        params3DR = dict(params.items() + params2.items())
        
        args = self._prepareCommand()
        # frealign program is already in the args script, that's why runJob('')
        self.runJob('', args % params3DR)
    
    def calculateOCCStep(self, iter, leaveDir, numberOfBlocks):
        imgSet = self.inputParticles.get()
        
        if iter == 1 and not leaveDir:
            ProtFrealignBase._createFilenameTemplates(self, iter)
            ProtFrealignBase._mergeAllParFiles(self, iter, numberOfBlocks)
            parFile = self._getFile('output_par', iter=iter)
            samplingRate = imgSet.getSamplingRate()
            numberOfClasses = self.numberOfRef
            
            self._createFilenameTemplates(iter)
            rootFn = self._getFile('output_par_class_tmp', iter=iter)
            
            args  = self._rsampleCommand()
            program = RSAMPLE_PATH
        else:
            args = self._occCommand()
            numberOfRef = self.numberOfRef
            for ref in range(1, self.numberOfRef + 1):
                if not leaveDir:
                    self._mergeAllParFiles(iter, ref, numberOfBlocks)
                args += '%s\n' % self._getFile('output_par_class', iter=iter, ref=ref)
                tmp += '%s\n' % self._getFile('output_par_class', iter=iter, ref=ref)
            args = args + tmp + 'eot'
            program = CALC_OCC_PATH
        
        self.runJob(program, args % locals())
        
        if leaveDir:
            self._leaveDir()
    
    def createOutputStep(self):
        lastIter = self.numberOfIterations.get()
        lastIterDir = self._iterWorkingDir(lastIter)
        volFn = join(lastIterDir, 'volume_iter_%03d.mrc' % lastIter)
        vol = Volume()
        vol.setSamplingRate(self.inputParticles.get().getSamplingRate())
        vol.setFileName(volFn)
        self._defineOutputs(outputVolume=vol)
    
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
        paramDics['outputParFn'] = self._getFile('output_par_block', iter=iter, ref=ref, block=block)
        paramDics['inputParFn'] = paramDics['outputParFn']
        paramDics['imgFnMatch'] = self._getFileName('match_block', block=block, iter=iter, ref=ref)
        paramDics['outputShiftFn'] = self._getFileName('shift_block', block=block, iter=iter, ref=ref)
        paramDics['3Dweigh'] = self._getFileName('weight_block', block=block, iter=iter, ref=ref)
        paramDics['FSC3DR1'] = self._getFileName('vol1_block', block=block, iter=iter, ref=ref)
        paramDics['FSC3DR2'] = self._getFileName('vol2_block', block=block, iter=iter, ref=ref)
        paramDics['VolPhResidual'] = self._getFileName('phase_block', block=block, iter=iter, ref=ref)
        paramDics['VolpointSpread'] = self._getFileName('spread_block', block=block, iter=iter, ref=ref)
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
                file1 = self._getFile('output_par_block', block=block, iter=iter, ref=ref)
                f1 = open(file1)
                
                for l in f1:
                    if not l.startswith('C'):
                        f2.write(l)
                f1.close()
            f2.close()
        else:
            file1 = self._getFile('output_par_block', block=1, iter=iter, ref=ref)
            copyFile(file1, file2)
    
    def _splitParFile(self, iter, ref, numberOfBlocks):
        """ This method split the parameter files that has been previosuly merged """
        
        prevIter = iter -1
        file1 = self._getFileName('input_par_class', iter=prevIter, ref=ref)
        if numberOfBlocks != 1:
            for block in range(1, numberOfBlocks + 1):
                f1 = open(file1)
                file2 = self._getFileName('output_par_block', block= block, iter=prevIter, ref=ref)
                f2 = open(file2, 'w+')
                initpart, finalPart = self._particlesInBlock(block)
                
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
            file2 = self._getFileName('output_par_block', block=1, iter=prevIter, ref=ref)
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
%(numberOfRef)d
1.0
"""