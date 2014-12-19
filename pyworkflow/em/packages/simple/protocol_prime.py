# **************************************************************************
# *
# * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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
This sub-package contains wrapper around Prime algorithm in Simple
"""

import os
from glob import glob

import pyworkflow.protocol.params as params
from pyworkflow.utils.path import cleanPath, cleanPattern

import pyworkflow.em as em  
import simple



class ProtPrime(em.ProtInitialVolume):
    """ Produces one or several initial volumes using reconstruct_significant """
    _label = 'prime'

    #--------------------------- DEFINE param functions --------------------------------------------
    
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputClasses', params.PointerParam, label="Input classes", important=True, 
                      pointerClass='SetOfClasses2D, SetOfAverages',
                      help='Select the input classes2D from the project.\n'
                           'It should be a SetOfClasses2D class with class representative')
        form.addParam('symmetryGroup', params.TextParam, default='c1',
                      label="Symmetry group",
                      help='cn or dn. For icosahedral viruses, use c5. \n'
                           'If no symmetry is present, give c1.')  
        form.addParam('Nvolumes', params.IntParam, default=1, 
                      label='Number of volumes', 
                      help="Number of volumes to reconstruct")
        form.addParam('maximumShift', params.IntParam, default=0, 
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Maximum shift (px):', 
                      help="Set to 0 for free shift search")
        form.addParam('shiftStep', params.IntParam, default=1, 
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Shift step (px):', 
                      help="Step for exhaustive shift search")
        form.addParam('outerMask', params.FloatParam, default=-1, 
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Outer mask radius (px):', 
                      help="Set to -1 for half the image size")
        form.addParam('dynamicFilter', params.BooleanParam, default=False, 
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Dynamic filtering:', 
                      help="Let the program estimate the maximum resolution")
        form.addParam('maxResolution', params.FloatParam, default=20, 
                      expertLevel=params.LEVEL_ADVANCED, condition="not dynamicFilter",
                      label='Max. Resolution (A):', 
                      help="The reconstructed volume will be limited to \n"
                           "this resolution in Angstroms.")
        form.addParam('fractionParticles', params.FloatParam, default=1, 
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Fraction of particles:', 
                      help="Fraction of particles to include in the refinement. \n"
                           "1=all particles, 0.8=80% of particles")
        form.addParam('molecularWeight', params.FloatParam, default=-1, 
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Molecular weight (kD):', 
                      help="Molecular weight in kilodaltons, \n"
                           "set to -1 for no constraint")
        form.addParam('keepIntermediate', params.BooleanParam, default=False, 
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Keep intermediate volumes',
                      help='Keep all volumes along iterations')

        form.addParallelSection(threads=8, mpi=0)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    
    def _insertAllSteps(self):
        """ Mainly prepare the command line for calling reconstruct_significant program"""
        
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('runPrime')
        if not self.keepIntermediate:
            self._insertFunctionStep('cleanPrime')
        self._insertFunctionStep('createOutputStep')        

    #--------------------------- STEPS functions --------------------------------------------        
    def convertInputStep(self):
        self.inputClasses.get().writeStack(self._getExtraPath("classes.spi:stk"))
            
    def runPrime(self):
        # simple_prime stk=stack.spi [vol1=invol.spi] [vol2=<refvol_2.spi> etc.] box=<image size(in pixels)> 
        #              smpd=<sampling distance(in A)> [ring2=<outer mask radius(in pixels){box/2}>] 
        #              [trs=<origin shift(in pixels){0}>] [trsstep=<origin shift stepsize{1}>] [lp=<low-pass limit{20}>]
        #              [dynlp=<yes|no{no}>] [nstates=nstates to reconstruct>] [frac=<fraction of ptcls to include{1}>]
        #              [mw=<molecular weight (in kD)>] [oritab=<previous rounds alignment doc>] [nthr=<nr of OpenMP threads{1}>]

        inputClasses = self.inputClasses.get()
        xdim, _, _ = inputClasses.getDimensions()
        args = "stk=classes.spi box=%d smpd=%f pgrp=%s" % (xdim, inputClasses.getSamplingRate(), self.symmetryGroup)
        
        if self.dynamicFilter:
            args += " dynlp=yes"
        else:
            args += " lp=%f" % self.maxResolution
        args += " trs=%d trsstep=%d" % (self.maximumShift, self.shiftStep)
        args += " nstates=%d" % self.Nvolumes
        args += " nthr=%d" % self.numberOfThreads
        
        if self.outerMask > 0:
            args += " ring2=%f" % self.outerMask
        args += " frac=%f" % self.fractionParticles
        
        if self.molecularWeight > 0:
            args += " mw=%f" % self.molecularWeight
        
        self.runJob("simple_prime", args, cwd=self._getExtraPath())

    def getLastIteration(self):
        lastIter = 1
        pattern = self._getExtraPath("recvol_state1_iter%d.spi")
        while os.path.exists(pattern % lastIter):
            lastIter += 1
        return lastIter - 1

    def cleanPrime(self):
        self._enterDir(self._getExtraPath())
        cleanPath("cmdline.txt")
        cleanPattern("*.txt")
        cleanPattern("startvol_state*.spi")
        # Get last iteration
        for i in range(1, self.getLastIteration()):
            cleanPattern("recvol_state*_iter%d.spi"%i)
        self._leaveDir()
    
    def createOutputStep(self):
        lastIter = self.getLastIteration()
        
        if lastIter <= 1:
            return
        
        if self.Nvolumes == 1:
            vol = em.Volume()
            vol.setLocation(self._getExtraPath('recvol_state1_iter%d.spi' % lastIter))
            vol.setSamplingRate(self.inputClasses.get().getSamplingRate())
            self._defineOutputs(outputVol=vol)
        else:
            vol = self._createSetOfVolumes()
            vol.setSamplingRate(self.inputClasses.get().getSamplingRate())
            fnVolumes=glob(self._getExtraPath('recvol_state*_iter%d.spi') % lastIter)
            fnVolumes.sort()
            for fnVolume in fnVolumes:
                aux=em.Volume()
                aux.setLocation(fnVolume)
                vol.append(aux)
            self._defineOutputs(outputVolumes=vol)

        self._defineSourceRelation(vol, self.inputClasses)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        summary.append("Input classes: %s" % self.inputClasses.get().getNameId())
        summary.append("Starting from: %d random volumes"%self.Nvolumes.get())
        return summary
    
    def _citations(self):
        return ['Elmlund2013']
    
    def _methods(self):
        if self.inputClasses.get() is not None:
            retval="We used *simple_prime* program [Elmlund2013] to produce an initial volume from the set of classes %s."
            return [retval % self.inputClasses.get().getNameId()]
        else:
            return []
