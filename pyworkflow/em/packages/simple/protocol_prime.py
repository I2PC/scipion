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
from pyworkflow.em.packages.xmipp3.convert import writeSetOfVolumes, volumeToRow
from pyworkflow.em.packages.xmipp3.xmipp3 import XmippMdRow
from protlib_filesystem import deleteFile
"""
This sub-package contains wrapper around reconstruct_significant Xmipp program
"""

from pyworkflow.em import *  



class ProtPrime(ProtInitialVolume):
    """ Produces one or several initial volumes using reconstruct_significant """
    _label = 'prime'

    #--------------------------- DEFINE param functions --------------------------------------------
    
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputClasses', PointerParam, label="Input classes", important=True, 
                      pointerClass='SetOfClasses2D, SetOfAverages',
                      help='Select the input classes2D from the project.\n'
                           'It should be a SetOfClasses2D class with class representative')
        form.addParam('symmetryGroup', TextParam, default='c1',
                      label="Symmetry group",
                      help='cn or dn. For icosahedral viruses, use c5. If no symmetry is present, give c1.')  
        form.addParam('Nvolumes', IntParam, label='Number of volumes', help="Number of volumes to reconstruct",
                      default=1)
        form.addParam('maximumShift', IntParam, default=0, expertLevel=LEVEL_ADVANCED,
                      label='Maximum shift (px):', help="Set to 0 for free shift search")
        form.addParam('shiftStep', IntParam, default=1, expertLevel=LEVEL_ADVANCED,
                      label='Shift step (px):', help="Step for exhaustive shift search")
        form.addParam('outerMask', FloatParam, default=-1, expertLevel=LEVEL_ADVANCED,
                      label='Outer mask radius (px):', help="Set to -1 for half the image size")
        form.addParam('dynamicFilter', BooleanParam, default=False, expertLevel=LEVEL_ADVANCED,
                      label='Dynamic filtering:', help="Let the program estimate the maximum resolution")
        form.addParam('maxResolution', FloatParam, default=20, expertLevel=LEVEL_ADVANCED, condition="not dynamicFilter",
                      label='Max. Resolution (A):', help="The reconstructed volume will be limited to this resolution in Angstroms.")
        form.addParam('fractionParticles', FloatParam, default=1, expertLevel=LEVEL_ADVANCED,
                      label='Fraction of particles:', help="Fraction of particles to include in the refinement. 1=all particles, 0.8=80% of particles")
        form.addParam('molecularWeight', FloatParam, default=-1, expertLevel=LEVEL_ADVANCED,
                      label='Molecular weight (kD):', help="Molecular weight in kilodaltons, set to -1 for no constraint")
        form.addParam('keepIntermediate', BooleanParam, default=False, expertLevel=LEVEL_ADVANCED,
                      label='Keep intermediate volumes',
                      help='Keep all volumes along iterations')

        form.addParallelSection(threads=4, mpi=0)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    
    def _insertAllSteps(self):
        """ Mainly prepare the command line for calling reconstruct_significant program"""
        
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('runPrime')
        self._insertFunctionStep('getLastIteration')
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

        inputClasses=self.inputClasses.get()
        xdim, _, _ = inputClasses.getDimensions()
        args="stk=classes.spi box=%d smpd=%f pgrp=%s"%(xdim,inputClasses.getSamplingRate(),self.symmetryGroup)
        if self.dynamicFilter:
            args+=" dynlp=yes"
        else:
            args+=" lp=%f"%self.maxResolution.get()
        args+=" trs=%d trsstep=%d"%(self.maximumShift.get(),self.shiftStep.get())
        args+=" nstates=%d"%self.Nvolumes.get()
        args+=" nthr=%d"%self.numberOfThreads.get()
        if self.outerMask>0:
            args+=" ring2=%f"%self.outerMask.get()
        args+=" frac=%f"%self.fractionParticles.get()
        if self.molecularWeight>0:
            args+=" mw=%f"%self.molecularWeight.get()
        
        
        self._enterDir(self._getExtraPath(''))
        self.runJob("simple_prime", args)
        self._leaveDir()

    def getLastIteration(self):
        self._enterDir(self._getExtraPath(''))
        self.lastIter=1
        while True:
            if os.path.exists("recvol_state1_iter%d.spi"%self.lastIter):
                self.lastIter+=1
            else:
                self.lastIter-=1
                break
        self._leaveDir()

    def cleanPrime(self):
        self._enterDir(self._getExtraPath(''))
        cleanPath("cmdline.txt")
        for fn in glob("*.txt"):
            cleanPath(fn)
        for fn in glob("startvol_state*.spi"):
            cleanPath(fn)

        # Get last iteration
        for i in range(1,self.lastIter):
            for fn in glob("recvol_state*_iter%d.spi"%i):
                cleanPath(fn)
        self._leaveDir()
    
    def createOutputStep(self):
        if self.lastIter<=1:
            return
        
        if self.Nvolumes==1:
            vol = Volume()
            print "Adding ",self._getExtraPath('recvol_state1_iter%d.spi'%self.lastIter)
            vol.setLocation(self._getExtraPath('recvol_state1_iter%d.spi'%self.lastIter))
            vol.setSamplingRate(self.inputClasses.get().getSamplingRate())
        else:
            vol = self._createSetOfVolumes()
            vol.setSamplingRate(self.inputClasses.get().getSamplingRate())
            fnVolumes=glob(self._getExtraPath('recvol_state*_iter%d.spi')%self.lastIter)
            print "Adding ",self._getExtraPath('recvol_state*_iter%d.spi'%self.lastIter)
            fnVolumes.sort()
            for fnVolume in fnVolumes:
                aux=Volume()
                aux.setLocation(fnVolume)
                vol.append(aux)

        self._defineOutputs(outputVol=vol)
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
        retval="We used prime [Elmlund2013] to produce an initial volume from the set of classes %s."%\
           self.inputClasses.get().getNameId()
        return [retval]
