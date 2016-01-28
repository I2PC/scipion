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
This sub-package contains wrapper around EMAN initialmodel program
"""

import os
import re
from os.path import exists
from glob import glob
import pyworkflow.em as em
from pyworkflow.em.packages.eman2.eman2 import getEmanProgram, validateVersion
from pyworkflow.em.packages.eman2.convert import createEmanProcess
from pyworkflow.protocol.params import (PointerParam, FloatParam, IntParam, EnumParam,
                                        StringParam, BooleanParam)
from pyworkflow.utils.path import cleanPattern, makePath, createLink
from pyworkflow.em.data import Volume
from pyworkflow.em.protocol import ProtRefine3D

from convert import rowToAlignment

# speed
SPEED_1 = 0
SPEED_2 = 1
SPEED_3 = 2
SPEED_4 = 3
SPEED_5 = 4
SPEED_6 = 5
SPEED_7 = 6

# modes to postprocess
FILTER_NONE = 0
FILTER_HP_1 = 1
FILTER_HP_2 = 2
FILTER_HP_3 = 3
FILTER_HP_4 = 4
FILTER_HP_5 = 5
FILTER_LP_1 = 6
FILTER_LP_2 = 7
FILTER_LP_3 = 8
FILTER_LP_4 = 9
FILTER_LP_5 = 10
FILTER_LP_6 = 11

                               
class EmanProtRefine(ProtRefine3D):
    """
    This Protocol wraps *e2refine_easy.py* Eman2 program.
This is the primary single particle refinement program in EMAN2.1+. Major
features of this program:
- While a range of command-line options still exist. You should not normally specify
more than a few basic requirements. The rest will be auto-selected for you.
- This program will split your data in half and automatically refine the halves
independently to produce a gold standard resolution curve for every step in the refinement.
- The gold standard FSC also permits us to automatically filter the structure
at each refinement step. The resolution you specify is a target, not the filter resolution.
    """
    _label = 'refine easy'

    def _createFilenameTemplates(self):
        """ Centralize the names of the files. """
        
        myDict = {
                  'partSet': 'sets/inputSet.lst',
                  'partFlipSet': 'sets/inputSet__ctf_flip.lst',
                  'data_scipion': self._getExtraPath('data_scipion_it%(iter)02d.sqlite'),
                  'projections': self._getExtraPath('projections_it%(iter)02d_%(half)s.sqlite'),
                  'volume': self._getExtraPath('refine_%(run)02d/threed_%(iter)02d.hdf'),
                  'classes': 'refine_%(run)02d/classes_%(iter)02d',
                  'classesEven': self._getExtraPath('refine_%(run)02d/classes_%(iter)02d_even.hdf'),
                  'classesOdd': self._getExtraPath('refine_%(run)02d/classes_%(iter)02d_odd.hdf'),
                  'cls': 'refine_%(run)02d/cls_result_%(iter)02d',
                  'clsEven': self._getExtraPath('refine_%(run)02d/cls_result_%(iter)02d_even.hdf'),
                  'clsOdd': self._getExtraPath('refine_%(run)02d/cls_result_%(iter)02d_odd.hdf'),
                  'angles': self._getExtraPath('projectionAngles_it%(iter)02d.txt'),
                  'mapEven': self._getExtraPath('refine_%(run)02d/threed_%(iter)02d_even.hdf'),
                  'mapOdd': self._getExtraPath('refine_%(run)02d/threed_%(iter)02d_odd.hdf'),
                  'mapFull': self._getExtraPath('refine_%(run)02d/threed_%(iter)02d.hdf'),
                  'fscUnmasked': self._getExtraPath('refine_%(run)02d/fsc_unmasked_%(iter)02d.txt'),
                  'fscMasked': self._getExtraPath('refine_%(run)02d/fsc_masked_%(iter)02d.txt'),
                  'fscMaskedTight': self._getExtraPath('refine_%(run)02d/fsc_maskedtight_%(iter)02d.txt'),
                  }
        
        self._updateFilenamesDict(myDict)
    
    def _createIterTemplates(self, currRun):
        """ Setup the regex on how to find iterations. """
        self._iterTemplate = self._getFileName('mapEven', run=currRun, iter=0).replace('00','??')
        # Iterations will be identify by threed_XX_ where XX is the iteration number
        # and is restricted to only 2 digits.
        self._iterRegex = re.compile('threed_(\d{2,2})_')
    
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('doContinue', BooleanParam, default=False,
              label='Continue from a previous run?',
              help='If you set to *Yes*, you should select a previous'
              'run of type *%s* class and most of the input parameters'
              'will be taken from it.' % self.getClassName())
        form.addParam('continueRun', PointerParam, pointerClass=self.getClassName(),
                      condition='doContinue', allowsNull=True,
                      label='Select previous run',
                      help='Select a previous run to continue from.')
        form.addParam('inputParticles', PointerParam, label="Input particles", important=True,
                      pointerClass='SetOfParticles', condition='not doContinue', allowsNull=True,
                      help='Select the input particles.\n')  
        form.addParam('input3DReference', PointerParam,
                      pointerClass='Volume', allowsNull=True,
                      label='Initial 3D reference volume:',
                      condition='not doContinue',
                      help='Input 3D reference reconstruction.\n')
        form.addParam('numberOfIterations', IntParam, default=2,
                      label='Number of iterations:',
                      help='Set the number of iterations. Iterative reconstruction'
                           'improves the overall normalization of the 2D images'
                           'as they are inserted into the reconstructed volume,'
                           'and allows for the exclusion of the poorer quality'
                           'images.')
        form.addParam('symmetry', StringParam, default='c1',
                      condition='not doContinue',
                      label='Symmetry group',
                      help='Set the symmetry; if no value is given then the model'
                           'is assumed to have no symmetry. Choices are: i(n),'
                           ' c(n), d(n), tet, icos, or oct.'
                           'See http://blake.bcm.edu/emanwiki/EMAN2/Symmetry'
                           'for a detailed descript of symmetry in Eman.')
        form.addParam('resol', FloatParam, default='10.0',
                      label='Resolution of this refinement run (A):',
                      help='Target resolution in A of this refinement run.'
                           'Usually works best in at least two steps'
                           '(low/medium) resolution, then final resolution)'
                           'when starting with a poor starting model.'
                           'Usually 3-4 iterations is sufficient.')
        form.addParam('molMass', FloatParam, default='500.0', 
                      label='Molecular mass of the specimen (kDa):',
                      help='Approximate molecular mass of the particle, in kDa.'
                           'This is used to runnormalize.bymass. Due to'
                           'resolution effects, not always the true mass.')
#         form.addParam('haveDataBeenPhaseFlipped', BooleanParam, default=False,
#                       label='Have data been phase-flipped?', condition='not doContinue',
#                       help='Set this to Yes if the images have been ctf-phase corrected during the '
#                            'pre-processing steps.')       
        form.addParam('doBreaksym', BooleanParam, default=False,
                       label='Do not impose symmetry?',
                       help='If set True, reconstruction will be asymmetric with'
                            '*Symmetry group* parameter specifying a known pseudosymmetry,'
                            'not an imposed symmetry.')
        form.addParam('useE2make3d', BooleanParam, default=False,
                       label='use e2make3d?',
                       help='Use the traditional e2make3d program instead of the'
                            'new e2make3dpar program.')
        form.addParam('speed', EnumParam,
                      choices=['1', '2', '3', '4', '5', '6', '7',],
                      label="Balance speed vs precision:", default=SPEED_5,
                      display=EnumParam.DISPLAY_COMBO,
                      help='Larger values sacrifice a bit of potential resolution for'
                           'significant speed increases. Set to 1 when really'
                           'pushing resolution. Set to 7 for initial refinements.')
        form.addParam('classKeep', FloatParam, default='0.9',
                      label='Fraction of particles to use in final average:',
                      help='The fraction of particles to keep in each class,'
                           'based on the similarity score.')
        form.addParam('m3dKeep', FloatParam, default='0.8',
                      label='Fraction of class-averages to use in 3-D map:',
                      help='The fraction of slices to keep in reconstruction.')
        form.addParam('useSetsfref', BooleanParam, default=True,
                       label='Use the setsfref option in class averaging?',
                       help='This matches the filtration of the class-averages to the'
                            'projections for easier comparison. May also improve'
                            'convergence.')
        form.addParam('doAutomask', BooleanParam, default=False,
                       label='Do automask to the class-average?',
                       help='This will apply an automask to the class-average'
                            'during iterative alignment for better accuracy. The'
                            'final class averages are unmasked.')
        form.addParam('doThreshold', BooleanParam, default=False,
                       label='Apply threshold before project the volume?',
                       help='Applies a threshold to the volume just before'
                            'generating projections. A sort of aggressive solvent'
                            'flattening for the reference.')
        form.addParam('m3dPostProcess', EnumParam,
                      choices=['None', 'filter.highpass.autopeak',
                               'filter.highpass.butterworth', 'filter.highpass.gauss',
                               'filter.highpass.tanh', 'filter.highpassl.tophat',
                               'filter.lowpass.autob', 'filter.lowpass.butterworth',
                               'filter.lowpass.gauss', 'filter.lowpass.randomphase',
                               'filter.lowpass.tanh', 'filter.lowpass.tophat'],
                      label="Mode to Fourier method:", default=FILTER_NONE,
                      display=EnumParam.DISPLAY_COMBO)

        form.addParallelSection(threads=4, mpi=0)
    
    #--------------------------- INSERT steps functions --------------------------------------------  
    def _insertAllSteps(self):        
        self._createFilenameTemplates()
        self._createIterTemplates(self._getRun())
        if self.doContinue:
            self.input3DReference.set(None)
            self.inputParticles.set(None)
            self._insertFunctionStep('createLinkSteps')
            args = self._prepareContinueParams()
        else:
            self._insertFunctionStep('convertImagesStep')
            args = self._prepareParams()
        self._insertFunctionStep('refineStep', args)
        self._insertFunctionStep('createOutputStep')
    
    #--------------------------- STEPS functions --------------------------------------------
    def createLinkSteps(self):
        continueRun = self.continueRun.get()
        prevPartDir = continueRun._getExtraPath("particles")
        currPartDir = self._getExtraPath("particles")
        runN = self._getRun() - 1
        prevRefDir = continueRun._getExtraPath("refine_%02d" % runN)
        currRefDir = self._getExtraPath("refine_%02d" % runN)
        prevSetsDir = continueRun._getExtraPath("sets")
        currSetsDir = self._getExtraPath("sets")

#         createLink(prevInfoDir, currInfoDir)
        createLink(prevPartDir, currPartDir)
        createLink(prevRefDir, currRefDir)
        createLink(prevSetsDir, currSetsDir)
    
    def convertImagesStep(self):
        from pyworkflow.em.packages.eman2.convert import writeSetOfParticles
        partSet = self._getInputParticles()
        partAlign = partSet.getAlignment()
        storePath = self._getExtraPath("particles")
        makePath(storePath)
        writeSetOfParticles(partSet, storePath, alignType=partAlign)
        if partSet.hasCTF():
            program = getEmanProgram('e2ctf.py')
            acq = partSet.getAcquisition()
            
            args = " --voltage %3d" % acq.getVoltage()
            args += " --cs %f" % acq.getSphericalAberration()
            args += " --ac %f" % (100 * acq.getAmplitudeContrast())
            if not partSet.isPhaseFlipped():
                args += " --phaseflip"
            args += " --computesf --apix %f --allparticles --autofit --curdefocusfix --storeparm -v 8" % (partSet.getSamplingRate())
            self.runJob(program, args, cwd=self._getExtraPath())
        
        program = getEmanProgram('e2buildsets.py')
        args = " --setname=inputSet --allparticles --minhisnr=-1"
        self.runJob(program, args, cwd=self._getExtraPath())
    
    def refineStep(self, args):
        """ Run the EMAN program to refine a volume. """
        if not self.doContinue:
            cleanPattern(self._getExtraPath('refine_01'))
        program = getEmanProgram('e2refine_easy.py')
        self.runJob(program, args, cwd=self._getExtraPath())
    
    def createOutputStep(self):
        iterN = self.numberOfIterations.get()
        partSet = self._getInputParticles()
        numRun = self._getRun()
        
        vol = Volume()
        vol.setFileName(self._getFileName("volume",run=numRun, iter=iterN))
        vol.copyInfo(partSet)
        
        self._execEmanProcess(numRun, iterN)
        newPartSet = self._createSetOfParticles()
        newPartSet.copyInfo(partSet)
        newPartSet.setAlignment(em.ALIGN_PROJ)
        newPartSet.copyItems(partSet,
                             updateItemCallback=self._createItemMatrix,
                             itemDataIterator=self._iterTextFile(iterN))
        
        self._defineOutputs(outputVolume=vol)
        self._defineSourceRelation(self._getInputParticlesPointer(), vol)
        self._defineOutputs(outputParticles=newPartSet)
        self._defineTransformRelation(self._getInputParticlesPointer(), newPartSet)
    
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        errors = []
        validateVersion(self, errors)

        particles = self._getInputParticles()
        samplingRate = particles.getSamplingRate()

        if self.resol <  2 * samplingRate:
            errors.append("\nTarget resolution is smaller than 2*samplingRate value. This is impossible.")
        
        if not self.doContinue:
            self._validateDim(particles, self.input3DReference.get(), errors,
                              'Input particles', 'Reference volume')

        return errors
    
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputVolumes'):
            summary.append("Output volumes not ready yet.")
        return summary
    
    #--------------------------- UTILS functions --------------------------------------------
    def _prepareParams(self):
        args1 = "--input=%(imgsFn)s --model=%(volume)s"
        args2 = self._commonParams()
        
        volume = os.path.relpath(self.input3DReference.get().getFileName(), self._getExtraPath()).replace(":mrc","")
        params = {'imgsFn': self._getParticlesStack(),
                  'volume': volume,
                  }
        
        args = args1 % params + args2
        return args
    
    def _prepareContinueParams(self):
        args1 = "--startfrom=refine_%02d" % (self._getRun() - 1)
        args2 = self._commonParams()
        args = args1 + args2
        return args
    
    def _commonParams(self):
        args = " --targetres=%(resol)f --speed=%(speed)d --sym=%(sym)s --iter=%(numberOfIterations)d"
        args += " --mass=%(molMass)f --apix=%(samplingRate)f --classkeep=%(classKeep)f"
        args += " --m3dkeep=%(m3dKeep)f --parallel=thread:%(threads)d --threads=%(threads)d"
        
        samplingRate = self._getInputParticles().getSamplingRate()
        params = {'resol': self.resol.get(),
                  'speed': int(self.getEnumText('speed')),
                  'numberOfIterations': self.numberOfIterations.get(),
                  'sym': self.symmetry.get(),
                  'molMass': self.molMass.get(),
                  'samplingRate': samplingRate,
                  'classKeep': self.classKeep.get(),
                  'm3dKeep': self.m3dKeep.get(),
                  'threads': self.numberOfThreads.get()
                  }
        args = args % params
         
        if self.doBreaksym:
            args += " --breaksym"
        if self.useE2make3d:
            args += " --m3dold"
        if self.useSetsfref:
            args += " --classrefsf"
        if self.doAutomask:
            args += " --classautomask"
        if self.doThreshold:
            args += " --prethreshold"
        if self.m3dPostProcess.get() > FILTER_NONE:
            args += " --m3dpostprocess=%s" % self.getEnumText('m3dPostProcess')
        return args
    
    def _getRun(self):
        
        if not self.doContinue:
            return 1
        else:
            files = sorted(glob(self.continueRun.get()._getExtraPath("refine*")))
            if files:
                f = files[-1]
                refineNumber = int(f.split("_")[-1]) + 1
            return refineNumber
    
    def _getBaseName(self, key, **args):
        """ Remove the folders and return the file from the filename. """
        return os.path.basename(self._getFileName(key, **args))
    
    def _getParticlesStack(self):
        if not self._getInputParticles().isPhaseFlipped() and self._getInputParticles().hasCTF():
            return self._getFileName("partFlipSet")
        else:
            return self._getFileName("partSet")
    
    def _iterTextFile(self, iterN):
        f = open(self._getFileName('angles', iter=iterN))
        
        for line in f:
            yield map(float, line.split())
            
        f.close()
    
    def _createItemMatrix(self, item, rowList):
        item.setTransform(rowToAlignment(rowList[1:], alignType=em.ALIGN_PROJ))
    
    def _getIterNumber(self, index):
        """ Return the list of iteration files, give the iterTemplate. """
        result = None
        files = sorted(glob(self._iterTemplate))
        if files:
            f = files[index]
            s = self._iterRegex.search(f)
            if s:
                result = int(s.group(1)) # group 1 is 3 digits iteration number
        return result
        
    def _lastIter(self):
        return self._getIterNumber(-1)

    def _firstIter(self):
        return self._getIterNumber(0) or 1

    def _getIterData(self, it):
        from convert import writeSqliteIterData
        data_sqlite = self._getFileName('data_scipion', iter=it)
        partSet = self._getInputParticles()
        if not exists(data_sqlite):
            clsFn = self._getFileName("cls", run=self._getRun(), iter=it)
            classesFn = self._getFileName("classes", run=self._getRun(), iter=it)
            angles = self._getBaseName('angles', iter=it)
            
            proc = createEmanProcess(args='read %s %s %s %s' 
                         % (self._getParticlesStack(), clsFn, classesFn, angles),
                         direc=self._getExtraPath())
            proc.wait()
            writeSqliteIterData(partSet, data_sqlite, self._createItemMatrix, self._iterTextFile(it))
        return data_sqlite
    
    def _getInputParticlesPointer(self):
        if self.doContinue:
            self.inputParticles.set(self.continueRun.get().inputParticles.get())
        return self.inputParticles
    
    def _getInputParticles(self):
        return self._getInputParticlesPointer().get()
    
    def _execEmanProcess(self, numRun, iterN):
        from pyworkflow.utils.path import cleanPath
        clsFn = self._getFileName("cls", run=numRun, iter=iterN)
        classesFn = self._getFileName("classes", run=numRun, iter=iterN)
        angles = self._getFileName('angles', iter=iterN)
        
        if not exists(angles) and exists(self._getFileName('clsEven', run=numRun, iter=iterN)):
            proc = createEmanProcess(args='read %s %s %s %s'
                                     % (self._getParticlesStack(), clsFn, classesFn,
                                        self._getBaseName('angles', iter=iterN)),
                                     direc=self._getExtraPath())
            proc.wait()

