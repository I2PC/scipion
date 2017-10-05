# **************************************************************************
# *
# * Authors:    C.O.S. Sorzano (coss@cnb.csic.es)
# *             David Maluenda (dmaluenda@cnb.csic.es)
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

from convert import rowToAlignment



                               
class EmanProtLocScale(em.ProtRefine3D):
    """
    This Protocol wraps *e2refine_easy.py* Eman2 program.
This is the primary single particle refinement program in EMAN2.1+. Major
features of this program:

- While a range of command-line options still exist. You should not normally
specify more than a few basic requirements. The rest will be auto-selected
for you.
- This program will split your data in half and automatically refine the halves
independently to produce a gold standard resolution curve for every step in the
refinement.
- The gold standard FSC also permits us to automatically filter the structure
at each refinement step. The resolution you specify is a target, not the filter
resolution.
    """
    _label = 'local scale'

    def _createFilenameTemplates(self):
        """ Centralize the names of the files. """
        
        myDict = {
                  'partSet': 'sets/inputSet.lst',
                  'partFlipSet': 'sets/inputSet__ctf_flip.lst',
                  'data_scipion': self._getExtraPath('data_scipion_it%(iter)02d.sqlite'),
                  'projections': self._getExtraPath('projections_it%(iter)02d_%(half)s.sqlite'),
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
                  'mapEvenUnmasked': self._getExtraPath('refine_%(run)02d/threed_even_unmasked.hdf'),
                  'mapOddUnmasked': self._getExtraPath('refine_%(run)02d/threed_odd_unmasked.hdf'),
                  'fscUnmasked': self._getExtraPath('refine_%(run)02d/fsc_unmasked_%(iter)02d.txt'),
                  'fscMasked': self._getExtraPath('refine_%(run)02d/fsc_masked_%(iter)02d.txt'),
                  'fscMaskedTight': self._getExtraPath('refine_%(run)02d/fsc_maskedtight_%(iter)02d.txt'),
                  }
        self._updateFilenamesDict(myDict)
    
    def _createIterTemplates(self, currRun):
        """ Setup the regex on how to find iterations. """
        self._iterTemplate = self._getFileName('mapFull', run=currRun,
                                               iter=1).replace('threed_01',
                                                               'threed_??')
        # Iterations will be identify by threed_XX_ where XX is the iteration
        #  number and is restricted to only 2 digits.
        self._iterRegex = re.compile('threed_(\d{2,2})')
    
    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

       
        form.addParam('inputVolume', PointerParam, pointerClass='Volume',
                    important=True, label='Input volume', help='Input EM volume')

        form.addParam('input3DReference', PointerParam,
                    pointerClass='Volume', label='3D reference volume:',
                    help='3D reference reconstruction.')

        form.addParam('binaryMask', PointerParam, pointerClass='VolumeMask',
                    label='3D mask', help='Binary mask: where 0 ignore voxels'
                          'and 1 let process it.')

        form.addParallelSection(mpi=0)
    
    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):        
        self._createFilenameTemplates()
        # self._createIterTemplates(self._getRun())
        # self._insertFunctionStep('convertImagesStep')
        args = self._prepareParams()
        self._insertFunctionStep('refineStep', args)
        self._insertFunctionStep('createOutputStep')
    
    #--------------------------- STEPS functions -------------------------------    
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
        cleanPattern(self._getExtraPath('refine_01'))
        program = getEmanProgram('e2refine_easy.py')
        self.runJob(program, args, cwd=self._getExtraPath())
    
    def createOutputStep(self):
        iterN = self.numberOfIterations.get()
        partSet = self._getInputParticles()
        numRun = self._getRun()
        
        vol = Volume()
        
        
        vol.setFileName(self._getFileName("mapFull",run=numRun, iter=iterN))
        halfMap1 = self._getFileName("mapEvenUnmasked", run=numRun)
        halfMap2 = self._getFileName("mapOddUnmasked", run=numRun)
        vol.setHalfMaps([halfMap1, halfMap2])
        vol.copyInfo(partSet)
        
        newPartSet = self._createSetOfParticles()
        newPartSet.copyInfo(partSet)
        self._fillDataFromIter(newPartSet, iterN)
        
        self._defineOutputs(outputVolume=vol)
        self._defineSourceRelation(self._getInputParticlesPointer(), vol)
        self._defineOutputs(outputParticles=newPartSet)
        self._defineTransformRelation(self._getInputParticlesPointer(), newPartSet)
    
    #--------------------------- INFO functions -------------------------------------------- 
    # def _validate(self):
    #     errors = []
    #     validateVersion(self, errors)

    #     particles = self._getInputParticles()
    #     samplingRate = particles.getSamplingRate()

    #     if self.resol <  2 * samplingRate:
    #         errors.append("\nTarget resolution is smaller than nyquist limit.")
        
    #     if not self.doContinue:
    #         self._validateDim(particles, self.input3DReference.get(), errors,
    #                           'Input particles', 'Reference volume')

    #     return errors
    
    # def _summary(self):
    #     summary = []
    #     if not hasattr(self, 'outputVolume'):
    #         summary.append("Output volumes not ready yet.")
    #     else:
    #         inputSize = self._getInputParticles().getSize()
    #         outputSize = self.outputParticles.getSize()
    #         diff = inputSize - outputSize
    #         if diff > 0:
    #             summary.append("Warning!!! There are %d particles "
    #                            "belonging to empty classes." % diff)
    #     return summary
    
    #--------------------------- UTILS functions --------------------------------------------
    def _prepareParams(self):
        # '-em', '--em_map', required=True, help='Input filename EM map')
        # '-mm', '--model_map', required=True, help='Input filename PDB map')
        # '-p', '--apix', type=float, required=True, help='pixel size in Angstrom')
        # '-ma', '--mask', help='Input filename mask')
        # '-w', '--window_size', type=int, help='window size in pixel')
        # '-o', '--outfile', required=True, help='Output filename')
        # '-mpi', '--mpi', action='store_true', default=False,
        #                  help='MPI version call by: \"{0}\"'.format(mpi_cmd))

        volumeFn = os.path.abspath(self.inputVolume.get().getFileName()).replace(":mrc","")

        args = '--em_map %s' % volumeFn

        modelFn = os.path.abspath(self.input3DReference.get().getFileName()).replace(":mrc","")
        
        args += ' --model_map %s' % modelFn
        args += ' --apix %s' % self.inputVolume.get().getSamplingRate()

        maskFn = '';
        if self.binaryMask.hasValue():
            maskFn = os.path.abspath(self.binaryMask.get().getFileName()).replace(":mrc","")

        args += ' --mask %s' % maskFn


        # maskPath = os.path.relpath(self.binaryMask.get().getFileName(), 
        #                          self._getExtraPath()).replace(":mrc","")

        program = getEmanProgram('locscale_mpi.py')
        self.runJob(program, args, cwd=self._getExtraPath())

        # print(args)



        args = "--em_map=%(volumeFn)s --model_map=%(refVolFn)s"
        args += "--apix=%(samplingRate)d"
        
        volume = os.path.relpath(self.input3DReference.get().getFileName(), 
                                 self._getExtraPath()).replace(":mrc","")
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
        if (not self._getInputParticles().isPhaseFlipped() and 
               self._getInputParticles().hasCTF()):
            return self._getFileName("partFlipSet")
        else:
            return self._getFileName("partSet")
    
    def _iterTextFile(self, iterN):
        f = open(self._getFileName('angles', iter=iterN))
        
        for line in f:
            yield map(float, line.split())
            
        f.close()
    
    def _createItemMatrix(self, item, rowList):
        if rowList[1] == 1:
            item.setTransform(rowToAlignment(rowList[2:], alignType=em.ALIGN_PROJ))
        else:
            setattr(item, "_appendItem", False)
    
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
        data_sqlite = self._getFileName('data_scipion', iter=it)
        if not exists(data_sqlite):
            iterImgSet = em.SetOfParticles(filename=data_sqlite)
            iterImgSet.copyInfo(self._getInputParticles())
            self._fillDataFromIter(iterImgSet, it)
            iterImgSet.write()
            iterImgSet.close()
        
        return data_sqlite
    
    def _getInputParticlesPointer(self):
        if self.doContinue:
            self.inputParticles.set(self.continueRun.get().inputParticles.get())
        return self.inputParticles
    
    def _getInputParticles(self):
        return self._getInputParticlesPointer().get()
    
    def _fillDataFromIter(self, imgSet, iterN):
        numRun = self._getRun()
        self._execEmanProcess(numRun, iterN)
        initPartSet = self._getInputParticles()
        imgSet.setAlignmentProj()
        partIter = iter(initPartSet.iterItems(orderBy=['_micId', 'id'],
                                              direction='ASC'))
        
        imgSet.copyItems(partIter,
                         updateItemCallback=self._createItemMatrix,
                         itemDataIterator=self._iterTextFile(iterN))
    
    def _execEmanProcess(self, numRun, iterN):
        clsFn = self._getFileName("cls", run=numRun, iter=iterN)
        classesFn = self._getFileName("classes", run=numRun, iter=iterN)
        angles = self._getFileName('angles', iter=iterN)
        
        if not exists(angles) and exists(self._getFileName('clsEven', run=numRun, iter=iterN)):
            proc = createEmanProcess(args='read %s %s %s %s'
                                     % (self._getParticlesStack(), clsFn, classesFn,
                                        self._getBaseName('angles', iter=iterN)),
                                        direc=self._getExtraPath())
            proc.wait()
