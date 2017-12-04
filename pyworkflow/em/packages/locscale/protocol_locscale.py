# **************************************************************************
# *
# * Authors:    David Maluenda (dmaluenda@cnb.csic.es)
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
from pyworkflow.em import ProtRefine3D
from pyworkflow.protocol.params import PointerParam
from pyworkflow.em.data import Volume
from locscale import *

                               
class ProtLocScale(ProtRefine3D):
    """ This Protocol refine volume using a PDB model previously prepared."""
    _label = 'local scale'

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

        form.addParallelSection(mpi=1)
    
    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):        
        args = self._prepareParams()
        self._insertFunctionStep('refineStep', args)
        self._insertFunctionStep('createOutputStep')
    
    #--------------------------- STEPS functions -------------------------------    
    def refineStep(self, args):
        """ Run the EMAN program to refine a volume. """
        program = getProgram('locscale_mpi.py')
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
    def _validate(self):
        errors = []
        validateEmanVersion(self, errors)

        return errors
    
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

    def _citations(self):
        return ['Jakobi2017']
    
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
        def getAbsPath(fileName):
            return os.path.abspath(fileName).replace(":mrc","")

        volumeFn = getAbsPath(self.inputVolume.get().getFileName())
        args  = ' --em_map %s' % volumeFn
        args += ' --apix %s' % self.inputVolume.get().getSamplingRate()
        
        modelFn = getAbsPath(self.input3DReference.get().getFileName())
        args += ' --model_map %s' % modelFn
        
        if self.binaryMask.hasValue():
            maskFn = getAbsPath(self.binaryMask.get().getFileName())
            args += ' --mask %s' % maskFn

        # if self.numberOfMpi>1:
        #     args += ' --mpi True'# % self.numberOfMpi

        args += ' -o %s' %self._getExtraPath('result.vol')

        return args