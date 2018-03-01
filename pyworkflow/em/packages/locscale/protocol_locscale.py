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
import pyworkflow.protocol.params as params
from pyworkflow.em.data import Volume
from pyworkflow.utils import replaceBaseExt, removeBaseExt, removeExt
from locscale import *

                               
class ProtLocScale(ProtRefine3D):
    """ This Protocol computes contrast-enhanced cryo-EM maps 
        by local amplitude scaling using a reference model.
    """
    
    _label = 'local scale'

    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputVolume', params.PointerParam, pointerClass='Volume',
                      important=True, label='Input volume',
                      help='Input EM volume')
        
        form.addParam('refObj', params.PointerParam, allowsNull=True,
                      label="Reference Volume", pointerClass='Volume',
                      help='Choose a model to take it as refference '
                           '(usually this volume should come from a PDB).')

        form.addParam('binaryMask', params.PointerParam,
                      pointerClass='VolumeMask',
                      label='3D mask', allowsNull=True,
                      help='Binary mask: 0 to ignore that voxel and 1 to process it.')

        form.addParam('patchSize', params.IntParam,
                      label='Patch size', help='Window size for local scale.\n'
                            'Recomended: 7 * average_map_resolution / pixel_size')

        form.addParallelSection(threads=0, mpi=4)
    
    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('refineStep')
        self._insertFunctionStep('createOutputStep')
    
    #--------------------------- STEPS functions -------------------------------
    def refineStep(self):
        """ Run the LocScale program (with EMAN enviroment)
            to refine a volume. 
        """
        self.info("Launching LocScale method")
        args = self.prepareParams()

        python, program_args = getProgram('locscale_mpi.py', args)
        self.runJob(python, program_args)
    
    def createOutputStep(self):
        volume = Volume()
        volume.setSamplingRate(self.getSampling())
        volume.setFileName(self.getOutputFn())
        self._defineOutputs(outputVolume=volume)
        self._defineTransformRelation(self.inputVolume, volume)
    
    #--------------------------- INFO functions --------------------------------
    def _validate(self):
        errors = []
        errors = validateEmanVersion(self, errors)

        inputSize = self.inputVolume.get().getDim()
        refSize = self.refObj.get().getDim()
        refSamp = self.refObj.get().getSamplingRate()

        if inputSize != refSize or self.getSampling() != refSamp:
            errors.append('Input volume and reference volume should be '
                          'of the same size and samplig rate')
        return errors

    def _warnings(self):
        warnings = []

        if self.binaryMask.hasValue() and self.binaryMask.get().getDim() != \
                                          self.inputVolume.get().getDim():
            warnings.append('Input volume and binary mask should be '
                          'of the same size')

        return warnings
    
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputVolume'):
            summary.append("Output volumes not ready yet.")
        else:
            summary.append('We obtained a sharpened volume of the %s '
                           'using %s as reference.'
                            % (self.getObjectTag('inputVolume'),
                               self.getObjectTag('refObj')))
        return summary

    def _methods(self):
        methods = []
        methods.append('LocScale has locally scaled the amplitude of the %s '
                       'using %s as reference with a window size of %d.'
                        % (self.getObjectTag('inputVolume'),
                           self.getObjectTag('refObj'),
                           self.patchSize))
        return methods


    def _citations(self):
        return ['Jakobi2017']
    
    #--------------------------- UTILS functions -------------------------------
    def prepareParams(self):
        """ The input params of the program are as follows:
        '-em', '--em_map', required=True, help='Input filename EM map')
        '-mm', '--model_map', required=True, help='Input filename PDB map')
        '-p', '--apix', type=float, required=True, help='pixel size in Angstrom')
        '-ma', '--mask', help='Input filename mask')
        '-w', '--window_size', type=int, help='window size in pixel')
        '-o', '--outfile', required=True, help='Output filename')
        '-mpi', '--mpi', action='store_true', default=False,
                         help='MPI version call by: \"{0}\"'.format(mpi_cmd)
        """
        args  =  "--em_map '%s'" % self.getInputFn()
        self.info("Input file: " + self.getInputFn())

        args += " --model_map '%s'" % self.getRefFn()
        self.info("Model file: " + self.getRefFn())

        args += " --apix %f" % self.getSampling()
        self.info("Samplig rate: %f" % self.getSampling())
        
        if self.binaryMask.hasValue():
            args += " --mask '%s'" % self.binaryMask.get().getFileName()
            self.info("Mask file: " + self.binaryMask.get().getFileName())

        args += " --window_size %d" % self.patchSize
        self.info("Window size: %d" % self.patchSize)

        if self.numberOfMpi>1:
            args += " -mpi"

        args += " -o '%s'" % self.getOutputFn()
        self.info("Output file: " + self.getOutputFn())

        return args

    def getSampling(self):
        return self.inputVolume.get().getSamplingRate()

    def getInputFn(self):
        return self.inputVolume.get().getFileName().replace(":mrc","")

    def getRefFn(self):
        return self.refObj.get().getFileName().replace(":mrc","")

    def getOutputFn(self):
        outputFnBase = removeBaseExt(self.getInputFn())
        return self._getExtraPath(outputFnBase) + '_scaled.mrc'
