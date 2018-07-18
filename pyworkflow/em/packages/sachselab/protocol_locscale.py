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

from locscale import *
from pyworkflow.em import Prot3D
from pyworkflow.protocol import params
from pyworkflow.em.data import Volume
from pyworkflow.utils import removeBaseExt
                               
class ProtLocScale(Prot3D):
    """ This Protocol computes contrast-enhanced cryo-EM maps 
        by local amplitude scaling using a reference model.
    """
    _label = 'locscale'

    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputVolume', params.PointerParam, pointerClass='Volume',
                      important=True, label='Input volume',
                      help='Input EM volume')
        
        form.addParam('refObj', params.PointerParam,
                      label="Reference Volume", pointerClass='Volume',
                      help='Choose a model to take it as refference '
                           '(usually this volume should come from a PDB).')

        form.addParam('binaryMask', params.PointerParam,
                      pointerClass='VolumeMask',
                      label='3D mask', allowsNull=True,
                      help='Binary mask (optional)')

        form.addParam('patchSize', params.IntParam, label='Patch size',
                      help='Window size for local scale.\n'
                           'Recommended: 7 * average_map_resolution / pixel_size')

        form.addParallelSection(threads=0, mpi=4)
    
    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertStep')
        self._insertFunctionStep('refineStep')
        self._insertFunctionStep('createOutputStep')
    
    #--------------------------- STEPS functions -------------------------------
    def convertStep(self):
        tmpFn = self._getTmpPath()
        self.inputVolFn = convertBinaryVol(self.inputVolume.get(), tmpFn)
        self.refVolFn = convertBinaryVol(self.refObj.get(), tmpFn)
        if self.binaryMask.hasValue():
            self.maskVolFn = convertBinaryVol(self.binaryMask.get(), tmpFn)

    def refineStep(self):
        """ Run the LocScale program (with EMAN enviroment)
            to refine a volume. 
        """
        self.info("Launching LocScale method")
        args = self.prepareParams()

        python, program = getEmanPythonProgram('locscale_mpi.py')
        program_args = program + ' ' + args
        self.runJob(python, program_args)
    
    def createOutputStep(self):
        """ Create the output volume 
        """
        outputVolume = Volume()
        outputVolume.setSamplingRate(self.getSampling())
        outputVolume.setFileName(self.getOutputFn())
        self._defineOutputs(outputVolume=outputVolume)
        self._defineTransformRelation(self.inputVolume, outputVolume)
    
    #--------------------------- INFO functions --------------------------------
    def _validate(self):
        """ We validate if eman is installed and if inputs make sense
        """
        errors = []
        errors = validateEmanVersion(errors)

        input = self.inputVolume.get()
        reference = self.refObj.get()
        if input is not None and reference is not None:
            inputSize = input.getDim()
            refSize = reference.getDim()
            refSamp = reference.getSamplingRate()

            if inputSize != refSize or self.getSampling() != refSamp:
                errors.append('Input volume and reference volume should be '
                              'of the same size and samplig rate')
        return errors

    def _warnings(self):
        """ The input volume and the mask should be of the same size
            but the progran can run.
        """
        warnings = []

        if self.binaryMask.hasValue() and \
            self.binaryMask.get().getDim() != self.inputVolume.get().getDim():
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
        """ The input params of the program are as follows (from source):
        '-em', '--em_map', required=True, help='Input filename EM map')
        '-mm', '--model_map', required=True, help='Input filename PDB map')
        '-p', '--apix', type=float, required=True, help='pixel size in Angstrom')
        '-ma', '--mask', help='Input filename mask')
        '-w', '--window_size', type=int, help='window size in pixel')
        '-o', '--outfile', required=True, help='Output filename')
        '-mpi', '--mpi', action='store_true', default=False,
                         help='MPI version call by: \"{0}\"'.format(mpi_cmd)
        """

        # Input volume
        args  =  "--em_map '%s'" % self.inputVolFn
        self.info("Input file: " + self.inputVolFn)

        # Reference volume
        args += " --model_map '%s'" % self.refVolFn
        self.info("Model file: " + self.refVolFn)

        # Samplig rate
        args += " --apix %f" % self.getSampling()
        self.info("Samplig rate: %f" % self.getSampling())
        
        # Mask
        if self.binaryMask.hasValue():
            args += " --mask '%s'" % self.maskVolFn
            self.info("Mask file: " + self.maskVolFn)

        # Windows size
        args += " --window_size %d" % self.patchSize
        self.info("Window size: %d" % self.patchSize)

        # MPI flag
        if self.numberOfMpi>1:
            args += " -mpi"

        # Output file
        args += " -o '%s'" % self.getOutputFn()
        self.info("Output file: " + self.getOutputFn())

        return args

    def getSampling(self):
        return self.inputVolume.get().getSamplingRate()

    def getOutputFn(self):
        """ Returns the scaled output file name
        """
        outputFnBase = removeBaseExt(self.inputVolFn)
        return self._getExtraPath(outputFnBase) + '_scaled.mrc'
