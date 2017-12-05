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
from pyworkflow.em import ProtRefine3D, downloadPdb
import pyworkflow.protocol.params as params
from pyworkflow.em.data import Volume
from pyworkflow.utils import replaceBaseExt, removeBaseExt, removeExt
from locscale import *

                               
class ProtLocScale(ProtRefine3D):
    """ This Protocol refine volume using a PDB model previously prepared."""
    _label = 'local scale'
    IMPORT_FROM_ID = 0
    IMPORT_OBJ = 1
    IMPORT_FROM_FILES = 2 

    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

       
        form.addParam('inputVolume', params.PointerParam, pointerClass='Volume',
                      important=True, label='Input volume',
                      help='Input EM volume')

        form.addParam('inputPdbData', params.EnumParam, 
                      choices=['id', 'object', 'file'],
                      label="Retrieve PDB model from", default=self.IMPORT_OBJ,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Model to apply the local scale.\n'
                           'The PDB can be retrieved from server, '
                           'using a pdb Object or from a local file.')
        form.addParam('pdbId', params.StringParam, 
                      label="Pdb Id ", allowsNull=True,
                      condition='inputPdbData == IMPORT_FROM_ID',
                      help='Type a pdb Id (four alphanumeric characters).')
        form.addParam('refObj', params.PointerParam, pointerClass='PdbFile',
                      label="Input pdb ", allowsNull=True,
                      condition='inputPdbData == IMPORT_OBJ',
                      help='Choose a pdb object.')
        form.addParam('pdbFile', params.FileParam,
                      label="File path", allowsNull=True,
                      condition='inputPdbData == IMPORT_FROM_FILES',
                      help='Specify a path to desired PDB structure.')

        form.addParam('binaryMask', params.PointerParam, pointerClass='VolumeMask',
                    label='3D mask',
                    help='Binary mask: 0 to ignore that voxel and 1 to process it.')

        form.addParam('doParalelize', params.BooleanParam, default=False,
                    label='Paralelize', help='Do it in paralel with 4 MPI.')
    
    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        if self.inputPdbData == self.IMPORT_FROM_ID:
            self._insertFunctionStep('pdbDownloadStep')
        args = self._prepareParams()
        self._insertFunctionStep('prepareInput')
        self._insertFunctionStep('refineStep', args)
        self._insertFunctionStep('createOutputStep')
    
    #--------------------------- STEPS functions -------------------------------
    def pdbDownloadStep(self):
        """Download all pdb files in file_list and unzip them."""
        downloadPdb(self.pdbId.get(), self._getPdbFileName(), self._log)

    def prepareInput(self):
        """Convert the pdb to a density map."""
        pdbFn = self._getPdbFileName()
        outFile = removeExt(self._getVolName())
        args = '-i %s --sampling %f -o %s' % (pdbFn, self.sampling, outFile)
        
        # args += ' --centerPDB'

        args += ' --size %d' % self.inputVolume.get().getDim()[0]

        self.info("Input file: " + pdbFn)
        self.info("Output file: " +outFile)
        
        program = "xmipp_volume_from_pdb"
        self.runJob(program, args)

    def refineStep(self, args):
        """ Run the LocScale program (with EMAN enviroment)
            to refine a volume. """
        program = getProgram('locscale_mpi.py')
        self.runJob(program, args, cwd=self._getExtraPath())
    
    def createOutputStep(self):
        volume = Volume()
        volume.setSamplingRate(self.sampling)
        volume.setFileName(self._getVolName())
        self._defineOutputs(outputVolume=volume)
        self._defineSourceRelation(self.inputVolume, volume)
    
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

        self.sampling = self.inputVolume.get().getSamplingRate()

        volumeFn = getAbsPath(self.inputVolume.get().getFileName())
        args  = ' --em_map %s' % volumeFn
        args += ' --apix %s' % self.sampling
        
        modelFn = getAbsPath(self._getVolName())
        args += ' --model_map %s' % modelFn

        if self.binaryMask.hasValue():
            maskFn = getAbsPath(self.binaryMask.get().getFileName())
            args += ' --mask %s' % maskFn

        if self.doParalelize:
            args += ' --mpi'

        inputFn = removeBaseExt(self.inputVolume.get().getFileName())
        outFn = inputFn + '_scaled.vol'
        args += ' -o %s' %self._getExtraPath(outFn)

        return args

    def _getPdbFileName(self):
        if self.inputPdbData == self.IMPORT_FROM_ID:
            return self._getExtraPath('%s.pdb' % self.pdbId.get())
        elif self.inputPdbData == self.IMPORT_OBJ:
            return self.refObj.get().getFileName()
        else:
            return self.pdbFile.get()
    
    def _getVolName(self):
        return self._getExtraPath(replaceBaseExt(self._getPdbFileName(), "vol"))