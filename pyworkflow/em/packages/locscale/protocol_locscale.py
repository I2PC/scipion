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
from pyworkflow.em.data import Volume, PdbFile
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
        form.addParam('refObj', params.PointerParam, allowsNull=True,
                      label="Input pdb", pointerClass='PdbFile, Volume',
                      condition='inputPdbData == IMPORT_OBJ',
                      help='Choose a pdb object.')
        form.addParam('inputPath', params.FileParam,
                      label="File path", allowsNull=True,
                      condition='inputPdbData == IMPORT_FROM_FILES',
                      help='Specify a path to desired PDB structure.')

        form.addParam('binaryMask', params.PointerParam, pointerClass='VolumeMask',
                    label='3D mask', allowsNull=True,
                    help='Binary mask: 0 to ignore that voxel and 1 to process it.')

        form.addParam('patchSize', params.IntParam,
                      label='Patch size', help='Window size for local scale.\n'
                            'Recomended: 7*average_map_resolution/pixel_size')

        form.addParam('doParalelize', params.BooleanParam, default=False,
                      label='Paralelize', help='Do it in paralel with 4 MPI.')
    
    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        if self.inputPdbData == self.IMPORT_FROM_ID:
            self._insertFunctionStep('pdbDownloadStep')
        if self.isRefInputPdb():
            self._insertFunctionStep('prepareInput')

        self._insertFunctionStep('refineStep')
        self._insertFunctionStep('createOutputStep')
    
    #--------------------------- STEPS functions -------------------------------
    def pdbDownloadStep(self):
        """Download the pdb file and unzip it."""
        downloadPdb(self.pdbId.get(), self.getPdbFn(), self._log)

    def prepareInput(self):
        """Convert the pdb to a density map."""
        pdbFn = self.getPdbFn()
        refVol = removeExt(self.getRefFn())

        args  = "-i '%s'" % pdbFn
        args += " -o '%s'" % refVol
        args += " --centerPDB"
        args += " --size %d" % self.inputVolume.get().getDim()[0]

        self.info("Input file: " + pdbFn)
        self.info("Output file: " + refVol)
        
        program = "xmipp_volume_from_pdb"

        import xmipp3
        xmipp3.runXmippProgram(program, args)


    def refineStep(self):
        """ Run the LocScale program (with EMAN enviroment)
            to refine a volume. """
        # self.info('Checking the origin of the map:')
        # check = getProgram('check_and_set_ori_zero.py')
        # checkArgs = '--map %s' % self.getAbsPath(self.getInputFn())
        # self.runJob(check, checkArgs, cwd=self._getExtraPath())

        # self.info('Checking the origin of the ref:')
        # check = getProgram('check_and_set_ori_zero.py')
        # checkArgs = '--map %s' % self.getAbsPath(self.getRefFn())
        # self.runJob(check, checkArgs, cwd=self._getExtraPath())

        args = self.prepareParams()

        self.info("Launching LocScale method")
        program = getProgram('locscale_mpi.py')
        self.runJob(program, args, cwd=self._getExtraPath())
    
    def createOutputStep(self):
        volume = Volume()
        volume.setSamplingRate(self.getSampling())
        volume.setFileName(self.getOutputFn())
        self._defineOutputs(outputVolume=volume)
        self._defineTransformRelation(self.inputVolume, volume)
    
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        errors = []
        validateEmanVersion(self, errors)

        return errors
    
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputVolume'):
            summary.append("Output volumes not ready yet.")
        else:
            methodsStr  = 'We obtained a sharpened volume of the '
            methodsStr += str(self.getObjectTag('inputVolume'))

            if self.inputPdbData == self.IMPORT_OBJ:
                methodsStr += ' using %s as reference.' % self.getObjectTag('refObj')
            elif self.inputPdbData == self.IMPORT_FROM_ID:
                methodsStr += ' using a dowloaded model of the %s pdb as reference.'\
                              % self.pdbId
            elif self.inputPdbData == self.IMPORT_FROM_FILES:
                methodsStr += ' using the %s as reference.' % self.inputPath

            methodsStr += '\nFor more information see https://doi.org/10.7554/eLife.27131'
            summary.append(methodsStr)
        return summary

    def _methods(self):
        methodsStr  = 'We obtained a sharpened volume of the '
        methodsStr += str(self.getObjectTag('inputVolume'))

        if self.inputPdbData == self.IMPORT_OBJ:
            methodsStr += ' using %s as reference.' % self.getObjectTag('refObj')
        elif self.inputPdbData == self.IMPORT_FROM_ID:
            methodsStr += ' using a dowloaded model of the %s pdb as reference.'\
                          % self.pdbId
        elif self.inputPdbData == self.IMPORT_FROM_FILES:
            methodsStr += ' using the %s as reference.' % self.inputPath

        methodsStr += '\nFor more information see https://doi.org/10.7554/eLife.27131'
        return [methodsStr]

    def _citations(self):
        return ['Jakobi2017']
    
    #--------------------------- UTILS functions --------------------------------------------
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

        args  = " --em_map '%s'" % self.getAbsPath(self.getInputFn())
        args += " --model_map '%s'" % self.getAbsPath(self.getRefFn())
        args += " --apix %f" % self.getSampling()
        
        if self.binaryMask.hasValue():
            maskFn = self.getAbsPath(self.binaryMask.get().getFileName())
            args += " --mask '%s'" % maskFn

        args += " --window_size %d" % self.patchSize
        args += " -o '%s'" % self.getAbsPath(self.getOutputFn())

        if self.doParalelize:
            args += " --mpi"

        return args

    def getSampling(self):
        return self.inputVolume.get().getSamplingRate()

    def getInputFn(self):
        return self.inputVolume.get().getFileName()

    def getPdbFn(self):
        if self.inputPdbData == self.IMPORT_FROM_ID:
            fn = self._getExtraPath('%s.pdb' % self.pdbId.get())
        elif self.inputPdbData == self.IMPORT_OBJ:
            fn = self.refObj.get().getFileName()
        elif self.inputPdbData == self.IMPORT_FROM_FILES:
            fn = self.inputPath.get()
        else:
            fn = None
        return fn

    def isRefInputPdb(self):
        """ The reference data can be either a pdb or a volume."""
        if self.inputPdbData == self.IMPORT_OBJ:
            result = isinstance(self.refObj.get(),PdbFile)
        elif self.inputPdbData == self.IMPORT_FROM_FILES:
            result = self.inputPath.endswith('.pdb')
        elif self.inputPdbData == self.IMPORT_FROM_ID:
            result = True
        else:
            result = False
        return result

    def getRefFn(self):
        if self.isRefInputPdb():
            fn = self._getExtraPath(replaceBaseExt(self.getPdbFn(),"vol"))
        else:
            fn = self.inputPath.get() if self.inputPdbData==self.IMPORT_FROM_FILES \
                 else self.refObj.get().getFileName()
        return fn

    def getOutputFn(self):
        inputFn = removeBaseExt(self.inputVolume.get().getFileName())
        return self._getExtraPath(inputFn) + '_scaled.mrc'

    def getAbsPath(self, fileName):
        return os.path.abspath(fileName).replace(":mrc","")
