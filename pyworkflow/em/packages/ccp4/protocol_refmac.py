# **************************************************************************
# *
# * Authors:     Marta Martinez (mmmtnez@cnb.csic.es)
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

from pyworkflow import VERSION_1_2
from pyworkflow.em.protocol import EMProtocol
import pyworkflow.protocol.constants as const
from pyworkflow.protocol.params import PointerParam, IntParam, FloatParam, BooleanParam
from pyworkflow.em import Volume, PdbFile
from pyworkflow.em.packages.ccp4.refmac_template import template
from pyworkflow.em.pdb_handler import fixCRYSrecordToPDBFile
import os
from convert import (adaptBinFileToCCP4, runCCP4Program)
import stat
from tempfile import mkdtemp


class CCP4ProtRunRefmac(EMProtocol):
    """ generates files for volumes and FSCs to submit structures to EMDB
    """
    _label = 'refmac'
    _program = ""
    _version = VERSION_1_2
    refmacScriptFileName = "refmac.sh"
    OutPdbFileName = "%s-refined.pdb"
    createMaskLogFileName="mask.log"
    refineLogFileName = "refine.log"
    fftLogFileName = "ifft.log"
    maskedMapFileName = "masked_fs"

    def __init__(self, **kwargs):
            EMProtocol.__init__(self, **kwargs)

    #--------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputVolume', PointerParam, label="Input Volume", important=True,
                      pointerClass='Volume',
                      help='This is the unit cell volume.')#que pasa si la extension no es mrc?
        form.addParam('inputStructure', PointerParam, label="Input PDB file", important=True,
                      pointerClass='PdbFile', help='Specify a PDB object.')
        form.addParam('maxResolution', FloatParam, default=5,
                      label='Max. Resolution (A):', help="Max resolution used in the refinement (Angstroms).")
        form.addParam('minResolution', FloatParam, default=200,
                      label='Min. Resolution (A):', help="Min resolution used in the refinement (Angstroms).")
        form.addParam('nRefCycle', IntParam, default=30, expertLevel=const.LEVEL_ADVANCED,
                      label='Number of refinement iterations:',
                      help='Specify the number of cycles of refinement.\n')
        form.addParam('weightMatrix', FloatParam, default=0.01, expertLevel=const.LEVEL_ADVANCED, label= 'Matrix refinement weight:',
                      help='Weight between the density map and the chemical constraints. Smaller means less weight for the EM map.\n')
        form.addParam('generateMaskedVolume', BooleanParam, default=True, label="Generate masked volume",
                      expertLevel=const.LEVEL_ADVANCED, important=True, help='If set to True, the masked volume will be generated')
        form.addParam('SFCALCmapradius', FloatParam, default=3, expertLevel=const.LEVEL_ADVANCED,
                      label='SFCALC mapradius:', help='Specify how much around molecule should be cut (Angstroms)')
        form.addParam('SFCALCmradius', FloatParam, default=3, expertLevel=const.LEVEL_ADVANCED,
                      label='SFCALC mradius:', help='Specify the radius (Angstroms)to calculate the mask around molecule')
        form.addParam('BFactorSet', FloatParam, default=0, expertLevel=const.LEVEL_ADVANCED,
                      label='B Factor:', help='Specify the B factor value prior to refinement')
        form.addParam('RefiSharpen', FloatParam, default=0, expertLevel=const.LEVEL_ADVANCED,
                      label='Map sharpening:', help='Specify the map sharpening to be used during refinement')


    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('fixCRYSrecordToPDBFileStep')
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('createScriptFileStep')
        self._insertFunctionStep('executeRefmacStep')
        self._insertFunctionStep('createRefmacOutputStep') #Llamada a Refmac y obtencion del output
        self._insertFunctionStep('writeFinalResultsTable') #Print output results

    def convertInputStep(self):
        """ convert 3Dmaps to MRC '.mrc' format
        """
        # get input 3D map filename
        inFileName  = self.inputVolume.get().getFileName()
        # create local copy of 3Dmap
        localInFileName = self._getVolumeFileName()
        adaptBinFileToCCP4(inFileName, localInFileName,
                           self.inputVolume.get().getOrigin().getShifts())

    # --------------------------- STEPS functions --------------------------------------------
    def createScriptFileStep(self):
        dict = {}
        dict['CCP4_HOME'] = os.environ['CCP4_HOME']
        dict['REFMAC_BIN'] = os.environ['REFMAC_BIN']
        #dict['PDBDIR'] =  os.path.dirname(self.inputStructure.get().getFileName())
        dict['PDBDIR'] =  os.path.dirname(self.fixedPDBFileName)
        pdfileName = os.path.splitext(os.path.basename(self.inputStructure.get().getFileName()))[0]
        dict['PDBFILE'] = pdfileName

        dict['MAPFILE'] = self.inputVolume.get().getFileName().replace(':mrc', '')

#       dict['CHIMERA_BIN'] = os.path.join(os.environ['CHIMERA_HOME'],"bin","chimera")
        dict['RESOMIN'] = self.minResolution.get()
        dict['RESOMAX'] = self.maxResolution.get()
        dict['NCYCLE'] = self.nRefCycle.get()
        dict['WEIGHT MATRIX'] = self.weightMatrix.get()
        dict['OUTPUTDIR'] = self._getExtraPath('')
        dict['MASKED_VOLUME'] = self.generateMaskedVolume.get()
        dict['SFCALC_mapradius'] = self.SFCALCmapradius.get()
        dict['SFCALC_mradius'] = self.SFCALCmradius.get()
        dict['XDIM'] = self.inputVolume.get().getDim()[0]
        dict['YDIM'] = self.inputVolume.get().getDim()[1]
        dict['ZDIM'] = self.inputVolume.get().getDim()[2]
        if self.BFactorSet.get() ==0:
            dict['BFACTOR_SET'] = "#BFACtor SET 0"
        else:
            dict['BFACTOR_SET'] = "BFACtor SET %f"%self.BFactorSet.get()
        if self.RefiSharpen.get() ==0:
            dict['REFI_SHARPEN'] = "#REFI sharpen SET 0"
        else:
            dict['REFI_SHARPEN'] = "REFI sharpen %f"%self.RefiSharpen.get()

        data = template%dict
        f = open(self._getScriptFileName(), "w")
        f.write(data)
        f.close()
        os.chmod(self._getScriptFileName(), stat.S_IEXEC | stat.S_IREAD | stat.S_IWRITE)


    def fixCRYSrecordToPDBFileStep(self):
        self.fixedPDBFileName = fixCRYSrecordToPDBFile(self.inputStructure.get().getFileName(),
                                     self._getTmpPath(),
                                     x= self.inputVolume.get().getDim()[0],
                                     y=self.inputVolume.get().getDim()[1],
                                     z= self.inputVolume.get().getDim()[2],
                                     alpha=90., beta=90.,gamma=90.
                                                       )
    def executeRefmacStep(self):
        # Generic is a env variable that coot uses as base dir for some
        # but not all files. "" force a trailing slash
        runCCP4Program(self._getScriptFileName(),"",{'GENERIC':self._getExtraPath("")})

    def createRefmacOutputStep(self):
        pdb = PdbFile()
        pdb.setFileName(self._getOutPdbFileName())
        self._defineOutputs(outputPdb=pdb)
        self._defineSourceRelation(self.inputStructure, self.outputPdb)
        self._defineSourceRelation(self.inputVolume, self.outputPdb)

    def writeFinalResultsTable(self):
        with open(self._getlogFileName()) as input_data:
            for line in input_data:
                if line.strip() == '$TEXT:Result: $$ Final results $$':
                    break
            for line in input_data:
                if line.strip() == '$$':
                    break
                print line

    # --------------------------- INFO functions --------------------------------------------


    # --------------------------- UTLIS functions --------------------------------------------
    def _getOutPdbFileName(self):
        pdfileName = os.path.splitext(os.path.basename(self.inputStructure.get().getFileName()))[0]
        return self._getExtraPath(self.OutPdbFileName%pdfileName)

    #def _getVolName(self):
        #return self.Volume.get()
        #pass
    def _getScriptFileName(self):
        return self._getTmpPath(self.refmacScriptFileName)

    def _getlogFileName(self):
        return self._getExtraPath(self.refineLogFileName)

    def _getVolumeFileName(self, baseFileName="tmp3DMapFile.mrc"):
        return self._getExtraPath(baseFileName)









