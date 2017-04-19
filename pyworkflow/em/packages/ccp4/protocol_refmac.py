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
from pyworkflow.protocol.params import PointerParam, IntParam, FloatParam
from pyworkflow.em import Volume, PdbFile
from pyworkflow.em.packages.ccp4.refmac_template import template
import os
from convert import (adapBinFileToCCP4,runCCP4Program)
import stat


class CCP4ProtRunRefmac(EMProtocol):
    """ generates files for volumes and FSCs to submit structures to EMDB
    """
    _label = 'refmac'
    _program = ""
    _version = VERSION_1_2
    refmacScriptFileName = "refmac.sh"
    refmacOutPDBFileName = "%s_final.pdb"

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
                      label='Number of refinement cycles:',
                      help='Specify the number of cycles of refinement.\n')
        form.addParam('weightMatrix', FloatParam, default=0.01, expertLevel=const.LEVEL_ADVANCED, label= 'Matrix refinement weight:',
                      help='Weight between density map and chemical constrain. Smaller means less weight for EM map\n')



    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('createScriptFile')
        self._insertFunctionStep('executeRefmac')
        #TODO: convert input file to mrc if needed
        #TODO: add CRYS record if neede
        #TODO: pass all parameters to script

    #    self._insertFunctionStep('createRefmacOutputStep') #Llamada a Refmac y obtencion del output

#    def convertInputStep(self):
#        """ convert 3Dmaps to MRC '.mrc' format
#        """
#        for vol in self.inputVolumes:
#            inFileName  = vol.get().getFileName()
#            if inFileName.endswith('.mrc'):
#                inFileName = inFileName + ":mrc"
#            outFileName = self._getVolumeFileName(inFileName)
#            if self.doNormalize:
#                img = ImageHandler()._img
#                img.read(inFileName)
#                mean, dev, min, max = img.computeStats()
#                img.inplaceMultiply(1./max)
#                mean2, dev2, min2, max2 = img.computeStats()
#                img.write(outFileName)
#            else:
#                adapBinFileToCCP4(inFileName, outFileName)

    # --------------------------- STEPS functions --------------------------------------------
    def createScriptFile(self):
        f = open(self._getScriptFileName(), "w")
        dict = {}
        dict['CCP4_HOME'] = os.environ['CCP4_HOME']
        dict['REFMAC_BIN'] = os.environ['REFMAC_BIN']
        dict['PDBDIR'] =  os.path.dirname(self.inputStructure.get().getFileName())
        pdfileName = os.path.splitext(os.path.basename(self.inputStructure.get().getFileName()))[0]
        dict['PDBFILE'] = pdfileName
        dict['MAPFILE'] = self.inputVolume.get().getFileName()
        dict['CHIMERA_BIN'] = os.path.join(os.environ['CHIMERA_HOME'],"bin","chimera")
        data = template%dict
        f.write(data)
        f.close()
        os.chmod(self._getScriptFileName(), stat.S_IEXEC | stat.S_IREAD)
        #TODO script needs to be modified so output goes to extra dir

    def fixAdCRYSrecordToPDBFile(self):
        #read input pdb file
        #search for CRYS RECORD
        #if available do nothing
        #else create a new PDB file
        #this new file should be the input
        pass

    def executeRefmac(self):
        runCCP4Program(self._getScriptFileName())

    def createRefmacOutputStep(self):
        pdb = PdbFile()
        pdb.setFileName(self._getOutPdbFileName())
        self._defineOutputs(outputPdb=pdb)
        self._defineSourceRelation(self.inputStructure, self.outputPdb)
        self._defineSourceRelation(self.inputVolume, self.outputPdb)

    # --------------------------- INFO functions --------------------------------------------


    # --------------------------- UTLIS functions --------------------------------------------
    def _getOutPdbFileName(self):
        pdfileName = os.path.splitext(os.path.basename(self.inputStructure.get().getFileName()))[0]
        return self._getExtraPath(self.refmacOutPDBFileName%pdfileName)

    #def _getVolName(self):
        #return self.Volume.get()
        #pass
    def _getScriptFileName(self):
        return self._getTmpPath(self.refmacScriptFileName)


