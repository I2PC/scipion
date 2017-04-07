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
from pyworkflow.protocol.params import StringParam, PointerParam, IntParam, FloatParam, EnumParam, FileParam
from pyworkflow.em import Volume

class ProtRunRefmac(EMProtocol):
    """ generates files for volumes and FSCs to submit structures to EMDB
    """
    _label = 'Run REFMAC5'
    _program = ""
    _version = VERSION_1_2
    IMPORT_FROM_ID = 0
    IMPORT_OBJ = 1
    IMPORT_FROM_FILES = 2

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    #--------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputVolumes', PointerParam, label="Input Volume", important=True,
                      pointerClass='Volume',
                      help='This is the unit cell volume')
        form.addParam('inputPdbData', EnumParam, choices=['id', 'object', 'file'],
                      label="Retrieve PDB from", default=self.IMPORT_FROM_ID,
                      display= EnumParam.DISPLAY_HLIST,
                      help='Retrieve PDB data from server, use a pdb Object, or a local file')
        form.addParam('pdbId', StringParam, condition='inputPdbData == IMPORT_FROM_ID',
                      label="Pdb Id ", allowsNull=True,
                      help='Type a pdb Id (four alphanumeric characters).')
        form.addParam('pdbObj', PointerParam, pointerClass='PdbFile',
                      label="Input pdb ", condition='inputPdbData == IMPORT_OBJ', allowsNull=True,
                      help='Specify a pdb object.')
        form.addParam('pdbFile', FileParam,
                      label="File path", condition='inputPdbData == IMPORT_FROM_FILES', allowsNull=True,
                      help='Specify a path to desired PDB structure.')
        form.addParam('maxResolution', FloatParam, default=5,
                      label='Max. Resolution (A):', help="The reconstructed volume will be limited to \n"
                           "this resolution in Angstroms.")
        form.addParam('minResolution', FloatParam, default=1,
                      label='Max. Resolution (A):', help="The reconstructed volume will be limited to \n"
                                                         "this resolution in Angstroms.")
        form.addParam('nRefCycle', IntParam, default=30,
                      label='Number of refinement cycles:',
                      help='Specify the number of cycles of refinement.\n')
        form.addParam('weightMatrix', FloatParam, default=0.01, label= 'Matrix refinement weight:',
                      help='Set the refinement weight.\n')



    #--------------------------- INSERT steps functions --------------------------------------------


        # --------------------------- UTLIS functions --------------------------------------------
        def _getPdbFileName(self):
            if self.inputPdbData == self.IMPORT_FROM_ID:
                return self._getExtraPath('%s.pdb' % self.pdbId.get())
            elif self.inputPdbData == self.IMPORT_OBJ:
                return self.pdbObj.get().getFileName()
            else:
                return self.pdbFile.get()

        def _getVolName(self):
            return self._getExtraPath(replaceBaseExt(self._getPdbFileName(), "vol"))

