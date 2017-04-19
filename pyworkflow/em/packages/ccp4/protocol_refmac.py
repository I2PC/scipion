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

class CCP4ProtRunRefmac(EMProtocol):
    """ generates files for volumes and FSCs to submit structures to EMDB
    """
    _label = 'refmac'
    _program = ""
    _version = VERSION_1_2

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    #--------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputVolumes', PointerParam, label="Input Volume", important=True,
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
    #def _insertAllSteps(self):
    #    self._insertFunctionStep(('modCoordFileStep') #Modificacion del fichero de coordenadas
    #    self._insertFunctionStep('createRefmacOutputStep') #Llamada a Refmac y obtencion del output

    # --------------------------- STEPS functions --------------------------------------------
    def modCoordFileStep(self):

        pass

    def createRefmacOutputStep(self):
        # vol = Volume()
        # vol.setLocation(self._getVolName())
        pass
            

    # --------------------------- INFO functions --------------------------------------------


        # --------------------------- UTLIS functions --------------------------------------------
    #def _getPdbFileName(self):
        #return self.pdbFile.get()
        #pass
    #def _getVolName(self):
        #return self.Volume.get()
        #pass
