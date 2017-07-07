# **************************************************************************
# *
# * Authors:  Carlos Oscar Sanchez Sorzano (coss@cnb.csic.es), May 2013
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

from pyworkflow.em import *  
from pyworkflow.utils import * 
from pyworkflow.em.convert import ImageHandler
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.em.packages.ccp4.convert import Ccp4Header

class XmippProt3DBionotes(ProtAnalysis3D):
    """ Protocol for checking annotations
    
    This is actually a wrapper to 3D Bionotes.
    See documentation at:
       http://3dbionotes.cnb.csic.es
    """
    _label = '3d bionotes'
    
    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputPDB', PointerParam, pointerClass='PdbFile',
                      label="Input PDB", important=True)
        form.addParam('inputVol', PointerParam, pointerClass='Volume',
                      label="Input volume", important=True)
    
    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('bionotesWrapper')
        
    #--------------------------- STEPS functions -------------------------------
    def bionotesWrapper(self):
        img = ImageHandler()
        fnVol = self._getExtraPath('volume.mrc')
        vol = self.inputVol.get()
        img.convert(vol,fnVol)

        ccp4header = Ccp4Header(fnVol, readHeader= True)
        ccp4header.setOffset(vol.getOrigin())
        ccp4header.setSampling(vol.getSamplingRate())
        ccp4header.writeHeader()

       
    #--------------------------- INFO functions --------------------------------
            