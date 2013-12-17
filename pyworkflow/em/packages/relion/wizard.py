# **************************************************************************
# *
# * Authors:     Jose Gutierrez (jose.gutierrez@cnb.csic.es)
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
This module implement some wizards
"""

import os
import Tkinter as tk
import ttk

from pyworkflow.viewer import Viewer, Wizard, DESKTOP_TKINTER, WEB_DJANGO

from pyworkflow.em.packages.xmipp3.wizard import XmippVolumeMaskRadiusWizard
from pyworkflow.em.packages.xmipp3.wizard import ListTreeProvider
from protocol_classify3d import Relion3DClassification

from pyworkflow import findResource

class RelionVolMaskRadiusWizard(XmippVolumeMaskRadiusWizard):
    
    _targets = [(Relion3DClassification, ['maskRadius'])]    

    def _getProvider(self, protocol):
        """ This should be implemented to return the list
        of object to be displayed in the tree.
        """
        if protocol.input3DReferences.hasValue():
            vols = [vol for vol in protocol.input3DReferences.get()]
            return ListTreeProvider(vols)
        return None
    
    @classmethod    
    def getView(self):
        return "wiz_volume_mask"
    
    
    
    
    
    