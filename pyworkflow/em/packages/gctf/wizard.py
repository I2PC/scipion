# **************************************************************************
# *
# * Authors:     Jose Gutierrez (jose.gutierrez@cnb.csic.es)
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
"""
This module implement some wizards
"""

import os
import Tkinter as tk
import ttk

from pyworkflow.em.constants import *

from protocol_gctf import ProtGctf
import pyworkflow.gui.dialog as dialog
from pyworkflow.em.wizard import *

from pyworkflow import findResource

#===============================================================================
# CTFs
#===============================================================================

class GctfCTFWizard(CtfWizard):
    _targets = [(ProtGctf, ['ctfDownFactor', 'lowRes', 'highRes'])]
    
    def _getParameters(self, protocol):
        
        label, value = self._getInputProtocol(self._targets, protocol)
        
        protParams = {}
        protParams['input']= protocol.inputMicrographs
        protParams['label']= label
        protParams['value']= value
        return protParams
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']
        return CtfWizard._getListProvider(self, _objs)
        
    def show(self, form):
        protocol = form.protocol
        params = self._getParameters(protocol)
        _value = params['value']
        _label = params['label']
        
#        form.setParamFromVar('inputMicrographs') # update selected input micrographs
        provider = self._getProvider(protocol)
        
        if provider is not None:
            args = {'unit': UNIT_PIXEL,
                    'downsample': _value[0],
                    'lf': _value[1],
                    'hf': _value[2],
                    'showInAngstroms': True
                    }
            d = CtfDownsampleDialog(form.root, provider, **args)

            if d.resultYes():
                form.setVar(_label[0], d.getDownsample())
                form.setVar(_label[1], d.getLowFreq())
                form.setVar(_label[2], d.getHighFreq())
        else:
            dialog.showWarning("Empty input", "Select elements first", form.root)    
    
    @classmethod    
    def getView(self):
        return "wiz_ctf_downsampling"  
