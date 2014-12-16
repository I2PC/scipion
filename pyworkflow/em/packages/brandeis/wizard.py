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

from pyworkflow.em.constants import *
from constants import *

from protocol_ctffind3 import ProtCTFFind
import pyworkflow.gui.dialog as dialog
from pyworkflow.em.wizard import *
from protocol_refinement import ProtFrealign
from protocol_ml_classification import ProtFrealignClassify

from pyworkflow import findResource

#===============================================================================
# CTFs
#===============================================================================

class BrandeisCTFWizard(CtfWizard):
    _targets = [(ProtCTFFind, ['ctfDownFactor', 'lowRes', 'highRes'])]
    
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
                    'hf': _value[2]
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

#===============================================================================
# MASKS
#===============================================================================

class FrealignVolRadiiWizard(VolumeMaskRadiiWizard):
    _targets = [(ProtFrealign, ['innerRadius', 'outerRadius']),
                (ProtFrealignClassify, ['innerRadius', 'outerRadius'])]
    
    def _getParameters(self, protocol):
        
        label, value = self._getInputProtocol(self._targets, protocol)
        
        protParams = {}
        protParams['input']= protocol.input3DReference
        protParams['label']= label
        protParams['value']= value
        return protParams  
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']
        return VolumeMaskRadiiWizard._getListProvider(self, _objs)
    
    def show(self, form):
        params = self._getParameters(form.protocol)
        _value = params['value']
        _label = params['label']
        VolumeMaskRadiiWizard.show(self, form, _value, _label, UNIT_ANGSTROM)


#===============================================================================
# FILTERS
#===============================================================================
 
class FrealignBandpassWizard(FilterParticlesWizard):
    _targets = [(ProtFrealign, ['lowResolRefine', 'highResolRefine']),
                (ProtFrealignClassify, ['lowResolRefine', 'highResolRefine'])]
    
    def _getParameters(self, protocol):
        
        label, value = self._getInputProtocol(self._targets, protocol)
        
        protParams = {}
        protParams['input']= protocol.inputParticles
        protParams['label']= label
        protParams['value']= value
        protParams['mode'] = FILTER_NO_DECAY
        return protParams  
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']
        return FilterParticlesWizard._getListProvider(self, _objs)
    
    def show(self, form):
        protocol = form.protocol
        provider = self._getProvider(protocol)
        params = self._getParameters(protocol)

        if provider is not None:
            
            args = {'mode': params['mode'],                   
                    'lowFreq': params['value'][0],
                    'highFreq': params['value'][1],
                    'unit': UNIT_ANGSTROM
                    }
            
            args['showDecay'] = False

            d = BandPassFilterDialog(form.root, provider, **args)
            
            if d.resultYes():
                form.setVar('highResolRefine', d.samplingRate/d.getHighFreq())
                form.setVar('lowResolRefine', d.samplingRate/d.getLowFreq())
                
        else:
            dialog.showWarning("Input particles", "Select particles first", form.root)  
            
 
class FrealignVolBandpassWizard(FilterVolumesWizard):
    _targets = [(ProtFrealign, ['resolution']),
                (ProtFrealignClassify, ['resolution'])]
    
    def _getParameters(self, protocol):
        
        label, value = self._getInputProtocol(self._targets, protocol)
        
        protParams = {}
        protParams['input']= protocol.input3DReference
        protParams['label']= label
        protParams['value']= value
        protParams['mode'] = FILTER_LOW_PASS_NO_DECAY
        return protParams
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']  
        return FilterVolumesWizard._getListProvider(self, _objs)
    
    def show(self, form):
        params = self._getParameters(form.protocol)
        protocol = form.protocol
        provider = self._getProvider(protocol)

        if provider is not None:
            
            args = {'mode': params['mode'],                   
                    'highFreq': params['value'],
                    'unit': UNIT_ANGSTROM
                    }
            
            args['showLowFreq'] = False
            args['showDecay'] = False

            d = BandPassFilterDialog(form.root, provider, **args)
            
            if d.resultYes():
                form.setVar('resolution', d.samplingRate/d.getHighFreq())
                
        else:
            dialog.showWarning("Input volumes", "Select volumes first", form.root)  

    