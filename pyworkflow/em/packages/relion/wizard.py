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

from pyworkflow.em.packages.xmipp3.constants import *

from constants import *
from pyworkflow.em import *
from pyworkflow.em.wizard import *
from protocol_classify3d import ProtRelionClassify3D<
from protocol_refine3d import ProtRelionRefine3D
from protocol_classify2d import ProtRelionClassify2D

#===============================================================================
# MASKS
#===============================================================================

class RelionPartMaskRadiusWizard(ParticleMaskRadiusWizard):
    _targets = [(ProtRelionClassify2D, ['maskRadiusA'])]
    
    def _getParameters(self, protocol):
        protParams = {}
        protParams['input']= protocol.inputParticles
        protParams['label']= "radius"
        protParams['value']= protocol.radius.get()
        return protParams  
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']    
        return ParticleMaskRadiusWizard._getListProvider(self, _objs)
    
    def show(self, form):
        params = self._getParameters(form.protocol)
        _value = params['value']
        _label = params['label']
        ParticleMaskRadiusWizard.show(self, form, _value, _label, UNIT_PIXEL)
    
    
class RelionVolMaskRadiusWizard(VolumeMaskRadiusWizard):
    _targets = [(ProtRelionClassify3D, ['maskRadiusA']),
                (ProtRelionRefine3D, ['maskRadiusA'])]
    
    def _getParameters(self, protocol):
        protParams = {}
        protParams['input']= protocol.inputVolumes
        protParams['label']= "radius"
        protParams['value']= protocol.radius.get()
        return protParams  
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']    
        return VolumeMaskRadiusWizard._getListProvider(self, _objs)
    
    def show(self, form):
        params = self._getParameters(form.protocol)
        _value = params['value']
        _label = params['label']
        VolumeMaskRadiusWizard.show(self, form, _value, _label, UNIT_PIXEL)
        

#===============================================================================
# FILTER
#===============================================================================

class RelionBandpassWizard(FilterParticlesWizard):
    _targets = [(ProtRelionClassify3D, ['iniLowPassFilter'])]
    
    def _getParameters(self, protocol):
        protParams = {}
        protParams['input']= protocol.inputParticles
        protParams['label']= "radius"
        protParams['value']= protocol.radius.get()
        return protParams  
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']
        return FilterParticlesWizard._getListProvider(self, _objs)
  
    def show(self, form):
        params = self._getParameters(form.protocol)
        protocol = form.protocol
        provider = self._getProvider(protocol)

        if provider is not None:
            self.mode = FILTER_LOW_PASS_NO_DECAY
            
            args = {'mode':  self.mode,                   
                    'highFreq': params['value'],
                    'unit': UNIT_ANGSTROM
                    }
            
            args['showLowFreq'] = False
            args['showDecay'] = False

            d = BandPassFilterDialog(form.root, provider, **args)
            
            if d.resultYes():
                form.setVar('iniLowPassFilter', 1/d.getHighFreq()*d.itemDim)
        else:
            dialog.showWarning("Input particles", "Select particles first", form.root)  

