# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (josue.gomez-blanco@mcgill.ca)
# *              Javier Vargas Balbuena (javier.vargasbalbuena@mcgill.ca)
# *
# * Department of Anatomy and Cell Biology, McGill University
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
from pyworkflow.em.constants import FILTER_LOW_PASS_NO_DECAY
from pyworkflow.em.wizard import (ParticleMaskRadiusWizard, UNIT_ANGSTROM,
                                  FilterVolumesWizard, BandPassFilterDialog,
                                  dialog)
from protocol_volume_selector import ProtInitialVolumeSelector

#===============================================================================
# MASKS
#===============================================================================

class BackRadiusWizard(ParticleMaskRadiusWizard):
    _unit = UNIT_ANGSTROM
    
    def _getProtocolImages(self, protocol):
        return protocol.inputParticles
    
    def _getParameters(self, protocol):
        
        label, value = self._getInputProtocol(self._targets, protocol)
        
        protParams = {}
        protParams['input']= self._getProtocolImages(protocol)
        protParams['label']= label
        protParams['value']= value
        return protParams  
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']    
        return ParticleMaskRadiusWizard._getListProvider(self, _objs)
    
    def show(self, form):
        params = self._getParameters(form.protocol)
        _value = params['value']
        _label = params['label']
        ParticleMaskRadiusWizard.show(self, form, _value,
                                      _label, units=self._unit)


class MaskDiameterWizard(BackRadiusWizard):
    _targets = [(ProtInitialVolumeSelector, ['maskDiameterA'])]

    def _getParameters(self, protocol):
        protParams = BackRadiusWizard._getParameters(self, protocol)
        # adjust to from diameter to radius
        protParams['value'] = protParams['value'] / 2

        return protParams

    def setVar(self, form, label, value):
        # adjust again from radius to diameter
        form.setVar(label, value * 2)

#===============================================================================
# FILTER
#===============================================================================
class BaseFilterWizard(FilterVolumesWizard):
    def _getParameters(self, protocol):
        
        label, value = self._getInputProtocol(self._targets, protocol)
        protParams = {}

        protParams['input'] = protocol.inputVolumes

        protParams['label'] = label
        protParams['value'] = value
        protParams['mode'] = FILTER_LOW_PASS_NO_DECAY
        return protParams  
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']    
        return FilterVolumesWizard._getListProvider(self, _objs)
    
    def show(self, form):
        params = self._getParameters(form.protocol)
        print(params['label'], self._targets)
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
                form.setVar(params['label'],
                            d.samplingRate/d.getHighFreq())

        else:
            dialog.showWarning("Input volumes",
                               "Select volumes first", form.root)

class TargetFilterWizard(BaseFilterWizard):
    _targets = [(ProtInitialVolumeSelector, ['targetResol'])]


class InitialPassFilterWizard(BaseFilterWizard):
    _targets = [(ProtInitialVolumeSelector, ['initialLowPassFilterA'])]
