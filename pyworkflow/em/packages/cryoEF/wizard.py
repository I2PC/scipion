# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
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

from pyworkflow.em.constants import UNIT_ANGSTROM
from pyworkflow.em.wizard import ParticleMaskRadiusWizard
from protocol_cryoef import ProtCryoEF


class cryoEFMaskDiameterWizard(ParticleMaskRadiusWizard):
    _targets = [(ProtCryoEF, ['diam'])]
    _unit = UNIT_ANGSTROM

    def _getParameters(self, protocol):
        label, value = self._getInputProtocol(self._targets, protocol)

        protParams = {}
        protParams['input'] = self._getProtocolImages(protocol)
        protParams['label'] = label
        protParams['value'] = value / 2
        return protParams

    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']
        return ParticleMaskRadiusWizard._getListProvider(self, _objs)

    def show(self, form):
        params = self._getParameters(form.protocol)
        _value = params['value']
        _label = params['label']
        ParticleMaskRadiusWizard.show(self, form, _value, _label, units=self._unit)

    def setVar(self, form, label, value):
        # adjust again from radius to diameter
        form.setVar(label, value * 2)

    def _getProtocolImages(self, protocol):
        return protocol.inputParticles
