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

from pyworkflow.em.constants import *
from constants import *

from pyworkflow.em.wizard import *
from protocol_classify3d import Relion3DClassification


class RelionVolMaskRadiusWizard(volumeMaskRadiusWizard):
    _targets = [(Relion3DClassification, ['maskRadius'])]


class RelionBandpassWizard(filterParticlesWizard):
    _targets = [(Relion3DClassification, ['iniLowPassFilter'])]

    def show(self, form):
        protocol = form.protocol
        provider = self._getProvider(protocol)

        if provider is not None:
            self.mode = FILTER_LOW_PASS_NO_DECAY
            
            args = {'mode':  self.mode,                   
                    'highFreq': protocol.iniLowPassFilter.get(),
                    'unit': UNIT_ANGSTROM
                    }
            
            args['showLowFreq'] = False
            args['showDecay'] = False

            d = bandPassFilterDialog(form.root, provider, **args)
            
            if d.resultYes():
                form.setVar('iniLowPassFilter', 1/d.getHighFreq()*d.itemDim)
        else:
            dialog.showWarning("Input particles", "Select particles first", form.root)  
    
    @classmethod    
    def getView(self):
        return "wiz_relion_bandpass"   

