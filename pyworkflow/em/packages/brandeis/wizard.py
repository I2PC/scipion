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
from pyworkflow.em.packages.xmipp3.constants import *
from pyworkflow.em.constants import *
from pyworkflow.viewer import Viewer, Wizard, DESKTOP_TKINTER, WEB_DJANGO
import pyworkflow.gui.dialog as dialog
from pyworkflow.em.packages.xmipp3.wizard import XmippRadiiWizard, XmippMaskRadiiPreviewDialog
from pyworkflow.em.packages.xmipp3.wizard import ListTreeProvider
from protocol_frealign import ProtFrealign

from pyworkflow import findResource

class FrealignRadiiWizard(XmippRadiiWizard):
    
    _targets = [(ProtFrealign, ['innerRadius', 'outerRadius'])]


#class FrealignRadiiWizard(XmippRadiiWizard):
#    
#    _targets = [(ProtFrealign, ['innerRadius', 'outerRadius'])]
#    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
#    
#    def show(self, form):
#       
#        protocol = form.protocol
#        provider = self._getProvider(protocol)
#
#        if provider is not None:
#            d = XmippMaskRadiiPreviewDialog(form.root, provider, innerRadius=protocol.innerRadius.get(), outerRadius=protocol.outerRadius.get())
#            if d.resultYes():
#                form.setVar('innerRadius', d.getRadius(d.radiusSliderIn))
#                form.setVar('outerRadius', d.getRadius(d.radiusSliderOut))
#        else:
#            dialog.showWarning("Empty input", "Select elements first", form.root)    
#    
#    @classmethod    
#    def getView(self):
#        return "wiz_volume_mask_radii"   

    
#class RelionBandpassWizard(XmippFilterParticlesWizard):
#    
#    _targets = [(Relion3DClassification, ['iniLowPassFilter'])]
#    _environments = [DESKTOP_TKINTER, WEB_DJANGO]    
#
#    def show(self, form):
#        protocol = form.protocol
#        provider = self._getProvider(protocol)
#
#        if provider is not None:
#            self.mode = FILTER_LOW_PASS_NO_DECAY
#            
#            args = {'mode':  self.mode,                   
#                    'highFreq': protocol.iniLowPassFilter.get(),
#                    'unit': UNIT_ANGSTROM
#                    }
#            
#            args['showLowFreq'] = False
#            args['showDecay'] = False
#
#            d = XmippBandPassFilterDialog(form.root, provider, **args)
#            
#            if d.resultYes():
#                print "SAMPLING RATE !!", d.samplingRate
#                print "ITEM DIM !!", d.itemDim
#               
##                1/self.hfSlider.get()*self.itemDim 
#                
#                form.setVar('iniLowPassFilter', 1/d.getHighFreq()*d.itemDim)
#                
#        else:
#            dialog.showWarning("Input particles", "Select particles first", form.root)  
#    
#    @classmethod    
#    def getView(self):
#        return "wiz_relion_bandpass"   

