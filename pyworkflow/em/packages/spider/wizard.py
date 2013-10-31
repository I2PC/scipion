# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
from pyworkflow.em import SetOfImages, SetOfMicrographs, Volume, ProtCTFMicrographs
from protocol_filters import *
from protocol_ca_pca import SpiderProtCAPCA
from protocol_align_apsr import SpiderProtAlignAPSR
from pyworkflow.em.packages.xmipp3.wizard import ListTreeProvider, XmippDownsampleDialog, XmippParticleMaskRadiusWizard, XmippRadiiWizard, XmippMaskPreviewDialog
import xmipp
import pyworkflow.gui.dialog as dialog
from pyworkflow.gui.widgets import LabelSlider
from pyworkflow import findResource

class SpiderProtMaskWizard(XmippParticleMaskRadiusWizard):
    
    _targets = [(SpiderProtCAPCA, ['maskRadius'])]
         
    
class SpiderProtMaskRadiiWizard(XmippRadiiWizard):
    
    _targets = [(SpiderProtAlignAPSR, ['innerRadius', 'outerRadius'])]    

    def _getText(self, obj):
        index = obj.getIndex()
        text = os.path.basename(obj.getFileName())
        if index:
            return "%03d@%s" % (index, text)
        return text
        
    def _getProvider(self, protocol):
        """ This should be implemented to return the list
        of object to be displayed in the tree.
        """
        provider = None
        if protocol.inputParticles.hasValue():
            particles = [] 
            for i, par in enumerate(protocol.inputParticles.get()):
                particles.append(par)
                if i == 100: # Limit the maximum number of particles to display
                    break
            provider = ListTreeProvider(particles)
            provider.getText = self._getText
            
        return provider

    @classmethod    
    def getView(self):
        return "wiz_particle_mask_radii"      
    
#class SpiderGaussianFilterWizard(Wizard):
#    _targets = [(SpiderProtFilter, ['filterRadius'])]
#    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
#
#    def _getText(self, obj):
#        index = obj.getIndex()
#        text = os.path.basename(obj.getFileName())
#        if index:
#            return "%03d@%s" % (index, text)
#        return text
#        
#    def _getProvider(self, protocol):
#        """ This should be implemented to return the list
#        of object to be displayed in the tree.
#        """
#        provider = None
#        if protocol.inputParticles.hasValue():
#            particles = [] 
#            for i, par in enumerate(protocol.inputParticles.get()):
#                particles.append(par)
#                if i == 100: # Limit the maximum number of particles to display
#                    break
#            provider = ListTreeProvider(particles)
#            provider.getText = self._getText
#            
#        return provider
#            
#    def show(self, form):
#        protocol = form.protocol
#        if protocol.inputParticles.hasValue():
#            pars = [par for par in protocol.inputParticles.get()]
#            d = SpiderGaussianFilterDialog(form.root, ListTreeProvider(pars), radius=protocol.filterRadius.get())
#            if d.resultYes():
#                form.setVar('filterRadius', d.getRadius())
#        else:
#            dialog.showWarning("Input particles", "Select particles first", form.root)
#    
#    @classmethod
#    def getView(self):
#        return "wiz_filter"
#
#
#class SpiderFermiFilterWizard(SpiderGaussianFilterWizard):
#    
#    _targets = [(SpiderProtFilter, ['lowFreq', 'highFreq'])]
#    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
#    
#    def show(self, form):
#        protocol = form.protocol
#        provider = self._getProvider(protocol)
#
#        if provider is not None:
#            d = SpiderFermiFilterDialog(form.root, provider, 
#                                          lowFreq=protocol.lowFreq.get(), highFreq=protocol.highFreq.get())
#            if d.resultYes():
#                form.setVar('lowFreq', d.getLowFreq())
#                form.setVar('highFreq', d.getHighFreq())
#        else:
#            dialog.showWarning("Input particles", "Select particles first", form.root)  
#    
#    @classmethod    
#    def getView(self):
#        return "wiz_bandpass"   
#
#
##--------------- Dialogs used by Wizards --------------------------        
#       
#        
#class SpiderGaussianFilterDialog(XmippDownsampleDialog):
#    
#    def _beforePreview(self):
#        XmippDownsampleDialog._beforePreview(self)
#        self.lastObj = None
#        self.rightPreviewLabel = "Filtered particle"
#        self.message = "Filtering particle..."
#        self.previewLabel = "Particle"
#        self.rightImage = xmipp.Image()
#        self.dim = 256           
#        self.dim_par = self.firstItem.getDim()[0]
#        
#    def _createControls(self, frame):
#        self.radiusSlider = LabelSlider(frame, 'Radius', from_=0, to=int(self.dim_par/2), value=self.radius, step=1, callback=lambda a, b, c:self.updateFilteredImage())
#        self.radiusSlider.grid(row=0, column=0, padx=5, pady=5) 
#        
#    def getRadius(self):
#        return int(self.radius.get())
#
#    def updateFilteredImage(self):
#        self.rightPreview.updateData(self.rightImage.getData())
#        
#    def _computeRightPreview(self):
#        """ This function should compute the right preview
#        using the self.lastObj that was selected
#        """
#        #xmipp.fastEstimateEnhancedPSD(self.rightImage, self.lastObj.getFileName(), self.getDownsample(), self.dim, 2)
#        #Do the spider filtering running spider command and set the self.rightImage
#        pass
#        
#
#class SpiderFermiFilterDialog(XmippDownsampleDialog):
#    
#    def _beforePreview(self):
#        XmippDownsampleDialog._beforePreview(self)
#        self.lastObj = None
#        self.rightPreviewLabel = "Filtered"
#        self.message = "Computing filtered image..."
#        self.previewLabel = "Image"
#        self.rightImage = xmipp.Image()
#
#    def _createControls(self, frame):
#        self.freqFrame = ttk.LabelFrame(frame, text="Frequencies", padding="5 5 5 5")
#        self.freqFrame.grid(row=0, column=0)
#        self.lfSlider = self.addFreqSlider('Low freq', self.lowFreq, col=0)
#        self.hfSlider = self.addFreqSlider('High freq', self.highFreq, col=1)
#
#    def addFreqSlider(self, label, value, col):
#        slider = LabelSlider(self.freqFrame, label, from_=0, to=0.5, value=value, callback=lambda a, b, c:self.updateFilteredImage())
#        slider.grid(row=0, column=col, padx=5, pady=5)
#        return slider
#    
#    def getLowFreq(self):
#        return self.lfSlider.get()
#        
#    def getHighFreq(self):
#        return self.hfSlider.get()
#
#    def updateFilteredImage(self):
#        self.rightPreview.updateData(self.rightImage.getData())
#                
#    def _computeRightPreview(self):
#        """ This function should compute the right preview
#        using the self.lastObj that was selected
#        """
#        #xmipp.bandPassFilter(self.rightImage, "%03d@%s" % (self.lastObj.getIndex(), self.lastObj.getFileName()), self.getLowFreq(), self.getHighFreq(), self.getFreqDecay(), self.dim)
#        #Do the spider filtering running spider command and set the self.rightImage
#        pass
#  


    