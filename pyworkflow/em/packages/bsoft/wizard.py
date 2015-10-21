# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Jose Gutierrez (jose.gutierrez@cnb.csic.es)
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

from pyworkflow.utils.path import cleanPath
import pyworkflow.gui.dialog as dialog
from pyworkflow.em.convert import ImageHandler
from pyworkflow.em.wizard import DownsampleDialog, ImagePreviewDialog, FilterParticlesWizard
    
from protocol_bfilter import BsoftProtBfilter



class BsoftFilterParticlesWizard(FilterParticlesWizard):    
    _targets = [(BsoftProtBfilter, ['filterType'])]
    
    def _getParameters(self, protocol):
        
        label, value = self._getInputProtocol(self._targets, protocol)
        
        protParams = {}
        protParams['input']= protocol.inputParticles
        protParams['label']= label
        protParams['value']= value

        return protParams
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']
        return FilterParticlesWizard._getListProvider(self, _objs)
    
    def show(self, form):
        protocol = form.protocol
        provider = self._getProvider(protocol)

        if provider is not None:
            d = BsoftFilterDialog(form.root, provider, protocolParent=protocol)
        else:
            dialog.showWarning("Input particles", "Select particles first", form.root)  
    
    
#--------------- Dialogs used by Wizards --------------------------        
       
class BsoftFilterDialog(DownsampleDialog):
    
    def _beforePreview(self):
        ImagePreviewDialog._beforePreview(self)
        self.lastObj = None
        self.rightPreviewLabel = "Filtered particle"
        self.message = "Filtering particle..."
        self.previewLabel = "Particle"
        self.rightImage = ImageHandler()._img
        
    def _createControls(self, frame):
        pass #FIXME
#         self.freqFrame = ttk.LabelFrame(frame, text="Frequencies", padding="5 5 5 5")
#         self.freqFrame.grid(row=0, column=0)
#         if self.protocolParent.filterType <= FILTER_SPACE_REAL:
#             self.radiusSlider = self.addFreqSlider('Radius', self.protocolParent.filterRadius.get(), col=0)
#         else:
#             self.lfSlider = self.addFreqSlider('Low freq', self.protocolParent.lowFreq.get(), col=0)
#             self.hfSlider = self.addFreqSlider('High freq', self.protocolParent.highFreq.get(), col=1)        
#             if self.protocolParent.filterType == FILTER_FERMI:
#                 self.tempSlider = self.addFreqSlider('Temperature', self.protocolParent.temperature.get(), col=2)
#         radiusButton = tk.Button(self.freqFrame, text='Preview', command=self._doPreview)
#         radiusButton.grid(row=0, column=3, padx=5, pady=5)
        
    def _doPreview(self, e=None):
        if self.lastObj is None:
            dialog.showError("Empty selection", "Select an item first before preview", self)
        else:
            self._computeRightPreview()
    
    def updateFilteredImage(self):
        self.rightPreview.updateData(self.rightImage.getData())
        
    def _computeRightPreview(self):
        """ This function should compute the right preview
        using the self.lastObj that was selected
        """
        from pyworkflow.em.packages.xmipp3 import locationToXmipp
        
        # Copy image to filter to Tmp project folder
        inputPath = os.path.join("Tmp", "bsoft_filter_input.spi")
        outputPath = os.path.join("Tmp", "bsoft_filter_output.spi")
        cleanPath(inputPath, outputPath)

        ih = ImageHandler()
        ih.convert(self.lastObj.getLocation(), inputPath) 
                
        self.protocolParent.runFilter(inputPath, outputPath)
        
        # Get output image and update filtered image
        img = ih._img
        img.read(outputPath)
        self.rightImage = img
        self.updateFilteredImage()

