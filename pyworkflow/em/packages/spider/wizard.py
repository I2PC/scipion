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
import Tkinter as tk
import ttk
from pyworkflow.em import SetOfImages, SetOfMicrographs, Volume, ProtCTFMicrographs

from constants import *
from pyworkflow.em.wizard import *
from pyworkflow.em.convert import ImageHandler

from protocol_filters import *
from protocol_ca_pca import SpiderProtCAPCA
from protocol_align_apsr import SpiderProtAlignAPSR

import pyworkflow.gui.dialog as dialog
from pyworkflow.gui.widgets import LabelSlider
from spider import SpiderShell
from convert import locationToSpider


class SpiderProtMaskWizard(ParticleMaskRadiusWizard):
    _targets = [(SpiderProtCAPCA, ['radius'])]
    
    def _getParameters(self, protocol):
        protParams = {}
        protParams['input']= protocol.inputParticles
        protParams['label']= "radius"
        protParams['value']= protocol.maskRadius.get()
        return protParams
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input'] 
        return ParticleMaskRadiusWizard._getListProvider(self, _objs)
    
    def show(self, form):
        params = self._getParameters(form.protocol)
        _value = params['value']
        _label = params['label']
        ParticleMaskRadiusWizard.show(self, form, _value, _label, UNIT_PIXEL)
        

    
class SpiderParticlesMaskRadiiWizard(ParticlesMaskRadiiWizard):
    _targets = [(SpiderProtAlignAPSR, ['innerRadius', 'outerRadius'])]        
    
    def _getParameters(self, protocol):
        protParams = {}
        protParams['input']= protocol.inputParticles
        protParams['label']= ["innerRadius", "outerRadius"]
        protParams['value']= [protocol.innerRadius.get(), protocol.outerRadius.get()]
        return protParams
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']
        return ParticlesMaskRadiiWizard._getListProvider(self, _objs)
    
    def show(self, form):
        params = self._getParameters(form.protocol)
        _value = params['value']
        _label = params['label']
        ParticlesMaskRadiiWizard.show(self, form, _value, _label, UNIT_PIXEL)
    



class SpiderFilterParticlesWizard(FilterParticlesWizard):    
    _targets = [(SpiderProtFilter, ['filterRadius', 'lowFreq', 'highFreq', 'temperature'])]
    
    def _getParameters(self, protocol):
        protParams = {}
        protParams['input']= protocol.inputParticles
        protParams['label']= ["lowFreq", "highFreq", "temperature"]
        protParams['value']= [protocol.innerRadius.get(), protocol.outerRadius.get()]
        return protParams
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']
        return FilterParticlesWizard._getListProvider(self, _objs)
        
    
    def show(self, form):
        protocol = form.protocol
        provider = self._getProvider(protocol)

        if provider is not None:
            d = SpiderFilterDialog(form.root, provider, 
                                          protocolParent=protocol)
            if d.resultYes():
                if protocol.filterType <= FILTER_GAUSSIAN:
                    form.setVar('filterRadius', d.getRadius())
                else:
                    form.setVar('lowFreq', d.getLowFreq())
                    form.setVar('highFreq', d.getHighFreq())
                    if protocol.filterType == FILTER_FERMI:
                        form.setVar('temperature', d.getTemperature())
        else:
            dialog.showWarning("Input particles", "Select particles first", form.root)  
    
    



#--------------- Dialogs used by Wizards --------------------------        
       
#class SpiderGaussianFilterDialog(XmippDownsampleDialog):
class SpiderFilterDialog(DownsampleDialog):
    
    def _beforePreview(self):
        ImagePreviewDialog._beforePreview(self)
        self.lastObj = None
        self.rightPreviewLabel = "Filtered particle"
        self.message = "Filtering particle..."
        self.previewLabel = "Particle"
        self.rightImage = ImageHandler()._img
        
    def _createControls(self, frame):
        self.freqFrame = ttk.LabelFrame(frame, text="Frequencies", padding="5 5 5 5")
        self.freqFrame.grid(row=0, column=0)
        if self.protocolParent.filterType <= FILTER_GAUSSIAN:
            self.radiusSlider = self.addFreqSlider('Radius', self.protocolParent.filterRadius.get(), col=0)
        else:
            self.lfSlider = self.addFreqSlider('Low freq', self.protocolParent.lowFreq.get(), col=0)
            self.hfSlider = self.addFreqSlider('High freq', self.protocolParent.highFreq.get(), col=1)        
            if self.protocolParent.filterType == FILTER_FERMI:
                self.tempSlider = self.addFreqSlider('Temperature', self.protocolParent.temperature.get(), col=2)
        radiusButton = tk.Button(self.freqFrame, text='Preview', command=self._doPreview)
        radiusButton.grid(row=0, column=3, padx=5, pady=5)
        
    def _doPreview(self, e=None):
        if self.lastObj is None:
            dialog.showError("Empty selection", "Select an item first before preview", self)
        else:
            self._computeRightPreview()
            
    def getRadius(self):
        return self.radiusSlider.get()
    
    def addFreqSlider(self, label, value, col):
        slider = LabelSlider(self.freqFrame, label, from_=0, to=0.5, value=value, callback=None)
        slider.grid(row=0, column=col, padx=5, pady=5)
        return slider
    
    def getLowFreq(self):
        return self.lfSlider.get()
        
    def getHighFreq(self):
        return self.hfSlider.get()

    def getTemperature(self):
        return self.tempSlider.get()
    
    def updateFilteredImage(self):
        self.rightPreview.updateData(self.rightImage.getData())
        
    def _computeRightPreview(self):
        """ This function should compute the right preview
        using the self.lastObj that was selected
        """
        from pyworkflow.em.packages.xmipp3 import locationToXmipp
        
        # Copy image to filter to Tmp project folder
        outputName = os.path.join("Tmp", "filtered_particle")
        outputPath = outputName + ".spi"

        outputLoc = (1, outputPath)
        ih = ImageHandler()
        ih.convert(self.lastObj.getLocation(), outputLoc) 
                
        outputLocSpiStr = locationToSpider(1, outputName)
        
        pars = {}
        pars["filterType"] = self.protocolParent.filterType.get()
        pars["filterMode"] = self.protocolParent.filterMode.get()
        pars["usePadding"] = self.protocolParent.usePadding.get()
        pars["op"] = "FQ"
        
        if self.protocolParent.filterType <= FILTER_GAUSSIAN:
            pars['filterRadius'] = self.getRadius()
        else:
            pars['lowFreq'] = self.getLowFreq()
            pars['highFreq'] = self.getHighFreq()
            
        if self.protocolParent.filterType == FILTER_FERMI:
            pars['temperature'] = self.getTemperature()

        filter_spider(outputLocSpiStr, outputLocSpiStr, **pars)
        
        # Get output image and update filtered image
        img = ImageHandler()._img
        locXmippStr = locationToXmipp(1, outputPath)
        img.read(locXmippStr)
        self.rightImage = img
        self.updateFilteredImage()

#TODO: Refactor this function to be used also by method filterParticles
def filter_spider(inputLocStr, outputLocStr, **pars):
    """ Function to filter an image located on inputLocStr and
    write it to outputLocStr. """
     
    spi = SpiderShell(ext='spi') # Create the Spider process to send commands         
    filterNumber = pars["filterType"] * 2 + 1
    # Consider low-pass or high-pass
    filterNumber += pars["filterMode"]
    OP = pars["op"]
    if not pars["usePadding"]:
        OP += ' NP'
        
    args = []
    
    if pars["filterType"] <= FILTER_GAUSSIAN:
        args.append(pars['filterRadius'])
    else:
        args.append('%f %f' % (pars['lowFreq'], pars['highFreq']))
        
    if pars["filterType"] == FILTER_FERMI:
        args.append(pars['temperature'])
        
    spi.runFunction(OP, inputLocStr, outputLocStr, filterNumber, *args)
    spi.close()
    