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
from pyworkflow.utils.path import removeExt
"""
This module implement some wizards
"""

import os
import Tkinter as tk
import ttk

from pyworkflow.em.constants import UNIT_PIXEL
from pyworkflow.em.convert import ImageHandler
from pyworkflow.em.wizard import (EmWizard, ParticleMaskRadiusWizard, ParticlesMaskRadiiWizard,
                                  FilterParticlesWizard, DownsampleDialog, ImagePreviewDialog,
                                  ListTreeProvider)
import pyworkflow.gui.dialog as dialog
from pyworkflow.gui.widgets import LabelSlider, HotButton

from spider import SpiderShell, runScript
from constants import FILTER_GAUSSIAN, FILTER_FERMI
from convert import locationToSpider
from protocol import (SpiderProtCAPCA, SpiderProtAlignAPSR, SpiderProtAlignPairwise, 
                      SpiderProtFilter, SpiderProtCustomMask)


#===============================================================================
# MASKS
#===============================================================================

class SpiderProtMaskWizard(ParticleMaskRadiusWizard):
    _targets = [(SpiderProtCAPCA, ['radius'])]
    
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
        
    
class SpiderParticlesMaskRadiiWizard(ParticlesMaskRadiiWizard):
    _targets = [(SpiderProtAlignAPSR, ['innerRadius', 'outerRadius']),
                (SpiderProtAlignPairwise, ['innerRadius', 'outerRadius'])]        
    
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
    

#===============================================================================
# FILTERS
#===============================================================================


class SpiderFilterParticlesWizard(FilterParticlesWizard):    
    _targets = [(SpiderProtFilter, ['filterRadius', 'lowFreq', 'highFreq', 'temperature'])]
    
    def _getParameters(self, protocol):
        protParams = {}
        protParams['input']= protocol.inputParticles
        protParams['label']= ["lowFreq", "highFreq", "temperature"]
        protParams['value']= [protocol.getAttributeValue(a) for a in protParams['label']]
        
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
    
#===============================================================================
# UTILS
#===============================================================================
    
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
    
    
#-------------- Custom mask Wizard -----------------------------

CUSTOMMASK_VARS = ['filterRadius1', 'sdFactor', 'filterRadius2', 'maskThreshold']
# Map between the index in the result stack
# and the label to the given image
MASKRESULT_LABELS = ['input image',
                     'filtered',
                     'thresholded',
                     'filtered mask',
                     'final mask',
                     'mask * image',
                     'inverted mask',
                     'inverted \nmask * image',
                     ]

class SpiderCustomMaskWizard(EmWizard):    
    _targets = [(SpiderProtCustomMask, CUSTOMMASK_VARS)]
    
    def _getParameters(self, protocol):
        protParams = {}
        protParams['input']= protocol.inputImage
        protParams['label']= CUSTOMMASK_VARS
        protParams['value']= [protocol.getAttributeValue(a) for a in protParams['label']]
        
        return protParams
    
    def _getProvider(self, protocol):
        return ListTreeProvider([self._getParameters(protocol)['input'].get()])
    
    def show(self, form):
        protocol = form.protocol
        provider = self._getProvider(protocol)

        if protocol.inputImage.get():
            d = CustomMaskDialog(form.root, provider, protocolParent=protocol)
            if d.resultYes():
                for varName in CUSTOMMASK_VARS:
                    form.setVar(varName, d.getVarValue(varName))
        else:
            dialog.showWarning("Input error", "Select the input image first", form.root)  


class CustomMaskDialog(ImagePreviewDialog):
        
    def _beforePreview(self):
        imgLocation = self.protocolParent.inputImage.get().getLocation()
        self.dim = ImageHandler().getDimensions(imgLocation)[0]
        self.lastObj = None
        self.rightPreviewLabel = "Final mask"
        self.message = "Generating mask..."
        self.rightImage = ImageHandler()._img
        
    def _createPreview(self, frame):
        """ Should be implemented by subclasses to 
        create the items preview. 
        """
        self._previews = []
        for i, label in enumerate(MASKRESULT_LABELS):
            self.previewLabel = label
            previewFrame = tk.Frame(frame)
            ImagePreviewDialog._createPreview(self, previewFrame)
            self._previews.append(self.preview) # store all previews created
            previewFrame.grid(row=i/4, column=i%4)
            
    def _itemSelected(self, obj):
        self.lastObj = obj
        dialog.FlashMessage(self, self.message, func=self._computeRightPreview)
           
    def _createVarWidgets(self, parent, varName, row, col):
        var = tk.StringVar()
        self._vars[varName] = var
        var.set(self.protocolParent.getAttributeValue(varName))
        varLabel = tk.Label(parent, text=varName)
        varLabel.grid(row=row, column=col*2, padx=5, pady=5)
        varEntry = tk.Entry(parent, width=10, textvariable=var)
        varEntry.grid(row=row, column=col*2+1, padx=5, pady=5)
    
    def _createControls(self, frame):
        self._vars = {}
        inputFrame = tk.Frame(frame)
        inputFrame.grid(row=0, column=0)
        
        for i, varName in enumerate(CUSTOMMASK_VARS):
            self._createVarWidgets(inputFrame, varName, i%2, i/2)
            
        previewBtn = HotButton(frame, text='Preview', command=self._computeRightPreview)
        previewBtn.grid(row=1, column=1, padx=5, pady=5)
            
    def getVarValue(self, varName):
        return self._vars[varName].get()
    
    def _computeRightPreview(self, e=None):
        """ This function should compute the right preview
        using the self.lastObj that was selected
        """
        
        #TODO, enter in project/Tmp folter
        params = {'[filter-radius1]': self.getVarValue('filterRadius1'),
                  '[sd-factor]': self.getVarValue('sdFactor'),
                  '[filter-radius2]': self.getVarValue('filterRadius2'),
                  '[mask-threshold2]': self.getVarValue('maskThreshold'),
                  '[input_image]': removeExt(self.lastObj.getFileName()),
                  '[output_mask]': 'stkmask',
                  }
        ext = self.protocolParent.getExt()
        runScript('mda/custommask.msa', ext, params)
        
        for i, preview in enumerate(self._previews):
            if i == 0:
                self.rightImage.read(self.lastObj.getFileName())
            else:
                self.rightImage.read('%d@stkmask.%s' % (i, ext))
            preview.updateData(self.rightImage.getData())
        
    
