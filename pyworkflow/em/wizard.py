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
from pyworkflow.em.constants import *
from pyworkflow.em.convert import ImageHandler

from pyworkflow.wizard import Wizard, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em import SetOfImages, SetOfMicrographs, Volume, ProtCTFMicrographs

import pyworkflow.gui.dialog as dialog
from pyworkflow.gui.widgets import LabelSlider
from pyworkflow.gui.tree import BoundTree, TreeProvider
from pyworkflow import findResource

import xmipp

class downsampleWizard(Wizard):
        
    def show(self, form):
        protocol = form.protocol
        if protocol.inputMicrographs.hasValue():
            mics = [mic.clone() for mic in protocol.inputMicrographs.get()]
            d = downsampleDialog(form.root, ListTreeProvider(mics), downsample=protocol.downFactor.get())
            if d.resultYes():
                form.setVar('downFactor', d.getDownsample())
        else:
            dialog.showWarning("Input micrographs", "Select micrographs first", form.root)
    
    @classmethod
    def getView(self):
        return "wiz_downsampling"
        
    
class ListTreeProvider(TreeProvider):
    """ Simple list tree provider. """
    def __init__(self, objList=None):
        self.objList = objList
        self.getColumns = lambda: [('Object', 150)]
        self.getObjects = lambda: self.objList
    
    def getObjectInfo(self, obj):
        info = {'key': obj.getObjId(), 'text': self.getText(obj), 'values': ()}
        return info
    
    def getText(self, obj):
        """ Get the text to display for an object. """
        return os.path.basename(obj.getFileName())
        
    
class ctfWizard(Wizard):
        
    def show(self, form):
        protocol = form.protocol
        
        form.setParamFromVar('inputMicrographs') # update selected input micrographs
        
        if protocol.inputMicrographs.hasValue():
            mics = [mic.clone() for mic in protocol.inputMicrographs.get()]
            #XmippDownsampleDialog(form.root, ListTreeProvider(mics))
            d = ctfDialog(form.root, ListTreeProvider(mics), 
                               lf=protocol.lowRes.get(), hf=protocol.highRes.get())
            if d.resultYes():
                form.setVar('lowRes', d.getLowFreq())
                form.setVar('highRes', d.getHighFreq())
            #dialog.showWarning("Input micrographs", "Select some micrographs first", form.root)
        else:
            dialog.showWarning("Input micrographs", "Select micrographs first", form.root)
    
    @classmethod    
    def getView(self):
        return "wiz_ctf"            
            
class maskRadiusWizard(Wizard):
    
    def _getProvider(self, protocol):
        """ This should be implemented to return the list
        """
        pass
        
    def show(self, form):
        protocol = form.protocol
        provider = self._getProvider(protocol)

        if provider is not None:
            d = maskPreviewDialog(form.root, 
                                       provider, 
                                       maskRadius=protocol.maskRadius.get(), 
                                       unit=UNIT_ANGSTROM)
            if d.resultYes():
                form.setVar('maskRadius', d.getRadius())
        else:
            dialog.showWarning("Empty input", "Select elements first", form.root)
                
                
class particleMaskRadiusWizard(maskRadiusWizard):
        
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
                particles.append(par.clone())
                if i == 100: # Limit the maximum number of particles to display
                    break
            provider = ListTreeProvider(particles)
            provider.getText = self._getText
            
        return provider
    
    @classmethod    
    def getView(self):
        return "wiz_particle_mask"       
    
class volumeMaskRadiusWizard(maskRadiusWizard):
        
    def _getProvider(self, protocol):
        """ This should be implemented to return the list
        of object to be displayed in the tree.
        """
        if protocol.input3DReferences.hasValue():
            vols = [vol.clone() for vol in protocol.input3DReferences.get()]
            return ListTreeProvider(vols)
        return None
    
    @classmethod    
    def getView(self):
        return "wiz_volume_mask"       
        
        
class radiiWizard(volumeMaskRadiusWizard):
    
    def show(self, form):
       
        protocol = form.protocol
        provider = self._getProvider(protocol)

        if provider is not None:
            d = maskRadiiPreviewDialog(form.root, 
                                            provider, 
                                            innerRadius=protocol.innerRadius.get(), 
                                            outerRadius=protocol.outerRadius.get(),
                                            unit=UNIT_ANGSTROM)
            if d.resultYes():
                form.setVar('innerRadius', d.getRadius(d.radiusSliderIn))
                form.setVar('outerRadius', d.getRadius(d.radiusSliderOut))
        else:
            dialog.showWarning("Empty input", "Select elements first", form.root)    
    
    @classmethod    
    def getView(self):
        return "wiz_volume_mask_radii"   


class filterParticlesWizard(Wizard):
    
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
                particles.append(par.clone())
                if i == 100: # Limit the maximum number of particles to display
                    break
            provider = ListTreeProvider(particles)
            provider.getText = self._getText
            
        return provider
    
class filterVolumesWizard(Wizard):
    
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
        if protocol.input3DReferences.hasValue():
            vols = [vol.clone() for vol in protocol.input3DReferences.get()]
            return ListTreeProvider(vols)
        return None


class bandpassParticleWizard(filterParticlesWizard):
    
    @classmethod    
    def getView(self):
        return "wiz_particle_bandpass"   
    
class bandpassVolumesWizard(filterVolumesWizard):
    
    @classmethod    
    def getView(self):
        return "wiz_volume_bandpass"  

    
class gaussianWizard(filterParticlesWizard):
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    @classmethod    
    def getView(self):
        return "wiz_gaussian"   
    

#--------------- Dialogs used by Wizards --------------------------
    
class previewDialog(dialog.Dialog):
    """ This will be the base class for several wizards.
    The layout of this wizard will be:
    1. Left panel(Items) that contains a list of items to preview
    2. Right-top panel (Preview) where some preview of the items will be displayed
    3. Right-bottom panel (Controls) where some controls can change the preview
    """
    def __init__(self, parent, provider, **args):
        """ 
        Params:
            parent: parent windows of the dialog.
            provider: the TreeProvider to populate items tree.
        """
        # Set the attributes in **args
        for k, v in args.iteritems():
            setattr(self, k, v)
            
        self.provider = provider
        self.firstItem = provider.getObjects()[0]
        buttons = [('Select', dialog.RESULT_YES), ('Cancel', dialog.RESULT_CANCEL)]
        dialog.Dialog.__init__(self, parent, "Wizard", buttons=buttons, default='Select', **args)

    def body(self, bodyFrame):
        bodyFrame.config()
        bodyFrame.columnconfigure(0, weight=1)
        bodyFrame.rowconfigure(0, weight=1)
        bodyFrame.columnconfigure(1, weight=1)
        # Create items frame
        itemsFrame = tk.Frame(bodyFrame, bg='white')
        itemsFrame.grid(row=0, column=0, padx=5, pady=5, sticky='news')
        itemsFrame.columnconfigure(0, weight=1)
        itemsFrame.rowconfigure(0, weight=1)
        itemsTree = BoundTree(itemsFrame, self.provider)
        itemsTree.grid(row=0, column=0, padx=5, pady=5, sticky='news')
        itemsTree.itemClick = self._itemSelected
        
        # Create preview frame
        previewFrame = tk.Frame(bodyFrame)
        previewFrame.grid(row=0, column=1, padx=5, pady=5)
        self._beforePreview()
        self._createPreview(previewFrame)
        
        # Create controls frame
        controlsFrame = tk.Frame(bodyFrame)
        controlsFrame.grid(row=1, column=1, padx=5, pady=5, sticky='news')
        self._createControls(controlsFrame)
        self._itemSelected(self.firstItem)
        itemsTree.selectItem(0) # Select the first item
    
    def _beforePreview(self):
        """ Called just before setting the preview.
        This is the place to set data values such as: labels, constants...
        """
        pass
    
    def _createPreview(self, frame):
        """ Should be implemented by subclasses to 
        create the items preview. 
        """
        pass
    
    def _createControls(self, frame):
        """ Create controls to be used. """
        pass
    
    def _itemSelected(self, obj):
        """ This will be call when an item in the tree is selected. """
        pass
        
    

class imagePreviewDialog(previewDialog):
    
    def _beforePreview(self):
        self.dim = 256
        self.previewLabel = ''
    
    def _createPreview(self, frame):
        """ Should be implemented by subclasses to 
        create the items preview. 
        """
        from pyworkflow.gui.matplotlib_image import ImagePreview
        self.preview = ImagePreview(frame, self.dim, label=self.previewLabel)
        self.preview.grid(row=0, column=0) 
        
    def _itemSelected(self, obj):
        
        index = obj.getIndex()
        filename = obj.getFileName()
        if index:
            filename = "%03d@%s" % (index, filename)
        
#        self.image = xmipp.Image()
        self.image = ImageHandler()._img
        
        try:
            self.image.readPreview(filename, self.dim)
            if filename.endswith('.psd'):
                self.image.convertPSD()
            self.Z = self.image.getData()
        except Exception, e:
            from pyworkflow.gui.matplotlib_image import getPngData
            self.Z = getPngData(findResource('no-image.png'))
            dialog.showError("Input particles", "Error reading image <%s>" % filename, self) 
        self.preview.updateData(self.Z)
       
        
class downsampleDialog(imagePreviewDialog):
    
    def _beforePreview(self):
        imagePreviewDialog._beforePreview(self)
        self.lastObj = None
        self.rightPreviewLabel = "PSD"
        self.message = "Computing PSD..."
        self.previewLabel = "Micrograph"
        self.rightImage = ImageHandler()._img
        
    def _createPreview(self, frame):
        """ Should be implemented by subclasses to 
        create the items preview. 
        """
        leftFrame = tk.Frame(frame)
        leftFrame.grid(row=0, column=0)
        
        rightFrame = tk.Frame(frame)
        rightFrame.grid(row=0, column=1)
        
        imagePreviewDialog._createPreview(self, leftFrame)
        self.rightPreview = self._createRightPreview(rightFrame)
        self.rightPreview.grid(row=0, column=0) 
                
    def _createRightPreview(self, rightFrame):
        from pyworkflow.gui.matplotlib_image import ImagePreview        
        return ImagePreview(rightFrame, self.dim, label=self.rightPreviewLabel)
        
    def _createControls(self, frame):
        self.downVar = tk.StringVar()
        self.downVar.set(getattr(self, 'downsample', 1))
        downFrame = tk.Frame(frame)
        downFrame.grid(row=0, column=0, sticky='nw')
        downLabel = tk.Label(downFrame, text='Downsample')
        downLabel.grid(row=0, column=0, padx=5, pady=5)
        downEntry = tk.Entry(downFrame, width=10, textvariable=self.downVar)
        downEntry.grid(row=0, column=1, padx=5, pady=5)
        downButton = tk.Button(downFrame, text='Preview', command=self._doPreview)
        downButton.grid(row=0, column=2, padx=5, pady=5)
        
    def getDownsample(self):
        return float(self.downVar.get())

    def _itemSelected(self, obj):
        self.lastObj = obj
        imagePreviewDialog._itemSelected(self, obj)
        
        dialog.FlashMessage(self, self.message, func=self._computeRightPreview)
        #self._computeRightPreview(obj)
        self.rightPreview.updateData(self.rightImage.getData())
        
    def _doPreview(self, e=None):
        if self.lastObj is None:
            dialog.showError("Empty selection", "Select an item first before preview", self)
        else:
            self._itemSelected(self.lastObj)
        
    def _computeRightPreview(self):
        """ This function should compute the right preview
        using the self.lastObj that was selected
        """
        xmipp.fastEstimateEnhancedPSD(self.rightImage, self.lastObj.getFileName(), self.getDownsample(), self.dim, 2)
        

class ctfDialog(downsampleDialog):
    
    def _createRightPreview(self, rightFrame):
        from pyworkflow.gui.matplotlib_image import PsdPreview  
        return  PsdPreview(rightFrame, 1.2*self.dim, self.lf, self.hf, label=self.rightPreviewLabel)

    def _createControls(self, frame):
        self.freqFrame = ttk.LabelFrame(frame, text="Frequencies", padding="5 5 5 5")
        self.freqFrame.grid(row=0, column=0)
        self.lfSlider = self.addFreqSlider('Low freq', self.lf, col=0)
        self.hfSlider = self.addFreqSlider('High freq', self.hf, col=1)
    
    def getDownsample(self):
        return 1.0 # Micrograph previously downsample, not taken into account here

    def addFreqSlider(self, label, value, col):
        slider = LabelSlider(self.freqFrame, label, from_=0, to=0.5, value=value, callback=lambda a, b, c:self.updateFreqRing())
        slider.grid(row=0, column=col, padx=5, pady=5)
        return slider
        
    def updateFreqRing(self):
        self.rightPreview.updateFreq(self.getLowFreq(), self.getHighFreq())
    
    def getLowFreq(self):
        return self.lfSlider.get()
        
    def getHighFreq(self):
        return self.hfSlider.get()

class bandPassFilterDialog(downsampleDialog):
    
    def _beforePreview(self):
        imagePreviewDialog._beforePreview(self)
        self.lastObj = None
        self.rightPreviewLabel = "Filtered"
        self.message = "Computing filtered image..."
        self.previewLabel = "Image"
        self.rightImage = ImageHandler()._img

    def _createControls(self, frame):
        self.freqFrame = ttk.LabelFrame(frame, text="Frequencies ("+self.unit+")", padding="5 5 5 5")
        self.freqFrame.grid(row=0, column=0)
        
        self.showLowFreq = getattr(self, 'showLowFreq', True)
        self.showHighFreq = getattr(self, 'showHighFreq', True)
        self.showDecay = getattr(self, 'showDecay', True)
            
        if (not self.showLowFreq) or (not self.showHighFreq):
            label_high = 'Freq'
            label_low = 'Freq'
        else:
            label_high = 'High freq'
            label_low = 'Low freq'
        
        self.samplingRate = 1.0
        self.sliTo = 0.5
        self.sliFrom = 0.
        if self.unit == UNIT_ANGSTROM:
            self.samplingRate = self.firstItem.getSamplingRate()
            self.itemDim,_,_,_ = self.firstItem.getDim()
            self.sliFrom = 2.*self.samplingRate
#            print "Sampling rate !!", self.samplingRate
#            print "from !!", self.sliFrom
            self.sliTo = 2.*self.itemDim*self.samplingRate
#            print " to !!", self.sliTo
            
        self.step = self.sliTo/1000
        
        if self.showLowFreq:
            self.lfSlider = self.addFreqSlider(label_low, self.lowFreq, col=0)
        if self.showHighFreq:
            self.hfSlider = self.addFreqSlider(label_high, self.highFreq, col=1)
        if self.showDecay:
            self.freqDecaySlider = self.addFreqSlider('Decay', self.freqDecay, col=2)
                                   
    def addFreqSlider(self, label, value, col):
        slider = LabelSlider(self.freqFrame, label, from_=self.sliFrom, to=self.sliTo, step=self.step , value=value, callback=lambda a, b, c:self.updateFilteredImage())
        slider.grid(row=0, column=col, padx=5, pady=5)
        return slider

    def updateFilteredImage(self):
        self._computeRightPreview()
        self.rightPreview.updateData(self.rightImage.getData())
        
    def _computeRightPreview(self):
        """ This function should compute the right preview
        using the self.lastObj that was selected
        """
        
#        print "LOW", self.getLowFreq()
#        print "HIGH", self.getHighFreq()
#        print "DECAY", self.getFreqDecay()
        
        xmipp.bandPassFilter(self.rightImage, "%03d@%s" % (self.lastObj.getIndex(), self.lastObj.getFileName()), self.getLowFreq(), self.getHighFreq(), self.getFreqDecay(), self.dim)

    def getLowFreq(self):
        if self.showLowFreq:
            if self.unit == UNIT_ANGSTROM:
                return 1/self.lfSlider.get()*self.itemDim
            else:    
                return self.lfSlider.get()
        return 0.
        
    def getHighFreq(self):
        if self.showHighFreq:
            if self.unit == UNIT_ANGSTROM:
                return 1/self.hfSlider.get()*self.itemDim
            else:    
                return self.hfSlider.get()
        return 1.0
        
    def getFreqDecay(self):
        if self.showDecay:
            if self.unit == UNIT_ANGSTROM:
                return 1/self.freqDecaySlider.get()*self.itemDim
            else:  
                return self.freqDecaySlider.get()
        return 0.    
    
class gaussianFilterDialog(bandPassFilterDialog):
    
    def _createControls(self, frame):
        self.freqVar = tk.StringVar()
        self.freqVar.set(getattr(self, 'freqSigma', 1))
        freqFrame = tk.Frame(frame)
        freqFrame.grid(row=0, column=0, sticky='nw')
        downEntry = tk.Entry(freqFrame, width=10, textvariable=self.freqVar)
        downEntry.grid(row=0, column=1, padx=5, pady=5)
        downButton = tk.Button(freqFrame, text='Preview', command=self._doPreview)
        downButton.grid(row=0, column=2, padx=5, pady=5)    
        
    def getFreqSigma(self):
        return float(self.freqVar.get())

    def _computeRightPreview(self):
        """ This function should compute the right preview
        using the self.lastObj that was selected
        """
        xmipp.gaussianFilter(self.rightImage, "%03d@%s" % (self.lastObj.getIndex(), self.lastObj.getFileName()), self.getFreqSigma(), self.dim)


class maskPreviewDialog(imagePreviewDialog):
    
    def _beforePreview(self):
        self.dim = 256
        self.samplingRate = 1
        self.unit = getattr(self, 'unit', UNIT_PIXEL)
        if self.unit != UNIT_PIXEL:
            self.samplingRate = self.firstItem.getSamplingRate()
            
        self.dim_par = self.firstItem.getDim()[0] * self.samplingRate
        
        self.ratio = self.dim / float(self.dim_par)
        
        self.previewLabel = 'Central slice'
    
    def _createPreview(self, frame):
        """ Should be implemented by subclasses to 
        create the items preview. 
        """
        from pyworkflow.gui.matplotlib_image import MaskPreview    
        if self.maskRadius == -1:
            self.iniRadius = self.dim_par/2 
        else:
            self.iniRadius = self.maskRadius
        self.preview = MaskPreview(frame, self.dim, label=self.previewLabel, outerRadius=self.iniRadius*self.ratio)
        self.preview.grid(row=0, column=0) 
    
    def _createControls(self, frame):
        self.addRadiusBox(frame) 
        
    def addRadiusBox(self, parent):
        self.radiusSlider = LabelSlider(parent, 'Outer radius ('+ self.unit +')', from_=0, to=int(self.dim_par/2), value=self.iniRadius, step=self.samplingRate, callback=lambda a, b, c:self.updateRadius())
        self.radiusSlider.grid(row=0, column=0, padx=5, pady=5) 
    
    def updateRadius(self):
        self.preview.updateMask(self.radiusSlider.get() * self.ratio)     
        
    def getRadius(self):
        return int(self.radiusSlider.get())
    

class maskRadiiPreviewDialog(maskPreviewDialog):

    def _createPreview(self, frame):
        """ Should be implemented by subclasses to 
        create the items preview. 
        """
        from pyworkflow.gui.matplotlib_image import MaskPreview    
        self.preview = MaskPreview(frame, self.dim, label=self.previewLabel, outerRadius=int(self.outerRadius)*self.ratio, innerRadius=0)
        self.preview.grid(row=0, column=0) 
    
    def _createControls(self, frame):
        self.radiusSliderOut = LabelSlider(frame, 'Outer radius', from_=0, to=int(self.dim_par/2), value=self.outerRadius, step=1, callback=lambda a, b, c:self.updateRadius(self.radiusSliderOut, self.radiusSliderIn))
        self.radiusSliderOut.grid(row=0, column=0, padx=5, pady=5) 

        self.radiusSliderIn = LabelSlider(frame, 'Inner radius', from_=0, to=int(self.dim_par/2), value=self.innerRadius, step=1, callback=lambda a, b, c:self.updateRadius(self.radiusSliderOut, self.radiusSliderIn))
        self.radiusSliderIn.grid(row=1, column=0, padx=5, pady=5) 
    
    def updateRadius(self, radiusSliderOut, radiusSliderIn):
        self.preview.updateMask(outerRadius = radiusSliderOut.get() * self.ratio, innerRadius = radiusSliderIn.get() * self.ratio)     
        
    def getRadius(self, radiusSlider):
        return int(radiusSlider.get())
    