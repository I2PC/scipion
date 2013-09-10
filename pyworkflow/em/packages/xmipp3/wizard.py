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
from pyworkflow.em import SetOfImages, SetOfMicrographs, Volume, DefCTFMicrographs
from protocol_projmatch import XmippDefProjMatch, XmippProtProjMatch 
from protocol_preprocess_micrographs import XmippDefPreprocessMicrograph
from protocol_filters import XmippDefMask, XmippProtMask
import pyworkflow.gui.dialog as dialog
from pyworkflow.gui.widgets import LabelSlider
from pyworkflow.gui.tree import BoundTree, TreeProvider
import xmipp


class XmippDownsampleWizard(Wizard):
    _targets = [(XmippDefPreprocessMicrograph, ['downFactor'])]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
        
    def show(self, form):
        protocol = form.protocol
        if protocol.inputMicrographs.hasValue():
            mics = [mic for mic in protocol.inputMicrographs.get()]
            d = XmippDownsampleDialog(form.root, ListTreeProvider(mics), downsample=protocol.downFactor.get())
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
        
    
class XmippCTFWizard(Wizard):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj
    """
    _targets = [(DefCTFMicrographs, ['lowRes', 'highRes'])]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
        
    def show(self, form):
        protocol = form.protocol
        
        form.setParamFromVar('inputMicrographs') # update selected input micrographs
        
        if protocol.inputMicrographs.hasValue():
            mics = [mic for mic in protocol.inputMicrographs.get()]
            #XmippDownsampleDialog(form.root, ListTreeProvider(mics))
            d = XmippCTFDialog(form.root, ListTreeProvider(mics), 
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
            
class XmippMaskRadiusWizard(Wizard):
    
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
        
    def _getProvider(self, protocol):
        """ This should be implemented to return the list
        of object to be displayed in the tree.
        """
        pass
        
    def show(self, form):
        protocol = form.protocol
        provider = self._getProvider(protocol)

        if provider is not None:
            d = XmippMaskPreviewDialog(form.root, provider, maskRadius=protocol.maskRadius.get())
            if d.resultYes():
                form.setVar('maskRadius', d.getRadius())
        else:
            dialog.showWarning("Empty input", "Select elements first", form.root)
                
                
class XmippParticleMaskRadiusWizard(XmippMaskRadiusWizard):

    _targets = [(XmippDefMask, ['maskRadius'])]
        
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
        return "wiz_particle_mask"       
    
class XmippVolumeMaskRadiusWizard(XmippMaskRadiusWizard):

    _targets = [(XmippDefProjMatch, ['maskRadius'])]
        
    def _getProvider(self, protocol):
        """ This should be implemented to return the list
        of object to be displayed in the tree.
        """
        if protocol.input3DReferences.hasValue():
            vols = [vol for vol in protocol.input3DReferences.get()]
            return ListTreeProvider(vols)
        return None
    
    @classmethod    
    def getView(self):
        return "wiz_volume_mask"       
        
        
class XmippRadiiWizard(Wizard):
    
    _targets = [(XmippDefProjMatch, ['innerRadius', 'outerRadius'])]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    def show(self, form):
        protocol = form.protocol
        dialog.showWarning("Correlation radii", "Not yet implemented the wizard to select radii", form.root)    
    
    @classmethod    
    def getView(self):
        return "wiz_volume_mask_radii"   

#--------------- Dialogs used by Wizards --------------------------
    
class XmippPreviewDialog(dialog.Dialog):
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
        
    

class XmippImagePreviewDialog(XmippPreviewDialog):
    
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
        filename = obj.getFileName()
        self.image = xmipp.Image()
        self.image.readPreview(filename, self.dim)
        if filename.endswith('.psd'):
            self.image.convertPSD()
        self.Z = self.image.getData()
        self.preview.updateData(self.Z)
       
        
class XmippDownsampleDialog(XmippImagePreviewDialog):
    
    def _beforePreview(self):
        XmippImagePreviewDialog._beforePreview(self)
        self.lastObj = None
        self.rightPreviewLabel = "PSD"
        self.message = "Computing PSD..."
        self.previewLabel = "Micrograph"
        self.rightImage = xmipp.Image()
        
    def _createPreview(self, frame):
        """ Should be implemented by subclasses to 
        create the items preview. 
        """
        leftFrame = tk.Frame(frame)
        leftFrame.grid(row=0, column=0)
        
        rightFrame = tk.Frame(frame)
        rightFrame.grid(row=0, column=1)
        
        XmippImagePreviewDialog._createPreview(self, leftFrame)
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
        XmippImagePreviewDialog._itemSelected(self, obj)
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
        

class XmippCTFDialog(XmippDownsampleDialog):
    
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

class XmippMaskPreviewDialog(XmippImagePreviewDialog):
    
    def _beforePreview(self):
        self.dim = 256
        self.previewLabel = 'Central slice'
    
    def _createPreview(self, frame):
        """ Should be implemented by subclasses to 
        create the items preview. 
        """
        from pyworkflow.gui.matplotlib_image import MaskPreview    
        self.preview = MaskPreview(frame, self.dim, label=self.previewLabel, outerRadius=self.maskRadius)
        self.preview.grid(row=0, column=0) 
    
    def _createControls(self, frame):
        self.addRadiusBox(frame) 
        
    def addRadiusBox(self, parent):
        print "maskRadius: %s" % self.maskRadius
        self.radiusSlider = LabelSlider(parent, 'Outer radius', from_=0, to=int(self.dim/2), value=self.maskRadius, step=1, callback=lambda a, b, c:self.updateRadius())
        self.radiusSlider.grid(row=0, column=0, padx=5, pady=5) 
    
    def updateRadius(self):
        self.preview.updateMask(self.radiusSlider.get())     
        
    def getRadius(self):
        return int(self.radiusSlider.get())
