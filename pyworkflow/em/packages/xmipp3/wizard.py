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
from pyworkflow.em.viewer import Viewer, Wizard
from pyworkflow.em import SetOfImages, SetOfMicrographs, DefCTFMicrographs
from protocol_projmatch import XmippDefProjMatch
import pyworkflow.gui.dialog as dialog
from pyworkflow.gui.tree import BoundTree, TreeProvider
import xmipp


class XmippWizardDownsample(Wizard):
    _targets = [(DefCTFMicrographs, ['lowRes', 'highRes'])]
        
    def show(self, form):
        protocol = form.protocol
        if not protocol.inputMicrographs.hasValue():
            dialog.showWarning("Input micrographs", "Select some micrographs first", form.root)
        else:
            dialog.showWarning("Micrographs", "OK", form.root)    
    
class ListTreeProvider(TreeProvider):
    """ Simple list tree provider. """
    def __init__(self, objList=None):
        self.objList = objList
        self.getColumns = lambda: [('Object', 300)]
        self.getObjects = lambda: self.objList
    
    def getObjectInfo(self, obj):
        info = {'key': obj.getObjId(), 'text': os.path.basename(obj.getFileName()), 'values': ()}
            
        return info
    
class XmippWizardCTF(Wizard):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj
    """
    _targets = [(DefCTFMicrographs, ['lowRes', 'highRes'])]
        
    def show(self, form):
        protocol = form.protocol
        
        if protocol.inputMicrographs.hasValue():
            mics = [mic for mic in protocol.inputMicrographs.get()]
            #XmippDownsampleDialog(form.root, ListTreeProvider(mics))
            XmippCTFDialog(form.root, ListTreeProvider(mics))
            #dialog.showWarning("Input micrographs", "Select some micrographs first", form.root)
        else:
            dialog.showWarning("Micrographs", "OK", form.root)
            
            
class XmippWizardMaskRadius(Wizard):

    _targets = [(XmippDefProjMatch, ['maskRadius'])]
        
    def show(self, form):
        protocol = form.protocol
        dialog.showWarning("Mask radius", "Not yet implemented the wizard to select mask radius", form.root)
        
        
class XmippWizardRadii(Wizard):
    
    _targets = [(XmippDefProjMatch, ['innerRadius', 'outerRadius'])]
    
    def show(self, form):
        protocol = form.protocol
        dialog.showWarning("Correlation radii", "Not yet implemented the wizard to select radii", form.root)    
    
    

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
        self.provider = provider
        dialog.Dialog.__init__(self, parent, "Title", **args)

    def body(self, bodyFrame):
        # Create items frame
        itemsFrame = tk.Frame(bodyFrame, bg='white')
        itemsFrame.grid(row=0, column=0, padx=5, pady=5, rowspan=2)
        itemsFrame.columnconfigure(0, weight=1)
        itemsFrame.rowconfigure(0, weight=1)
        itemsTree = BoundTree(itemsFrame, self.provider)
        itemsTree.grid(row=0, column=0, padx=5, pady=5, sticky='news')
        itemsTree.itemClick = self._itemSelected
        
        # Create preview frame
        previewFrame = tk.Frame(bodyFrame)
        previewFrame.grid(row=0, column=1, padx=5, pady=5)
        self._createPreview(previewFrame)
        
        # Create controls frame
        controlsFrame = tk.Frame(bodyFrame)
        controlsFrame.grid(row=1, column=1, padx=5, pady=5)
        self._createControls(controlsFrame)
    
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
        
    

class XmippDownsampleDialog(XmippPreviewDialog):
    
    def _createPreview(self, frame):
        """ Should be implemented by subclasses to 
        create the items preview. 
        """
        from protlib_gui_figure import ImagePreview
        self.dim = 256
        label = 'kk'
        self.preview = ImagePreview(frame, self.dim, label=label)
        
    def _itemSelected(self, obj):
        filename = obj.getFileName()
        
        self.image = xmipp.Image()        
        self.image.readPreview(filename, self.dim)
        if filename.endswith('.psd'):
            self.image.convertPSD()
        self.Z = self.image.getData()
        self.preview.updateData(self.Z)
        
        
class XmippCTFDialog(XmippDownsampleDialog):
    
    def _createPreview(self, frame):
        """ Should be implemented by subclasses to 
        create the items preview. 
        """
        leftFrame = tk.Frame(frame)
        leftFrame.grid(row=0, column=0)
        
        rightFrame = tk.Frame(frame)
        rightFrame.grid(row=0, column=1)
        
        XmippDownsampleDialog._createPreview(self, leftFrame)
        from protlib_gui_figure import ImagePreview
        self.rightPreview = ImagePreview(rightFrame, self.dim, label='pp') 
        
    def _itemSelected(self, obj):
        XmippDownsampleDialog._itemSelected(self, obj)
        self.rightPreview.updateData(self.Z)        
        
    
