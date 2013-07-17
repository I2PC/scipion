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
from pyworkflow.em.viewer import Viewer, Wizard
from pyworkflow.em import SetOfImages, SetOfMicrographs, DefCTFMicrographs
import pyworkflow.gui.dialog as dialog
import xmipp

class XmippDialogPreview(dialog.Dialog):
    """ This will be the base class for several wizards.
    The layout of this wizard will be:
    1. Left panel(Items) that contains a list of items to preview
    2. Right-top panel (Preview) where some preview of the items will be displayed
    3. Right-bottom panel (Controls) where some controls can change the preview
    """
    def __init__(self, parent, **args):
        Viewer.__init__(self)
        dialog.Dialog.__init__(self, parent, **args) 

    def body(self, bodyFrame):
        #bodyFrame.config(bg='white')
        frame = tk.Frame(bodyFrame, bg='white')
        frame.grid(row=0, column=0, padx=20, pady=20)
        label = tk.Label(bodyFrame, text=self.entryLabel, bg='white', bd=0)
        label.grid(row=0, column=0, sticky='nw', padx=(15, 10), pady=15)
        self.entry = tk.Entry(bodyFrame, bg=gui.cfgEntryBgColor, width=self.entryWidth)
        self.entry.grid(row=0, column=1, sticky='new', padx=(0,15), pady=15)
        self.initial_focus = self.entry
        
    def apply(self):
        self.value = self.entry.get()
        
    def validate(self):
        if len(self.entry.get().strip()) == 0:
            showError("Validation error", "Value is empty", self)
            return False
        return True
    
class XmippWizardCTF(Wizard):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj
    """
    _targets = [(DefCTFMicrographs, ['lowRes', 'highRes'])]
        
    def show(self, form):
        protocol = form.protocol
        if not protocol.inputMicrographs.hasValue():
            dialog.showWarning("Input micrographs", "Select some micrographs first", form.root)
        else:
            dialog.showWarning("Micrographs", "OK", form.root)
