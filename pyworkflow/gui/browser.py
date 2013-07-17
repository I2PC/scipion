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
In this module a simple ObjectBrowser is implemented.
This class can be subclasses to extend its functionality.
A concrete use of ObjectBrowser is FileBrowser, where the
elements to inspect and preview are files.
"""

import os
        
import Tkinter as tk
import ttk

import gui
from tree import BoundTree
from text import TaggedText


class ObjectBrowser(tk.Frame):
    """ This class will implement a simple object browser.
    Basically, it will display a list of elements at the left
    panel and can display a preview and description on the
    right panel for the selected element.
    An ObjectView will be used to grab information for
    each element such as: icon, preview and description.
    A TreeProvider will be used to populate the list (Tree).
    """
    def __init__(self, parent, treeProvider, 
                 showPreview=True, **args):
        tk.Frame.__init__(self, parent, **args)
        self.treeProvider = treeProvider
        gui.configureWeigths(self)
        # The main layout will be two panes, 
        # At the left containing the elements list
        # and the right containing the preview and description
        p = tk.PanedWindow(self, orient=tk.HORIZONTAL)
        p.grid(row=0, column=0, sticky='news')
        
        leftPanel = tk.Frame(p)
        gui.configureWeigths(leftPanel)
        self._fillLeftPanel(leftPanel)
        p.add(leftPanel, padx=5, pady=5)
        p.paneconfig(leftPanel, minsize=300)
        
        if showPreview:
            rightPanel = tk.Frame(p, bg='blue')
            gui.configureWeigths(rightPanel)
            self._fillRightPanel(rightPanel)
            p.add(rightPanel, padx=5, pady=5)    
            p.paneconfig(rightPanel, minsize=200)    
        
            # Register a callback when the item is clicked
            self.tree.itemClick = self._itemClicked
        
    def _fillLeftPanel(self, frame):
        self.tree = BoundTree(frame, self.treeProvider)
        self.tree.grid(row=0, column=0, sticky='news')
        self.itemConfig = self.tree.itemConfig
        self.getImage = self.tree.getImage
    
    def _fillRightPanel(self, frame):
        top = tk.Frame(frame)
        top.grid(row=0, column=0, sticky='news')
        gui.configureWeigths(top)
        top.rowconfigure(0, minsize=200)
        self._fillRightTop(top)
        
        bottom = tk.Frame(frame)
        bottom.grid(row=1, column=0, sticky='news')
        gui.configureWeigths(bottom)
        bottom.rowconfigure(1, weight=1)
        self._fillRightBottom(bottom)
        
    def _fillRightTop(self, top):
        self.noImage = self.getImage('no-image128.png')
        self.label = tk.Label(top, image=self.noImage)
        self.label.grid(row=0, column=0, sticky='news')
        
    def _fillRightBottom(self, bottom):
        self.text = TaggedText(bottom, width=40, height=15, bg='white')
        self.text.grid(row=0, column=0, sticky='news')
        
    def _itemClicked(self, obj):
        img, desc = self.treeProvider.getObjectPreview(obj)
        self.text.clear()
        img = self.getImage(img)
        if img is None:
            img = self.noImage
        self.label.config(image=img)        
        if desc is not None:
            self.text.addText(desc)
    
    
if __name__ == '__main__':
    showBrowseDialog(path='.', main=True, seltype="none")
