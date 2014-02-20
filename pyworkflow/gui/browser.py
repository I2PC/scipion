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
In this module a simple ObjectBrowser is implemented.
This class can be subclasses to extend its functionality.
A concrete use of ObjectBrowser is FileBrowser, where the
elements to inspect and preview are files.
"""

import os
import os.path

import Tkinter as tk
import ttk

from tree import Tree
import gui
import tooltip
import widgets
import matplotlib_image
from pyworkflow.utils.path import findResource, dirname, getHomePath
from tree import BoundTree, FileTreeProvider
from text import TaggedText
from pyworkflow.utils.properties import Message, Icon


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
        gui.configureWeigths(frame)
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
      
# Some constants for the type of selection
# when the file browser is opened
SELECT_NONE = 0 # No selection, just browse files                  
SELECT_FILE = 1
SELECT_FOLDER = 2
SELECT_PATH = 3 # Can be either file or folder


class FileBrowser(ObjectBrowser):
    """ The FileBrowser is a particular class of ObjectBrowser
    where the "objects" are just files and directories.
    """
    def __init__(self, parent, initialDir='.', 
                 selectionType=SELECT_FILE, selectionSingle=True, 
                 allowFilter=True, filterFunction=None, previewDim=144,
                 showHidden=False):
        """ 
        """
        tp = FileTreeProvider(initialDir, showHidden)
        ObjectBrowser.__init__(self, parent, tp)
        
        if selectionType != SELECT_NONE:
            buttonsFrame = tk.Frame(self)
            self._fillButtonsFrame(buttonsFrame)
            buttonsFrame.grid(row=1, column=0)

    def _fillLeftPanel(self, frame):
        """ Redefine this method to include a buttons toolbar and
        also include a filter bar at the bottom of the Tree.
        """
        # Tree with files
        frame.columnconfigure(0, weight=1)
        
        treeFrame = tk.Frame(frame)
        ObjectBrowser._fillLeftPanel(self, treeFrame)
        # Register the double-click event
        self.tree.itemDoubleClick = self._itemDoubleClick
        treeFrame.grid(row=1, column=0, sticky='news')
        # Toolbar frame
        toolbarFrame = tk.Frame(frame)
        self._fillToolbar(toolbarFrame)
        toolbarFrame.grid(row=0, column=0, sticky='new')
        # Filter frame
        tk.Label(frame, text="Filter").grid(row=2, column=0)
        
        frame.rowconfigure(1, weight=1)

    def _fillToolbar(self, frame):
        """ Fill the toolbar frame with some buttons. """
        self._col = 0
        
        def addButton(text, image, command):
            btn = tk.Label(frame, text=text, image=self.getImage(image), 
                       compound=tk.LEFT, cursor='hand2')
            btn.bind('<Button-1>', command)
            btn.grid(row=0, column=self._col, sticky='nw',
                     padx=(0, 5), pady=5)
            self._col += 1
            
        addButton('Refresh', Icon.ACTION_REFRESH, self._actionRefresh)
        addButton('Home', Icon.HOME, self._actionHome)
        addButton('Back', Icon.ARROW_LEFT, self._actionUp)
        addButton('Up', Icon.ARROW_UP, self._actionUp)
        
    def _fillButtonsFrame(self, frame):
        """ Add button to the bottom frame if the selectMode
        is distinct from SELECT_NONE.
        """
        tk.Button(frame, text="Select", image=self.getImage(Icon.BUTTON_SELECT),
                        compound=tk.LEFT).grid(row=0, column=0, padx=(0,5))
        tk.Button(frame, text="Close", image=self.getImage(Icon.BUTTON_CLOSE),
                        compound=tk.LEFT).grid(row=0, column=1)                        
                
    def _actionRefresh(self, e=None):
        self.tree.update()
        
    def _goDir(self, newDir):
        self.treeProvider.setDir(newDir)
        self.tree.update()
        
    def _actionUp(self, e=None):
        self._goDir(dirname(self.treeProvider.getDir()))
        
    def _actionHome(self, e=None):
        self._goDir(getHomePath())
        
    def _itemDoubleClick(self, obj):
        if obj.isDir():
            self._goDir(obj.getPath())
