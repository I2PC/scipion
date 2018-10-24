# **************************************************************************
# *
# * Authors:    Yunior C. Fonseca Reyna (cfonseca@cnb.csic.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os.path
import stat

from Tkinter import *
import Tkinter as tk
import ttk
import Tix

import gui
from tree import BoundTree, TreeProvider, Tree
from text import TaggedText, openTextFileEditor
from widgets import Button, HotButton
from pyworkflow.config import MenuConfig
from pyworkflow.utils.properties import Icon
from pyworkflow.gui.form import *


class CheckboxTreeview(ttk.Treeview):
    """
        Treeview widget with checkboxes left of each item.
        The checkboxes are done via the image attribute of the item, so to keep
        the checkbox, you cannot add an image to the item.
    """

    def __init__(self, master=None, **kw):
        ttk.Treeview.__init__(self, master, **kw)
        # checkboxes are implemented with pictures
        self.im_checked = tk.PhotoImage(file='/home/yunior/Escritorio/icon/fa-checked.png')
        self.im_unchecked = tk.PhotoImage(file='/home/yunior/Escritorio/icon/fa-unchecked.png')
        self.im_tristate = tk.PhotoImage(file='/home/yunior/Escritorio/icon/fa-checkmark.png')
        self.tag_configure("unchecked", image=self.im_unchecked)
        self.tag_configure("tristate", image=self.im_tristate)
        self.tag_configure("checked", image=self.im_checked)
        # check / uncheck boxes on click
        self.bind("<Button-3>", self.box_click, True)


    def insert(self, parent, index, iid=None, **kw):
        """ same method as for standard treeview but add the tag 'unchecked'
            automatically if no tag among ('checked', 'unchecked', 'tristate')
            is given """
        if not "tags" in kw:
            kw["tags"] = ("unchecked",)
        elif not ("unchecked" in kw["tags"] or "checked" in kw["tags"]
                  or "tristate" in kw["tags"]):
            kw["tags"] = ("unchecked",)
        ttk.Treeview.insert(self, parent, index, iid, **kw)

    def check_descendant(self, item):
        """ check the boxes of item's descendants """
        children = self.get_children(item)
        for iid in children:
            self.item(iid, tags=("checked",))
            self.check_descendant(iid)

    def check_ancestor(self, item):
        """ check the box of item and change the state of the boxes of item's
            ancestors accordingly """
        self.item(item, tags=("checked",))
        parent = self.parent(item)
        if parent:
            self.item(parent, tags=("checked",))


    def tristate_parent(self, item):
        """ put the box of item in tristate and change the state of the boxes of
            item's ancestors accordingly """
        self.item(item, tags=("tristate",))
        parent = self.parent(item)
        if parent:
            self.tristate_parent(parent)

    def uncheck_descendant(self, item):
        """ uncheck the boxes of item's descendant """
        children = self.get_children(item)
        for iid in children:
            self.item(iid, tags=("unchecked",))
            self.uncheck_descendant(iid)
        self.item(item, tags=("unchecked",))

    def uncheck_ancestor(self, item):
        """ uncheck the box of item and change the state of the boxes of item's
            ancestors accordingly """
        self.item(item, tags=("unchecked",))
        parent = self.parent(item)
        if parent:
            children = self.get_children(parent)
            b = ["unchecked" in self.item(c, "tags") for c in children]
            if False in b:
                # at least one box is checked and item's box is unchecked
                self.tristate_parent(parent)
            else:
                # no box is checked
                self.uncheck_ancestor(parent)

    def box_click(self, event):
        """ check or uncheck box when clicked """
        x, y, widget = event.x, event.y, event.widget
        elem = widget.identify("element", x, y)
        if "image" in elem:
            # a box was clicked
            item = self.identify_row(y)
            tags = self.item(item, "tags")
            if "unchecked" in tags:
                self.check_descendant(item)
                self.check_ancestor(item)
            else:
                if "checked" in tags:
                    self.uncheck_descendant(item)

class PluginBrowser(tk.Frame):
    """ This class will implement a frame.
        It will display a list of plugin at the left
        panel. A TreeProvider will be used to populate the list (Tree).
        """
    def __init__(self, parent,  **args):
        tk.Frame.__init__(self, parent, **args)
        self._lastSelected = None
        gui.configureWeigths(self)
        # The main layout will be two panes,
        # At the left containing the plugin list
        # and the right containing the description
        mainFrame = tk.PanedWindow(parent, orient=tk.HORIZONTAL)
        mainFrame.grid(row=0, column=0, sticky='news')

        leftPanel = tk.Frame(mainFrame)
        leftPanel.grid(row=0, column=0, padx=0, pady=0)

        self._fillLeftPanel(leftPanel)


        # Add the Plugin list at left
        mainFrame.add(leftPanel, padx=0, pady=0)
        mainFrame.paneconfig(leftPanel, minsize=200)


    def _fillLeftPanel(self, frame):
        gui.configureWeigths(frame)
        self.tree = CheckboxTreeview(frame, show="tree")
        self.tree.grid(row=0, column=0, sticky='news')

        self.yscrollbar = ttk.Scrollbar(frame, orient='vertical',
                                        command=self.tree.yview)
        self.yscrollbar.grid(row=0, column=1, sticky='news')
        self.tree.configure(yscrollcommand=self.yscrollbar.set)
        self.yscrollbar.configure(command=self.tree.yview)

        self.tree.insert("", 0, "Appion", text="Appion")
        self.tree.insert("Appion", "end", "dogpicker", text="dogpicker")

        self.tree.insert("", 0, "Xmipp", text="Xmipp")
        self.tree.insert("Xmipp", "end", "XmippBin", text="XmippBin")
        self.tree.insert("Xmipp", "end", "XmippSource", text="XmippSource")




class PluginManagerWindow(gui.Window):
    """ Windows to hold a plugin manager frame inside. """

    def __init__(self, title, master=None, **kwargs):
        if 'minsize' not in kwargs:
            kwargs['minsize'] = (300, 300)
        gui.Window.__init__(self, title, master, **kwargs)

        menu = MenuConfig()

        fileMenu = menu.addSubMenu('File')
        fileMenu.addSubMenu('Browse files', 'browse', icon='fa-folder-open.png')
        fileMenu.addSubMenu('Exit', 'exit', icon='fa-sign-out.png')

        self.menuCfg = menu
        gui.Window.createMainMenu(self, self.menuCfg)

    def onExit(self):
        self.close()


class PluginManager(PluginManagerWindow):
    """ Windows to hold a frame inside. """

    def __init__(self, title, master=None, path=None,
                 onSelect=None, shortCuts=None, **kwargs):
        PluginManagerWindow.__init__(self, title, master, **kwargs)

        browser = PluginBrowser(self.root, **kwargs)
