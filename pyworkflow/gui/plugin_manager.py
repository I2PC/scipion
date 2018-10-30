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
import tkFont
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
from install.plugin_funcs import PluginRepository, PluginInfo

PLUGIN = 'plugin'
BINARY = 'binary'
UNCHECKED = 'unchecked'
CHECKED = 'checked'
INSTALL = 'install'
UNINSTALL = 'uninstall'
pluginRepo = PluginRepository()
pluginDict = pluginRepo.getPlugins(getPipData=True)


class PluginTreeview(ttk.Treeview):
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
        self.im_install = tk.PhotoImage(file='/home/yunior/Escritorio/icon/fa-install.png')
        self.im_uninstall = tk.PhotoImage(file='/home/yunior/Escritorio/icon/fa-uninstall.png')
        self.tag_configure(UNCHECKED, image=self.im_unchecked)
        self.tag_configure(CHECKED, image=self.im_checked)
        self.tag_configure(INSTALL, image=self.im_install)
        self.tag_configure(UNINSTALL, image=self.im_uninstall)
        self.selectedItem = None

    def insert(self, parent, index, iid=None, **kw):
        """ same method as for standard treeview but add the tag 'unchecked'
            automatically if no tag among ('checked', 'unchecked')
            is given """
        if not "tags" in kw:
            kw["tags"] = (UNCHECKED,)
        elif not (UNCHECKED in kw["tags"] or CHECKED in kw["tags"]):
            kw["tags"] = (UNCHECKED,)
        ttk.Treeview.insert(self, parent, index, iid, **kw)

    def check_item(self, item):
        """ check the box of item and change the state of the boxes of item's
            ancestors accordingly """
        if UNCHECKED in self.item(item, 'tags'):
            self.item(item, tags=(INSTALL,))
        else:
            self.item(item, tags=(CHECKED,))

    def uncheck_item(self, item):
        """ uncheck the boxes of item's descendant """
        if CHECKED in self.item(item, 'tags'):
            self.item(item, tags=(UNINSTALL,))
            children = self.get_children(item)
            for iid in children:
                self.delete(iid)
        else:

            self.item(item, tags=(UNCHECKED,))


class Operation:
    """
    This class contain the operation details
    """
    def __init__(self, objName, objType=PLUGIN, objStatus=INSTALL):
        self.objName = objName
        self.objType = objType
        self.status = objStatus

    def getObjName(self):
        return self.objName

    def getObjType(self):
        return self.objType

    def getObjStatus(self):
        return self.status

    def setObjStatus(self):
        if self.status == INSTALL:
            self.status = UNINSTALL
        else:
            self.status = INSTALL


class OperationList:
    """
    This class contain a plugins/binaries operation list and allow execute it
    """
    def __init__(self):
        self.operationList = []

    def insertOperation(self, operation):
        index = self.operationIndex(operation)
        if index is not None:
            self.removeOperation(index)
        else:
            self.operationList.append(operation)

    def removeOperation(self, index):
        self.operationList.pop(index)

    def operationIndex(self, operation):
        index = None
        for i in range(0, len(self.operationList)):
            if self.operationList[i].getObjName() == operation.getObjName():
                return i
        return index


class PluginBrowser(tk.Frame):
    """ This class will implement a frame.
        It will display a list of plugin at the left
        panel. A TreeProvider will be used to populate the list (Tree).
        """
    def __init__(self, parent,  **args):
        tk.Frame.__init__(self, parent, **args)
        self._lastSelected = None
        self.operationList = OperationList()
        gui.configureWeigths(self)
        # The main layout will be two panes,
        # At the left containing the plugin list
        # and the right containing the description
        mainFrame = tk.PanedWindow(parent, orient=tk.HORIZONTAL)
        mainFrame.grid(row=0, column=0, sticky='news')
        # ---------------------------------------------------------------
        # Left Panel
        leftPanel = tk.Frame(mainFrame)  # Create a left panel to put the tree
        leftPanel.grid(row=0, column=0, padx=0, pady=0)
        self._fillLeftPanel(leftPanel)  # Fill the left panel
        # ---------------------------------------------------------------
        # Right Panel: will be two vertical panes
        # At the Top contain the plugin or binary information
        # At the Bottom contain a system terminal that show the operation steps
        rightPanel = tk.PanedWindow(mainFrame, orient=tk.VERTICAL)
        rightPanel.grid(row=0, column=1, padx=0, pady=0)

        # Top Panel
        topPanel = ttk.Frame(rightPanel)  # Panel to put the plugin information
        topPanel.pack(side=TOP, fill=BOTH, expand=Y)

        self.dataCols = ('Name                      ',
                         'Version            ',
                         'Description               ',
                         'Url                       ',
                         'Author                    ')
        self.topPanelTree = ttk.Treeview(topPanel, columns=self.dataCols,
                                         show='headings')
        self.topPanelTree.grid(row=0, column=0, sticky='news')

        # configure column headings
        for c in self.dataCols:
            self.topPanelTree.heading(c, text=c.title())
            self.topPanelTree.column(c, width=tkFont.Font().measure(c.title()))

        # configure horizontal scroollbar
        xsb = ttk.Scrollbar(topPanel, orient='horizontal',
                                        command=self.topPanelTree.xview)
        xsb.grid(row=1, column=0, sticky='news')
        self.topPanelTree.configure(yscrollcommand=xsb.set)
        xsb.configure(command=self.topPanelTree.xview)
        topPanel.rowconfigure(0, weight=1)
        topPanel.columnconfigure(0, weight=1)

        # Bottom Panel
        # This section show the plugin operation and a console
        bottomPanel = ttk.Frame(rightPanel)
        tabControl = ttk.Notebook(bottomPanel)  # Create Tab Control
        operationTab = ttk.Frame(tabControl)    # Create a operation tab
        consoleTab = ttk.Frame(tabControl)    # Create a console
        tabControl.add(operationTab, text='Operations')  # Add the Operation tab
        tabControl.add(consoleTab, text='Console')
        tabControl.pack(expand=1, fill="both")    # Pack to make visible


        # Add the right panels to Right Panel
        rightPanel.add(topPanel, padx=0, pady=0)
        rightPanel.add(bottomPanel, padx=0, pady=0)

        self._fillTopRightPanel(topPanel)
        self._fillBottomRightPanel(bottomPanel)

        # Add the Plugin list at left
        mainFrame.add(leftPanel, padx=0, pady=0)
        mainFrame.paneconfig(leftPanel, minsize=200)

        # Add the Plugins or Binaries information
        mainFrame.add(rightPanel, padx=0, pady=0)
        mainFrame.paneconfig(rightPanel, minsize=200)

    def _fillTopRightPanel(self, topFrame):
        """
        Fill the top Panel with the plugins list
        """
        pass

    def _fillBottomRightPanel(self, topFrame):
        """
        Fill the top Panel with the plugins list
        """
        pass

    def _fillLeftPanel(self, leftFrame):
        """
        Fill the left Panel with the plugins list
        """
        gui.configureWeigths(leftFrame)
        self.tree = PluginTreeview(leftFrame, show="tree")
        self.tree.grid(row=0, column=0, sticky='news')

        self.yscrollbar = ttk.Scrollbar(leftFrame, orient='vertical',
                                        command=self.tree.yview)
        self.yscrollbar.grid(row=0, column=1, sticky='news')
        self.tree.configure(yscrollcommand=self.yscrollbar.set)
        self.yscrollbar.configure(command=self.tree.yview)

        # check / uncheck boxes(plugin or binary) on right click
        self.tree.bind("<Button-3>", self.box_rightClick, True)
        # show the plugin or binary information on click
        self.tree.bind("<Button-1>", self.objectInformation, True)

        # Load all plugins and fill the tree view
        self.loadPlugins()

    def objectInformation(self, event):
        """Show the plugin or binary information"""
        x, y, widget = event.x, event.y, event.widget
        item = self.tree.selectedItem = self.tree.identify_row(y)
        if self.tree.selectedItem is not None and \
                self.isPlugin(self.tree.item(self.tree.selectedItem,
                                             "values")[0]):
            self.showPluginInformation(item)

    def box_rightClick(self, event):
        """ check or uncheck a plugin or binary box when clicked """
        x, y, widget = event.x, event.y, event.widget
        elem = widget.identify("element", x, y)
        if "image" in elem:
            # a box was clicked
            self.tree.selectedItem = self.tree.identify_row(y)
            tags = self.tree.item(self.tree.selectedItem, "tags")
            type = self.tree.item(self.tree.selectedItem, "value")
            if tags[0] in [UNCHECKED, UNINSTALL]:
                self.tree.check_item(self.tree.selectedItem)
                self.reloadInstalledPlugin(self.tree.selectedItem)
            else:
                self.tree.uncheck_item(self.tree.selectedItem)
            self.operationList.insertOperation(Operation(self.tree.selectedItem,
                                                         type))
            x = 10

    def isPlugin(self, value):
        return value == 'plugin'

    def showPluginInformation(self, pluginName):
        """Shows the information associated with a given plugin"""
        #dataCols = ('Name', 'Version', 'Description', 'Home Page', 'Author')
        plugin = pluginDict.get(pluginName, None)
        if plugin is not None:
            pluginInfo = [(plugin.getPipName(), plugin.getPipVersion(),
                           plugin.getSummary(), plugin.getHomePage(),
                           plugin.getAuthor())]
            self.topPanelTree.delete(*self.topPanelTree.get_children())
            self.topPanelTree.insert('', 'end', values=pluginInfo[0])

            for idx, val in enumerate(pluginInfo):
                iwidth = tkFont.Font().measure(self.dataCols[idx])
                if self.topPanelTree.column(self.dataCols[idx], 'width') <= iwidth:
                    self.topPanelTree.column(self.dataCols[idx], width=iwidth)

    def reloadInstalledPlugin(self, pluginName):
        """
        Reload a given installed plugin and update a tree view
        """
        plugin = pluginDict.get(pluginName, None)
        if plugin is not None:
            # Insert all binaries of plugin on the tree
            pluginBinaryList = plugin.getInstallenv()
            if pluginBinaryList is not None:
                binaryList = pluginBinaryList.getPackages()
                keys = sorted(binaryList.keys())
                for k in keys:
                    pVersions = binaryList[k]
                    for binary, version in pVersions:
                        installed = pluginBinaryList._isInstalled(binary,
                                                                  version)
                        tag = UNCHECKED
                        if installed:
                            tag = CHECKED
                        binaryName = str(binary + '-' + version)
                        self.tree.insert(pluginName, "end", binaryName,
                                         text=binaryName, tags=tag,
                                         values='binary')

    def loadPlugins(self):
        """
        Load all plugins and fill the tree view widget
        """
        for pluginObj in pluginDict:
            plugin = pluginDict.get(pluginObj, None)
            if plugin is not None:
                tag = UNCHECKED
                if plugin.isInstalled():
                    # Insert the plugin name in the tree
                    tag = CHECKED
                    self.tree.insert("", 0, pluginObj, text=pluginObj, tags=tag,
                                     values='plugin')
                    # Insert all binaries of plugin on the tree
                    pluginBinaryList = plugin.getInstallenv()
                    if pluginBinaryList is not None:
                        binaryList = pluginBinaryList.getPackages()
                        keys = sorted(binaryList.keys())
                        for k in keys:
                            pVersions = binaryList[k]
                            for binary, version in pVersions:
                                installed = pluginBinaryList._isInstalled(binary,
                                                                     version)
                                tag = UNCHECKED
                                if installed:
                                    tag = CHECKED
                                binaryName = str(binary + '-' + version)
                                self.tree.insert(pluginObj, "end", binaryName,
                                                 text=binaryName, tags=tag,
                                                 values='binary')
                else:
                    self.tree.insert("", 0, pluginObj, text=pluginObj, tags=tag,
                                     values='plugin')


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
        import time
        time.sleep(5)
        browser = PluginBrowser(self.root, **kwargs)
