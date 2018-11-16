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

from Tkinter import *
import webbrowser
import threading
from pyworkflow.config import MenuConfig
from pyworkflow.utils.log import ScipionLogger
from pyworkflow.gui.text import TextFileViewer
from pyworkflow.gui import *
from install.plugin_funcs import PluginRepository, PluginInfo

from pyworkflow.utils.properties import *
from pyworkflow.utils import redStr
PLUGIN_LOG_NAME = 'Plugin.log'
PLUGIN_ERRORS_LOG_NAME = 'Plugin.err'

pluginRepo = PluginRepository()
pluginDict = None


class PluginTree(ttk.Treeview):
    """
        Treeview widget with checkboxes left of each item.
        The checkboxes are done via the image attribute of the item, so to keep
        the checkbox, you cannot add an image to the item. """
    def __init__(self, master=None, **kw):
        ttk.Treeview.__init__(self, master, **kw)

        self.im_checked = gui.getImage(Icon.CHECKED)
        self.im_unchecked = gui.getImage(Icon.UNCHECKED)
        self.im_install = gui.getImage(Icon.INSTALL)
        self.im_uninstall = gui.getImage(Icon.UNINSTALL)
        self.im_to_install = gui.getImage(Icon.TO_INSTALL)
        self.im_installed = gui.getImage(Icon.INSTALLED)
        self.im_processing = gui.getImage(Icon.PROCESSING)
        self.im_failure = gui.getImage(Icon.FAILURE)
        self.tag_configure(PluginStates.UNCHECKED, image=self.im_unchecked)
        self.tag_configure(PluginStates.CHECKED, image=self.im_checked)
        self.tag_configure(PluginStates.INSTALL, image=self.im_install)
        self.tag_configure(PluginStates.UNINSTALL, image=self.im_uninstall)
        self.tag_configure(PluginStates.TO_INSTALL, image=self.im_to_install)
        self.tag_configure(PluginStates.INSTALLED, image=self.im_installed)
        self.tag_configure(PluginStates.PRECESSING, image=self.im_processing)
        self.tag_configure(PluginStates.FAILURE, image=self.im_failure)
        self.selectedItem = None

    def insert(self, parent, index, iid=None, **kw):
        """ same method as for standard treeview but add the tag 'unchecked'
            automatically if no tag among ('checked', 'unchecked')
            is given """
        if not "tags" in kw:
            kw["tags"] = (PluginStates.UNCHECKED,)
        elif not (PluginStates.UNCHECKED in kw["tags"] or
                  PluginStates.CHECKED in kw["tags"] or
                  PluginStates.TO_INSTALL in kw["tags"] or
                  PluginStates.INSTALLED in kw["tags"]):
            kw["tags"] = (PluginStates.UNCHECKED,)
        ttk.Treeview.insert(self, parent, index, iid, **kw)

    def check_item(self, item):
        """ check the box of item and change the state of the boxes of item's
            ancestors accordingly """
        if PluginStates.UNCHECKED in self.item(item, 'tags'):
            self.item(item, tags=(PluginStates.INSTALL,))
        else:
            self.item(item, tags=(PluginStates.CHECKED,))

    def uncheck_item(self, item):
        """ uncheck the boxes of item's descendant """
        if PluginStates.CHECKED in self.item(item, 'tags'):
            self.item(item, tags=(PluginStates.UNINSTALL,))
            children = self.get_children(item)
            for iid in children:
                self.delete(iid)
        else:
            self.item(item, tags=(PluginStates.UNCHECKED,))

    def processing_item(self, item):
        """change the box item to processing item"""
        self.item(item, tags=(PluginStates.PRECESSING,))

    def installed_item(self, item):
        """change the box item to processing item"""
        self.item(item, tags=(PluginStates.INSTALLED,))

    def failure_item(self, item):
        """change the box item to failure item"""
        self.item(item, tags=(PluginStates.FAILURE,))

    def disable(self):
        self.state(('disabled',))

    def enable(self):
        self.state(('!disabled',))

    def is_disabled(self):
        return 'disabled' in self.state()

    def is_enabled(self):
        return not self.is_disabled()


class Operation:
    """
    This class contain the object(plugin/binary) operation details
    """
    def __init__(self, objName, objType=PluginStates.PLUGIN,
                 objStatus=PluginStates.INSTALL, objParent=None):
        self.objName = objName
        self.objType = objType
        self.objStatus = objStatus
        self.objParent = objParent

    def getObjName(self):
        """
        Get the object(plugin/binary) name
        """
        return self.objName

    def getObjType(self):
        """
        Get the object type (plugin or binary)
        """
        return self.objType

    def getObjStatus(self):
        """
        Get the object status (installed, uninstalled, to install,...)
        """
        return self.objStatus

    def setObjStatus(self, status):
        """
        Set the object status
        """
        self.objStatus = status

    def getObjParent(self):
        """
        Get the object parent in the tree. If the object is a bynary, this
        method return None
        """
        return self.objParent

    def runOperation(self):
        """
        This method install or uninstall a plugin/binary operation
        """
        if self.objType == PluginStates.PLUGIN:
            if self.objStatus == PluginStates.INSTALL:
                plugin = pluginDict.get(self.objName, None)
                if plugin is not None:
                    installed = plugin.installPipModule()
                    if installed:
                        plugin.installBin()
            else:
                plugin = PluginInfo(self.objName, self.objName, remote=False)
                if plugin is not None:
                    plugin.uninstallBins()
                    plugin.uninstallPip()
        else:
            plugin = PluginInfo(self.objParent, self.objParent, remote=False)
            if self.objStatus == PluginStates.INSTALL:
                if plugin is not None:
                    plugin.installBin([self.objName])
            else:
                plugin.uninstallBins([self.objName])


class OperationList:
    """
    This class contain a plugins/binaries operations list and allow execute it
    """
    def __init__(self):
        self.operationList = []

    def insertOperation(self, operation):
        """
        This method insert into the list a given operation. If the operation was
        inserted previously is eliminated
        """
        index = self.operationIndex(operation)
        if index is not None:
            self.removeOperation(index)
        else:
            tag = PluginStates.UNINSTALL
            if operation.getObjStatus() == PluginStates.UNCHECKED:
                tag = PluginStates.INSTALL
            operation.setObjStatus(tag)
            self.operationList.append(operation)

    def removeOperation(self, index):
        """
        Remove an operation from the list based on an index
        """
        self.operationList.pop(index)

    def operationIndex(self, operation):
        """
        Returns the index of an operation within the list
        """
        index = None
        for i in range(0, len(self.operationList)):
            if self.operationList[i].getObjName() == operation.getObjName():
                return i
        return index

    def getOperations(self, op):
        """
        Return the operation List. If the operation is not None return a list
        with only the operation op
        """
        if op is None:
            return self.operationList
        else:
            return [self.getOperationByName(op.getObjName())]

    def getOperationByName(self, opName):
        """
        Return an operation that match with a given name
        """
        operation = [op for op in self.operationList
                     if op.getObjName() == opName]
        if len(operation):
            return operation[0]
        return None

    def applyOperations(self):
        """
        Execute a operation list
        """
        for op in self.operationList:
            op.runOperation()

    def clearOperations(self):
        """
        Clear the operation List
        """
        del self.operationList[:]


class PluginBrowser(tk.Frame):
    """ This class will implement a frame.
        It will display a list of plugin at the left
        panel. A TreeProvider will be used to populate the list (Tree).
        At the right panel provide a plugin/binary information(top panel) and
        a list of operation (bottom panel)
        """
    def __init__(self, master,  **args):
        tk.Frame.__init__(self, master, **args)
        self._lastSelected = None
        self.operationList = OperationList()
        gui.configureWeigths(self)

        # Creating the layout where all application elements will be placed
        parentFrame = tk.Frame(master)
        parentFrame.grid(row=0, column=0, sticky='news')
        gui.configureWeigths(parentFrame, 1)

        # The main layout will be two panes:
        # At the left containing the plugin list
        # and the right containing a description and the operation list
        mainFrame = tk.PanedWindow(parentFrame, orient=tk.HORIZONTAL)
        mainFrame.grid(row=1, column=0, sticky='news')

        # ---------------------------------------------------------------
        # Left Panel
        leftPanel = tk.Frame(mainFrame)  # Create a left panel to put the tree
        leftPanel.grid(row=0, column=0, padx=0, pady=0, sticky='news')
        self._fillLeftPanel(leftPanel)  # Fill the left panel

        # ---------------------------------------------------------------
        # Right Panel: will be two vertical panes
        # At the Top contain the plugin or binary information
        # At the Bottom contain a tab widget with an operation list and
        # a system terminal that show the operation steps
        rightPanel = tk.PanedWindow(mainFrame, orient=tk.VERTICAL)
        rightPanel.grid(row=0, column=1, padx=0, pady=0, sticky='news')

        # Top Panel
        # Panel to put the plugin information
        topPanel = ttk.Frame(rightPanel)
        topPanel.pack(side=TOP, fill=BOTH, expand=Y)
        topPanel.configure(cursor='hand1')
        self._createRightTopPanel(topPanel)

        # Bottom Panel
        # This section show the plugin operation and a console
        bottomPanel = ttk.Frame(rightPanel)
        tabControl = ttk.Notebook(bottomPanel)  # Create Tab Control
        tabControl.grid(row=1, column=0, sticky='news')

        operationTab = ttk.Frame(tabControl)    # Create a operation tab
        operationTab.grid(row=0, column=0, padx=0, pady=0)
        self._fillRightBottomOperationsPanel(operationTab)
        consoleTab = ttk.Frame(tabControl)    # Create a console
        self._fillRightBottomOutputLogPanel(consoleTab)

        tabControl.add(operationTab, text='Operations')  # Add the Operation tab
        tabControl.add(consoleTab, text='Output Log')
        tabControl.pack(expand=1, fill="both")    # Pack to make visible

        # Add the widgets to Right Panel
        rightPanel.add(topPanel, padx=0, pady=0)
        rightPanel.add(bottomPanel, padx=0, pady=0)

        # Add the Plugin list at left
        mainFrame.add(leftPanel, padx=0, pady=0)
        mainFrame.paneconfig(leftPanel, minsize=200)

        # Add the Plugins or Binaries information and Operation list at right
        mainFrame.add(rightPanel, padx=0, pady=0)
        mainFrame.paneconfig(rightPanel, minsize=200)

    def _fillToolbar(self, frame):
        """ Fill the toolbar frame with some buttons. """
        self._col = 0
        self.executeOpsBtn = self._addButton(frame, '', Icon.TO_INSTALL,
                                       Message.EXECUTE_PLUGINS_MANAGER_OPERATION,
                                          'disable', self._applyAllOperations)
        self.cancelOpsBtn = self._addButton(frame, '', Icon.DELETE_OPERATION,
                              Message.CANCEL_SELECTED_OPERATION, 'disable',
                                         self._deleteSelectedOperation)

    def _addButton(self, frame, text, image, tooltip, state, command):
        btn = IconButton(frame, text, image, command=command,
                         tooltip=tooltip, bg=None)
        btn.config(relief="flat", activebackground=None, compound='left',
                   fg='black', overrelief="raised",
                   state=state)
        btn.bind('<Button-1>', command)
        btn.grid(row=0, column=self._col, sticky='nw',
                 padx=3, pady=7)
        self._col += 1
        return btn

    def _fillLeftPanel(self, leftFrame):
        """
        Fill the left Panel with the plugins list
        """
        gui.configureWeigths(leftFrame)
        self.tree = PluginTree(leftFrame, show="tree")
        self.tree.grid(row=0, column=0, sticky='news')

        self.yscrollbar = ttk.Scrollbar(leftFrame, orient='vertical',
                                        command=self.tree.yview)
        self.yscrollbar.grid(row=0, column=1, sticky='news')
        self.tree.configure(yscrollcommand=self.yscrollbar.set)
        self.yscrollbar.configure(command=self.tree.yview)

        # check / uncheck boxes(plugin or binary) on right click
        self.tree.bind("<Button-1>", self._onPluginTreeClick, True)

        # Load all plugins and fill the tree view
        self.loadPlugins()

    def _createRightTopPanel(self, topPanel):
        """
        Create a right top panel
        """
        self.topPanelTree = ttk.Treeview(topPanel, show='tree', cursor='hand2')
        self.topPanelTree.grid(row=0, column=0, sticky='news')

        # configure vertical scroollbar
        ysb = ttk.Scrollbar(topPanel, orient='vertical',
                            command=self.topPanelTree.yview)
        ysb.grid(row=0, column=1, sticky='news')
        self.topPanelTree.configure(yscrollcommand=ysb.set)
        ysb.configure(command=self.topPanelTree.yview)
        xsb = ttk.Scrollbar(topPanel, orient='horizontal',
                            command=self.topPanelTree.yview)
        xsb.grid(row=1, column=0, sticky='news')
        self.topPanelTree.configure(xscrollcommand=xsb.set)
        xsb.configure(command=self.topPanelTree.xview)
        topPanel.rowconfigure(0, weight=1)
        topPanel.columnconfigure(0, weight=1)
        self.topPanelTree.bind("<Button-1>", self.linkToWebSite, True)

    def _fillRightBottomOperationsPanel(self, panel):
        """
        Create the Operations Tab
        """
        gui.configureWeigths(panel)
        # Define a Tool Bar
        opPanel = tk.Frame(panel)
        opPanel.grid(row=0, column=0, sticky='news')
        gui.configureWeigths(opPanel, 1)

        toolBarFrame = tk.Frame(opPanel)
        toolBarFrame.grid(row=0, column=0, sticky=W)
        self._fillToolbar(toolBarFrame)
        gui.configureWeigths(toolBarFrame)

        # Fill the operation tab
        self.operationTree = PluginTree(opPanel, show="tree")
        self.operationTree.grid(row=1, column=0, sticky='news')
        yscrollbar = ttk.Scrollbar(opPanel, orient='vertical',
                                        command=self.operationTree.yview)
        yscrollbar.grid(row=1, column=1, sticky='news')
        self.operationTree.configure(yscrollcommand=yscrollbar.set)
        yscrollbar.configure(command=self.operationTree.yview)
        self.operationTree.bind("<Button-1>", self.operationInformation, True)

    def _fillRightBottomOutputLogPanel(self, panel):
        """ Create and fill the output log with two tabs(plugin.log and
        plugin.err)"""
        # Fill the Output Log
        gui.configureWeigths(panel)
        self.terminal = tk.Frame(panel)
        self.terminal.grid(row=0, column=0, sticky='news')
        gui.configureWeigths(self.terminal)

        self.Textlog = TextFileViewer(self.terminal, font='black')
        self.Textlog.grid(row=0, column=0, sticky='news')

        self.file_log_path = os.path.join(os.environ['SCIPION_LOGS'],
                                     PLUGIN_LOG_NAME)
        self.file_errors_path = os.path.join(os.environ['SCIPION_LOGS'],
                                        PLUGIN_ERRORS_LOG_NAME)
        self.fileLog = open(self.file_log_path, 'w', 0)
        self.fileLogErr = open(self.file_errors_path, 'w', 0)
        self.plug_log = ScipionLogger(self.file_log_path)
        self.plug_errors_log = ScipionLogger(self.file_errors_path)

    def _onPluginTreeClick(self, event):
        """ check or uncheck a plugin or binary box when clicked """
        if self.tree.is_enabled():
            x, y, widget = event.x, event.y, event.widget
            elem = widget.identify("element", x, y)
            self.tree.selectedItem = self.tree.identify_row(y)
            if "image" in elem:
                # a box was clicked
                tags = self.tree.item(self.tree.selectedItem, "tags")
                objType = self.tree.item(self.tree.selectedItem, "value")
                parent = self.tree.parent(self.tree.selectedItem)
                operation = Operation(self.tree.selectedItem, objType[0], tags[0],
                                      parent)
                self.operationList.insertOperation(operation)
                if tags[0] == PluginStates.UNCHECKED:
                    self.tree.check_item(self.tree.selectedItem)
                elif tags[0] == PluginStates.UNINSTALL:
                    if objType[0] == PluginStates.PLUGIN:
                        self.reloadInstalledPlugin(self.tree.selectedItem)
                else:
                    children = self.tree.get_children(self.tree.selectedItem)
                    for iid in children:
                        self.deleteOperation(iid)
                    self.tree.uncheck_item(self.tree.selectedItem)
                self.showPluginInformation(self.tree.selectedItem)
                self.showOperationList()
            else:
                if self.tree.selectedItem is not None:
                    if self.isPlugin(self.tree.item(self.tree.selectedItem,
                                                    "values")[0]):
                        self.showPluginInformation(self.tree.selectedItem)
                    else:
                        parent = self.tree.parent(self.tree.selectedItem)
                        self.showPluginInformation(parent)
            if len(self.operationList.getOperations(None)):
                self.executeOpsBtn.config(state='normal')
            else:
                self.executeOpsBtn.config(state='disable')

    def _deleteSelectedOperation(self, e=None):
        """
        Delete a selected operation
        """
        if self.operationTree.selectedItem:
            item = self.operationTree.selectedItem
            operation = self.operationList.getOperationByName(item)
            index = self.operationList.operationIndex(operation)
            if index is not None:
                self.operationList.removeOperation(index)
                self.showOperationList()
                if PluginStates.INSTALL in self.tree.item(item, 'tags'):
                    self.tree.item(self.operationTree.selectedItem,
                                   tags=(PluginStates.UNCHECKED,))
                else:
                    if operation.getObjType() == PluginStates.PLUGIN:
                        self.reloadInstalledPlugin(item)
                    else:
                        self.reloadInstalledPlugin(self.tree.parent(item))
                self.operationTree.selectedItem = None
                self.cancelOpsBtn.config(state='disable')
                if not len(self.operationList.getOperations(None)):
                    self.executeOpsBtn.config(state='disable')

    def _applyAllOperations(self, event=None):
        """
        Execute the operation list
        """
        # Disable the execute and cancel button
        self.executeOpsBtn.config(state='disable')
        self.cancelOpsBtn.config(state='disable')
        # Disable the TreeView
        self.tree.disable()
        # Create two tabs where the log and errors will appears
        self.Textlog.createWidgets([self.file_log_path, self.file_errors_path])
        if event is not None:
            threadOp = threading.Thread(target=self._applyOperations,
                                        args=(None,))
            threadOp.start()

    def _applyOperations(self, operation=None):
        """
        Execute one operation. If operation is None, then execute the operation
        list
        """
        # Take the standard system out and errors
        oldstdout = sys.stdout
        oldstderr = sys.stderr
        sys.stdout = self.fileLog
        sys.stderr = self.fileLogErr
        for op in self.operationList.getOperations(operation):
            item = op.getObjName()
            try:
                self.operationTree.processing_item(item)
                self.operationTree.update()
                op.runOperation()
                self.operationTree.installed_item(item)
                self.operationTree.update()
                self.Textlog.refreshAll(goEnd=True)
                self.Textlog.update()
                if op.getObjStatus() == PluginStates.INSTALL:
                    if op.getObjType() == PluginStates.PLUGIN:
                        self.reloadInstalledPlugin(item)
                    else:
                        self.reloadInstalledPlugin(self.tree.parent(item))
                else:
                    self.tree.uncheck_item(item)
            except AssertionError as err:
                self.operationTree.failure_item(item)
                self.tree.uncheck_item(item)
                self.operationTree.update()
                strErr = str('Error executing the operation: ' +
                             op.getObjStatus() + ' ' +
                             op.getObjName())
                self.plug_log.info(redStr(strErr), False)
                self.plug_errors_log.error(redStr(strErr), False)
                self.Textlog.refreshAll(goEnd=True)
                self.Textlog.update()
        self.operationList.clearOperations()
        sys.stdout.flush()
        sys.stderr.flush()
        sys.stdout = oldstdout
        sys.stderr = oldstderr
        # Enable the treeview
        self.tree.enable()

    def linkToWebSite(self, event):
        """
        Load the plugin url
        """
        x, y, widget = event.x, event.y, event.widget
        item = self.topPanelTree.selectedItem = self.topPanelTree.identify_row(y)
        if len(self.topPanelTree.selectedItem) and \
                self.topPanelTree.item(item, 'value')[0] == 'pluginUrl':
            browser = webbrowser.get()
            browser.open(item)

    def operationInformation(self, event):
        """Update the operationTree selected item"""
        x, y, widget = event.x, event.y, event.widget
        elem = widget.identify("element", x, y)
        item = self.operationTree.selectedItem = self.operationTree.identify_row(y)
        if len(item) and len(self.operationList.getOperations(None)):
            self.cancelOpsBtn.config(state='normal')

    def deleteOperation(self, operationName):
        """
        Delete an operation given the object name
        """
        for op in self.operationList.getOperations(None):
            if operationName == op.getObjName():
                self.operationList.insertOperation(op)

    def showOperationList(self):
        """
        Shows the operation list at left bottom panel
        :return:
        """
        self.operationTree.delete(*self.operationTree.get_children())
        for op in self.operationList.getOperations(None):
            self.operationTree.insert("", 'end', op.getObjName(),
                                      text=str(op.getObjStatus().upper() +
                                               ' --> ' + op.getObjName()),
                                      tags=PluginStates.TO_INSTALL)
        self.operationTree.update()

    def isPlugin(self, value):
        return value == PluginStates.PLUGIN

    def showPluginInformation(self, pluginName):
        """Shows the information associated with a given plugin"""
        plugin =  pluginDict.get(pluginName, None)
        if plugin is not None:
            pluginName = plugin.getPipName()
            pluginVersion = plugin.getPipVersion()
            pluginDescription = plugin.getSummary()
            pluginUrl = plugin.getHomePage()
            pluginAuthor = plugin.getAuthor()
            self.topPanelTree.delete(*self.topPanelTree.get_children())

            self.topPanelTree.tag_configure('pluginUrl', foreground='blue')

            self.topPanelTree.insert('', 'end', pluginName,
                                     text='Name:          ' + pluginName,
                                     values='pluginName')
            self.topPanelTree.insert('', 'end', pluginVersion,
                                     text='Version:        ' + pluginVersion,
                                     values='pluginVersion')
            self.topPanelTree.insert('', 'end', pluginDescription,
                                     text='Description:  ' + pluginDescription,
                                     values='pluginDescription')
            self.topPanelTree.insert('', 'end', pluginUrl,
                                     text='URL:             ' + pluginUrl,
                                     values='pluginUrl', tags=('pluginUrl',))
            self.topPanelTree.insert('', 'end', pluginAuthor,
                                     text='Author:         ' + pluginAuthor,
                                     values='pluginAuthor')

    def reloadInstalledPlugin(self, pluginName):
        """
        Reload a given plugin and update the tree view
        """
        plugin = PluginInfo(pluginName, pluginName, remote=False)
        if plugin is not None:
            # Insert all binaries of plugin on the tree
            if plugin.isInstalled():
                pluginBinaryList = plugin.getInstallenv()
                if pluginBinaryList is not None:
                    binaryList = pluginBinaryList.getPackages()
                    keys = sorted(binaryList.keys())
                    for k in keys:
                        pVersions = binaryList[k]
                        for binary, version in pVersions:
                            installed = pluginBinaryList._isInstalled(binary,
                                                                      version)
                            tag = PluginStates.UNCHECKED
                            if installed:
                                tag = PluginStates.CHECKED
                            binaryName = str(binary + '-' + version)
                            if binaryName in self.tree.get_children(pluginName):
                                self.tree.item(binaryName, tags=(tag,))
                            else:
                                self.tree.insert(pluginName, "end", binaryName,
                                                 text=binaryName, tags=tag,
                                                 values='binary')
                self.tree.item(pluginName, tags=(PluginStates.CHECKED,))
            else:
                if PluginStates.UNINSTALL in self.tree.item(pluginName, 'tags'):
                    self.tree.item(pluginName, tags=(PluginStates.UNCHECKED,))
        else:
            self.tree.item(pluginName, tags=(PluginStates.UNCHECKED,))

    def loadPlugins(self):
        """
        Load all plugins and fill the tree view widget
        """
        global pluginDict
        pluginDict = pluginRepo.getPlugins(getPipData=True)
        pluginList = sorted(pluginDict.keys(), reverse=True)
        self.tree.delete(*self.tree.get_children())
        for pluginName in pluginList:
            plugin = PluginInfo(pluginName, pluginName, remote=False)
            if plugin is not None:
                tag = PluginStates.UNCHECKED
                if plugin._getPlugin():
                    # Insert the plugin name in the tree
                    tag = PluginStates.CHECKED
                    self.tree.insert("", 0, pluginName, text=pluginName, tags=tag,
                                     values=PluginStates.PLUGIN)
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
                                tag = PluginStates.UNCHECKED
                                if installed:
                                    tag = PluginStates.CHECKED
                                binaryName = str(binary + '-' + version)
                                self.tree.insert(pluginName, "end", binaryName,
                                                 text=binaryName, tags=tag,
                                                 values=PluginStates.BINARY)
                else:
                    self.tree.insert("", 0, pluginName, text=pluginName, tags=tag,
                                     values=PluginStates.PLUGIN)


class PluginManagerWindow(gui.Window):
    """
     Windows to hold a plugin manager frame inside.
    """
    def __init__(self, title, master=None, **kwargs):
        if 'minsize' not in kwargs:
            kwargs['minsize'] = (900, 300)
        gui.Window.__init__(self, title, master, **kwargs)

        menu = MenuConfig()
        fileMenu = menu.addSubMenu('File')
        fileMenu.addSubMenu('Exit', 'exit', icon='fa-sign-out.png')

        helpMenu = menu.addSubMenu('Help')
        helpMenu.addSubMenu('Help', 'help', icon='fa-question-circle.png')
        self.menuCfg = menu
        gui.Window.createMainMenu(self, self.menuCfg)

    def onExit(self):
        self.close()

    def onBrowsePlugin(self):
        pass

    def onHelp(self):
        PluginHelp('Plugin Manager Glossary', self).show()


class PluginHelp(gui.Window):
    """
    Windows to hold a plugin manager help
    """
    def __init__(self, title, master=None, **kwargs):
        if 'minsize' not in kwargs:
            kwargs['minsize'] = (500, 300)
            gui.Window.__init__(self, title, master, **kwargs)
        self.root.resizable(0, 0)
        self.createHelp()

    def createHelp(self):
        helpFrame = tk.Frame(self.root)
        helpFrame.grid(row=0, column=0, sticky='news')
        photo = PhotoImage(file=gui.findResource(Icon.CHECKED))
        btn = Label(helpFrame, image=photo)
        btn.photo = photo
        btn.grid(row=0, column=0, sticky='sw', padx=10, pady=5)
        btn = Label(helpFrame, text='INSTALLED Plugin/Binary')
        btn.grid(row=0, column=1, sticky='sw', padx=0, pady=5)

        photo = PhotoImage(file=gui.findResource(Icon.UNCHECKED))
        btn = Label(helpFrame, image=photo)
        btn.photo = photo
        btn.grid(row=1, column=0, sticky='sw', padx=10, pady=5)
        btn = Label(helpFrame, text='UNINSTALLED Plugin/Binary')
        btn.grid(row=1, column=1, sticky='sw', padx=0, pady=5)

        photo = PhotoImage(file=gui.findResource(Icon.INSTALL))
        btn = Label(helpFrame, image=photo)
        btn.photo = photo
        btn.grid(row=2, column=0, sticky='sw', padx=10, pady=5)
        btn = Label(helpFrame, text='Plugin/Binary TO INSTALL')
        btn.grid(row=2, column=1, sticky='sw', padx=0, pady=5)

        photo = PhotoImage(file=gui.findResource(Icon.UNINSTALL))
        btn = Label(helpFrame, image=photo)
        btn.photo = photo
        btn.grid(row=3, column=0, sticky='sw', padx=10, pady=5)
        btn = Label(helpFrame, text='Plugin/Binary TO UNINSTALL')
        btn.grid(row=3, column=1, sticky='sw', padx=0, pady=5)

        photo = PhotoImage(file=gui.findResource(Icon.TO_INSTALL))
        btn = Label(helpFrame, image=photo)
        btn.photo = photo
        btn.grid(row=4, column=0, sticky='sw', padx=10, pady=5)
        btn = Label(helpFrame, text='Execute the selected operations')
        btn.grid(row=4, column=1, sticky='sw', padx=0, pady=5)

        photo = PhotoImage(file=gui.findResource(Icon.DELETE_OPERATION))
        btn = Label(helpFrame, image=photo)
        btn.photo = photo
        btn.grid(row=5, column=0, sticky='sw', padx=10, pady=5)
        btn = Label(helpFrame, text='Cancel a selected operation')
        btn.grid(row=5, column=1, sticky='sw', padx=0, pady=5)


class PluginManager(PluginManagerWindow):
    """
    Windows to hold a frame inside.
    """
    def __init__(self, title, master=None, path=None,
                 onSelect=None, shortCuts=None, **kwargs):
        PluginManagerWindow.__init__(self, title, master, **kwargs)
        browser = PluginBrowser(self.root, **kwargs)
