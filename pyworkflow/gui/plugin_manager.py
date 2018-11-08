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
from pyworkflow.config import MenuConfig
from pyworkflow.utils.log import ScipionLogger
from pyworkflow.gui.text import TextFileViewer
from pyworkflow.gui import *
from pyworkflow.gui.form import *
from install.plugin_funcs import PluginRepository, PluginInfo
from pyworkflow.utils import *

PLUGIN = 'plugin'
BINARY = 'binary'
UNCHECKED = 'unchecked'
CHECKED = 'checked'
INSTALL = 'install'
UNINSTALL = 'uninstall'
TO_INSTALL = 'to_install'
INSTALLED = 'installed'
PRECESSING = 'processing'
FAILURE = 'failure'
PLUGIN_LOG_NAME = 'Plugin.log'
PLUGIN_ERRORS_LOG_NAME = 'Plugin.err'

pluginRepo = PluginRepository()
pluginDict = pluginRepo.getPlugins(getPipData=True)


class PluginTree(ttk.Treeview):
    """
        Treeview widget with checkboxes left of each item.
        The checkboxes are done via the image attribute of the item, so to keep
        the checkbox, you cannot add an image to the item.
    """

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
        self.tag_configure(UNCHECKED, image=self.im_unchecked)
        self.tag_configure(CHECKED, image=self.im_checked)
        self.tag_configure(INSTALL, image=self.im_install)
        self.tag_configure(UNINSTALL, image=self.im_uninstall)
        self.tag_configure(TO_INSTALL, image=self.im_to_install)
        self.tag_configure(INSTALLED, image=self.im_installed)
        self.tag_configure(PRECESSING, image=self.im_processing)
        self.tag_configure(FAILURE, image=self.im_failure)
        self.selectedItem = None

    def insert(self, parent, index, iid=None, **kw):
        """ same method as for standard treeview but add the tag 'unchecked'
            automatically if no tag among ('checked', 'unchecked')
            is given """
        if not "tags" in kw:
            kw["tags"] = (UNCHECKED,)
        elif not (UNCHECKED in kw["tags"] or CHECKED in kw["tags"] or
                  TO_INSTALL in kw["tags"] or INSTALLED in kw["tags"]):
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

    def processing_item(self, item):
        """change the box item to processing item"""
        self.item(item, tags=(PRECESSING,))

    def installed_item(self, item):
        """change the box item to processing item"""
        self.item(item, tags=(INSTALLED,))

    def failure_item(self, item):
        """change the box item to failure item"""
        self.item(item, tags=(FAILURE,))


class Operation:
    """
    This class contain the object(plugin/binary) operation details
    """
    def __init__(self, objName, objType=PLUGIN, objStatus=INSTALL,
                 objParent=None):
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
        if self.objType == PLUGIN:
            if self.objStatus == INSTALL:
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
            if self.objStatus == INSTALL:
                if plugin is not None:
                    plugin.installBin([self.objName])
            else:
                plugin.uninstallBins([self.objName])


class OperationList:
    """
    This class contain a plugins/binaries operation list and allow execute it
    """
    def __init__(self):
        self.operationList = []

    def insertOperation(self, operation):
        """
        This method insert into the list a given operation
        """
        index = self.operationIndex(operation)
        if index is not None:
            self.removeOperation(index)
        else:
            tag = UNINSTALL
            if operation.getObjStatus() == UNCHECKED:
                tag = INSTALL
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

    def getOperations(self):
        """
        Return the operation List
        """
        return self.operationList

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
        self.operationList = []


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

        # Define a Tool Bar
        toolBarFrame = tk.Frame(parentFrame)
        toolBarFrame.grid(row=0, column=0, sticky=W)
        self._fillToolbar(toolBarFrame)
        gui.configureWeigths(toolBarFrame)
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
        self._addButton(frame, 'Apply', Icon.ACTION_EXECUTE,
                        Message.EXECUTE_PLUGINS_MANAGER_OPERATION, 'normal',
                        self._applyOperations)

    def _addButton(self, frame, text, image, tooltip, state, command):
        btn = IconButton(frame, text, image, command=command,
                         tooltip=tooltip, bg=None)
        btn.config(relief="flat", activebackground=None, compound='left',
                   fg='black', overrelief="raised",
                   state=state)
        btn.bind('<Button-1>', command)
        btn.grid(row=0, column=self._col, sticky='nw',
                 padx=(0, 5), pady=0)
        self._col += 1

    def _fillLeftPanel(self, leftFrame):
        """
        Fill the left Panel with the plugins list
        """
        gui.configureWeigths(leftFrame)
        pluginColumn = ('Value')
        self.tree = PluginTree(leftFrame, show="tree", columns=pluginColumn)
        self.tree.grid(row=0, column=0, sticky='news')

        self.yscrollbar = ttk.Scrollbar(leftFrame, orient='vertical',
                                        command=self.tree.yview)
        self.yscrollbar.grid(row=0, column=1, sticky='news')
        self.tree.configure(yscrollcommand=self.yscrollbar.set)
        self.yscrollbar.configure(command=self.tree.yview)

        # check / uncheck boxes(plugin or binary) on right click
        self.tree.bind("<Button-3>", self._box_rightClick, True)
        # show the plugin or binary information on click
        self.tree.bind("<Button-1>", self.objectInformation, True)

        # Load all plugins and fill the tree view
        self.loadPlugins()

    def _fillRightBottomOperationsPanel(self, panel):
        # Fill the operation tab
        gui.configureWeigths(panel)
        self.operationTree = PluginTree(panel, show="tree")
        self.operationTree.grid(row=0, column=0, sticky='news')

        yscrollbar = ttk.Scrollbar(panel, orient='vertical',
                                        command=self.operationTree.yview)
        yscrollbar.grid(row=0, column=1, sticky='news')
        self.operationTree.configure(yscrollcommand=yscrollbar.set)
        yscrollbar.configure(command=self.operationTree.yview)

    def _fillRightBottomOutputLogPanel(self, panel):
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

    def _box_rightClick(self, event):
        """ check or uncheck a plugin or binary box when clicked """
        x, y, widget = event.x, event.y, event.widget
        elem = widget.identify("element", x, y)
        if "image" in elem:
            # a box was clicked
            self.tree.selectedItem = self.tree.identify_row(y)
            tags = self.tree.item(self.tree.selectedItem, "tags")
            objType = self.tree.item(self.tree.selectedItem, "value")
            parent = self.tree.parent(self.tree.selectedItem)
            operation = Operation(self.tree.selectedItem, objType[0], tags[0],
                                  parent)
            self.operationList.insertOperation(operation)
            if tags[0] in [UNCHECKED, UNINSTALL]:
                self.tree.check_item(self.tree.selectedItem)
                if objType[0] == PLUGIN:
                    self.reloadInstalledPlugin(self.tree.selectedItem)
            else:
                children = self.tree.get_children(self.tree.selectedItem)
                for iid in children:
                    self.deleteOperation(iid)
                self.tree.uncheck_item(self.tree.selectedItem)
            self.showPluginInformation(self.tree.selectedItem)
            self.showOperationList()

    def _applyOperations(self, e=None):
        """
        Execute all operation
        """
        # Take the standard system out and errors
        oldstdout = sys.stdout
        oldstderr = sys.stderr
        sys.stdout = self.fileLog
        sys.stderr = self.fileLogErr
        # Create two tabs where the log and errors will appears
        self.Textlog.createWidgets([self.file_log_path, self.file_errors_path])
        for op in self.operationList.getOperations():
            item = op.getObjName()
            try:
                self.operationTree.processing_item(item)
                self.operationTree.update()
                op.runOperation()
                self.Textlog.refreshAll(goEnd=True)
                self.Textlog.update()
                self.operationTree.installed_item(item)
                self.operationTree.update()
                if op.getObjType() == PLUGIN:
                    self.reloadInstalledPlugin(item)
                else:
                    self.reloadInstalledPlugin(self.tree.parent(item))
            except AssertionError as err:
                self.operationTree.failure_item(item)
                self.tree.uncheck_item(item)
                self.operationTree.update()
                strErr = str('Error executing the operation: ' +
                             op.getObjStatus() + ' ' +
                             op.getObjName() + '\n')
                sys.stdout.write(redStr(strErr))
                sys.stderr.write(redStr(strErr))
                self.Textlog.refreshAll(goEnd=True)
                self.Textlog.update()
        self.operationList.clearOperations()
        sys.stdout.flush()
        sys.stderr.flush()
        sys.stdout = oldstdout
        sys.stderr = oldstderr

    def objectInformation(self, event):
        """Show the plugin or binary information"""
        x, y, widget = event.x, event.y, event.widget
        item = self.tree.selectedItem = self.tree.identify_row(y)
        if self.tree.selectedItem is not None:
            if self.isPlugin(self.tree.item(self.tree.selectedItem,
                                            "values")[0]):
                self.showPluginInformation(item)
            else:
                parent = self.tree.parent(item)
                self.showPluginInformation(parent)

    def deleteOperation(self, operationName):
        """
        Delete an operation given the object name
        """
        for op in self.operationList.getOperations():
            if operationName == op.getObjName():
                self.operationList.insertOperation(op)

    def showOperationList(self):
        """
        Shows the operation list at left bottom panel
        :return:
        """
        self.operationTree.delete(*self.operationTree.get_children())
        for op in self.operationList.getOperations():
            self.operationTree.insert("", 'end', op.getObjName(),
                                      text=str(op.getObjStatus().upper() +
                                               ' --> ' + op.getObjName()),
                                      tags=TO_INSTALL)
        self.operationTree.update()

    def isPlugin(self, value):
        return value == PLUGIN

    def showPluginInformation(self, pluginName):
        """Shows the information associated with a given plugin"""
        plugin =  pluginDict.get(pluginName, None)
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
        Reload a given plugin and update the tree view
        """
        plugin = PluginInfo(pluginName, pluginName, remote=False)
        if plugin is not None:
            # Insert all binaries of plugin on the tree
            if plugin._getPlugin():
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
                            if binaryName in self.tree.get_children(pluginName):
                                self.tree.item(binaryName, tags=(tag,))
                            else:
                                self.tree.insert(pluginName, "end", binaryName,
                                                 text=binaryName, tags=tag,
                                                 values='binary')
                self.tree.item(pluginName, tags=(CHECKED,))

            else:
                if UNINSTALL in self.tree.item(pluginName, 'tags'):
                    self.tree.item(pluginName, tags=(UNCHECKED,))

    def loadPlugins(self):
        """
        Load all plugins and fill the tree view widget
        """
        self.tree.delete(*self.tree.get_children())
        for pluginObj in pluginDict:
            plugin = pluginDict.get(pluginObj, None)
            if plugin is not None:
                tag = UNCHECKED
                if plugin.isInstalled():
                    # Insert the plugin name in the tree
                    tag = CHECKED
                    self.tree.insert("", 0, pluginObj, text=pluginObj, tags=tag,
                                     values=PLUGIN)
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
                                                 values=BINARY)
                else:
                    self.tree.insert("", 0, pluginObj, text=pluginObj, tags=tag,
                                     values=PLUGIN)


class PluginManagerWindow(gui.Window):
    """
     Windows to hold a plugin manager frame inside.
    """
    def __init__(self, title, master=None, **kwargs):
        if 'minsize' not in kwargs:
            kwargs['minsize'] = (300, 300)
        gui.Window.__init__(self, title, master, **kwargs)

        menu = MenuConfig()

        fileMenu = menu.addSubMenu('File')
        # #fileMenu.addSubMenu('Browse Plugin', 'browse', icon='fa-folder-open.png')
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
        btn.grid(row=0, column=1, sticky='sw', padx=0, pady=0)

        photo = PhotoImage(file=gui.findResource(Icon.UNCHECKED))
        btn = Label(helpFrame, image=photo)
        btn.photo = photo
        btn.grid(row=1, column=0, sticky='sw', padx=10, pady=5)
        btn = Label(helpFrame, text='UNINSTALLED Plugin/Binary')
        btn.grid(row=1, column=1, sticky='sw', padx=0, pady=0)

        photo = PhotoImage(file=gui.findResource(Icon.INSTALL))
        btn = Label(helpFrame, image=photo)
        btn.photo = photo
        btn.grid(row=2, column=0, sticky='sw', padx=10, pady=5)
        btn = Label(helpFrame, text='Plugin/Binary TO INSTALL')
        btn.grid(row=2, column=1, sticky='sw', padx=0, pady=0)

        photo = PhotoImage(file=gui.findResource(Icon.UNINSTALL))
        btn = Label(helpFrame, image=photo)
        btn.photo = photo
        btn.grid(row=3, column=0, sticky='sw', padx=10, pady=5)
        btn = Label(helpFrame, text='Plugin/Binary TO UNINSTALL')
        btn.grid(row=3, column=1, sticky='sw', padx=0, pady=0)

        photo = PhotoImage(file=gui.findResource(Icon.ACTION_EXECUTE))
        btn = Label(helpFrame, image=photo)
        btn.photo = photo
        btn.grid(row=4, column=0, sticky='sw', padx=10, pady=5)
        btn = Label(helpFrame, text='Apply the selected operations')
        btn.grid(row=4, column=1, sticky='sw', padx=0, pady=0)

        photo = PhotoImage(file=gui.findResource(Icon.TO_INSTALL))
        btn = Label(helpFrame, image=photo)
        btn.photo = photo
        btn.grid(row=5, column=0, sticky='sw', padx=10, pady=5)
        btn = Label(helpFrame, text='Operation waiting to be executed')
        btn.grid(row=5, column=1, sticky='sw', padx=0, pady=0)

        btn = Label(helpFrame, text='Right-Click')
        btn.photo = photo
        btn.grid(row=6, column=0, sticky='sw', padx=10, pady=5)
        btn = Label(helpFrame, text='Change or Undo the status of the '
                                    'plugins/binaries in the tree view')
        btn.grid(row=6, column=1, sticky='sw', padx=0, pady=0)


class PluginManager(PluginManagerWindow):
    """
    Windows to hold a frame inside.
    """
    def __init__(self, title, master=None, path=None,
                 onSelect=None, shortCuts=None, **kwargs):
        PluginManagerWindow.__init__(self, title, master, **kwargs)
        browser = PluginBrowser(self.root, **kwargs)
