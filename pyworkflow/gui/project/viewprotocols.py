#!/usr/bin/env python
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
View with the protocols inside the main project window.
"""

import os
import pickle
from collections import OrderedDict
import Tkinter as tk
import ttk
import datetime as dt
import pyworkflow.object as pwobj
import pyworkflow.utils as pwutils
import pyworkflow.protocol as pwprot
import pyworkflow.gui as pwgui
import pyworkflow.em as em
from pyworkflow.config import isAFinalProtocol
from pyworkflow.em.wizard import ListTreeProvider
from pyworkflow.gui.dialog import askColor, ListDialog
from pyworkflow.viewer import DESKTOP_TKINTER, ProtocolViewer
from pyworkflow.utils.properties import Message, Icon, Color

from constants import STATUS_COLORS

DEFAULT_BOX_COLOR = '#f8f8f8'

ACTION_EDIT = Message.LABEL_EDIT
ACTION_COPY = Message.LABEL_COPY
ACTION_DELETE = Message.LABEL_DELETE
ACTION_REFRESH = Message.LABEL_REFRESH
ACTION_STEPS = Message.LABEL_STEPS
ACTION_BROWSE = Message.LABEL_BROWSE
ACTION_DB = Message.LABEL_DB
ACTION_TREE = Message.LABEL_TREE
ACTION_LIST = Message.LABEL_LIST
ACTION_STOP = Message.LABEL_STOP
ACTION_DEFAULT = Message.LABEL_DEFAULT
ACTION_CONTINUE = Message.LABEL_CONTINUE
ACTION_RESULTS = Message.LABEL_ANALYZE
ACTION_EXPORT = Message.LABEL_EXPORT
ACTION_SWITCH_VIEW = 'Switch_View'
ACTION_COLLAPSE = 'Collapse'
ACTION_EXPAND = 'Expand'
ACTION_LABELS = 'Labels'

RUNS_TREE = Icon.RUNS_TREE
RUNS_LIST = Icon.RUNS_LIST

VIEW_LIST = 0
VIEW_TREE = 1
VIEW_TREE_SMALL = 2
 
ActionIcons = {
    ACTION_EDIT: Icon.ACTION_EDIT , 
    ACTION_COPY: Icon.ACTION_COPY ,
    ACTION_DELETE:  Icon.ACTION_DELETE,
    ACTION_REFRESH:  Icon.ACTION_REFRESH,
    ACTION_STEPS:  Icon.ACTION_STEPS,
    ACTION_BROWSE:  Icon.ACTION_BROWSE,
    ACTION_DB: Icon.ACTION_DB,
    ACTION_TREE:  None, # should be set
    ACTION_LIST:  Icon.ACTION_LIST,
    ACTION_STOP: Icon.ACTION_STOP,
    ACTION_CONTINUE: Icon.ACTION_CONTINUE,
    ACTION_RESULTS: Icon.ACTION_RESULTS,
    ACTION_EXPORT: Icon.ACTION_EXPORT,
    ACTION_COLLAPSE: 'fa-minus-square.png',
    ACTION_EXPAND: 'fa-plus-square.png',
    ACTION_LABELS: Icon.TAGS
              }



def populateTree(self, tree, treeItems, prefix, obj, subclassedDict, level=0):
    text = obj.text.get()
    if text:
        value = obj.value.get(text)
        key = '%s.%s' % (prefix, value)
        
        img = obj.icon.get('')
        tag = obj.tag.get('')
            
        if len(img):
            img = self.getImage(img)
        
        protClassName = value.split('.')[-1] # Take last part
        emProtocolsDict = em.getProtocols()
        prot = emProtocolsDict.get(protClassName, None)
        
        if tag == 'protocol' and text == 'default':
            if prot is None:
                raise Exception("Protocol className '%s' not found!!!. \n"
                                "Fix your config/scipion.conf configuration."
                                % protClassName)
                 
            text = prot.getClassLabel()
             
        item = tree.insert(prefix, 'end', key, text=text, image=img, tags=(tag))
        treeItems[item] = obj
        # Check if the attribute should be open or close
        openItem = getattr(obj, 'openItem', level < 2)
        if openItem:
            tree.item(item, open=True)
            
        if obj.value.hasValue() and tag == 'protocol_base':
            if prot is not None:
                tree.item(item, image=self.getImage('class_obj.gif'))
                
                for k, v in emProtocolsDict.iteritems():
                    if (not k in subclassedDict and
                        not v is prot and issubclass(v, prot)):
                        key = '%s.%s' % (item, k)
                        t = v.getClassLabel()
                        tree.insert(item, 'end', key, text=t, tags=('protocol'))
                        
            else:
                raise Exception("Class '%s' not found" % obj.value.get())
    else:
        key = prefix
    
    for sub in obj:
        populateTree(self, tree, treeItems, key, sub, subclassedDict, level+1)
    

class RunsTreeProvider(pwgui.tree.ProjectRunsTreeProvider):
    """Provide runs info to populate tree"""

    def __init__(self, project, actionFunc):
        pwgui.tree.ProjectRunsTreeProvider.__init__(self, project)
        self.actionFunc = actionFunc
        self._selection = project.getSettings().runSelection
    
    def getActionsFromSelection(self):
        """ Return the list of options available for selection. """
        n = len(self._selection)
        single = n == 1
        if n:
            prot = self.project.getProtocol(self._selection[0]) 
            status = prot.getStatus()
            nodeInfo = self.project.getSettings().getNodeById(prot.getObjId())
            expanded = nodeInfo.isExpanded() if nodeInfo else True
        else:
            status = None
        
        return [   (ACTION_EDIT, single),
                   (ACTION_COPY, True), 
                   (ACTION_DELETE, status != pwprot.STATUS_RUNNING),
                   (ACTION_STEPS, single),
                   (ACTION_BROWSE, single),
                   (ACTION_DB, single),
                   (ACTION_STOP, status == pwprot.STATUS_RUNNING and single),
                   (ACTION_EXPORT, not single),
                   (ACTION_COLLAPSE, single and status and expanded),
                   (ACTION_EXPAND, single and status and not expanded),
                   (ACTION_LABELS, True)
                   ]
        
    def getObjectActions(self, obj):
        
        def addAction(a):
            if a:
                text = a
                action = a
                a = (text, lambda: self.actionFunc(action),
                     ActionIcons.get(action, None))
            return a
            
        actions = [addAction(a)
                   for a, cond in self.getActionsFromSelection() if cond]
        
        return actions 
    
    
class ProtocolTreeProvider(pwgui.tree.ObjectTreeProvider):
    """Create the tree elements for a Protocol run"""
    def __init__(self, protocol):
        self.protocol = protocol
        # This list is create to group the protocol parameters
        # in the tree display
        self.status = pwobj.List(objName='_status')
        self.params = pwobj.List(objName='_params')
        self.statusList = ['status', 'initTime', 'endTime', 'error',
                           'interactive', 'mode']

        objList = [] if protocol is None else [protocol]
        pwgui.tree.ObjectTreeProvider.__init__(self, objList)
        
        
class StepsTreeProvider(pwgui.tree.TreeProvider):
    """Create the tree elements for a Protocol run"""
    def __init__(self, stepsList):
        for i, s in enumerate(stepsList):
            if not s._index:
                s._index = i + 1
            
        self._stepsList = stepsList
        self.getColumns = lambda: [('Index', 50), ('Step', 200),
                                   ('Time', 150), ('Class', 100)]
        self._parentDict = {}
    
    def getObjects(self):
        return self._stepsList
        
    def getObjectInfo(self, obj):
        info = {'key': obj._index, 
                'values': (str(obj), pwutils.prettyDelta(obj.getElapsedTime()),
                           obj.getClassName())}
            
        return info
    
    def getObjectPreview(self, obj):
        args = pickle.loads(obj.argsStr.get())
        msg = "*Prerequisites*: %s \n" % str(obj._prerequisites)
        msg += "*Arguments*: " + '\n  '.join([str(a) for a in args])
        if hasattr(obj, 'resultFiles'):
            results = pickle.loads(obj.resultFiles.get())
            if len(results):
                msg += "\n*Result files:* " + '\n  '.join(results)

        return None, msg
    

class StepsWindow(pwgui.browser.BrowserWindow):
    def __init__(self, title, parentWindow, protocol, **args):
        self._protocol = protocol
        provider = StepsTreeProvider(protocol.loadSteps())
        pwgui.browser.BrowserWindow.__init__(self, title, parentWindow,
                                             weight=False, **args)
        # Create buttons toolbar
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(1, weight=1)
        
        toolbar = tk.Frame(self.root)
        toolbar.grid(row=0, column=0, sticky='nw', padx=5, pady=5)
        btn = tk.Label(toolbar, text="Tree", image=self.getImage(Icon.RUNS_TREE), 
                       compound=tk.LEFT, cursor='hand2')
        btn.bind('<Button-1>', self._showTree)
        btn.grid(row=0, column=0, sticky='nw')
        # Create and set browser
        browser = pwgui.browser.ObjectBrowser(self.root, provider,
                                              showPreviewTop=False)
        self.setBrowser(browser, row=1, column=0)
        
    def _showTree(self, e=None):
        g = self._protocol.getStepsGraph()
        w = pwgui.Window("Protocol steps", self, minsize=(800, 600))
        root = w.root
        canvas = pwgui.Canvas(root, width=600, height=500)
        canvas.grid(row=0, column=0, sticky='nsew')
        canvas.drawGraph(g, pwgui.graph.LevelTreeLayout())
        w.show()
    

class SearchProtocolWindow(pwgui.Window):
    def __init__(self, parentWindow, **kwargs):
        pwgui.Window.__init__(self, title="Search Protocol",
                              masterWindow=parentWindow)
        content = tk.Frame(self.root, bg='white')
        self._createContent(content)
        content.grid(row=0, column=0, sticky='news')
        content.columnconfigure(0, weight=1)
        content.rowconfigure(1, weight=1)
        
    def _createContent(self, content):
        self._createSearchBox(content)
        self._createResultsBox(content)
        
    def _createSearchBox(self, content):
        """ Create the Frame with Search widgets """
        frame = tk.Frame(content, bg='white')
                
        label = tk.Label(frame, text="Search", bg='white')
        label.grid(row=0, column=0, sticky='nw')
        self._searchVar = tk.StringVar()
        entry = tk.Entry(frame, bg='white', textvariable=self._searchVar)
        entry.bind('<Return>', self._onSearchClick)
        entry.focus_set()
        entry.grid(row=0, column=1, sticky='nw')
        btn = pwgui.widgets.IconButton(frame, "Search",
                                       imagePath=Icon.ACTION_SEARCH,
                         command=self._onSearchClick)
        btn.grid(row=0, column=2, sticky='nw')
        
        frame.grid(row=0, column=0, sticky='new', padx=5, pady=(10, 5))

    def _createResultsBox(self, content):
        frame = tk.Frame(content, bg=Color.LIGHT_GREY_COLOR, padx=5, pady=5)
        pwgui.configureWeigths(frame)
        self._resultsTree = self.master.getViewWidget()._createProtocolsTree(frame)
        self._resultsTree.grid(row=0, column=0, sticky='news')
        frame.grid(row=1, column=0, sticky='news', padx=5, pady=5)
        
    def _onSearchClick(self, e=None):
        self._resultsTree.clear()
        keyword = self._searchVar.get().lower()
        emProtocolsDict = em.getProtocols()
        protList = []
        
        for key, prot in emProtocolsDict.iteritems():
            if isAFinalProtocol(prot,key):
                label = prot.getClassLabel().lower()
                if keyword in label:
                    protList.append((key, label))
        
        protList.sort(key=lambda x: x[1]) # sort by label
        for key, label in protList:             
            self._resultsTree.insert('', 'end', key,
                                     text=label, tags=('protocol'))
        

class RunIOTreeProvider(pwgui.tree.TreeProvider):
    """Create the tree elements from a Protocol Run input/output childs"""
    def __init__(self, parent, protocol, mapper):
        #TreeProvider.__init__(self)
        self.parent = parent
        self.protocol = protocol
        self.mapper = mapper

    def getColumns(self):
        return [('Attribute', 200), ('Info', 100)]
    
    def getObjects(self):
        objs = []
        if self.protocol:
            # Store a dict with input parents (input, PointerList)
            self.inputParentDict = OrderedDict()
            inputs = []
            inputObj = pwobj.String(Message.LABEL_INPUT)
            inputObj._icon = Icon.ACTION_IN
            self.inputParentDict['_input'] = inputObj
            inputParents = [inputObj]
            
            for key, attr in self.protocol.iterInputAttributes():
                attr._parentKey = key
                # Repeated keys means there are inside a pointerList
                # since the same key is yielded for all items inside
                # so update the parent dict with a new object
                if key in self.inputParentDict:
                    if self.inputParentDict[key] == inputObj:
                        parentObj = pwobj.String(key)
                        parentObj._icon = Icon.ACTION_IN
                        parentObj._parentKey = '_input'
                        inputParents.append(parentObj)
                        self.inputParentDict[key] = parentObj
                else:
                    self.inputParentDict[key] = inputObj
                inputs.append(attr) 
                    
            outputs = [attr for _, attr in self.protocol.iterOutputAttributes(em.EMObject)]
            self.outputStr = pwobj.String(Message.LABEL_OUTPUT)
            objs = inputParents + inputs + [self.outputStr] + outputs                
        return objs
    
    def _visualizeObject(self, ViewerClass, obj):
        viewer = ViewerClass(project=self.protocol.getProject(),
                             protocol=self.protocol,
                             parent=self.parent.windows)
        viewer.visualize(obj)
        
    def _editObject(self, obj):
        """Open the Edit GUI Form given an instance"""
        pwgui.dialog.EditObjectDialog(self.parent, Message.TITLE_EDIT_OBJECT,
                                      obj, self.mapper)
        
    def _deleteObject(self, obj):
        """ Remove unnecesary output, specially for Coordinates. """
        prot = self.protocol
        try:
            objLabel = self.getObjectLabel(obj, prot)
            if self.parent.windows.askYesNo("Delete object",
                                            "Are you sure to delete *%s* object?"
                                            % objLabel):
                prot.getProject().deleteProtocolOutput(prot, obj)
                self.parent._fillSummary()
                self.parent.windows.showInfo("Object *%s* successfuly deleted."
                                             % objLabel)
        except Exception, ex:
            self.parent.windows.showError(str(ex))
        
    def getObjectPreview(self, obj):
        desc = "<name>: " + obj.getName()
        
        return (None, desc)
    
    def getObjectActions(self, obj):
        if isinstance(obj, pwobj.Pointer):
            obj = obj.get()
            isPointer = True
        else:
            isPointer = False
        actions = []    
        
        viewers = em.findViewers(obj.getClassName(), DESKTOP_TKINTER)
        def yatevalejosemi(v):
            return lambda: self._visualizeObject(v, obj)
        for v in viewers:
            actions.append(('Open with %s' % v.__name__,
                            yatevalejosemi(v), 
                            Icon.ACTION_VISUALIZE))
            
        # EDIT 
        actions.append((Message.LABEL_EDIT, 
                        lambda : self._editObject(obj),
                        Icon.ACTION_EDIT))
        
        # DELETE
        # Special case to allow delete outputCoordinates
        # since we can end up with several outputs and
        # we may want to clean up
        if self.protocol.allowsDelete(obj) and not isPointer:
            actions.append((Message.LABEL_DELETE_ACTION,
                           lambda: self._deleteObject(obj),
                           Icon.ACTION_DELETE))
        
        return actions
    
    def getObjectLabel(self, obj, parent):
        """ We will try to show in the list the string representation
        that is more readable for the user to pick the desired object.
        """
        label = 'None'
        if obj:
            label = obj.getObjLabel()
            if not len(label.strip()):
                if parent:
                    parentLabel = parent.getObjLabel()
                else:
                    parentLabel = 'None'
                label = "%s -> %s" % (parentLabel, obj.getLastName())
        return label
        
    def getObjectInfo(self, obj):
        if obj is None or not obj.hasValue():
            return None
        
        if isinstance(obj, pwobj.String):
            value = obj.get()
            info = {'key': value, 'text': value, 'values': (''), 'open': True}
            if hasattr(obj, '_parentKey'):
                info['parent'] = self.inputParentDict[obj._parentKey]
        else:
            # All attributes are considered output, unless they are pointers
            image = Icon.ACTION_OUT
            parent = self.outputStr
            
            if isinstance(obj, pwobj.Pointer):
                name = obj.getLastName()
                # Remove ugly item notations inside lists
                name = name.replace('__item__000', '')
                # Consider Pointer as inputs
                image = getattr(obj, '_icon', '')
                parent = self.inputParentDict[obj._parentKey]
                
                suffix = ''
                if obj.hasExtended():
                    extendedValue = obj.getExtended()
                    if obj.hasExtended():
                        suffix = '[%s]' % extendedValue
                    elif obj.hasExtended():
                        suffix = '[Item %s]' % extendedValue
                    if obj.get() is None:
                        labelObj = obj.getObjValue()
                        suffix = ''
                    else:
                        labelObj = obj.get()
                else:
                    labelObj = obj.get()
                    
                objKey = obj._parentKey + str(labelObj.getObjId())
                label = self.getObjectLabel(labelObj,
                                            self.mapper.getParent(labelObj))
                name += '   (from %s %s)' % (label, suffix)
            else:
                name = self.getObjectLabel(obj, self.protocol)
                objKey = str(obj.getObjId())
                labelObj = obj
                
            info = {'key': objKey, 'parent': parent, 'image': image,
                    'text': name, 'values': (str(labelObj),)}
        return info     


class ProtocolsView(tk.Frame):

    """ What you see when the "Protocols" tab is selected.

    In the main project window there are three tabs: "Protocols | Data | Hosts".
    This extended tk.Frame is what will appear when Protocols is on.
    """

    RUNS_CANVAS_NAME = "runs_canvas"
    def __init__(self, parent, windows, **args):
        tk.Frame.__init__(self, parent, **args)
        # Load global configuration
        self.windows = windows
        self.project = windows.project
        self.root = windows.root
        self.getImage = windows.getImage
        self.protCfg = windows.protCfg
        self.icon = windows.icon
        self.settings = windows.getSettings()
        self.runsView = self.settings.getRunsView()
        self._loadSelection()        
        self._items = {}
        self._lastSelectedProtId = None
        self._lastStatus = None
        self.selectingArea = False

        self.style = ttk.Style()
        self.root.bind("<F5>", self.refreshRuns)
        self.root.bind("<Control-f>", self._findProtocol)
        self.root.bind("<Control-a>", self._selectAllProtocols)
        self.root.bind("<Control-t>", self._toggleColorScheme)

        # To bind key press to mehtods
        # Listen to any key: send event to keyPressed method
        self.root.bind("<Key>", self.keyPressed)
        self.keybinds = dict()

        # Register key binds
        self._bindKeyPress('Delete', self._onDelPressed)

        self.__autoRefresh = None
        self.__autoRefreshCounter = 3 # start by 3 secs  

        c = self.createContent()
        pwgui.configureWeigths(self)
        c.grid(row=0, column=0, sticky='news')

    def _bindKeyPress(self, key, method):

        self.keybinds[key] = method

    def keyPressed(self, event):

        if self.keybinds.has_key(event.keysym):

            method = self.keybinds[event.keysym]

            method()


    def createContent(self):
        """ Create the Protocols View for the Project.
        It has two panes:
            Left: containing the Protocol classes tree
            Right: containing the Runs list
        """
        p = tk.PanedWindow(self, orient=tk.HORIZONTAL, bg='white')

        bgColor = Color.LIGHT_GREY_COLOR
        # Left pane, contains Protocols Pane
        leftFrame = tk.Frame(p, bg=bgColor)
        leftFrame.columnconfigure(0, weight=1)
        leftFrame.rowconfigure(1, weight=1)


        # Protocols Tree Pane
        protFrame = tk.Frame(leftFrame, width=300, height=500, bg=bgColor)
        protFrame.grid(row=1, column=0, sticky='news', padx=5, pady=5)
        protFrame.columnconfigure(0, weight=1)
        protFrame.rowconfigure(1, weight=1)
        self._createProtocolsPanel(protFrame, bgColor)
        self.updateProtocolsTree(self.protCfg)
        # Create the right Pane that will be composed by:
        # a Action Buttons TOOLBAR in the top
        # and another vertical Pane with:
        # Runs History (at Top)
        # Sectected run info (at Bottom)
        rightFrame = tk.Frame(p, bg='white')
        rightFrame.columnconfigure(0, weight=1)
        rightFrame.rowconfigure(1, weight=1)
        #rightFrame.rowconfigure(0, minsize=label.winfo_reqheight())

        # Create the Action Buttons TOOLBAR
        toolbar = tk.Frame(rightFrame, bg='white')
        toolbar.grid(row=0, column=0, sticky='news')
        pwgui.configureWeigths(toolbar)
        #toolbar.columnconfigure(0, weight=1)
        toolbar.columnconfigure(1, weight=1)

        self.runsToolbar = tk.Frame(toolbar, bg='white')
        self.runsToolbar.grid(row=0, column=0, sticky='sw')
        # On the left of the toolbar will be other
        # actions that can be applied to all runs (refresh, graph view...)
        self.allToolbar = tk.Frame(toolbar, bg='white')
        self.allToolbar.grid(row=0, column=10, sticky='se')
        self.createActionToolbar()

        # Create the Run History tree
        v = ttk.PanedWindow(rightFrame, orient=tk.VERTICAL)
        #runsFrame = ttk.Labelframe(v, text=' History ', width=500, height=500)
        runsFrame = tk.Frame(v, bg='white')
        #runsFrame.grid(row=1, column=0, sticky='news', pady=5)
        self.runsTree = self.createRunsTree(runsFrame)
        pwgui.configureWeigths(runsFrame)

        self.createRunsGraph(runsFrame)

        if self.runsView == VIEW_LIST:
            treeWidget = self.runsTree
        else:
            treeWidget = self.runsGraphCanvas

        treeWidget.grid(row=0, column=0, sticky='news')

        # Create the Selected Run Info
        infoFrame = tk.Frame(v)
        infoFrame.columnconfigure(0, weight=1)
        infoFrame.rowconfigure(1, weight=1)
        # Create the Analyze results button
        btnAnalyze = pwgui.Button(infoFrame, text=Message.LABEL_ANALYZE,
                                  fg='white', bg=Color.RED_COLOR,
                                  image=self.getImage(Icon.ACTION_VISUALIZE),
                                  compound=tk.LEFT,
                                  activeforeground='white',
                                  activebackground='#A60C0C',
                                  command=self._analyzeResultsClicked)
        btnAnalyze.grid(row=0, column=0, sticky='ne', padx=15)
        #self.style.configure("W.TNotebook")#, background='white')
        tab = ttk.Notebook(infoFrame)#, style='W.TNotebook')

        # Summary tab
        dframe = tk.Frame(tab, bg='white')
        pwgui.configureWeigths(dframe, row=0)
        pwgui.configureWeigths(dframe, row=2)
        provider = RunIOTreeProvider(self, self.getSelectedProtocol(),
                                     self.project.mapper)
        self.style.configure("NoBorder.Treeview", background='white',
                             borderwidth=0, font=self.windows.font)
        self.infoTree = pwgui.browser.BoundTree(dframe, provider, height=6,
                                                show='tree',
                                                style="NoBorder.Treeview")
        self.infoTree.grid(row=0, column=0, sticky='news')
        label = tk.Label(dframe, text='SUMMARY', bg='white',
                         font=self.windows.fontBold)
        label.grid(row=1, column=0, sticky='nw', padx=(15, 0))

        hView = {'sci-open': self._viewObject}
        self.summaryText = pwgui.text.TaggedText(dframe, width=40, height=5,
                                                 bg='white', bd=0,
                                                 font=self.windows.font,
                                                 handlers=hView)
        self.summaryText.grid(row=2, column=0, sticky='news', padx=(30, 0))

        # Method tab
        mframe = tk.Frame(tab)
        pwgui.configureWeigths(mframe)
        self.methodText = pwgui.text.TaggedText(mframe, width=40, height=15,
                                                bg='white', handlers=hView)
        self.methodText.grid(row=0, column=0, sticky='news')

        #Logs
        ologframe = tk.Frame(tab)
        pwgui.configureWeigths(ologframe)
        self.outputViewer = pwgui.text.TextFileViewer(ologframe, allowOpen=True,
                                                      font=self.windows.font)
        self.outputViewer.grid(row=0, column=0, sticky='news')
        self.outputViewer.windows = self.windows

        self._updateSelection()

        # Add all tabs
        tab.add(dframe, text=Message.LABEL_SUMMARY)
        tab.add(mframe, text=Message.LABEL_METHODS)
        tab.add(ologframe, text=Message.LABEL_LOGS_OUTPUT)
#         tab.add(elogframe, text=Message.LABEL_LOGS_ERROR)
#         tab.add(slogframe, text=Message.LABEL_LOGS_SCIPION)
        tab.grid(row=1, column=0, sticky='news')

        v.add(runsFrame, weight=3)
        v.add(infoFrame, weight=1)
        v.grid(row=1, column=0, sticky='news')

        # Add sub-windows to PanedWindows
        p.add(leftFrame, padx=5, pady=5, sticky='news')
        p.add(rightFrame, padx=5, pady=5)
        p.paneconfig(leftFrame, minsize=300)
        p.paneconfig(rightFrame, minsize=400)

        return p

    def _viewObject(self, objId):
        """ Call appropriate viewer for objId. """
        obj = self.project.getObject(int(objId))
        viewerClasses = em.findViewers(obj.getClassName(), DESKTOP_TKINTER)
        if not viewerClasses:
            return  # TODO: protest nicely
        viewer = viewerClasses[0](project=self.project, parent=self.windows)
        viewer.visualize(obj)

    def _loadSelection(self):
        """ Load selected items, remove if some do not exists. """
        self._selection = self.settings.runSelection
        for protId in list(self._selection):
            try:
                self.project.getProtocol(protId)
            except Exception:
                self._selection.remove(protId)

    def refreshRuns(self, e=None, initRefreshCounter=True, checkPids=False):
        """ Refresh the status of displayed runs.
         Params:
            e: Tk event input
            initRefreshCounter: if True the refresh counter will be set to 3 secs
             then only case when False is from _automaticRefreshRuns where the
             refresh time is doubled each time to avoid refreshing too often.
        """
        if pwutils.envVarOn('SCIPION_DEBUG'):
            import psutil
            proc = psutil.Process(os.getpid())
            mem = psutil.virtual_memory()
            print "------------- refreshing ---------- "
            files = proc.get_open_files()
            print "  open files: ", len(files)
            for f in files:
                print "    - %s, %s" % (f.path, f.fd)
            print "  memory percent: ", proc.get_memory_percent()

        self.updateRunsGraph(True, checkPids=checkPids)
        self.updateRunsTree(False)


        if initRefreshCounter:
            self.__autoRefreshCounter = 3 # start by 3 secs
            if self.__autoRefresh:
                self.runsTree.after_cancel(self.__autoRefresh)
                self.__autoRefresh = self.runsTree.after(self.__autoRefreshCounter*1000,
                                                         self._automaticRefreshRuns)

        
    def _automaticRefreshRuns(self, e=None):
        """ Schedule automatic refresh increasing the time between refreshes. """
        self.refreshRuns(initRefreshCounter=False)
        secs = self.__autoRefreshCounter
        # double the number of seconds up to 30 min
        self.__autoRefreshCounter = min(2*secs, 1800)
        self.__autoRefresh = self.runsTree.after(secs*1000,
                                                 self._automaticRefreshRuns)
                
    def _findProtocol(self, e=None):
        """ Find a desired protocol by typing some keyword. """
        window = SearchProtocolWindow(self.windows)
        window.show()         
        
    def createActionToolbar(self):
        """ Prepare the buttons that will be available for protocol actions. """
       
        self.actionList = [ACTION_EDIT, ACTION_COPY, ACTION_DELETE, 
                           ACTION_STEPS, ACTION_BROWSE, ACTION_DB,
                           ACTION_STOP, ACTION_CONTINUE, ACTION_RESULTS, 
                           ACTION_EXPORT, ACTION_COLLAPSE, ACTION_EXPAND,
                           ACTION_LABELS]
        self.actionButtons = {}
        
        def addButton(action, text, toolbar):
            btn = tk.Label(toolbar, text=text,
                           image=self.getImage(ActionIcons.get(action, None)),
                           compound=tk.LEFT, cursor='hand2', bg='white')
            btn.bind('<Button-1>', lambda e: self._runActionClicked(action))
            return btn
        
        for action in self.actionList:
            self.actionButtons[action] = addButton(action, action,
                                                   self.runsToolbar)
            
        ActionIcons[ACTION_TREE] = RUNS_TREE
            
        self.viewButtons = {}
        
        # Add combo for switch between views
        viewFrame = tk.Frame(self.allToolbar)
        viewFrame.grid(row=0, column=0)
        self._createViewCombo(viewFrame)
        
        # Add refresh Tree button
        btn = addButton(ACTION_TREE, "  ", self.allToolbar)
        pwgui.tooltip.ToolTip(btn, "Re-organize the node positions.", 1500)
        self.viewButtons[ACTION_TREE] = btn       
        if self.runsView != VIEW_LIST:
            btn.grid(row=0, column=1)
        
        # Add refresh button
        btn = addButton(ACTION_REFRESH, ACTION_REFRESH, self.allToolbar)
        btn.grid(row=0, column=2)
        self.viewButtons[ACTION_REFRESH] = btn
            
    def _createViewCombo(self, parent):
        """ Create the select-view combobox. """
        label = tk.Label(parent, text='View:', bg='white')
        label.grid(row=0, column=0)
        viewChoices = ['List', 'Tree', 'Tree - small']
        self.switchCombo = pwgui.widgets.ComboBox(parent, width=10, 
                                    choices=viewChoices, 
                                    values=[VIEW_LIST, VIEW_TREE, VIEW_TREE_SMALL],
                                    initial=viewChoices[self.runsView],
                                    onChange=lambda e: self._runActionClicked(ACTION_SWITCH_VIEW))
        self.switchCombo.grid(row=0, column=1)
            
    def _updateActionToolbar(self):
        """ Update which action buttons should be visible. """
        
        def displayAction(action, i, cond=True):
            """ Show/hide the action button if the condition is met. """

            # If action present (set color is not in the toolbar but in the context menu)
            if self.actionButtons.has_key(action):
                if cond:
                    self.actionButtons[action].grid(row=0, column=i, sticky='sw',
                                                    padx=(0, 5), ipadx=0)
                else:
                    self.actionButtons[action].grid_remove()
            else:
                # print action + " not in toolbar."
                pass

        for i, actionTuple in enumerate(self.provider.getActionsFromSelection()):
            action, cond = actionTuple
            displayAction(action, i, cond)          
        
    def _createProtocolsTree(self, parent, background=Color.LIGHT_GREY_COLOR):
        self.style.configure("W.Treeview", background=background, borderwidth=0,
                             fieldbackground=background)
        t = pwgui.tree.Tree(parent, show='tree', style='W.Treeview')
        t.column('#0', minwidth=300)
        t.tag_configure('protocol', image=self.getImage('python_file.gif'))
        t.tag_bind('protocol', '<Double-1>', self._protocolItemClick)
        t.tag_bind('protocol', '<Return>', self._protocolItemClick)
        t.tag_configure('protocol_base', image=self.getImage('class_obj.gif'))
        t.tag_configure('protocol_group', image=self.getImage('class_obj.gif'))
        t.tag_configure('section', font=self.windows.fontBold)
        return t
        
    def _createProtocolsPanel(self, parent, bgColor):
        """Create the protocols Tree displayed in left panel"""
        comboFrame = tk.Frame(parent, bg=bgColor)
        tk.Label(comboFrame, text='View', bg=bgColor).grid(row=0, column=0,
                                                           padx=(0, 5), pady=5)
        choices = self.project.getProtocolViews() 
        initialChoice = self.settings.getProtocolView()
        combo = pwgui.widgets.ComboBox(comboFrame, choices=choices,
                                       initial=initialChoice)
        combo.setChangeCallback(self._onSelectProtocols)
        combo.grid(row=0, column=1)
        comboFrame.grid(row=0, column=0, padx=5, pady=5, sticky='nw')
        
        t = self._createProtocolsTree(parent)
        t.grid(row=1, column=0, sticky='news')
        # Program automatic refresh
        t.after(3000, self._automaticRefreshRuns)
        self.protTree = t

    def _onSelectProtocols(self, combo):
        """ This function will be called when a protocol menu
        is selected. The index of the new menu is passed. 
        """
        protView = combo.getText()
        self.settings.setProtocolView(protView)
        self.protCfg = self.project.getCurrentProtocolView()
        self.updateProtocolsTree(self.protCfg)
                
    def updateProtocolsTree(self, protCfg):
        self.protCfg = protCfg
        self.protTree.clear()
        self.protTree.unbind('<<TreeviewOpen>>')
        self.protTree.unbind('<<TreeviewClose>>')
        self.protTreeItems = {}
        subclassedDict = {} # Check which classes serve as base to not show them
        emProtocolsDict = em.getProtocols()
        for _, v1 in emProtocolsDict.iteritems():
            for k2, v2 in emProtocolsDict.iteritems():
                if v1 is not v2 and issubclass(v1, v2):
                    subclassedDict[k2] = True
        populateTree(self, self.protTree, self.protTreeItems, '', self.protCfg,
                     subclassedDict)
        self.protTree.bind('<<TreeviewOpen>>',
                           lambda e: self._treeViewItemChange(True))
        self.protTree.bind('<<TreeviewClose>>',
                           lambda e: self._treeViewItemChange(False))
        
    def _treeViewItemChange(self, openItem):
        item = self.protTree.focus()
        if item in self.protTreeItems:
            self.protTreeItems[item].openItem.set(openItem)
        
    def createRunsTree(self, parent):
        self.provider = RunsTreeProvider(self.project, self._runActionClicked)

        # This line triggers the getRuns for the first time.
        # Ne need to force the check pids here, temporary
        self.provider._checkPids = True
        t = pwgui.tree.BoundTree(parent, self.provider)
        self.provider._checkPids = False

        t.itemDoubleClick = self._runItemDoubleClick
        t.itemClick = self._runTreeItemClick
            
        return t
   
    def updateRunsTree(self, refresh=False):
        self.provider.setRefresh(refresh)
        self.runsTree.update()
        self.updateRunsTreeSelection()

    def updateRunsTreeSelection(self):
        for prot in self._iterSelectedProtocols():
            treeId = self.provider.getObjectFromId(prot.getObjId())._treeId
            self.runsTree.selection_add(treeId)
    
    def createRunsGraph(self, parent):
        self.runsGraphCanvas = pwgui.Canvas(parent, width=400, height=400, 
                                tooltipCallback=self._runItemTooltip,
                                tooltipDelay=1000, name=ProtocolsView.RUNS_CANVAS_NAME, takefocus=True, highlightthickness=0)

        self.runsGraphCanvas.onClickCallback = self._runItemClick
        self.runsGraphCanvas.onDoubleClickCallback = self._runItemDoubleClick
        self.runsGraphCanvas.onRightClickCallback = self._runItemRightClick
        self.runsGraphCanvas.onControlClickCallback = self._runItemControlClick
        self.runsGraphCanvas.onAreaSelected = self._selectItemsWithinArea

        parent.grid_columnconfigure(0, weight=1)
        parent.grid_rowconfigure(0, weight=1)
        
        self.settings.getNodes().updateDict()
        self.settings.getLabels().updateDict()

        self.updateRunsGraph()

    def updateRunsGraph(self, refresh=False, reorganize=False, checkPids=False):

        self.runsGraph = self.project.getRunsGraph(refresh=refresh,
                                                   checkPids=checkPids)
        self.runsGraphCanvas.clear()
        
        # Check if there are positions stored
        if reorganize or len(self.settings.getNodes()) == 0:
            # Create layout to arrange nodes as a level tree
            layout = pwgui.graph.LevelTreeLayout()
        else:
            layout = pwgui.graph.BasicLayout()
            
        # Create empty nodeInfo for new runs
        for node in self.runsGraph.getNodes():
            nodeId = node.run.getObjId() if node.run else 0
            nodeInfo = self.settings.getNodeById(nodeId)
            if nodeInfo is None:
                self.settings.addNode(nodeId, x=0, y=0, expanded=True) 
            
        self.runsGraphCanvas.drawGraph(self.runsGraph, layout,
                                       drawNode=self.createRunItem)
        
    def createRunItem(self, canvas, node):

        nodeId = node.run.getObjId() if node.run else 0
        nodeInfo = self.settings.getNodeById(nodeId)

        # Extend attributes: use some from nodeInfo
        node.expanded = nodeInfo.isExpanded()
        node.x, node.y = nodeInfo.getPosition()
        nodeText = self._getNodeText(node)

        # Get status color
        statusColor = self._getStatusColor(node)

        # Get the box color (depends on color mode: label or status)
        boxColor = self._getBoxColor(nodeInfo, statusColor, node)

        # Draw the box
        item = RunBox(nodeInfo, self.runsGraphCanvas,
                      nodeText, node.x, node.y,
                      bgColor=boxColor, textColor='black')
        # No border
        item.margin = 0

        # Paint the oval..if apply.
        self._paintOval(item, statusColor)

        # Paint the bottom line (for now only labels are painted there).
        self._paintBottomLine(item)

        if nodeId in self._selection:
            item.setSelected(True)
        return item

    def _getBoxColor(self, nodeInfo, statusColor, node):

        # If the color has to go to the box
        if self.settings.statusColorMode():
            boxColor = statusColor

        elif self.settings.ageColorMode():

            if node.run:

                # Get the latest activity timestamp
                ts = node.run.getObjCreation().datetime(fs=False)

                if node.run.initTime.hasValue():
                    ts = node.run.initTime.datetime()

                if node.run.endTime.hasValue():
                    ts = node.run.endTime.datetime()

                age = dt.datetime.now() - ts

                boxColor = self._ageColor('#6666ff', age.days * 24)
            else:
                boxColor = DEFAULT_BOX_COLOR

        # ... box is for the labels.
        elif self.settings.labelsColorMode():
            # If there is only one label use the box for the color.
            if self._getLabelsCount(nodeInfo) == 1:

                labelId = nodeInfo.getLabels()[0]
                label = self.settings.getLabels().getLabel(labelId)

                # If there is no label (it has been deleted)
                if label is None:
                    nodeInfo.getLabels().remove(labelId)
                    boxColor = DEFAULT_BOX_COLOR
                else:
                    boxColor = label.getColor()

            else:
                boxColor = DEFAULT_BOX_COLOR
        else:
            boxColor = DEFAULT_BOX_COLOR

        return boxColor

    def _ageColor(self, rgbColor, ageInDays):

        if ageInDays > 255: ageInDays = 255

        # Mek it a percentage: 1 = 100% white, 0 = same rgbColor
        ageInDays /= 255.

        return pwutils.rgb_to_hex(pwutils.lighter(pwutils.hex_to_rgb(rgbColor),
                                                  ageInDays))


    @staticmethod
    def _getLabelsCount(nodeInfo):

        return 0 if nodeInfo.getLabels() is None else len(nodeInfo.getLabels())

    def _paintBottomLine(self, item):

        if self.settings.labelsColorMode():

            self._addLabels(item)

    def _paintOval(self, item, statusColor):
        # Show the status as a circle in the top right corner
        if not self.settings.statusColorMode():
            # Option: Status item.
            (topLeftX, topLeftY, bottomRightX, bottomRightY) = self.runsGraphCanvas.bbox(item.id)
            statusSize = 10
            statusX = bottomRightX - (statusSize + 3)
            statusY = topLeftY + 3

            pwgui.Oval(self.runsGraphCanvas, statusX, statusY, statusSize,
                       color=statusColor, anchor=item)

        # in statusColorMode
        else:
            # Show a black circle if there is any label
            if self._getLabelsCount(item.nodeInfo) > 0:

                (topLeftX, topLeftY, bottomRightX, bottomRightY) = self.runsGraphCanvas.bbox(item.id)
                statusSize = 10
                statusX = bottomRightX - (statusSize + 3)
                statusY = topLeftY + 3

                pwgui.Oval(self.runsGraphCanvas, statusX, statusY, statusSize,
                           color='black', anchor=item)

    @staticmethod
    def _getStatusColor(node):

        # If it is a run node (not PROJECT)
        if node.run:
            status = node.run.status.get(pwprot.STATUS_FAILED)
            return STATUS_COLORS[status]
        else:
            return '#ADD8E6'  # Lightblue

    def _getNodeText(self, node):

        nodeText = node.getLabel()

        if node.run:

            if node.expanded:
                expandedStr = ''
            else:
                expandedStr = ' (+)'
            if self.runsView == VIEW_TREE_SMALL:
                nodeText = node.getName() + expandedStr
            else:
                nodeText = nodeText + expandedStr + '\n' + node.run.getStatusMessage()

        return nodeText

    def _addLabels(self, item):

        # If there is only one label it should be already used in the box color.
        if self._getLabelsCount(item.nodeInfo) < 2: return

        # Get the positions of the box
        (topLeftX, topLeftY, bottomRightX, bottomRightY) = self.runsGraphCanvas.bbox(item.id)

        # Get the width of the box
        boxWidth = bottomRightX - topLeftX

        # Set the size
        marginV = 3
        marginH = 2
        labelWidth = (boxWidth - (2 * marginH)) / len(item.nodeInfo.getLabels())
        labelHeight = 6

        # Leave some margin on the right and bottom
        labelX = bottomRightX - marginH
        labelY = bottomRightY - (labelHeight + marginV)

        for index, labelId in enumerate(item.nodeInfo.getLabels()):

            # Get the label
            label = self.settings.getLabels().getLabel(labelId)

            # If not none
            if label is not None:
                # Move X one label to the left
                if index == len(item.nodeInfo.getLabels()) - 1:
                    labelX = topLeftX + marginH
                else:
                    labelX -= labelWidth

                pwgui.Rectangle(self.runsGraphCanvas, labelX, labelY,
                                labelWidth, labelHeight, color=label.getColor(),
                                anchor=item)
            else:

                item.nodeInfo.getLabels().remove(labelId)

    def switchRunsView(self):
        previousView = self.runsView
        viewValue = self.switchCombo.getValue()
        self.runsView = viewValue
        self.settings.setRunsView(viewValue)
        
        if viewValue == VIEW_LIST:
            self.runsTree.grid(row=0, column=0, sticky='news')
            self.runsGraphCanvas.frame.grid_remove()
            self.updateRunsTreeSelection()
            self.viewButtons[ACTION_TREE].grid_remove()
        else:            
            self.runsTree.grid_remove()
            self.updateRunsGraph(reorganize=(previousView!=VIEW_LIST))
            self.runsGraphCanvas.frame.grid(row=0, column=0, sticky='news')
            self.viewButtons[ACTION_TREE].grid(row=0, column=1)

    def _protocolItemClick(self, e=None):
        # Get the tree widget that originated the event
        # it could be the left panel protocols tree or just
        # the search protocol dialog tree
        tree = e.widget
        protClassName = tree.getFirst().split('.')[-1]
        protClass = em.getProtocols().get(protClassName)
        prot = self.project.newProtocol(protClass)
        self._openProtocolForm(prot)

    def _toggleColorScheme(self, e=None):

        currentMode = self.settings.getColorMode()

        if currentMode >= len(self.settings.COLOR_MODES) - 1:
            currentMode = -1

        nextColorMode = currentMode + 1

        self.settings.setColorMode(nextColorMode)
        self._updateActionToolbar()
        self.updateRunsGraph()

    def _selectAllProtocols(self, e=None):
        self._selection.clear()
        for prot in self.project.getRuns():
            self._selection.append(prot.getObjId())
        self._updateSelection()
        self.updateRunsGraph()

    def _deleteSelectedProtocols(self, e=None):

        for selection in self._selection:
            self.project.getProtocol(self._selection[0])


        self._updateSelection()
        self.updateRunsGraph()


    def _updateSelection(self):
        self._fillSummary()
        self._fillMethod()
        self._fillLogs()

        last = self.getSelectedProtocol()
        self._lastSelectedProtId = last.getObjId() if last else None

        self._updateActionToolbar()
        
    def _runTreeItemClick(self, item=None):
        self._selection.clear()
        for prot in self.runsTree.iterSelectedObjects():
            self._selection.append(prot.getObjId())
        self._updateSelection()
                         
    def _selectItemProtocol(self, prot):
        """ Call this function when a new box (item) of a protocol
        is selected. It should be called either from itemClick
        or itemRightClick
        """
        self._selection.clear()
        self.settings.dataSelection.clear()
        self._selection.append(prot.getObjId())

        # Select output data too
        self.toggleDataSelection(prot,True)

        self._updateSelection()
        self.runsGraphCanvas.update_idletasks()
        
    def _deselectItems(self, item):
        """ Deselect all items except the item one
        """
        g = self.project.getRunsGraph(refresh=False)
        
        for node in g.getNodes():
            if node.run and node.run.getObjId() in self._selection:
                # This option is only for compatibility with all projects
                if hasattr(node, 'item'):
                    node.item.setSelected(False)
        item.setSelected(True)
        
    def _runItemClick(self, item=None):

        self.runsGraphCanvas.focus_set()

        # Get last selected item for tree or graph
        if self.runsView == VIEW_LIST:
            prot = self.project.mapper.selectById(int(self.runsTree.getFirst()))
        else:
            prot = item.node.run
            if prot is None:  # in case it is the main "Project" node
                return
            self._deselectItems(item)
        self._selectItemProtocol(prot)
        
    def _runItemDoubleClick(self, e=None):
        self._runActionClicked(ACTION_EDIT)
        
    def _runItemRightClick(self, item=None):
        prot = item.node.run
        if prot is None:  # in case it is the main "Project" node
            return
        n = len(self._selection)
        # Only select item with right-click if there is a single
        # item selection, not for multiple selection
        if n <= 1:
            self._deselectItems(item)
            self._selectItemProtocol(prot)
        return self.provider.getObjectActions(prot)
        
    def _runItemControlClick(self, item=None):
        # Get last selected item for tree or graph
        if self.runsView == VIEW_LIST:
            prot = self.project.mapper.selectById(int(self.runsTree.getFirst()))        
        else:
            prot = item.node.run
            protId = prot.getObjId()
            if protId in self._selection:
                item.setSelected(False)
                self._selection.remove(protId)

                # Remove data selected
                self.toggleDataSelection(prot, False)
            else:

                item.setSelected(True)
                self._selection.append(prot.getObjId())

                # Select output data too
                self.toggleDataSelection(prot, True)

        self._updateSelection()


    def toggleDataSelection(self, prot, append):

        # Go through the data selection
        for paramName, output in prot.iterOutputEM():
            if append:
                self.settings.dataSelection.append(output.getObjId())
            else:
                self.settings.dataSelection.remove(output.getObjId())

    def _runItemTooltip(self, tw, item):
        """ Create the contents of the tooltip to be displayed
        for the given item.
        Params:
            tw: a tk.TopLevel instance (ToolTipWindow)
            item: the selected item.
        """
        prot = item.node.run
        
        if prot:
            tm = '*%s*\n' % prot.getRunName()
            tm += '   Id: %s\n' % prot.getObjId()
            tm += 'State: %s\n' % prot.getStatusMessage()
            tm += ' Time: %s\n' % pwutils.prettyDelta(prot.getElapsedTime()) 
            if not hasattr(tw, 'tooltipText'):
                frame = tk.Frame(tw)
                frame.grid(row=0, column=0)
                tw.tooltipText = pwgui.dialog.createMessageBody(frame, tm, None,
                                                textPad=0,
                                                textBg=Color.LIGHT_GREY_COLOR_2)
                tw.tooltipText.config(bd=1, relief=tk.RAISED)
            else:
                pwgui.dialog.fillMessageText(tw.tooltipText, tm)

    def _selectItemsWithinArea(self, x1, y1, x2, y2, enclosed=False):
        """
        Parameters
        ----------
        x1: x coordinate of first corner of the area
        y1: y coordinate of first corner of the area
        x2: x coordinate of second corner of the area
        y2: y coordinate of second corner of the area
        enclosed: Default True. Returns enclosed items,
                  overlapping items otherwise.
        Returns
        -------
        Nothing

        """

        # NOT working properly: Commented for the moment.
        return

        if enclosed:
            items = self.runsGraphCanvas.find_enclosed(x1, y1, x2, y2)
        else:
            items = self.runsGraphCanvas.find_overlapping(x1, y1, x2, y2)

        update = False

        for itemId in items:
            if itemId in self.runsGraphCanvas.items:

                item = self.runsGraphCanvas.items[itemId]
                if not item.node.isRoot():
                    item.setSelected(True)
                    self._selection.append(itemId)
                    update = True

        if update is not None: self._updateSelection()
        
    def _openProtocolForm(self, prot):
        """Open the Protocol GUI Form given a Protocol instance"""
        
        w = pwgui.form.FormWindow(Message.TITLE_NAME_RUN + prot.getClassName(),
                                  prot, self._executeSaveProtocol, self.windows,
                                  hostList=self.project.getHostNames(),
                                  updateProtocolCallback=self._updateProtocol(prot))
        w.adjustSize()
        w.show(center=True)

    def _browseSteps(self):
        """ Open a new window with the steps list. """
        window = StepsWindow(Message.TITLE_BROWSE_DATA, self.windows, 
                             self.getSelectedProtocol(), icon=self.icon)
        window.show()        
    
    def _browseRunData(self):
        provider = ProtocolTreeProvider(self.getSelectedProtocol())
        window = pwgui.browser.BrowserWindow(Message.TITLE_BROWSE_DATA,
                                             self.windows, icon=self.icon)
        window.setBrowser(pwgui.browser.ObjectBrowser(window.root, provider))
        window.itemConfig(self.getSelectedProtocol(), open=True)  
        window.show()
        
    def _browseRunDirectory(self):
        """ Open a file browser to inspect the files generated by the run. """
        protocol = self.getSelectedProtocol()
        workingDir = protocol.getWorkingDir()
        if os.path.exists(workingDir):
            window = pwgui.browser.FileBrowserWindow("Browsing: " + workingDir, 
                                       master=self.windows, 
                                       path=workingDir)
            window.show()
        else:
            self.windows.showInfo("Protocol working dir does not exists: \n %s"
                                  % workingDir)
        
    def _iterSelectedProtocols(self):
        for protId in sorted(self._selection):
            prot = self.project.getProtocol(protId)
            if prot:
                yield prot
            
    def _getSelectedProtocols(self):
        return [prot for prot in self._iterSelectedProtocols()]


    def _iterSelectedNodes(self):

        for protId in sorted(self._selection):

            node = self.settings.getNodeById(protId)

            yield node

    def _getSelectedNodes(self):
        return [node for node in self._iterSelectedNodes()]


    def getSelectedProtocol(self):
        if self._selection:
            return self.project.getProtocol(self._selection[0])
        return None
            
    def _fillSummary(self):
        self.summaryText.setReadOnly(False)
        self.summaryText.clear()
        self.infoTree.clear()
        n = len(self._selection)
        
        if n == 1:
            prot = self.getSelectedProtocol()

            if prot:
                provider = RunIOTreeProvider(self, prot, self.project.mapper)
                self.infoTree.setProvider(provider)
                self.infoTree.grid(row=0, column=0, sticky='news')
                self.infoTree.update_idletasks()
                # Update summary
                self.summaryText.addText(prot.summary())
            else:
                self.infoTree.clear()
        
        elif n > 1:
            self.infoTree.clear()
            for prot in self._iterSelectedProtocols():
                self.summaryText.addLine('> _%s_' % prot.getRunName())
                for line in prot.summary():
                    self.summaryText.addLine(line)
                self.summaryText.addLine('')
        self.summaryText.setReadOnly(True)
        
    def _fillMethod(self):
        self.methodText.setReadOnly(False)
        self.methodText.clear()
        self.methodText.addLine("*METHODS:*")
        cites = OrderedDict()
        
        for prot in self._iterSelectedProtocols():
            self.methodText.addLine('> _%s_' % prot.getRunName())
            for line in prot.getParsedMethods():
                self.methodText.addLine(line)
            cites.update(prot.getCitations())
            cites.update(prot.getPackageCitations())
            self.methodText.addLine('')
        
        if cites:
            self.methodText.addLine('*REFERENCES:*')
            for cite in cites.values():
                self.methodText.addLine(cite)
        
        self.methodText.setReadOnly(True)
        
    def _fillLogs(self):
        prot = self.getSelectedProtocol()

        if len(self._selection) != 1 or not prot:
            self.outputViewer.clear()
            self._lastStatus = None
        elif prot.getObjId() != self._lastSelectedProtId:
            self._lastStatus = prot.getStatus()
            i = self.outputViewer.getIndex()
            self.outputViewer.clear()
            # Right now skip the err tab since we are redirecting
            # stderr to stdout
            out, _, log = prot.getLogPaths()
            self.outputViewer.addFile(out)
            self.outputViewer.addFile(log)
            self.outputViewer.setIndex(i) # Preserve the last selected tab
            self.outputViewer.selectedText().goEnd()
            # when there are not logs, force re-load next time
            if (not os.path.exists(out) or
                not os.path.exists(log)):
                self._lastStatus = None
                
        elif  prot.isActive() or prot.getStatus() != self._lastStatus:
            doClear = self._lastStatus is None                
            self._lastStatus = prot.getStatus()
            self.outputViewer.refreshAll(clear=doClear, goEnd=doClear)

    def _scheduleRunsUpdate(self, secs=1):
        #self.runsTree.after(secs*1000, self.refreshRuns)
        self.windows.enqueue(self.refreshRuns)
        
    def executeProtocol(self, prot):
        """ Function to execute a protocol called not
        directly from the Form "Execute" button.
        """
        # We need to equeue the execute action
        # to be executed in the same thread
        self.windows.enqueue(lambda: self._executeSaveProtocol(prot))
        
    def _executeSaveProtocol(self, prot, onlySave=False):
        if onlySave:
            self.project.saveProtocol(prot)
            msg = Message.LABEL_SAVED_FORM
#            msg = "Protocol successfully saved."
        else:
            self.project.launchProtocol(prot)
            # Select the launched protocol to display its summary, methods..etc
            self._selection.clear()
            self._selection.append(prot.getObjId())
            self._updateSelection()
            self._lastStatus = None # clear lastStatus to force re-load the logs
            msg = ""
            
        # Update runs list display, even in save we
        # need to get the updated copy of the protocol
        self._scheduleRunsUpdate()

        return msg
    
    def _updateProtocol(self, prot):
        """ Callback to notify about the change of a protocol
        label or comment. 
        """
        self._scheduleRunsUpdate()
        
    def _continueProtocol(self, prot):
        self.project.continueProtocol(prot)
        self._scheduleRunsUpdate()

    def _onDelPressed(self):
        # This function will be connected to the key 'Del' press event
        # We need to check if the canvas have the focus and then
        # proceed with the delete action

        # get the widget with the focus
        widget = self.focus_get()

        # Call the delete action only if the widget is the canvas
        if str(widget).endswith(ProtocolsView.RUNS_CANVAS_NAME):
            self._deleteProtocol()

    def _deleteProtocol(self):
        protocols = self._getSelectedProtocols()

        if len(protocols) == 0:
            return

        protStr = '\n  - '.join(['*%s*' % p.getRunName() for p in protocols])

        if pwgui.dialog.askYesNo(Message.TITLE_DELETE_FORM,
                                 Message.LABEL_DELETE_FORM % protStr,
                                 self.root):
            self.project.deleteProtocol(*protocols)
            self._selection.clear()
            self._updateSelection()
            self._scheduleRunsUpdate()
            
    def _copyProtocols(self):
        protocols = self._getSelectedProtocols()
        if len(protocols) == 1:
            newProt = self.project.copyProtocol(protocols[0])
            if newProt is None:
                self.windows.showError("Error copying protocol.!!!")
            else:
                self._openProtocolForm(newProt)
        else:
            self.project.copyProtocol(protocols)
            self.refreshRuns()

    def _selectLabels(self):
        selectedNodes = self._getSelectedNodes()

        if selectedNodes:
            dlg = self.windows.manageLabels()

            if dlg.resultYes():
                for node in selectedNodes:
                    node.setLabels([label.getName() for label in dlg.values])

                self.updateRunsGraph()

    def _exportProtocols(self):
        protocols = self._getSelectedProtocols()
        
        def _export(obj):
            filename = os.path.join(browser.getCurrentDir(), 
                                    browser.getEntryValue())
            try:
                self.project.exportProtocols(protocols, filename)
                self.windows.showInfo("Workflow successfully saved to '%s' "
                                      % filename)
            except Exception, ex:
                self.windows.showError(str(ex))
            
        browser = pwgui.browser.FileBrowserWindow("Choose .json file to save workflow", 
                                    master=self.windows, 
                                    path=self.project.getPath(''), 
                                    onSelect=_export,
                                    entryLabel='File', entryValue='workflow.json')
        browser.show()
            
    def _stopProtocol(self, prot):
        if pwgui.dialog.askYesNo(Message.TITLE_STOP_FORM, Message.LABEL_STOP_FORM, self.root):
            self.project.stopProtocol(prot)
            self._lastStatus = None # force logs to re-load
            self._scheduleRunsUpdate()

    def _analyzeResults(self, prot):        

        viewers = em.findViewers(prot.getClassName(), DESKTOP_TKINTER)
        if len(viewers):
            # Instanciate the first available viewer
            # TODO: If there are more than one viewer we should display
            # TODO: a selection menu
            firstViewer = viewers[0](project=self.project, protocol=prot,
                                     parent=self.windows)

            if isinstance(firstViewer, ProtocolViewer):
                firstViewer.visualize(prot, windows=self.windows)
            else:
                firstViewer.visualize(prot)
        else:
            for _, output in prot.iterOutputAttributes(em.EMObject):
                viewers = em.findViewers(output.getClassName(), DESKTOP_TKINTER)
                if len(viewers):
                    # Instanciate the first available viewer
                    # TODO: If there are more than one viewer we should display
                    # TODO: a selection menu
                    viewerclass = viewers[0]
                    firstViewer = viewerclass(project=self.project, 
                                              protocol=prot,
                                              parent=self.windows)
                    # FIXME:Probably o longer needed protocol on args, already provided on init
                    firstViewer.visualize(output, windows=self.windows,
                                          protocol=prot)
            
    def _analyzeResultsClicked(self, e=None):
        """ Function called when button "Analyze results" is called. """
        prot = self.getSelectedProtocol()
        if os.path.exists(prot._getPath()):
            self._analyzeResults(prot)
        else:
            self.windows.showInfo("Selected protocol hasn't been run yet.")
                
    def _runActionClicked(self, action):
        prot = self.getSelectedProtocol()
        if prot:
            try:
                if action == ACTION_DEFAULT:
                    pass
                elif action == ACTION_EDIT:
                    self._openProtocolForm(prot)
                elif action == ACTION_COPY:
                    self._copyProtocols()
                elif action == ACTION_DELETE:
                    self._deleteProtocol()
                elif action == ACTION_STEPS:
                    self._browseSteps()                    
                elif action == ACTION_BROWSE:
                    self._browseRunDirectory()
                elif action == ACTION_DB:
                    self._browseRunData()
                elif action == ACTION_STOP:
                    self._stopProtocol(prot)
                elif action == ACTION_CONTINUE:
                    self._continueProtocol(prot)
                elif action == ACTION_RESULTS:
                    self._analyzeResults(prot)
                elif action == ACTION_EXPORT: 
                    self._exportProtocols()
                elif action == ACTION_COLLAPSE:
                    nodeInfo = self.settings.getNodeById(prot.getObjId())
                    nodeInfo.setExpanded(False)
                    self.updateRunsGraph(True, reorganize=True)
                    self._updateActionToolbar()
                elif action == ACTION_EXPAND:
                    nodeInfo = self.settings.getNodeById(prot.getObjId())
                    nodeInfo.setExpanded(True)
                    self.updateRunsGraph(True, reorganize=True)
                    self._updateActionToolbar()
                elif action == ACTION_LABELS:
                    self._selectLabels()

            except Exception, ex:
                self.windows.showError(str(ex))
                if pwutils.envVarOn('SCIPION_DEBUG'):
                    import traceback
                    traceback.print_exc()
 
        # Following actions do not need a select run
        if action == ACTION_TREE:
            self.updateRunsGraph(True, reorganize=True)
        elif action == ACTION_REFRESH:
            self.refreshRuns(checkPids=True)
        elif action == ACTION_SWITCH_VIEW:
            self.switchRunsView()
    

class RunBox(pwgui.TextBox):
    """ Just override TextBox move method to keep track of 
    position changes in the graph.
    """
    def __init__(self, nodeInfo, canvas, text, x, y, bgColor, textColor):
        pwgui.TextBox.__init__(self, canvas, text, x, y, bgColor, textColor)
        self.nodeInfo = nodeInfo
        canvas.addItem(self)
        
    def move(self, dx, dy):
        pwgui.TextBox.move(self, dx, dy)
        self.nodeInfo.setPosition(self.x, self.y)

    def moveTo(self, x, y):
        pwgui.TextBox.moveTo(self, x, y)
        self.nodeInfo.setPosition(self.x, self.y)
