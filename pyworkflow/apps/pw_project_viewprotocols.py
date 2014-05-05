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
Main project window application
"""
import Tkinter as tk
import ttk
import tkFont

from pyworkflow.gui.tree import TreeProvider, BoundTree
from pyworkflow.gui.browser import ObjectBrowser, FileBrowser
from pyworkflow.protocol.protocol import *

import pyworkflow as pw
from pyworkflow.object import *
from pyworkflow.em import *
from pyworkflow.protocol import *
from pyworkflow.protocol.params import *
from pyworkflow.mapper import SqliteMapper, XmlMapper
from pyworkflow.project import Project
from pyworkflow.utils import prettyDelta
from pyworkflow.utils.properties import Message, Icon, Color

import pyworkflow.gui as gui
from pyworkflow.gui import getImage
from pyworkflow.gui.tree import Tree, ObjectTreeProvider, DbTreeProvider, ProjectRunsTreeProvider
from pyworkflow.gui.form import FormWindow, editObject
from pyworkflow.gui.dialog import askYesNo, EditObjectDialog
from pyworkflow.gui.text import TaggedText, TextFileViewer
from pyworkflow.gui import Canvas
from pyworkflow.gui.graph import LevelTree
from pyworkflow.gui.widgets import ComboBox

from config import *
from pw_browser import BrowserWindow

ACTION_EDIT = Message.LABEL_EDIT
ACTION_COPY = Message.LABEL_COPY
ACTION_DELETE = Message.LABEL_DELETE
ACTION_REFRESH = Message.LABEL_REFRESH
ACTION_STEPS = Message.LABEL_STEPS
ACTION_BROWSE = Message.LABEL_BROWSE
ACTION_TREE = Message.LABEL_TREE
ACTION_LIST = Message.LABEL_LIST
ACTION_STOP = Message.LABEL_STOP
ACTION_DEFAULT = Message.LABEL_DEFAULT
ACTION_CONTINUE = Message.LABEL_CONTINUE
ACTION_RESULTS = Message.LABEL_ANALYZE

RUNS_TREE = Icon.RUNS_TREE
RUNS_LIST = Icon.RUNS_LIST
 
ActionIcons = {
    ACTION_EDIT: Icon.ACTION_EDIT , 
    ACTION_COPY: Icon.ACTION_COPY ,
    ACTION_DELETE:  Icon.ACTION_DELETE,
    ACTION_REFRESH:  Icon.ACTION_REFRESH,
    ACTION_STEPS:  Icon.ACTION_STEPS,
    ACTION_BROWSE:  Icon.ACTION_BROWSE,
    ACTION_TREE:  None, # should be set
    ACTION_LIST:  Icon.ACTION_LIST,
    ACTION_STOP: Icon.ACTION_STOP,
    ACTION_CONTINUE: Icon.ACTION_CONTINUE,
    ACTION_RESULTS: Icon.ACTION_RESULTS,
               }

STATUS_COLORS = {
               STATUS_SAVED: '#D9F1FA', 
               STATUS_LAUNCHED: '#D9F1FA', 
               STATUS_RUNNING: '#FCCE62', 
               STATUS_FINISHED: '#D2F5CB', 
               STATUS_FAILED: '#F5CCCB', 
               STATUS_WAITING_APPROVAL: '#F3F5CB',
               STATUS_ABORTED: '#F5CCCB',
               #STATUS_SAVED: '#124EB0',
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
        item = tree.insert(prefix, 'end', key, text=text, image=img, tags=(tag))
        treeItems[item] = obj
        # Check if the attribute should be open or close
        openItem = obj.getAttributeValue('openItem', level < 2)
        if openItem:
            tree.item(item, open=True)
            
        if obj.value.hasValue() and tag == 'protocol_base':
            protClassName = value.split('.')[-1] # Take last part
            prot = emProtocolsDict.get(protClassName, None)
            if prot is not None:
                tree.item(item, image=self.getImage('class_obj.gif'))
                
                for k, v in emProtocolsDict.iteritems():
                    if not k in subclassedDict and not v is prot and issubclass(v, prot):# and Protocol.hasDefinition(v):
                        key = '%s.%s' % (item, k)
                        t = v.getClassLabel()
                        tree.insert(item, 'end', key, text=t, tags=('protocol'))
                        
            else:
                raise Exception("Class '%s' not found" % obj.value.get())
    else:
        key = prefix
    
    for sub in obj:
        populateTree(self, tree, treeItems, key, sub, subclassedDict, level+1)
    

class RunsTreeProvider(ProjectRunsTreeProvider):
    """Provide runs info to populate tree"""
    def __init__(self, project, actionFunc):
        ProjectRunsTreeProvider.__init__(self, project)
        self.actionFunc = actionFunc
        self._selection = project.getSettings().runSelection
    
    def getActionsFromSelection(self):
        """ Return the list of options available for selection. """
        n = len(self._selection)
        single = n == 1
        if n:
            prot = self.project.getProtocol(self._selection[0]) 
            status = prot.getStatus()
        else:
            status = None
        
        return [   (ACTION_EDIT, single),
                   (ACTION_COPY, True), 
                   (ACTION_DELETE, status != STATUS_RUNNING),
                   (ACTION_STEPS, single),
                   (ACTION_BROWSE, single),
                   (ACTION_STOP, status == STATUS_RUNNING and single),
                   (ACTION_CONTINUE, status == STATUS_WAITING_APPROVAL and single)
                   ]
        
    def getObjectActions(self, obj):
        
        def addAction(a):
            if a:
                text = a
                action = a
                a = (text, lambda: self.actionFunc(action), ActionIcons[action])
            return a
            
        actions = [addAction(a) for a, cond in self.getActionsFromSelection() if cond]
        
        return actions 
    
    
class ProtocolTreeProvider(ObjectTreeProvider):
    """Create the tree elements for a Protocol run"""
    def __init__(self, protocol):
        self.protocol = protocol
        # This list is create to group the protocol parameters
        # in the tree display
        self.status = List(objName='_status')
        self.params = List(objName='_params')
        self.statusList = ['status', 'initTime', 'endTime', 'error', 'isInteractive', 'mode']
        if protocol is None:
            objList = []
        else:
            objList = [protocol]
        ObjectTreeProvider.__init__(self, objList)
        
        
class StepsTreeProvider(TreeProvider):
    """Create the tree elements for a Protocol run"""
    def __init__(self, stepsList):
        for i, s in enumerate(stepsList):
            if not s._index:
                s._index = i + 1
            
        self._stepsList = stepsList
        self.getColumns = lambda: [('Index', 50), ('Step', 200), ('Time', 150), ('Class', 100)]
        self._parentDict = {}
    
    def getObjects(self):
        return self._stepsList
        
    def getObjectInfo(self, obj):
        info = {'key': obj._index, 
                'values': (str(obj), prettyDelta(obj.getElapsedTime()), obj.getClassName())}
            
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
    

class StepsWindow(BrowserWindow):
    def __init__(self, title, parentWindow, protocol, **args):
        self._protocol = protocol
        provider = StepsTreeProvider(protocol._steps)
        BrowserWindow.__init__(self, title, parentWindow, weight=False, **args)
        # Create buttons toolbar
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(1, weight=1)
        
        toolbar = tk.Frame(self.root)
        toolbar.grid(row=0, column=0, sticky='nw', padx=5, pady=5)
        btn = tk.Label(toolbar, text="Tree", image=self.getImage(Icon.ACTION_STEPS), 
                       compound=tk.LEFT, cursor='hand2')
        btn.bind('<Button-1>', self._showTree)
        btn.grid(row=0, column=0, sticky='nw')
        # Create and set browser
        browser = ObjectBrowser(self.root, provider, showPreviewTop=False)
        self.setBrowser(browser, row=1, column=0)
        
    def _showTree(self, e=None):
        g = self._protocol.getStepsGraph()
        w = gui.Window("Protocol steps", self, minsize=(800, 600))
        root = w.root
        canvas = Canvas(root, width=600, height=500)
        canvas.grid(row=0, column=0, sticky='nsew')
        lt = LevelTree(g)
        lt.setCanvas(canvas)
        lt.paint()
        canvas.updateScrollRegion()
        w.show()
    

class RunIOTreeProvider(TreeProvider):
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
            inputs = [attr for n, attr in self.protocol.iterInputAttributes()]
            outputs = [attr for n, attr in self.protocol.iterOutputAttributes(EMObject)]
            self.inputStr = String(Message.LABEL_INPUT)
            self.outputStr = String(Message.LABEL_OUTPUT)
            objs = [self.inputStr, self.outputStr] + inputs + outputs                
        return objs
    
    def _visualizeObject(self, ViewerClass, obj):
        viewer = ViewerClass(project=self.protocol.getProject(),
                             protocol=self.protocol,
                             parent=self.parent.windows)
        viewer.visualize(obj)
        
    def _editObject(self, obj):
        """Open the Edit GUI Form given an instance"""
        EditObjectDialog(self.parent, Message.TITLE_EDIT_OBJECT, obj, self.mapper)
        
    def getObjectPreview(self, obj):
        desc = "<name>: " + obj.getName()
        
        return (None, desc)
    
    def getObjectActions(self, obj):
        if isinstance(obj, Pointer):
            obj = obj.get()
        actions = []    
        
        viewers = findViewers(obj.getClassName(), DESKTOP_TKINTER)
        for v in viewers:
            actions.append(('Open with %s' % v.__name__, 
                            lambda : self._visualizeObject(v, obj), 
                            Icon.ACTION_VISUALIZE))
            
        # EDIT 
        actions.append((Message.LABEL_EDIT, 
                        lambda : self._editObject(obj),
                        Icon.ACTION_EDIT))
            
        return actions
    
    def getObjectInfo(self, obj):
        if obj is None or not obj.hasValue():
            return None
        
        if isinstance(obj, String):
            value = obj.get()
            info = {'key': value, 'text': value, 'values': (''), 'open': True}
        else:
            image = Icon.ACTION_OUT
            parent = self.outputStr
            name = obj.getLastName()
            
            if isinstance(obj, Pointer):
                obj = obj.get()
                if obj is None:
                    return None
                image = Icon.ACTION_IN
                parent = self.inputStr
#                objName = self.mapper.getFullName(obj)
                objName = obj.getNameId()
                name += '   (from %s)' % objName
            info = {'key': obj.getObjId(), 'parent': parent, 'image': image,
                    'text': name, 'values': (str(obj),)}
        return info     
    
   
class ProtocolsView(tk.Frame):
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
        self.showGraph = self.settings.graphView.get()
        
        self._selection = self.settings.runSelection
        self._items = {}
        
        self.style = ttk.Style()
        self.root.bind("<F5>", self.refreshRuns)
        self.__autoRefreshCounter = 3 # start by 3 secs  

        c = self.createContent()
        gui.configureWeigths(self)
        c.grid(row=0, column=0, sticky='news')
        
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
        self.protTree = self.createProtocolsTree(protFrame, bgColor)
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
        gui.configureWeigths(toolbar)
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
        gui.configureWeigths(runsFrame)
        
        self.createRunsGraph(runsFrame)
        
        if self.showGraph:
            treeWidget = self.runsGraph
        else:
            treeWidget = self.runsTree
            
        treeWidget.grid(row=0, column=0, sticky='news')
        
        # Create the Selected Run Info
        infoFrame = tk.Frame(v)
        infoFrame.columnconfigure(0, weight=1)
        infoFrame.rowconfigure(1, weight=1)
        # Create the Analyze results button
        btnAnalyze = gui.Button(infoFrame, text=Message.LABEL_ANALYZE, fg='white', bg=Color.RED_COLOR, # font=self.font, 
                          image=self.getImage(Icon.ACTION_VISUALIZE), compound=tk.LEFT, 
                        activeforeground='white', activebackground='#A60C0C', command=self._analyzeResultsClicked)
        btnAnalyze.grid(row=0, column=0, sticky='ne', padx=15)
        #self.style.configure("W.TNotebook")#, background='white')
        tab = ttk.Notebook(infoFrame)#, style='W.TNotebook')

        # Summary tab
        dframe = tk.Frame(tab, bg='white')
        gui.configureWeigths(dframe, row=0)
        gui.configureWeigths(dframe, row=2)
        provider = RunIOTreeProvider(self, self.getSelectedProtocol(), self.project.mapper)
        self.style.configure("NoBorder.Treeview", background='white', borderwidth=0, font=self.windows.font)
        self.infoTree = BoundTree(dframe, provider, height=6, show='tree', style="NoBorder.Treeview") 
        self.infoTree.grid(row=0, column=0, sticky='news')
        label = tk.Label(dframe, text='SUMMARY', bg='white', font=self.windows.fontBold)
        label.grid(row=1, column=0, sticky='nw', padx=(15, 0))  
        self.summaryText = TaggedText(dframe, width=40, height=5, bg='white', bd=0)
        self.summaryText.grid(row=2, column=0, sticky='news', padx=(30, 0))        
        
        # Method tab
        mframe = tk.Frame(tab)
        gui.configureWeigths(mframe)
        self.methodText = TaggedText(mframe, width=40, height=15, bg='white')
        self.methodText.grid(row=0, column=0, sticky='news')   
        
        #Logs 
        #TODO: join 3 logs in just one tab
        ologframe = tk.Frame(tab)
        gui.configureWeigths(ologframe)
        self.outputViewer = TextFileViewer(ologframe, allowOpen=True)
        self.outputViewer.grid(row=0, column=0, sticky='news')
        self.outputViewer.windows = self.windows
        
        self._updateSelection()
#         self.outputLogText = TaggedText(ologframe, width=40, height=15, 
#                                         bg='black', foreground='white')
#         self.outputLogText.grid(row=0, column=0, sticky='news')
#         elogframe = tk.Frame(tab)
#         gui.configureWeigths(elogframe)
#         self.errorLogText = TaggedText(elogframe, width=40, height=15, 
#                                        bg='black', foreground='white')
#         self.errorLogText.grid(row=0, column=0, sticky='news')  
#         slogframe = tk.Frame(tab)
#         gui.configureWeigths(slogframe)
#         self.scipionLogText = TaggedText(slogframe, width=40, height=15, 
#                                          bg='black', foreground='white')
#         self.scipionLogText.grid(row=0, column=0, sticky='news')
        
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
        
        
    def refreshRuns(self, e=None, autoRefresh=False):
        """ Refresh the status of diplayed runs. """
        self.updateRunsTree(True)
        self.updateRunsGraph(True)
        if not autoRefresh:
            self.__autoRefreshCounter = 3 # start by 3 secs  
            if self.__autoRefresh:
                self.runsTree.after_cancel(self.__autoRefresh)
                self.__autoRefresh = self.runsTree.after(self.__autoRefreshCounter*1000, self._automaticRefreshRuns)
        
    def _automaticRefreshRuns(self, e=None):
        # Schedule an automatic refresh after 1 sec
        self.refreshRuns(autoRefresh=True)
        secs = self.__autoRefreshCounter
        # double the number of seconds up to 30 min
        self.__autoRefreshCounter = min(2*secs, 1800)
        self.__autoRefresh = self.runsTree.after(secs*1000, self._automaticRefreshRuns)
                
    def createActionToolbar(self):
        """ Prepare the buttons that will be available for protocol actions. """
       
        self.actionList = [ACTION_EDIT, ACTION_COPY, ACTION_DELETE, ACTION_STEPS, ACTION_BROWSE, 
                           ACTION_STOP, ACTION_CONTINUE, ACTION_RESULTS]
        self.actionButtons = {}
        
        def addButton(action, text, toolbar):
            btn = tk.Label(toolbar, text=text, image=self.getImage(ActionIcons[action]), 
                       compound=tk.LEFT, cursor='hand2', bg='white')
            btn.bind('<Button-1>', lambda e: self._runActionClicked(action))
            return btn
        
        for action in self.actionList:
            self.actionButtons[action] = addButton(action, action, self.runsToolbar)
            
        if self.showGraph:
            ActionIcons[ACTION_TREE] = RUNS_LIST
        else:
            ActionIcons[ACTION_TREE] = RUNS_TREE
            
        self.viewButtons = {}
        
        for i, action in enumerate([ACTION_TREE, ACTION_REFRESH]):
            btn = addButton(action, action, self.allToolbar)
            btn.grid(row=0, column=i)
            self.viewButtons[action] = btn
        
            
    def _updateActionToolbar(self):
        """ Update which action buttons should be visible. """
        
        def displayAction(action, i, cond=True):
            """ Show/hide the action button if the condition is met. """
            if cond:
                self.actionButtons[action].grid(row=0, column=i, sticky='sw', padx=(0, 5), ipadx=0)
            else:
                self.actionButtons[action].grid_remove()  
                
        for i, actionTuple in enumerate(self.provider.getActionsFromSelection()):
            action, cond = actionTuple
            displayAction(action, i, cond)          
        
    def createProtocolsTree(self, parent, bgColor):
        """Create the protocols Tree displayed in left panel"""
        comboFrame = tk.Frame(parent, bg=bgColor)
        tk.Label(comboFrame, text='View', bg=bgColor).grid(row=0, column=0, padx=(0, 5), pady=5)
        choices = [pm.text.get() for pm in self.settings.protMenuList]
        initialChoice = choices[self.settings.protMenuList.getIndex()]
        combo = ComboBox(comboFrame, choices=choices, initial=initialChoice)
        combo.setChangeCallback(self._onSelectProtocols)
        combo.grid(row=0, column=1)
        comboFrame.grid(row=0, column=0, padx=5, pady=5, sticky='nw')
        
        self.style.configure("W.Treeview", background=Color.LIGHT_GREY_COLOR, borderwidth=0)
        tree = Tree(parent, show='tree', style='W.Treeview')
        tree.column('#0', minwidth=300)
        tree.tag_configure('protocol', image=self.getImage('python_file.gif'))
        tree.tag_bind('protocol', '<Double-1>', self._protocolItemClick)
        tree.tag_configure('protocol_base', image=self.getImage('class_obj.gif'))
        tree.tag_configure('section', font=self.windows.fontBold)
        tree.grid(row=1, column=0, sticky='news')
        # Program automatic refresh
        tree.after(3000, self._automaticRefreshRuns)
        return tree

    def _onSelectProtocols(self, combo):
        """ This function will be called when a protocol menu
        is selected. The index of the new menu is passed. 
        """
        protIndex = combo.getIndex()
        self.protCfg = self.settings.setCurrentProtocolMenu(protIndex)
        self.updateProtocolsTree(self.protCfg)
                
    def updateProtocolsTree(self, protCfg):
        self.protCfg = protCfg
        self.protTree.clear()
        self.protTree.unbind('<<TreeviewOpen>>')
        self.protTree.unbind('<<TreeviewClose>>')
        self.protTreeItems = {}
        subclassedDict = {} # Check which classes serve as base to not show them
        for k1, v1 in emProtocolsDict.iteritems():
            for k2, v2 in emProtocolsDict.iteritems():
                if v1 is not v2 and issubclass(v1, v2):
                    subclassedDict[k2] = True
        populateTree(self, self.protTree, self.protTreeItems, '', self.protCfg, subclassedDict)
        self.protTree.bind('<<TreeviewOpen>>', lambda e: self._treeViewItemChange(True))
        self.protTree.bind('<<TreeviewClose>>', lambda e: self._treeViewItemChange(False))
        
    def _treeViewItemChange(self, openItem):
        item = self.protTree.focus()
        if item in self.protTreeItems:
            self.protTreeItems[item].openItem.set(openItem)
        
    def createRunsTree(self, parent):
        self.provider = RunsTreeProvider(self.project, self._runActionClicked)
        tree = BoundTree(parent, self.provider)        
        tree.itemDoubleClick = self._runItemDoubleClick
        tree.itemClick = self._runTreeItemClick
            
        return tree
   
    def updateRunsTree(self, refresh=False):
        self.provider.setRefresh(refresh)
        self.runsTree.update()
        self.updateRunsTreeSelection()

    def updateRunsTreeSelection(self):
        for prot in self._iterSelectedProtocols():
            treeId = self.provider.getObjectFromId(prot.getObjId())._treeId
            self.runsTree.selection_add(treeId)
    
    def createRunsGraph(self, parent):
        self.runsGraph = Canvas(parent, width=400, height=400)
        self.runsGraph.onClickCallback = self._runItemClick
        self.runsGraph.onDoubleClickCallback = self._runItemDoubleClick
        self.runsGraph.onRightClickCallback = lambda e: self.provider.getObjectActions(e.node.run)
        self.runsGraph.onControlClickCallback = self._runItemControlClick
        parent.grid_columnconfigure(0, weight=1)
        parent.grid_rowconfigure(0, weight=1)
        
        self.updateRunsGraph()        

    def updateRunsGraph(self, refresh=False):      
        g = self.project.getRunsGraph(refresh=refresh)
        lt = LevelTree(g)
        self.runsGraph.clear()
        lt.setCanvas(self.runsGraph)
        lt.paint(self.createRunItem)
        #self.updateRunsGraphSelection()
        self.runsGraph.updateScrollRegion()

    def createRunItem(self, canvas, node, y):
        nodeText = node.label
        textColor = 'black'
        color = '#ADD8E6' #Lightblue
            
        if node.run:
            status = node.run.status.get(STATUS_FAILED)
            nodeText = nodeText + '\n' + node.run.getStatusMessage()
            color = STATUS_COLORS[status]
        
        item = self.runsGraph.createTextbox(nodeText, 100, y, bgColor=color, textColor=textColor)
        item.run = node.run
        node.item = item
        
        if node.run and node.run.getObjId() in self._selection:
            item.setSelected(True)
        return item
        
    def switchRunsView(self):
        self.showGraph = not self.showGraph
        
        if self.showGraph:
            ActionIcons[ACTION_TREE] = RUNS_LIST
            show = self.runsGraph.frame
            hide = self.runsTree
            #self.updateRunsGraphSelection()
        else:
            ActionIcons[ACTION_TREE] = RUNS_TREE
            show = self.runsTree
            hide = self.runsGraph.frame
            self.updateRunsTreeSelection()
            
        self.viewButtons[ACTION_TREE].config(image=self.getImage(ActionIcons[ACTION_TREE]))
        hide.grid_remove()
        show.grid(row=0, column=0, sticky='news')
        self.settings.graphView.set(self.showGraph)
    
    def _protocolItemClick(self, e=None):
        protClassName = self.protTree.getFirst().split('.')[-1]
        protClass = emProtocolsDict.get(protClassName)
        prot = self.project.newProtocol(protClass)
        self._openProtocolForm(prot)
        
    def _updateSelection(self, prot=None):
        self._fillSummary()
        self._fillMethod()
        self._fillLogs()
        
        #FIXME:
        self._updateActionToolbar()
        
    def _runTreeItemClick(self, item=None):
        self._selection.clear()
        for prot in self.runsTree.iterSelectedObjects():
            self._selection.append(prot.getObjId())
        self._updateSelection()
                         
    def _runItemClick(self, item=None):
        # Get last selected item for tree or graph
        if self.showGraph:
            prot = item.node.run
            g = self.project.getRunsGraph(refresh=False)
            for node in g.getNodes():
                if node.run and node.run.getObjId() in self._selection:
                    node.item.setSelected(False)
            item.setSelected(True)
        else:
            prot = self.project.mapper.selectById(int(self.runsTree.getFirst()))
        
        self._selection.clear()
        self._selection.append(prot.getObjId())
        self._updateSelection(prot)
        self.runsGraph.update_idletasks()
        
    def _runItemDoubleClick(self, e=None):
        self._runActionClicked(ACTION_EDIT)
        
    def _runItemControlClick(self, item=None):
        # Get last selected item for tree or graph
        if self.showGraph:
            prot = item.node.run
            protId = prot.getObjId()
            if protId in self._selection:
                item.setSelected(False)
                self._selection.remove(protId)
            else:
                item.setSelected(True)
                self._selection.append(prot.getObjId())
        else:
            prot = self.project.mapper.selectById(int(self.runsTree.getFirst()))        
        
        self._updateSelection(prot)
        
    def _openProtocolForm(self, prot):
        """Open the Protocol GUI Form given a Protocol instance"""
        hosts = [host.getLabel() for host in self.settings.getHosts()]
        w = FormWindow(Message.TITLE_NAME_RUN + prot.getClassName(), prot, 
                       self._executeSaveProtocol, self.windows,
                       hostList=hosts)
        w.adjustSize()
        w.show(center=True)
        
    def _browseSteps(self):
        """ Open a new window with the steps list. """
        window = StepsWindow(Message.TITLE_BROWSE_DATA, self.windows, 
                             self.getSelectedProtocol(), icon=self.icon)
        window.show()        
    
    def _browseRunData(self):
        provider = ProtocolTreeProvider(self.getSelectedProtocol())
        window = BrowserWindow(Message.TITLE_BROWSE_DATA, self.windows, icon=self.icon)
        window.setBrowser(ObjectBrowser(window.root, provider))
        window.itemConfig(self.getSelectedProtocol(), open=True)  
        window.show()
        
    def _iterSelectedProtocols(self):
        for protId in sorted(self._selection):
            prot = self.project.getProtocol(protId)
            if prot:
                yield prot
            
    def _getSelectedProtocols(self):
        return [prot for prot in self._iterSelectedProtocols()]
            
    def getSelectedProtocol(self):
        if self._selection:
            return self.project.getProtocol(self._selection[0])
        return None
            
    def _fillSummary(self):
        self.summaryText.clear()
        n = len(self._selection)
        
        if n == 1:
            prot = self.getSelectedProtocol()
            provider = RunIOTreeProvider(self, prot, self.project.mapper)
            self.infoTree.setProvider(provider)
            self.infoTree.grid(row=0, column=0, sticky='news')
            self.infoTree.update_idletasks()
            # Update summary
            self.summaryText.addText(prot.summary())
        
        elif n > 1:
            self.infoTree.clear()
            for prot in self._iterSelectedProtocols():
                self.summaryText.addLine('> _%s_' % prot.getRunName())
                for line in prot.summary():
                    self.summaryText.addLine(line)
                self.summaryText.addLine('')
        
    def _fillMethod(self):
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
        n = len(self._selection)
        
        if n == 1:
            prot = self.project.getProtocol(self._selection[0])
            i = self.outputViewer.getIndex()
            self.outputViewer.clear()
            for f in prot.getLogPaths():
                self.outputViewer.addFile(f)
            self.outputViewer.setIndex(i) # Preserve the last selected tab
        
    def _scheduleRunsUpdate(self, secs=1):
        self.runsTree.after(secs*1000, self.refreshRuns)
        
    def _executeSaveProtocol(self, prot, onlySave=False):
        if onlySave:
            self.project.saveProtocol(prot)
            msg = Message.LABEL_SAVED_FORM
#            msg = "Protocol sucessfully saved."
        else:
            self.project.launchProtocol(prot)
            self._scheduleRunsUpdate()
            msg = ""
            
        
        return msg
        
    def _continueProtocol(self, prot):
        self.project.continueProtocol(prot)
        self._scheduleRunsUpdate()
        
    def _deleteProtocol(self):
        if askYesNo(Message.TITLE_DELETE_FORM, Message.LABEL_DELETE_FORM, self.root):
            self.project.deleteProtocol(*self._getSelectedProtocols())
            self._selection.clear()
            self._scheduleRunsUpdate()
            
    def _stopProtocol(self, prot):
        if askYesNo(Message.TITLE_STOP_FORM, Message.LABEL_STOP_FORM, self.root):
            self.project.stopProtocol(prot)
            self._scheduleRunsUpdate()

    def _analyzeResults(self, prot):        

        viewers = findViewers(prot.getClassName(), DESKTOP_TKINTER)
        if len(viewers):
            #TODO: If there are more than one viewer we should display a selection menu
            firstViewer = viewers[0](project=self.project) # Instanciate the first available viewer
            firstViewer.visualize(prot, windows=self.windows)
        else:
            for _, output in prot.iterOutputAttributes(EMObject):
                viewers = findViewers(output.getClassName(), DESKTOP_TKINTER)
                if len(viewers):
                    #TODO: If there are more than one viewer we should display a selection menu
                    viewerclass = viewers[0]
                    firstViewer = viewerclass(project=self.project) # Instanciate the first available viewer
                    firstViewer.visualize(output, windows=self.windows, protocol=prot)
            
        
    def _analyzeResultsClicked(self, e=None):
        """ this method should be called when button "Analyze results" is called. """
        self._analyzeResults(self.getSelectedProtocol())
                
    def _runActionClicked(self, action):
        prot = self.getSelectedProtocol()
        if prot:
            try:
                if action == ACTION_DEFAULT:
                    pass
                elif action == ACTION_EDIT:
                    self._openProtocolForm(prot)
                elif action == ACTION_COPY:
                    newProt = self.project.copyProtocol(prot)
                    self._openProtocolForm(newProt)
                elif action == ACTION_DELETE:
                    self._deleteProtocol()
                elif action == ACTION_STEPS:
                    self._browseSteps()                    
                elif action == ACTION_BROWSE:
                    self._browseRunData()
                elif action == ACTION_STOP:
                    self._stopProtocol(prot)
                elif action == ACTION_CONTINUE:
                    self._continueProtocol(prot)
                elif action == ACTION_RESULTS:
                    self._analyzeResults(prot)
            except Exception, ex:
                self.windows.showError(str(ex))
 
        # Following actions do not need a select run
        if action == ACTION_TREE:
            self.switchRunsView()
        elif action == ACTION_REFRESH:
            self.refreshRuns()
    
