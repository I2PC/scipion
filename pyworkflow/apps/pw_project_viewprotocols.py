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
import os, sys
from os.path import join, exists, basename

import Tkinter as tk
import ttk
import tkFont

from pyworkflow.gui.tree import TreeProvider, BoundTree
from pyworkflow.protocol.protocol import *

import pyworkflow as pw
from pyworkflow.object import *
from pyworkflow.em import *
from pyworkflow.protocol import *
from pyworkflow.protocol.params import *
from pyworkflow.mapper import SqliteMapper, XmlMapper
from pyworkflow.project import Project
from pyworkflow.utils.properties import Message, Icon, Color

import pyworkflow.gui as gui
from pyworkflow.gui import getImage
from pyworkflow.gui.tree import Tree, ObjectTreeProvider, DbTreeProvider, ProjectRunsTreeProvider
from pyworkflow.gui.form import FormWindow, editObject
from pyworkflow.gui.dialog import askYesNo
from pyworkflow.gui.text import TaggedText
from pyworkflow.gui import Canvas
from pyworkflow.gui.graph import LevelTree
from pyworkflow.gui.widgets import ComboBox

from config import *
from pw_browser import BrowserWindow

ACTION_EDIT = Message.LABEL_EDIT
ACTION_COPY = Message.LABEL_COPY
ACTION_DELETE = Message.LABEL_DELETE
ACTION_REFRESH = Message.LABEL_REFRESH
ACTION_STEPS = Message.LABEL_BROWSE
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
    
    def getObjectActions(self, obj):
        prot = obj # Object should be a protocol
        actionsList = [(ACTION_EDIT, Message.LABEL_EDIT_ACTION),
                       (ACTION_COPY, Message.LABEL_COPY_ACTION),
                       (ACTION_DELETE, Message.LABEL_DELETE_ACTION),
                       #(None, None),
                       #(ACTION_STOP, 'Stop'),
                       (ACTION_STEPS, Message.LABEL_BROWSE_ACTION)
                       ]
        status = prot.status.get()
        if status == STATUS_RUNNING:
            actionsList.insert(0, (ACTION_STOP, Message.LABEL_STOP_ACTION))
            actionsList.insert(1, None)
        elif status == STATUS_WAITING_APPROVAL:
            actionsList.insert(0, (ACTION_CONTINUE, Message.LABEL_CONTINUE_ACTION))
            actionsList.insert(1, None)
        
        actions = []
        def appendAction(a):
            v = a
            if v is not None:
                action = a[0]
                text = a[1]
                v = (text, lambda: self.actionFunc(action), ActionIcons[action])
            actions.append(v)
            
        for a in actionsList:
            appendAction(a)
            
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
    

class RunIOTreeProvider(TreeProvider):
    """Create the tree elements from a Protocol Run input/output childs"""
    def __init__(self, protocol, mapper):
        #TreeProvider.__init__(self)
        self.protocol = protocol
        self.mapper = mapper
        self.viewer = XmippViewer()

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
    
    def visualize(self, Viewer, obj):
        viewer = Viewer(project=self.protocol.getProject())
        viewer.visualize(obj)
        
    def edit(self, obj):
        """Open the Edit GUI Form given an instance"""
        root = tk.Frame()
        d = editObject(self, "Object Editor", root, obj, self.mapper)
        
        if d.resultYes():
            pass
#            print "Object edited"
        
    def getObjectPreview(self, obj):
        desc = "<name>: " + obj.getName()
        
        return (None, desc)
    
    def getObjectActions(self, obj):
        from pyworkflow.viewer import DESKTOP_TKINTER
        
        if isinstance(obj, Pointer):
            obj = obj.get()
        actions = []    
        
        
        viewers = findViewers(obj.getClassName(), DESKTOP_TKINTER)
        for v in viewers:
            actions.append(('Open with %s' % v.__name__, lambda : self.visualize(v, obj)))
            
        # EDIT 
        actions.append((Message.LABEL_EDIT, lambda : self.edit(obj)))
            
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
        self.selectedProtocol = None
        self.windows = windows
        self.project = windows.project
        self.root = windows.root
        self.getImage = windows.getImage
        self.protCfg = windows.protCfg
        self.icon = windows.icon
        self.settings = windows.getSettings()
        self.showGraph = self.settings.graphView.get()
        
        self.style = ttk.Style()
        self.root.bind("<F5>", self.refreshRuns)
        self.__autoRefreshCounter = 3 # start by 3 secs  
        # Hide the right-click menu
        #self.root.bind('<FocusOut>', self._unpostMenu)
        #self.root.bind("<Key>", self._unpostMenu)
        #self.root.bind('<Button-1>', self._unpostMenu)
        
        #self.menuRun = tk.Menu(self.root, tearoff=0)
        c = self.createContent()
        gui.configureWeigths(self)
        c.grid(row=0, column=0, sticky='news')
        
        #self.viewer = XmippViewer()
        
    def createContent(self):
        """ Create the Protocols View for the Project.
        It has two panes:
            Left: containing the Protocol classes tree
            Right: containing the Runs list
        """
        p = tk.PanedWindow(self, orient=tk.HORIZONTAL, bg='white')
        
        # Left pane, contains Protocols Pane
        leftFrame = tk.Frame(p, bg='white')
        leftFrame.columnconfigure(0, weight=1)
        leftFrame.rowconfigure(1, weight=1)

        
        # Protocols Tree Pane        
        bgColor = Color.LIGHT_GREY_COLOR
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
        # Data tab
        dframe = tk.Frame(tab)
        gui.configureWeigths(dframe)
        provider = RunIOTreeProvider(self.selectedProtocol, self.project.mapper)
        self.infoTree = BoundTree(dframe, provider) 
        TaggedText(dframe, width=40, height=15, bg='white')
        self.infoTree.grid(row=0, column=0, sticky='news')  
        # Summary tab
        sframe = tk.Frame(tab)
        gui.configureWeigths(sframe)
        self.summaryText = TaggedText(sframe, width=40, height=15, bg='white')
        self.summaryText.grid(row=0, column=0, sticky='news')        
        # Method tab
        mframe = tk.Frame(tab)
        gui.configureWeigths(mframe)
        self.methodText = TaggedText(mframe, width=40, height=15, bg='white')
        self.methodText.grid(row=0, column=0, sticky='news')   
        #Logs 
        #TODO: join 3 logs in just one tab
        ologframe = tk.Frame(tab)
        gui.configureWeigths(ologframe)
        self.outputLogText = TaggedText(ologframe, width=40, height=15, 
                                        bg='black', foreground='white')
        self.outputLogText.grid(row=0, column=0, sticky='news')
        elogframe = tk.Frame(tab)
        gui.configureWeigths(elogframe)
        self.errorLogText = TaggedText(elogframe, width=40, height=15, 
                                       bg='black', foreground='white')
        self.errorLogText.grid(row=0, column=0, sticky='news')  
        slogframe = tk.Frame(tab)
        gui.configureWeigths(slogframe)
        self.scipionLogText = TaggedText(slogframe, width=40, height=15, 
                                         bg='black', foreground='white')
        self.scipionLogText.grid(row=0, column=0, sticky='news')
        
        tab.add(dframe, text=Message.LABEL_DATA)
        tab.add(sframe, text=Message.LABEL_SUMMARY)   
        tab.add(mframe, text=Message.LABEL_METHODS)
        tab.add(ologframe, text=Message.LABEL_LOGS_OUTPUT)
        tab.add(elogframe, text=Message.LABEL_LOGS_ERROR)
        tab.add(slogframe, text=Message.LABEL_LOGS_SCIPION)
        tab.grid(row=1, column=0, sticky='news')
        
        v.add(runsFrame, weight=3)
        v.add(infoFrame, weight=1)
        v.grid(row=1, column=0, sticky='news')
        
        # Add sub-windows to PanedWindows
        p.add(leftFrame, padx=5, pady=5)
        p.add(rightFrame, padx=5, pady=5)
        p.paneconfig(leftFrame, minsize=300)
        p.paneconfig(rightFrame, minsize=400)        
        
        return p
        
        
    def refreshRuns(self, e=None, autoRefresh=False):
        """ Refresh the status of diplayed runs. """
        self.updateRunsTree()
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
       
        self.actionList = [ACTION_EDIT, ACTION_COPY, ACTION_DELETE, ACTION_STEPS, 
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
        
            
    def updateActionToolbar(self):
        """ Update which action buttons should be visible. """
        status = self.selectedProtocol.status.get()
        
        def displayAction(action, i, cond=True):
            """ Show/hide the action button if the condition is met. """
            if cond:
                self.actionButtons[action].grid(row=0, column=i, sticky='sw', padx=(0, 5), ipadx=0)
            else:
                self.actionButtons[action].grid_remove()            
                
        for i, action in enumerate(self.actionList[:-2]):
            displayAction(action, i)            
            
        displayAction(ACTION_DELETE, 2, status != STATUS_RUNNING)     
        displayAction(ACTION_STOP, 4, status == STATUS_RUNNING)
        displayAction(ACTION_CONTINUE, 5, status == STATUS_WAITING_APPROVAL)
#        displayAction(ACTION_RESULTS, 6, status != STATUS_RUNNING)
        
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
        f = tkFont.Font(family='helvetica', size='10', weight='bold')
        tree.tag_configure('section', font=f)
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
        tree.itemClick = self._runItemClick
            
        return tree
   
    def updateRunsTree(self):
        self.runsTree.update()
        self.updateRunsTreeSelection()

    def updateRunsTreeSelection(self):
        if not self.showGraph:
            if self.selectedProtocol is not None:
                for i, obj in enumerate(self.runsTree._objects):
                    if self.selectedProtocol.getObjId() == obj.getObjId():
                        self.runsTree.selectChildByIndex(i)
                        break
    
    def createRunsGraph(self, parent):
        self.runsGraph = Canvas(parent, width=400, height=400)
        self.runsGraph.onClickCallback = self._runItemClick
        self.runsGraph.onDoubleClickCallback = self._runItemDoubleClick
        parent.grid_columnconfigure(0, weight=1)
        parent.grid_rowconfigure(0, weight=1)
        
        self.updateRunsGraph()

    def updateRunsGraph(self, refresh=False):      
        g = self.project.getRunsGraph(refresh=refresh)
        #g.printDot()
        lt = LevelTree(g)
        self.runsGraph.clear()
        lt.setCanvas(self.runsGraph)
        lt.paint(self.createRunItem)
        self.updateRunsGraphSelection()
        self.runsGraph.updateScrollRegion()

    def updateRunsGraphSelection(self):
        if self.showGraph:
            if self.selectedProtocol is not None:
                for item in self.runsGraph.items.values():
                    #item.setSelected(True)
                    run = item.run
                    if run and self.selectedProtocol.getObjId() == run.getObjId():
                        self.runsGraph.selectItem(item)
                        break
                            
    def createRunItem(self, canvas, node, y):
        """ If not nodeBuildFunc is specified, this one will be used by default."""
        nodeText = node.label
        textColor = 'black'
        color = '#ADD8E6' #Lightblue
            
        if node.run:
            status = node.run.status.get(STATUS_FAILED)
            nodeText = nodeText + '\n' + node.run.getStatusMessage()
            color = STATUS_COLORS[status]
        
        item = self.runsGraph.createTextbox(nodeText, 100, y, bgColor=color, textColor=textColor)
        item.run = node.run
        return item
        
    
    def switchRunsView(self):
        self.showGraph = not self.showGraph
        
        if self.showGraph:
            ActionIcons[ACTION_TREE] = RUNS_LIST
            show = self.runsGraph.frame
            hide = self.runsTree
            self.updateRunsGraphSelection()
        else:
            ActionIcons[ACTION_TREE] = RUNS_TREE
            show = self.runsTree
            hide = self.runsGraph.frame
            self.updateRunsTreeSelection()
            
        self.viewButtons[ACTION_TREE].config(image=self.getImage(ActionIcons[ACTION_TREE]))
        hide.grid_remove()
        show.grid(row=0, column=0, sticky='news')
        self.settings.graphView.set(self.showGraph)
        #self.settings.write()
        
    
    def _protocolItemClick(self, e=None):
        protClassName = self.protTree.getFirst().split('.')[-1]
        protClass = emProtocolsDict.get(protClassName)
        prot = self.project.newProtocol(protClass)
        self._openProtocolForm(prot)
        
    def _selectProtocol(self, prot):
        if prot is not None:
            prot.mapper = self.project.mapper
            del self.selectedProtocol
            self.selectedProtocol = prot
            # TODO self.settings.selectedProtocol.set(prot)
            self.updateActionToolbar()
            self._fillData()
            self._fillSummary()
            self._fillMethod()
            self._fillLogs()
        else:
            pass #TODO: implement what to do
                    
    def _runItemClick(self, e=None):
        # Get last selected item for tree or graph
        if self.showGraph:
            prot = e.node.run
        else:
            prot = self.project.mapper.selectById(int(self.runsTree.getFirst()))
        self._selectProtocol(prot)
        
    def _runItemDoubleClick(self, e=None):
        self._runActionClicked(ACTION_EDIT)
        
    def _openProtocolForm(self, prot):
        """Open the Protocol GUI Form given a Protocol instance"""
        hosts = [host.getLabel() for host in self.settings.getHosts()]
        w = FormWindow(Message.TITLE_NAME_RUN + prot.getClassName(), prot, 
                       self._executeSaveProtocol, self.windows,
                       hostList=hosts)
        w.show(center=True)
        
    def _browseRunData(self):
        provider = ProtocolTreeProvider(self.selectedProtocol)
        window = BrowserWindow(Message.TITLE_BROWSE_DATA, provider, self.windows, icon=self.icon)
        window.itemConfig(self.selectedProtocol, open=True)  
        window.show()
        
    def _fillData(self):
        provider = RunIOTreeProvider(self.selectedProtocol, self.project.mapper)
        self.infoTree.setProvider(provider)
        #self.infoTree.itemConfig(self.selectedProtocol, open=True)  
        
    def _fillSummary(self):
        self.summaryText.clear()
        self.summaryText.addText(self.selectedProtocol.summary())
        
    def _fillMethod(self):
        self.methodText.clear()
        self.methodText.addText(self.selectedProtocol.methods())
        
    def _fillLogs(self):
        #TODO: REMOVE THIS...READ LOGS DIRECTLY FROM FILE
        fOutString, fErrString, fScipionString = self.selectedProtocol.getLogsAsStrings()
        
        self.outputLogText.clear()
        self.outputLogText.addText(fOutString)
        self.outputLogText.goEnd()
        
        self.errorLogText.clear()
        self.errorLogText.addText(fErrString)
        self.errorLogText.goEnd()
        
        self.scipionLogText.clear()
        self.scipionLogText.addText(fScipionString)
        self.scipionLogText.goEnd()
        
    def _scheduleRunsUpdate(self, secs=1):
        self.runsTree.after(secs*1000, self.refreshRuns)
        
    def _executeSaveProtocol(self, prot, onlySave=False):
        if onlySave:
            self.project.saveProtocol(prot)
            msg = Message.LABEL_SAVED_FORM
#            msg = "Protocol sucessfully saved."
        else:
            self.project.launchProtocol(prot)
            msg = ""
            
        self._scheduleRunsUpdate()
        
        return msg
        
    def _continueProtocol(self, prot):
        self.project.continueProtocol(prot)
        self._scheduleRunsUpdate()
        
    def _deleteProtocol(self, prot):
        if askYesNo(Message.TITLE_DELETE_FORM, Message.LABEL_DELETE_FORM, self.root):
            self.project.deleteProtocol(prot)
            self._scheduleRunsUpdate()
            
    def _stopProtocol(self, prot):
        if askYesNo(Message.TITLE_STOP_FORM, Message.LABEL_STOP_FORM, self.root):
            self.project.stopProtocol(prot)
            self._scheduleRunsUpdate()

    def _analyzeResults(self, prot):
        from pyworkflow.em import findViewers
        from pyworkflow.viewer import DESKTOP_TKINTER

        viewers = findViewers(prot.getClassName(), DESKTOP_TKINTER)
        if len(viewers):
            #TODO: If there are more than one viewer we should display a selection menu
            firstViewer = viewers[0](project=self.project) # Instanciate the first available viewer
            firstViewer.visualize(prot, windows=self.windows)
        else:
            for key, output in prot.iterOutputAttributes(EMObject):
                viewers = findViewers(output.getClassName(), DESKTOP_TKINTER)
                if len(viewers):
                    #TODO: If there are more than one viewer we should display a selection menu
                    viewerclass = viewers[0]
                    print viewerclass
                    firstViewer = viewerclass(project=self.project) # Instanciate the first available viewer
                    firstViewer.visualize(output, windows=self.windows, protocol=prot)
            
        
    def _analyzeResultsClicked(self, e=None):
        """ this method should be called when button "Analyze results" is called. """
        self._analyzeResults(self.selectedProtocol)
                
    def _runActionClicked(self, action):
        prot = self.selectedProtocol
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
                    self._deleteProtocol(prot)
                elif action == ACTION_STEPS:
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
    
