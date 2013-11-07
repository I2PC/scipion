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

import pyworkflow.gui as gui
from pyworkflow.gui import getImage
from pyworkflow.gui.tree import Tree, ObjectTreeProvider, DbTreeProvider, ProjectRunsTreeProvider
from pyworkflow.gui.form import FormWindow
from pyworkflow.gui.dialog import askYesNo
from pyworkflow.gui.text import TaggedText
from pyworkflow.gui import Canvas
from pyworkflow.gui.graph import LevelTree

from config import *
from pw_browser import BrowserWindow

ACTION_EDIT = 'Edit'
ACTION_COPY = 'Copy'
ACTION_DELETE = 'Delete'
ACTION_REFRESH = 'Refresh'
ACTION_STEPS = 'Browse'
ACTION_TREE = 'Tree'
ACTION_STOP = 'Stop'
ACTION_DEFAULT = 'Default'
ACTION_CONTINUE = 'Continue'
ACTION_RESULTS = 'Analyze results'

ActionIcons = {
    ACTION_EDIT:  'edit.gif',
    ACTION_COPY:  'copy.gif',
    ACTION_DELETE:  'delete.gif',
    ACTION_REFRESH:  'refresh.gif',
    ACTION_STEPS:  'run_steps.gif',
    ACTION_TREE:  'tree2.gif',
    ACTION_STOP: 'stop.gif',
    ACTION_CONTINUE: 'play.png',
    ACTION_RESULTS: 'visualize.gif'
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

def populateTree(self, tree, treeItems, prefix, obj, level=0):
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
                    if not v is prot and issubclass(v, prot):# and Protocol.hasDefinition(v):
                        key = '%s.%s' % (item, k)
                        tree.insert(item, 'end', key, text=k, tags=('protocol'))
                        
            else:
                raise Exception("Class '%s' not found" % obj.value.get())
    else:
        key = prefix
    
    for sub in obj:
        populateTree(self, tree, treeItems, key, sub, level+1)
    

class RunsTreeProvider(ProjectRunsTreeProvider):
    """Provide runs info to populate tree"""
    def __init__(self, project, actionFunc):
        ProjectRunsTreeProvider.__init__(self, project)
        self.actionFunc = actionFunc
    
    def getObjectActions(self, obj):
        prot = obj # Object should be a protocol
        actionsList = [(ACTION_EDIT, 'Edit     '),
                       (ACTION_COPY, 'Copy   '),
                       (ACTION_DELETE, 'Delete    '),
                       #(None, None),
                       #(ACTION_STOP, 'Stop'),
                       (ACTION_STEPS, 'Browse ')
                       ]
        status = prot.status.get()
        if status == STATUS_RUNNING:
            actionsList.insert(0, (ACTION_STOP, 'Stop execution'))
            actionsList.insert(1, None)
        elif status == STATUS_WAITING_APPROVAL:
            actionsList.insert(0, (ACTION_CONTINUE, 'Approve continue'))
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
#        self.viewer = XmippViewer()
#        
#    def show(self, obj):
#        self.viewer.visualize(obj)
#        
#    def getObjectPreview(self, obj):
#        desc = "<name>: " + obj.getName()
#        
#        return (None, desc)
#    
#    def getObjectActions(self, obj):
#        if isinstance(obj, Pointer):
#            obj = obj.get()
#            
#        if isinstance(obj, SetOfMicrographs):
#            return [('Open Micrographs with Xmipp', lambda: self.viewer.visualize(obj))]
#        if isinstance(obj, SetOfImages):
#            return [('Open Images with Xmipp', lambda: self.viewer.visualize(obj))]
#        if isinstance(obj, XmippSetOfClasses2D):
#            return [('Open SetOfClasses2D with Xmipp', lambda: self.viewer.visualize(obj))]
#        return []   
#    
#    def getObjectInfo(self, obj):
#        info = ObjectTreeProvider.getObjectInfo(self, obj)
#        attrName = obj.getLastName()
#        if hasattr(self.protocol, attrName):
#            if isinstance(obj, Pointer) and obj.hasValue():
#                info['image'] = 'db_input.gif'
#            else:
#                if (self.protocol._definition.hasParam(attrName) or
#                    attrName in ['numberOfMpi', 'numberOfThreads']):
#                    info['parent'] = self.params
#                elif attrName in self.statusList:
#                    if info['parent'] is self.protocol:
#                        info['parent'] = self.status
#                    
#            if attrName.startswith('output'):# in self.protocol._outputs:
#                info['image'] = 'db_output.gif'
#        if obj is self.params or obj is self.status:
#            info['parent'] = self.protocol
#        return info     
#    
#    def _getChilds(self, obj):
#        childs = ObjectTreeProvider._getChilds(self, obj)
#        if obj is self.protocol:
#            childs.insert(0, self.status)
#            childs.insert(1, self.params)
#        return childs
    

class RunIOTreeProvider(TreeProvider):
    """Create the tree elements from a Protocol Run input/output childs"""
    def __init__(self, protocol, mapper):
        #TreeProvider.__init__(self)
        self.protocol = protocol
        self.mapper = mapper
        self.viewer = XmippViewer()

    def getColumns(self):
        return [('Attribute', 200), ('Class', 100)]
    
    def getObjects(self):
        objs = []
        if self.protocol:
            inputs = [attr for n, attr in self.protocol.iterInputAttributes()]
            outputs = [attr for n, attr in self.protocol.iterOutputAttributes(EMObject)]
            self.inputStr = String('Input')
            self.outputStr = String('Output')
            objs = [self.inputStr, self.outputStr] + inputs + outputs                
        return objs
    
    def visualize(self, Viewer, obj):
        Viewer().visualize(obj)
        
    def getObjectPreview(self, obj):
        desc = "<name>: " + obj.getName()
        
        return (None, desc)
    
    def getObjectActions(self, obj):
        from pyworkflow.em import findViewers
        from pyworkflow.viewer import DESKTOP_TKINTER
        
        if isinstance(obj, Pointer):
            obj = obj.get()
            
        actions = []    
        viewers = findViewers(obj.getClassName(), DESKTOP_TKINTER)
        for v in viewers:
            actions.append(('Open with %s' % v.__name__, lambda : self.visualize(v, obj)))
            
#        if isinstance(obj, SetOfMicrographs):
#            return [('Open Micrographs with Xmipp', lambda: self.viewer.visualize(obj))]
#        if isinstance(obj, XmippSetOfCoordinates):
#            return [('Open Coordinates with Xmipp', lambda: self.viewer.visualize(obj))]
#        if isinstance(obj, SetOfImages):
#            return [('Open Images with Xmipp', lambda: self.viewer.visualize(obj))]
#        if isinstance(obj, XmippSetOfClasses2D):
#            return [('Open SetOfClasses2D with Xmipp', lambda: self.viewer.visualize(obj))]
        return actions
    
    def getObjectInfo(self, obj):
        if obj is None or not obj.hasValue():
            return None
        
        if isinstance(obj, String):
            value = obj.get()
            info = {'key': value, 'text': value, 'values': (''), 'open': True}
        else:
            image = 'db_output.gif'
            parent = self.outputStr
            name = obj.getLastName()
            
            if isinstance(obj, Pointer):
                obj = obj.get()
                image = 'db_input.gif'
                parent = self.inputStr
                parentObj = self.mapper.getParent(obj)
                name += '   (from %s.%s)' % (parentObj.getLastName(), obj.getLastName())
            info = {'key': obj.getObjId(), 'parent': parent, 'image': image,
                    'text': name, 'values': (obj.getClassName())}
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
        
        self.root.bind("<F5>", self.refreshRuns)
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
        p = tk.PanedWindow(self, orient=tk.HORIZONTAL)
        
        # Left pane, contains Protocols Pane
        leftFrame = tk.Frame(p)
        leftFrame.columnconfigure(0, weight=1)
        leftFrame.rowconfigure(1, weight=1)

        # Protocols Tree Pane        
        protFrame = ttk.Labelframe(leftFrame, text=self.protCfg.text.get(), width=300, height=500)
        protFrame.grid(row=1, column=0, sticky='news', padx=5, pady=5)
        gui.configureWeigths(protFrame)
        self.protTree = self.createProtocolsTree(protFrame)
        self.updateProtocolsTree(self.protCfg)
        
        # Create the right Pane that will be composed by:
        # a Action Buttons TOOLBAR in the top
        # and another vertical Pane with:
        # Runs History (at Top)
        # Sectected run info (at Bottom)
        rightFrame = tk.Frame(p)
        rightFrame.columnconfigure(0, weight=1)
        rightFrame.rowconfigure(1, weight=1)
        #rightFrame.rowconfigure(0, minsize=label.winfo_reqheight())
        
        # Create the Action Buttons TOOLBAR
        toolbar = tk.Frame(rightFrame)
        toolbar.grid(row=0, column=0, sticky='news')
        gui.configureWeigths(toolbar)
        #toolbar.columnconfigure(0, weight=1)
        toolbar.columnconfigure(1, weight=1)
        
        self.runsToolbar = tk.Frame(toolbar)
        self.runsToolbar.grid(row=0, column=0, sticky='sw')
        # On the left of the toolbar will be other
        # actions that can be applied to all runs (refresh, graph view...)
        self.allToolbar = tk.Frame(toolbar)
        self.allToolbar.grid(row=0, column=10, sticky='se')
        self.createActionToolbar()

        # Create the Run History tree
        v = ttk.PanedWindow(rightFrame, orient=tk.VERTICAL)
        runsFrame = ttk.Labelframe(v, text=' History ', width=500, height=500)
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
        #infoFrame.columnconfigure(0, weight=1)
        gui.configureWeigths(infoFrame)
        
        tab = ttk.Notebook(infoFrame)
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
        #self.summaryText.addText("\nSummary should go <HERE!!!>\n More info here.")
        
        tab.add(dframe, text="Data")
        tab.add(sframe, text="Summary")     
        tab.grid(row=0, column=0, sticky='news')
        
        v.add(runsFrame, weight=3)
        v.add(infoFrame, weight=1)
        v.grid(row=1, column=0, sticky='news')
        
        # Add sub-windows to PanedWindows
        p.add(leftFrame, padx=5, pady=5)
        p.add(rightFrame, padx=5, pady=5)
        p.paneconfig(leftFrame, minsize=300)
        p.paneconfig(rightFrame, minsize=400)        
        
        return p
        
        
    def refreshRuns(self, e=None):
        """ Refresh the status of diplayed runs. """
        self.runsTree.update()
        self.updateRunsGraph()
        
    def createActionToolbar(self):
        """ Prepare the buttons that will be available for protocol actions. """
       
        self.actionList = [ACTION_EDIT, ACTION_COPY, ACTION_DELETE, ACTION_STEPS, 
                           ACTION_STOP, ACTION_CONTINUE, ACTION_RESULTS]
        self.actionButtons = {}
        
        def addButton(action, text, toolbar):
            btn = tk.Label(toolbar, text=text, image=self.getImage(ActionIcons[action]), 
                       compound=tk.LEFT, cursor='hand2')
            btn.bind('<Button-1>', lambda e: self._runActionClicked(action))
            return btn
        
        for action in self.actionList:
            self.actionButtons[action] = addButton(action, action, self.runsToolbar)
            
        for i, action in enumerate([ACTION_TREE, ACTION_REFRESH]):
            btn = addButton(action, action, self.allToolbar)
            btn.grid(row=0, column=i)
        
            
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
        displayAction(ACTION_RESULTS, 6, status != STATUS_RUNNING)
        
    def createProtocolsTree(self, parent):
        """Create the protocols Tree displayed in left panel"""
        tree = Tree(parent, show='tree')
        tree.column('#0', minwidth=300)
        tree.tag_configure('protocol', image=self.getImage('python_file.gif'))
        tree.tag_bind('protocol', '<Double-1>', self._protocolItemClick)
        tree.tag_configure('protocol_base', image=self.getImage('class_obj.gif'))
        f = tkFont.Font(family='verdana', size='10', weight='bold')
        tree.tag_configure('section', font=f)
        tree.grid(row=0, column=0, sticky='news')
        return tree
    
    def updateProtocolsTree(self, protCfg):
        self.protCfg = protCfg
        self.protTree.clear()
        self.protTree.unbind('<<TreeviewOpen>>')
        self.protTree.unbind('<<TreeviewClose>>')
        self.protTreeItems = {}
        populateTree(self, self.protTree, self.protTreeItems, '', self.protCfg)
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
    
    def createRunsGraph(self, parent):
        self.runsGraph = Canvas(parent, width=400, height=400)
        self.runsGraph.onClickCallback = self._runItemClick
        self.runsGraph.onDoubleClickCallback = self._runItemDoubleClick
        parent.grid_columnconfigure(0, weight=1)
        parent.grid_rowconfigure(0, weight=1)
        
        self.updateRunsGraph()

    def updateRunsGraph(self):      
        g = self.project.getRunsGraph(refresh=False)
        #g.printDot()
        lt = LevelTree(g)
        self.runsGraph.clear()
        lt.setCanvas(self.runsGraph)
        lt.paint(self.createRunItem)
        self.runsGraph.updateScrollRegion()
        
    def createRunItem(self, canvas, node, y):
        """ If not nodeBuildFunc is specified, this one will be used by default."""
        nodeText = node.label
        textColor = 'black'
        color = '#ADD8E6' #Lightblue
            
        if node.run:
            status = node.run.status.get(STATUS_FAILED)
            nodeText = nodeText + '\n' + status
            color = STATUS_COLORS[status]
        
        return self.runsGraph.createTextbox(nodeText, 100, y, bgColor=color, textColor=textColor)
        
    
    def switchRunsView(self):
        self.showGraph = not self.showGraph
        
        if self.showGraph:
            show = self.runsGraph.frame
            hide = self.runsTree
        else:
            show = self.runsTree
            hide = self.runsGraph.frame
            
        hide.grid_remove()
        show.grid(row=0, column=0, sticky='news')
        self.settings.graphView.set(self.showGraph)
        self.settings.write()
        
    
    def _protocolItemClick(self, e=None):
        protClassName = self.protTree.getFirst().split('.')[-1]
        protClass = emProtocolsDict.get(protClassName)
        prot = protClass()
        prot.mapper = self.project.mapper
        self._openProtocolForm(prot)
        
    def _runItemClick(self, e=None):
        # Get last selected item for tree or graph
        if self.showGraph:
            prot = e.node.run
        else:
            prot = self.project.mapper.selectById(int(self.runsTree.getFirst()))
        
        if prot is not None:
            prot.mapper = self.project.mapper
            self.selectedProtocol = prot
            self.updateActionToolbar()
            self._fillData()
            self._fillSummary()
        else:
            pass #TODO: implement what to do
        
    def _runItemDoubleClick(self, e=None):
        self._runActionClicked(ACTION_EDIT)
        
    def _openProtocolForm(self, prot):
        """Open the Protocol GUI Form given a Protocol instance"""
        hosts = [host.getLabel() for host in self.settings.getHosts()]
        w = FormWindow("Protocol Run: " + prot.getClassName(), prot, 
                       self._executeSaveProtocol, self.windows,
                       hostList=hosts)
        w.show(center=True)
        
    def _browseRunData(self):
        provider = ProtocolTreeProvider(self.selectedProtocol)
        window = BrowserWindow("Protocol data", provider, self.windows, icon=self.icon)
        window.itemConfig(self.selectedProtocol, open=True)  
        window.show()
        
    def _fillData(self):
        provider = RunIOTreeProvider(self.selectedProtocol, self.project.mapper)
        self.infoTree.setProvider(provider)
        #self.infoTree.itemConfig(self.selectedProtocol, open=True)  
        
    def _fillSummary(self):
        self.summaryText.clear()
        self.summaryText.addText(self.selectedProtocol.summary())
        
    def _executeSaveProtocol(self, prot, onlySave=False):
        if onlySave:
            self.project.saveProtocol(prot)
        else:
            self.project.launchProtocol(prot)
        self.runsTree.after(1000, self.runsTree.update)
        
    def _continueProtocol(self, prot):
        self.project.continueProtocol(prot)
        self.runsTree.after(1000, self.runsTree.update)        
        
    def _deleteProtocol(self, prot):
        if askYesNo("Confirm DELETE", "<ALL DATA> related to this <protocol run> will be <DELETED>. \n"
                    "Do you really want to continue?", self.root):
            self.project.deleteProtocol(prot)
            self.runsTree.update()
            
    def _stopProtocol(self, prot):
        if askYesNo("Confirm STOP", "Do you really want to <STOP> this run?", self.root):
            self.project.stopProtocol(prot)
            self.runsTree.update()        

    def _analyzeResults(self, prot):
        from pyworkflow.em import findViewers
        from pyworkflow.viewer import DESKTOP_TKINTER

        self.viewers = findViewers(prot.getClassName(), DESKTOP_TKINTER)
        self.viewer = self.viewers[0]()
        self.viewer.project = self.project
        self.viewer.visualize(prot, windows=self.windows)
        
                
    def _runActionClicked(self, action):
        prot = self.selectedProtocol
        if prot:
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
 
        # Following actions do not need a select run
        if action == ACTION_TREE:
            self.switchRunsView()
        elif action == ACTION_REFRESH:
            self.refreshRuns()
    
