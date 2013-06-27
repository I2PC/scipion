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
from pyworkflow.gui.tree import Tree, ObjectTreeProvider, DbTreeProvider
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

ActionIcons = {
    ACTION_EDIT:  'edit.gif',
    ACTION_COPY:  'copy.gif',
    ACTION_DELETE:  'delete.gif',
    ACTION_REFRESH:  'refresh.gif',
    ACTION_STEPS:  'run_steps.gif',
    ACTION_TREE:  'tree2.gif',
    ACTION_STOP: 'stop.gif',
    ACTION_CONTINUE: 'play.png'
               }


def populateTree(self, tree, prefix, obj, level=0):
    text = obj.text.get()
    if text:
        value = obj.value.get(text)
        key = '%s.%s' % (prefix, value)
        img = obj.icon.get('')
        tag = obj.tag.get('')
            
        if len(img):
            img = self.getImage(img)
        item = tree.insert(prefix, 'end', key, text=text, image=img, tags=(tag))
        
        if level < 2:
            tree.item(item, open=True)
        if obj.value.hasValue() and tag == 'protocol_base':
            protClassName = value.split('.')[-1] # Take last part
            prot = emProtocolsDict.get(protClassName, None)
            if prot is not None:
                tree.item(item, image=self.getImage('class_obj.gif'))
                for k, v in emProtocolsDict.iteritems():
                    if not v is prot and issubclass(v, prot):
                        key = '%s.%s' % (item, k)
                        tree.insert(item, 'end', key, text=k, tags=('protocol'))
                        
            else:
                raise Exception("Class '%s' not found" % obj.value.get())
    else:
        key = prefix
    
    for sub in obj:
        populateTree(self, tree, key, sub, level+1)
    
def getMapper(fn, classesDict):
    """Select what Mapper to use depending on
    the filename extension"""
    if fn.endswith('.xml'):
        return XmlMapper(fn, classesDict)
    elif fn.endswith('.sqlite'):
        return SqliteMapper(fn, classesDict)
    return None

    
def loadConfig(config, name):
    c = getattr(config, name) 
    fn = getConfigPath(c.get())
    if not os.path.exists(fn):
        raise Exception('loadMenuConfig: menu file "%s" not found' % fn )
    mapper = ConfigMapper(getConfigPath(fn), globals())
    menuConfig = mapper.getConfig()
    return menuConfig


class RunsTreeProvider(TreeProvider):
    """Provide runs info to populate tree"""
    def __init__(self, mapper, actionFunc):
        self.actionFunc = actionFunc
        self.getObjects = lambda: mapper.selectByClass('Protocol')
        
    def getColumns(self):
        return [('Run', 250), ('State', 100), ('Time', 100)]
    
    def getObjectInfo(self, obj):
        return {'key': obj.getObjId(),
                'text': obj.getRunName(),
                'values': (obj.status.get(), obj.getElapsedTime())}
      
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
        self.viewer = XmippViewer()
        
    def show(self, obj):
        self.viewer.visualize(obj)
        
    def getObjectPreview(self, obj):
        desc = "<name>: " + obj.getName()
        
        return (None, desc)
    
    def getObjectActions(self, obj):
        if isinstance(obj, Pointer):
            obj = obj.get()
            
        if isinstance(obj, SetOfMicrographs):
            return [('Open Micrographs with Xmipp', lambda: self.viewer.visualize(obj))]
        if isinstance(obj, SetOfImages):
            return [('Open Images with Xmipp', lambda: self.viewer.visualize(obj))]
        if isinstance(obj, XmippClassification2D):
            return [('Open Classification2D with Xmipp', lambda: self.viewer.visualize(obj))]
        return []   
    
    def getObjectInfo(self, obj):
        info = ObjectTreeProvider.getObjectInfo(self, obj)
        attrName = obj.getLastName()
        if hasattr(self.protocol, attrName):
            if isinstance(obj, Pointer) and obj.hasValue():
                info['image'] = 'db_input.gif'
            else:
                if (self.protocol._definition.hasParam(attrName) or
                    attrName in ['numberOfMpi', 'numberOfThreads']):
                    info['parent'] = self.params
                elif attrName in self.statusList:
                    if info['parent'] is self.protocol:
                        info['parent'] = self.status
                    
            if attrName.startswith('output'):# in self.protocol._outputs:
                info['image'] = 'db_output.gif'
        if obj is self.params or obj is self.status:
            info['parent'] = self.protocol
        return info     
    
    def _getChilds(self, obj):
        childs = ObjectTreeProvider._getChilds(self, obj)
        if obj is self.protocol:
            childs.insert(0, self.status)
            childs.insert(1, self.params)
        return childs
    

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
    
    def show(self, obj):
        self.viewer.visualize(obj)
        
    def getObjectPreview(self, obj):
        desc = "<name>: " + obj.getName()
        
        return (None, desc)
    
    def getObjectActions(self, obj):
        if isinstance(obj, Pointer):
            obj = obj.get()
            
        if isinstance(obj, SetOfMicrographs):
            return [('Open Micrographs with Xmipp', lambda: self.viewer.visualize(obj))]
        if isinstance(obj, SetOfImages):
            return [('Open Images with Xmipp', lambda: self.viewer.visualize(obj))]
        if isinstance(obj, XmippClassification2D):
            return [('Open Classification2D with Xmipp', lambda: self.viewer.visualize(obj))]
        return []  
    
    def getObjectInfo(self, obj):
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
    
VIEW_PROTOCOLS = 'Protocols'
VIEW_DATA = 'Data'
VIEW_HOSTS = 'Hosts'
   
class ProjectWindow(gui.Window):
    def __init__(self, path, master=None):
        # Load global configuration
        self.projName = 'Project: ' + basename(path)
        self.projPath = path
        self.loadProject()
        self.loadProjectConfig()
        self.icon = self.generalCfg.icon.get()
        self.selectedProtocol = None
        self.showGraph = False
        
        gui.Window.__init__(self, self.projName, master, icon=self.icon, minsize=(900,500))
        
        content = tk.Frame(self.root)
        content.columnconfigure(0, weight=1)
        content.rowconfigure(1, weight=1)
        content.grid(row=0, column=0, sticky='news')
        self.content = content
        
        self.createMainMenu(self.menuCfg)
        
        header = self.createHeaderFrame(content)
        header.grid(row=0, column=0, sticky='new')
        
        self.view, self.viewWidget = None, None
        self.viewFuncs = {VIEW_PROTOCOLS: self.createProtocolsView,
                          VIEW_DATA: self.createDataView,
                          VIEW_HOSTS: self.createHostsView
                          }
        self.switchView(VIEW_PROTOCOLS)

    def createHeaderFrame(self, parent):
        """ Create the Header frame at the top of the windows.
        It has (from left to right):
            - Main application Logo
            - Project Name
            - View selection combobox
        """
        header = tk.Frame(parent, bg='white')        
        header.columnconfigure(1, weight=1)
        header.columnconfigure(2, weight=1)
        # Create the SCIPION logo label
        logoImg = self.getImage(self.generalCfg.logo.get())
        logoLabel = tk.Label(header, image=logoImg, 
                             borderwidth=0, anchor='nw', bg='white')
        logoLabel.grid(row=0, column=0, sticky='nw', padx=5)
        # Create the Project Name label
        self.projNameFont = tkFont.Font(size=12, family='verdana', weight='bold')
        projLabel = tk.Label(header, text=self.projName, font=self.projNameFont,
                             borderwidth=0, anchor='nw', bg='white')
        projLabel.grid(row=0, column=1, sticky='sw', padx=(20, 5), pady=10)
        # Create view selection frame
        viewFrame = tk.Frame(header, bg='white')
        viewFrame.grid(row=0, column=2, sticky='se', padx=5, pady=10)
        viewLabel = tk.Label(viewFrame, text='View:', bg='white')
        viewLabel.grid(row=0, column=0, padx=5)
        self.viewVar = tk.StringVar()
        self.viewVar.set(VIEW_PROTOCOLS)
        viewCombo = ttk.Combobox(viewFrame, textvariable=self.viewVar, state='readonly')
        viewCombo['values'] = [VIEW_PROTOCOLS, VIEW_DATA, VIEW_HOSTS]
        viewCombo.grid(row=0, column=1)
        viewCombo.bind('<<ComboboxSelected>>', self._viewComboSelected)
        
        return header
    
    def getHostsMapper(self):
        return self.project.hostsMapper
    
    def _viewComboSelected(self, e=None):
        if self.viewVar.get() != self.view:
            self.switchView(self.viewVar.get())
        
    def switchView(self, newView):
        # Destroy the previous view if existing:
        if self.viewWidget:
            self.viewWidget.grid_forget()
            self.viewWidget.destroy()
        # Create the new view
        self.viewWidget = self.viewFuncs[newView](self.content)
        # Grid in the second row (1)
        self.viewWidget.grid(row=1, column=0, sticky='news')
        self.view = newView
        
    def createProtocolsView(self, parent):
        """ Create the Protocols View for the Project.
        It has two panes:
            Left: containing the Protocol classes tree
            Right: containing the Runs list
        """
        from pw_project_viewprotocols import ProtocolsView
        p = ProtocolsView(parent, self)
        return p
        
    def createDataView(self, parent):
        dataFrame = tk.Frame(parent)
        dataLabel = tk.Label(dataFrame, text='DATA VIEW not implemented.',
                             font=self.projNameFont)
        dataLabel.grid(row=0, column=0, padx=50, pady=50)
        return dataFrame

    def createHostsView(self, parent):
        from pw_project_viewhosts import HostsView
        return HostsView(parent, self)
    
    def loadProject(self):
        self.project = Project(self.projPath)
        self.project.load()
                
    def loadProjectConfig(self):
        self.configMapper = ConfigMapper(getConfigPath('configuration.xml'), globals())
        self.generalCfg = self.configMapper.getConfig()
        self.menuCfg = loadConfig(self.generalCfg, 'menu')
        self.protCfg = loadConfig(self.generalCfg, 'protocols')


if __name__ == '__main__':
    from pyworkflow.manager import Manager
    if len(sys.argv) > 1:
        manager = Manager()
        projName = os.path.basename(sys.argv[1])
        projPath = manager.getProjectPath(projName)
        projWindow = ProjectWindow(projPath)
        projWindow.show()
    else:
        print "usage: pw_project.py PROJECT_NAME"
