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

import pyworkflow as pw
from pyworkflow.object import *
from pyworkflow.em import *
from pyworkflow.protocol import *
from pyworkflow.protocol.params import *
from pyworkflow.mapper import SqliteMapper, XmlMapper
from pyworkflow.project import Project

import pyworkflow.gui as gui
from pyworkflow.gui import getImage
from pyworkflow.gui.tree import Tree
from pyworkflow.gui.form import FormWindow
from config import *


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
        
        if level < 3:
            tree.item(item, open=True)
        if obj.value.hasValue() and tag == 'protocol_base':
            protClassName = value.split('.')[-1] # Take last part
            prot = emProtocolsDict.get(protClassName, None)
            print "found protocol_base: ", protClassName
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
    mapper = ConfigXmlMapper(getConfigPath(fn), globals())
    menuConfig = mapper.getConfig()
    return menuConfig



class ProjectWindow(gui.Window):
    def __init__(self, path, master=None):
        # Load global configuration
        self.projName = 'Project: ' + basename(path)
        self.projPath = path
    
        gui.Window.__init__(self, self.projName, master, icon='scipion_bn.xbm', minsize=(900,500))
        
        self.configMapper = ConfigXmlMapper(getConfigPath('configuration.xml'), globals())
        parent = self.root

        self.loadProjectConfig()
        self.createMainMenu(self.menuCfg)
        
        # The main layout will be two panes, 
        # At the left containing the Protocols
        # and the right containing the Runs
        p = tk.PanedWindow(parent, orient=tk.HORIZONTAL)
        
        # Left pane, contains SCIPION Logo and Protocols Pane
        leftFrame = tk.Frame(p)
        leftFrame.columnconfigure(0, weight=1)
        leftFrame.rowconfigure(1, weight=1)
        logo = self.getImage('scipion_logo.gif')
        label = tk.Label(leftFrame, image=logo, borderwidth=0, 
                         anchor='nw', bg='white')
        label.grid(row=0, column=0, sticky='new', padx=5)

        # Protocols Tree Pane        
        protFrame = ttk.Labelframe(leftFrame, text=' Protocols ', width=300, height=500)
        protFrame.grid(row=1, column=0, sticky='news', padx=5, pady=5)
        gui.configureWeigths(protFrame)
        self.protTree = self.createProtocolsTree(protFrame)
        
        # Runs history Pane
        runsFrame = ttk.Labelframe(p, text=' History ', width=500, height=500)
        self.runsTree = self.createRunsTree(runsFrame)
        
        gui.configureWeigths(runsFrame)
        
        # Add sub-windows to PanedWindows
        p.add(leftFrame, padx=5, pady=5)
        p.add(runsFrame, padx=5, pady=5)
        p.paneconfig(leftFrame, minsize=300)
        p.paneconfig(runsFrame, minsize=300)        
        p.grid(row=0, column=0, sticky='news')
        
        # Event bindings
        self.root.bind("<F5>", self.updateRunsTree)
        
    def loadProjectConfig(self):
        self.project = Project(self.projPath)
        self.project.load()
        self.generalCfg = self.configMapper.getConfig()
        self.menuCfg = loadConfig(self.generalCfg, 'menu')
        self.protCfg = loadConfig(self.generalCfg, 'protocols')
        
    def createProtocolsTree(self, parent):
        """Create the protocols Tree displayed in left panel"""
        tree = Tree(parent, show='tree')
        tree.column('#0', minwidth=300)
        tree.tag_configure('protocol', image=self.getImage('python_file.gif'))
        tree.tag_bind('protocol', '<Double-1>', self.protocolItemClick)
        tree.tag_configure('protocol_base', image=self.getImage('class_obj.gif'))
        f = tkFont.Font(family='verdana', size='10', weight='bold')
        tree.tag_configure('section', font=f)
        populateTree(self, tree, '', self.protCfg)
        tree.grid(row=0, column=0, sticky='news')
        return tree
        
    def createRunsTree(self, parent):
        columns = ('State', 'Modified')
        tree = Tree(parent, columns=columns)
        for c in columns:
            tree.column(c, anchor='e', width=100)
            tree.heading(c, text=c) 
        tree.column('#0', width=250)
        tree.heading('#0', text='Run')
        #tree.bind('<<TreeviewSelect>>', self.selectTreeRun)
        tree.grid(row=0, column=0, sticky='news')
        tree.bind('<Double-1>', self.runItemClick)
        self.updateRunsTree(tree)
        #tree.bind("<Button-3>", self.onRightClick)
        return tree
    
    def updateRunsTree(self, tree):
        tree.clear()
        objList = self.project.mapper.selectAll()
        for obj in objList:
            t = '%s.%02d' % (obj.getClassName(), obj.getId())
            tree.insert('',  'end', obj.getId(), text=t, values=(obj.status.get(), obj.endTime.get()))
        
        #tree.after(1000, self.updateRunsTree, tree)
    
    def protocolItemClick(self, e=None):
        print self.protTree.getFirst()
        protClassName = self.protTree.getFirst().split('.')[-1]
        protClass = emProtocolsDict.get(protClassName)
        prot = protClass()
        prot.mapper = self.project.mapper
        self._openProtocolForm(prot)
        
    def runItemClick(self, e=None):
        prot = self.project.mapper.selectById(int(self.runsTree.getFirst()))
        prot.mapper = self.project.mapper
        prot.setId(None) # Force insert a new run
        #prot.printAll()
        self._openProtocolForm(prot)
        
    def _openProtocolForm(self, prot):
        """Open the Protocol GUI Form given a Protocol instance"""
        w = FormWindow("Protocol Run: " + prot.getClassName(), prot, self)
        w.show(center=True)
        
    
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
