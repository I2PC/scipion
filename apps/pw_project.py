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
from os.path import join, exists
from inspect import isclass
        
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
from pyworkflow.gui.widgets import Tree
from pyworkflow.gui.form import FormWindow
from config import *


def loadSubclasses():
    """This function will try to find
    all sub-classes of className located in
    the pyworkflow/protocol/packages/* path and made them
    available in the menu"""
    path = join(pw.HOME, 'em', 'packages')
    sys.path.append(path)
    folders = os.listdir(path)
    subclasses = {}
    for f in folders:
        if exists(join(path, f, '__init__.py')):
            m = __import__(f)
            print "imported: ", f
            for k, v in m.__dict__.iteritems():
                if isclass(v) and issubclass(v, Protocol):
                    subclasses[k] = v
    for k, v in globals().iteritems():
        if isclass(v) and issubclass(v, Protocol):
                    subclasses[k] = v
    return subclasses
    
# All protocol subclasses
subclasses =  loadSubclasses()

    
def populateTree(tree, prefix, obj, level=0):
    text = obj.text.get()
    if text:
        value = obj.value.get(text)
        key = '%s.%s' % (prefix, value)
        img = obj.icon.get('')
        tag = obj.tag.get('')
            
        if len(img):
            img = gui.getImage(img)
        item = tree.insert(prefix, 'end', key, text=text, image=img, tags=(tag))
        
        if level < 3:
            tree.item(item, open=True)
        if obj.value.hasValue() and tag == 'protocol_base':
            protName = value.split('.')[-1] # Take last part
            prot = subclasses.get(value, None)
            print "found protocol_base: ", value
            if not prot is None:
                tree.item(item, image=gui.getImage('class_obj.gif'))
                for k, v in subclasses.iteritems():
                    if not v is prot and issubclass(v, prot):
                        tree.insert(item, 'end', item+k, text=k, image=gui.getImage('python_file.gif'))
                        
            else:
                raise Exception("Class '%s' not found" % obj.value.get())
    else:
        key = prefix
    
    for sub in obj:
        populateTree(tree, key, sub, level+1)
    
def getMapper(fn, classesDict):
    """Select what Mapper to use depending on
    the filename extension"""
    if fn.endswith('.xml'):
        return XmlMapper(fn, classesDict)
    elif fn.endswith('.sqlite'):
        return SqliteMapper(fn, classesDict)
    return None

def addMenuChilds(root, menu, menuConfig):
    for sub in menuConfig:
        if len(sub):
            submenu = tk.Menu(root, tearoff=0)
            menu.add_cascade(label=sub.text.get(), menu=submenu)
            addMenuChilds(root, submenu, sub)
        else:
            menu.add_command(label=sub.text.get(), compound=tk.LEFT,
                             image=gui.getImage(sub.icon.get()))
 
def createMainMenu(root, menuConfig):
    menu = tk.Menu(root)
    addMenuChilds(root, menu, menuConfig)
    root.config(menu=menu)
    
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
        self.configMapper = ConfigXmlMapper(getConfigPath('configuration.xml'), globals())
        self.projName = 'Project: ' + basename(path)
        self.projPath = path
    
        gui.Window.__init__(self, self.projName, master, icon='scipion_bn.xbm', minsize=(900,500))
        
        parent = self.root

        self.loadProjectConfig()
        createMainMenu(parent, self.menuCfg)
        
        # The main layout will be two panes, 
        # At the left containing the Protocols
        # and the right containing the Runs
        p = tk.PanedWindow(parent, orient=tk.HORIZONTAL)
        
        # Left pane, contains SCIPION Logo and Protocols Pane
        leftFrame = tk.Frame(p)
        leftFrame.columnconfigure(0, weight=1)
        leftFrame.rowconfigure(1, weight=1)
        logo = gui.getImage('scipion_logo.gif')
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
        self.runsTree = self.createRunsTree(runsFrame, path)
        
        gui.configureWeigths(runsFrame)
        
        # Add sub-windows to PanedWindows
        p.add(leftFrame, padx=5, pady=5)
        p.add(runsFrame, padx=5, pady=5)
        p.paneconfig(leftFrame, minsize=300)
        p.paneconfig(runsFrame, minsize=300)        
        p.grid(row=0, column=0, sticky='news')
        
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
        tree.tag_configure('protocol', image=gui.getImage('python_file.gif'))
        tree.tag_bind('protocol', '<Double-1>', self.protocolItemClick)
        tree.tag_configure('protocol_base', image=gui.getImage('class_obj.gif'))
        f = tkFont.Font(family='verdana', size='10', weight='bold')
        tree.tag_configure('section', font=f)
        populateTree(tree, '', self.protCfg)
        tree.grid(row=0, column=0, sticky='news')
        return tree
        
    def createRunsTree(self, parent, projPath):
        columns = ('State', 'Modified')
        tree = Tree(parent, columns=columns)
        for c in columns:
            tree.column(c, anchor='e', width=100)
            tree.heading(c, text=c) 
        tree.column('#0', width=250)
        tree.heading('#0', text='Run')
        #tree.bind('<<TreeviewSelect>>', self.selectTreeRun)
        tree.bind('<Double-1>', self.runItemClick)
        #tree.bind("<Button-3>", self.onRightClick)
        objList = self.project.mapper.getAll()
        for obj in objList:
            t = '%s.%02d' % (obj.getClassName(), obj.getId())
            tree.insert('',  'end', obj.getId(), text=t, values=(obj.status.get(), obj.endTime.get()))
        tree.grid(row=0, column=0, sticky='news')
        return tree
    
    def protocolItemClick(self, e=None):
        print "protocol clicked"
        print self.protTree.getFirst()
        protName = self.protTree.getFirst().split('.')[-1]
        protClass = subclasses.get(protName)
        self._openProtocolForm(protClass())
        
    def runItemClick(self, e=None):
        print "run clicked"
        print self.runsTree.getFirst()
        prot = self.project.mapper.get(int(self.runsTree.getFirst()))
        #prot.printAll()
        self._openProtocolForm(prot)
        
    def _openProtocolForm(self, prot):
        """Open the Protocol GUI Form given a Protocol instance"""
        w = FormWindow(prot)
        w.show()
        
        
    
if __name__ == '__main__':
    import sys
    path = '.'
    if len(sys.argv) > 1:
        path = sys.argv[1]
    window = ProjectWindow(path)
    window.show()
