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
from gui.gui import getImage
from inspect import isclass
"""
Main project window application
"""
import os
from os.path import join, exists
import sys
        
import Tkinter as tk
import ttk
import tkFont

import pyworkflow as pw
from pyworkflow.object import *
from pyworkflow.em import *
from pyworkflow.mapper import SqliteMapper, XmlMapper
from pyworkflow.project import Project
import gui
from gui.widgets import Tree
from protocol import *
from protocol.params import *
from config import *


def loadSubclasses():
    """This function will try to find
    all sub-classes of className located in
    the pyworkflow/protocol/packages/* path and made them
    available in the menu"""
    path = join(pw.HOME, 'em', 'packages')
    sys.path.append(path)
    folders = os.listdir(path)
    #print "path: ", path
    subclasses = {}
    for f in folders:
        if exists(join(path, f, '__init__.py')):
            m = __import__(f)
            for k, v in m.__dict__.iteritems():
                if isclass(v) and issubclass(v, Protocol):
                    subclasses[k] = v
    return subclasses
    
subclasses =  loadSubclasses()
    
def populateTree(tree, prefix, obj, level=0):
    text = obj.text.get()
    if text:
        key = '%s.%s' % (prefix, text)
        img = obj.icon.get()
        if img is None:
            img = ''
            pass
        else:
            img = gui.getImage(img)
        item = tree.insert(prefix, 'end', key, text=text, image=img)
        
        if level < 3:
            tree.item(item, open=True)
        if not obj.isEmpty() and obj.action.hasValue():
            prot = globals().get(obj.action.get(), None)
            if not prot is None:
                tree.item(item, image=gui.getImage('class_obj.gif'))
                for k, v in subclasses.iteritems():
                    if not v is prot and issubclass(v, prot):
                        tree.insert(item, 'end', item+k, text=k, image=gui.getImage('python_file.gif'))
                        
            else:
                raise Exception("Class '%s' not found" % obj.action.get())
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

def createHistoryTree(parent, path):
    columns = ('State', 'Modified')
    tree = Tree(parent, columns=columns)
    for c in columns:
        tree.column(c, anchor='e', width=100)
        tree.heading(c, text=c) 
    tree.column('#0', width=250)
    tree.heading('#0', text='Run')
    #tree.bind('<<TreeviewSelect>>', self.selectTreeRun)
    #tree.bind('<Double-1>', lambda e:self.runButtonClick('ACTION_DEFAULT'))
    #tree.bind("<Button-3>", self.onRightClick)
    proj = Project(path)
    proj.load()
    #mapper = SqliteMapper('kk.sqlite', globals())
    objList = proj.mapper.getAll()
    
    c = 0
    from pyworkflow.utils import prettyDate
    for obj in objList:
        t = '%s.%02d' % (obj.getClassName(), obj.getId())
        tree.insert('',  'end', str(c), text=t, values=(obj.status.get(), obj.endTime.get()))
        c += 1
    return tree

def createProjectGUI(path):
    # Load global configuration
    mapper = ConfigXmlMapper(getConfigPath('configuration.xml'), globals())
    config = mapper.getConfig()

    window = gui.Window("Project window", icon='scipion_bn.xbm', minsize=(900,500))
    
    parent = window.root
    menuConfig = loadConfig(config, 'menu')
    createMainMenu(parent, menuConfig)
    p = tk.PanedWindow(parent, orient=tk.HORIZONTAL)
    # first pane, which would get widgets gridded into it:
    f1 = tk.Frame(p)
    f1.columnconfigure(0, weight=1)
    f1.rowconfigure(1, weight=1)
    logo = gui.getImage('scipion_logo.gif')
    label = tk.Label(f1, image=logo, borderwidth=0, 
                     anchor='nw', bg='white')
    label.grid(row=0, column=0, sticky='new', padx=5)
    
    lf = ttk.Labelframe(f1, text=' Protocols ', width=300, height=500)
    lf.grid(row=1, column=0, sticky='news', padx=5, pady=5)
    # second pane
    #f2 = tk.Frame(p)
    #f2.grid(row=0, column=0, sticky='news')
    
    
    f2 = ttk.Labelframe(p, text=' History ', width=500, height=500)
    #lf.grid(row=1, column=0, sticky='news')
    
    tree = createHistoryTree(f2, path)
    tree.grid(row=0, column=0, sticky='news')
    gui.configureWeigths(f2)
    
    p.add(f1, padx=5, pady=5)
    p.add(f2, padx=5, pady=5)
    p.paneconfig(f1, minsize=300)
    p.paneconfig(f2, minsize=300)
    
    p.grid(row=0, column=0, sticky='news')
    
    #gui.configureWeigths(f1)
    gui.configureWeigths(lf)
    tree = Tree(lf, show='tree')
    #tree.heading(c, text=c)
    #tree.heading('#0', text='Object')
    tree.column('#0', minwidth=300)
    objList = loadConfig(config, 'protocols')
    populateTree(tree, '', objList)
    tree.grid(row=0, column=0, sticky='news')   
    
    
    window.show()            
    
if __name__ == '__main__':
    import sys
    path = '.'
    if len(sys.argv) > 1:
        path = sys.argv[1]
    createProjectGUI(path)
