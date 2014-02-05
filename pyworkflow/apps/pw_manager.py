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
from inspect import isclass
from pyworkflow.utils.utils import prettyDate
from pyworkflow.apps.config import loadSettings
from pyworkflow.utils.properties import Message
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
from pyworkflow.manager import Manager
from pyworkflow.mapper import SqliteMapper, XmlMapper
from pyworkflow.protocol import *
from pyworkflow.protocol.params import *
import config
from pyworkflow.em import *
import pyworkflow.gui as gui
from pyworkflow.gui.widgets import Button
from pyworkflow.gui.text import TaggedText
from pyworkflow.gui.dialog import askString, askYesNo, showError


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
            img = self.getImage(img)
        item = tree.insert(prefix, 'end', key, text=text, image=img)
        
        if level < 3:
            tree.item(item, open=True)
        if not obj.isEmpty() and obj.action.hasValue():
            prot = globals().get(obj.action.get(), None)
            if not prot is None:
                tree.item(item, image=self.getImage('step.gif'))
                for k, v in subclasses.iteritems():
                    if not v is prot and issubclass(v, prot):
                        tree.insert(item, 'end', item+k, text=k, image=self.getImage('python_file.gif'))
                        
            else:
                raise Exception("Class '%s' not found" % obj.action.get())
    else:
        key = prefix
    
    for sub in obj:
        populateTree(tree, key, sub, level+1)
    
#def getMapper(fn, classesDict):
#    """Select what Mapper to use depending on
#    the filename extension"""
#    if fn.endswith('.xml'):
#        return XmlMapper(fn, classesDict)
#    elif fn.endswith('.sqlite'):
#        return SqliteMapper(fn, classesDict)
#    return None
#
#
#    
#def loadConfig(config, name):
#    c = getattr(config, name) 
#    fn = getConfigPath(c.get())
#    if not os.path.exists(fn):
#        raise Exception('loadMenuConfig: menu file "%s" not found' % fn )
#    mapper = ConfigMapper(getConfigPath(fn), globals())
#    menuConfig = mapper.getConfig()
#    return menuConfig


class ManagerWindow(gui.WindowBase):
    """Windows to manage projects"""
    def __init__(self, **args):
        # Load global configuration
        settings = loadSettings(pw.SETTINGS)
        self.menuCfg = settings.getCurrentMenu()
        self.generalCfg = settings.getConfig()
        
        gui.WindowBase.__init__(self, Message.LABEL_PROJECTS, minsize=(750, 500), **args)
        self.manager = Manager()
        
        self.switchView(gui.VIEW_PROJECTS)
         

    def createProjectsView(self, parent):
        """ Create the Projects View.
        It has one panel and a menu:
            Top: menu
            Bottom: panel containing the projects list
        """
        return ProjectsView(parent, self)

            
class ProjectsView(tk.Frame):    
    def __init__(self, parent, windows, **args): 
        tk.Frame.__init__(self, parent, bg='white', **args)
        self.windows = windows
        self.manager = windows.manager
        self.root = windows.root
        
        #tkFont.Font(size=12, family='verdana', weight='bold')
        self.projNameFont = tkFont.Font(size=12, family='helvetica', weight='bold')
        self.projDateFont = tkFont.Font(size=8, family='helvetica')
        self.projDelFont = tkFont.Font(size=8, family='helvetica', weight='bold')
        self.manager = Manager()
        btn = Button(self, text=Message.LABEL_CREATE_PROJECT, font=self.projNameFont, 
                     command=self.createNewProject)
        btn.grid(row=0, column=0, sticky='nw', padx=10, pady=10)
        
        self.columnconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)
        text = TaggedText(self, width=40, height=15, bd=0, bg='white')
        text.grid(row=1, column=0, sticky='news')
      
        self.createProjectList(text)
        text.setReadOnly(True)
        self.text = text
    
    def createProjectList(self, text):
        """Load the list of projects"""
        r = 0
        text.clear()
        parent = tk.Frame(text, bg='white')    
        parent.columnconfigure(0, weight=1)
        colors = ['white', '#EAEBFF']
        for i, p in enumerate(self.manager.listProjects()):
            frame = self.createProjectLabel(parent, p, color=colors[i%2])
            frame.grid(row=r, column=0, padx=10, pady=5, sticky='new')
            r += 1
        text.window_create(tk.INSERT, window=parent)
        text.bindWidget(parent)
      
    def createProjectLabel(self, parent, projInfo, color):
        frame = tk.Frame(parent, bg=color)
        label = tk.Label(frame, text=projInfo.projName, anchor='nw', bg=color, 
                         justify=tk.LEFT, font=self.projNameFont, cursor='hand1')
        label.grid(row=0, column=0, padx=2, pady=2, sticky='nw')
        label.bind('<Button-1>', lambda e: self.openProject(projInfo.projName))
        dateLabel = tk.Label(frame, text='   ' + Message.LABEL_MODIFIED + ' ' + prettyDate(projInfo.mTime), 
                             font=self.projDateFont, bg=color)
        dateLabel.grid(row=1, column=0)
        delLabel = tk.Label(frame, text=Message.LABEL_DELETE_PROJECT, font=self.projDelFont, bg=color, cursor='hand1')
        delLabel.grid(row=1, column=1, padx=10)
        delLabel.bind('<Button-1>', lambda e: self.deleteProject(projInfo.projName))
        
        return frame
    
    def createNewProject(self, e=None):
        projName =  askString(Message.LABEL_CREATE_PROJECT, Message.TITLE_CREATE_PROJECT_NAME, self.root, 30)
        if not projName is None:
            self.manager.createProject(projName)
            self.createProjectList(self.text)
    
    def openProject(self, projName):
        projPath = self.manager.getProjectPath(projName)
        from pw_project import ProjectWindow
        projWindow = ProjectWindow(projPath, self.windows)
        projWindow.show()
          
    def deleteProject(self, projName):
        if askYesNo(Message.TITLE_DELETE_PROJECT, 
                     "Project *%s*" % projName + Message.MESSAGE_CREATE_PROJECT , self.root):
            self.manager.deleteProject(projName)
            self.createProjectList(self.text)
          
if __name__ == '__main__':
    ManagerWindow().show()
    
