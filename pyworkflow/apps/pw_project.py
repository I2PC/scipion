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
import pyworkflow.apps.config as config

from config import *
from pw_browser import BrowserWindow

from pyworkflow.gui.plotter import Plotter
Plotter.setInteractive(True)   
    
VIEW_PROTOCOLS = 'Protocols'
VIEW_DATA = 'Data'
VIEW_HOSTS = 'Hosts'
   
   
class ProjectWindow(gui.Window):
    def __init__(self, path, master=None):
        # Load global configuration
        self.projName = 'Project: ' + basename(path)
        self.projPath = path
        self.loadProject()
        self.icon = self.generalCfg.icon.get()
        self.selectedProtocol = None
        self.showGraph = False
        
        gui.Window.__init__(self, self.projName, master, icon=self.icon, minsize=(900,500))
        
        content = tk.Frame(self.root)
        content.columnconfigure(0, weight=1)
        content.rowconfigure(1, weight=1)
        content.grid(row=0, column=0, sticky='news')
        self.content = content
        
        self.createMainMenu()
        
        header = self.createHeaderFrame(content)
        header.grid(row=0, column=0, sticky='new')
        
        self.view, self.viewWidget = None, None
        self.viewFuncs = {VIEW_PROTOCOLS: self.createProtocolsView,
                          VIEW_DATA: self.createDataView,
                          VIEW_HOSTS: self.createHostsView
                          }
        self.switchView(VIEW_PROTOCOLS)

    def createMainMenu(self):
        gui.Window.createMainMenu(self, self.menuCfg) 
        
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
        logoLabel.grid(row=0, column=0, sticky='nw', padx=5, pady=5)
        # Create the Project Name label
        self.projNameFont = tkFont.Font(size=-28, family='helvetica')
        projLabel = tk.Label(header, text=self.projName, font=self.projNameFont,
                             borderwidth=0, anchor='nw', bg='white', fg='#707070')
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
        
        # Create header line
#        headerLine = tk.Frame(header,bg='red', height=10, width=100%)
        headerLine = Canvas(header, height=-10, bg='red')
        headerLine.grid(row=1, column=0, columnspan=3, sticky='nw')
        
        return header
    
    def getSettings(self):
        return self.settings
    
    def saveSettings(self):
        self.settings.write()
        
    def _onClosing(self):
        self.saveSettings() 
        gui.Window._onClosing(self)
    
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
        
#    def handleResize(self):
#        print self._w, self._h
#        
#    def handleMove(self):
#        print self._x, self._y
        
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
        self.settings = self.project.getSettings()
        self.generalCfg = self.settings.getConfig()
        self.menuCfg = self.settings.getCurrentMenu()
        self.protCfg = self.settings.getCurrentProtocolMenu()
                


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
