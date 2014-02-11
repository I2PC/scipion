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
from pyworkflow.utils.properties import Message

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
   
   
class ProjectWindow(gui.WindowBase):
    def __init__(self, path, master=None):
        # Load global configuration
        self.projName = Message.LABEL_PROJECT + basename(path)
        self.projPath = path
        self.loadProject()
        self.icon = self.generalCfg.icon.get()
        self.selectedProtocol = None
        self.showGraph = False
        
        gui.WindowBase.__init__(self, self.projName, master, icon=self.icon, minsize=(900,500))
        
        self.switchView(gui.VIEW_PROTOCOLS)
    
    def getSettings(self):
        return self.settings
    
    def saveSettings(self):
        self.settings.write()
        
    def _onClosing(self):
        try:
            self.saveSettings()
        except Exception, ex:
            print Message.NO_SAVE_SETTINGS + str(ex) 
        gui.Window._onClosing(self)
        
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
        return DataView(parent, self)
#        dataFrame = tk.Frame(parent)
#        dataLabel = tk.Label(dataFrame, text='DATA VIEW not implemented.',
#                             font=self.projNameFont)
#        dataLabel.grid(row=0, column=0, padx=50, pady=50)
#        return dataFrame

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
                
                
class DataView(tk.Frame):
    def __init__(self, parent, windows, **args):
        tk.Frame.__init__(self, parent, **args)
        dataLabel = tk.Label(self, text='DATA VIEW not implemented.',)
        dataLabel.grid(row=0, column=0, padx=50, pady=50)

if __name__ == '__main__':
    from pyworkflow.manager import Manager
    if len(sys.argv) > 1:
        manager = Manager()
        projName = os.path.basename(sys.argv[1])
        if projName == 'last': # Get last project
            projName = manager.listProjects()[0].projName
            
        projPath = manager.getProjectPath(projName)
        projWindow = ProjectWindow(projPath)
        projWindow.show()
    else:
        print "usage: pw_project.py PROJECT_NAME"
