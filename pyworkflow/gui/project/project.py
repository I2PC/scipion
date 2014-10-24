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
Main Project window implementation.
It is composed by three panels:
1. Left: protocol tree.
2. Right upper: VIEWS (Data/Protocols)
3. Summary/Details
"""

from os.path import basename

import pyworkflow as pw
from pyworkflow.manager import Manager
from pyworkflow.config import * # We need this to retrieve object from mapper
from pyworkflow.project import Project
from pyworkflow.gui import Message
from pyworkflow.gui.plotter import Plotter

from base import ProjectBaseWindow, VIEW_PROTOCOLS, VIEW_PROJECTS
   
   
   
class ProjectWindow(ProjectBaseWindow):
    """ Main window for working in a Project. """
    def __init__(self, path, master=None):
        # Load global configuration
        self.projName = Message.LABEL_PROJECT + basename(path)
        self.projPath = path
        self.loadProject()
        self.icon = self.generalCfg.icon.get()
        self.selectedProtocol = None
        self.showGraph = False
        
        Plotter.setInteractive(True)   
        
        ProjectBaseWindow.__init__(self, self.projName, master, icon=self.icon, minsize=(900,500))
        
        self.switchView(VIEW_PROTOCOLS)
    
    def createHeaderFrame(self, parent):
        """Create the header and add the view selection frame at the right."""
        header = ProjectBaseWindow.createHeaderFrame(self, parent)
        self.addViewList(header)
        return header

    def getSettings(self):
        return self.settings
    
    def saveSettings(self):
        self.settings.write()
        
    def _onClosing(self):
        try:
            self.saveSettings()
        except Exception, ex:
            print Message.NO_SAVE_SETTINGS + str(ex) 
        ProjectBaseWindow._onClosing(self)
     
    def loadProject(self):
        self.project = Project(self.projPath)
        self.project.load()
        self.settings = self.project.getSettings()
        self.generalCfg = self.settings.getConfig()
        self.menuCfg = self.settings.getCurrentMenu()
        self.protCfg = self.settings.getCurrentProtocolMenu()
        
        
class ProjectManagerWindow(ProjectBaseWindow):
    """ Windows to manage all projects. """
    def __init__(self, **args):
        # Load global configuration
        settings = ProjectSettings()
        self.menuCfg = settings.getCurrentMenu()
        self.generalCfg = settings.getConfig()
        
        ProjectBaseWindow.__init__(self, Message.LABEL_PROJECTS, minsize=(750, 500), **args)
        self.manager = Manager()
        
        self.switchView(VIEW_PROJECTS)

