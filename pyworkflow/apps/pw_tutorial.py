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
Launch main project window 
"""

import sys

import pyworkflow.tests as tests
import pyworkflow.em as em
from pyworkflow.manager import Manager
from pyworkflow.gui.project import ProjectWindow


class TutorialIntro(tests.BaseTest):
    
    def __init__(self):
        # Create a new project
        tests.setupTestProject(self.__class__)
        self.ds = tests.DataSet.getDataSet('xmipp_tutorial')
        
    def getProjectName(self):
        return self.proj.getName()
    
    def setup(self):
        """ Run an Import particles protocol. """
        self.proj.loadProtocols(self.ds.getFile('workflow.json'))
        settings = self.proj.getSettings()
        settings.setGraphView(True)
        settings.write()
        
        # Update the path of imports
        protImportMics = self.proj.getProtocolsByClass('ProtImportMicrographs')[0]
        protImportMics.filesPath.set('kkk.mrc')
        self.proj.saveProtocol(protImportMics)
        
    
if __name__ == '__main__':

    # Add callback for remote debugging if available.
    try:
        from rpdb2 import start_embedded_debugger
        from signal import signal, SIGUSR2
        signal(SIGUSR2, lambda sig, frame: start_embedded_debugger('a'))
    except ImportError:
        pass

    if len(sys.argv) > 1:
        manager = Manager()
        tutorialName = sys.argv[1]
        
        tutorial = TutorialIntro()
        tutorial.setup()
        
        projPath = manager.getProjectPath(tutorial.getProjectName())
        projWindow = ProjectWindow(projPath)
        projWindow.show()
    else:
        pass #TODO: print list of all tutorials
