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
        projName = self.__class__.__name__
        manager = Manager()
        if manager.hasProject(projName):
            project = manager.loadProject(projName)
        else:
            project = manager.createProject(projName)
            
            # Create a new project
            self.outputPath = project.path
            self.ds = tests.DataSet.getDataSet('xmipp_tutorial')
            
            project.loadProtocols(self.ds.getFile('workflow.json'))
            
            # Use graph view as default
            settings = project.getSettings()
            settings.setRunsView(1) # graph view
            settings.write()
            
            # Update the path of imports
            protImportMics = project.getProtocolsByClass('ProtImportMicrographs')[0]
            protImportMics.filesPath.set(self.ds.getFile('allMics'))
            project.saveProtocol(protImportMics)
            
            protImportVol = project.getProtocolsByClass('ProtImportVolumes')[0]
            protImportVol.filesPath.set(self.ds.getFile('vol110'))
            project.saveProtocol(protImportVol)

        self.project = project
        
        
class Tutorial2D(tests.BaseTest):
    
    def __init__(self):
        projName = self.__class__.__name__
        manager = Manager()
        if manager.hasProject(projName):
            project = manager.loadProject(projName)
        else:
            project = manager.createProject(projName)
            
            # Create a new project
            self.outputPath = project.path
            self.ds = tests.DataSet.getDataSet('2d_analysis')
            
            project.loadProtocols(self.ds.getFile('workflow.json'))
            
            # Use graph view as default
            #settings = project.getSettings()
            #settings.setRunsView(1) # graph view
            #settings.write()
            
            # Update the path of imports
            protImportParticles = project.getProtocolsByClass('ProtImportParticles')[0]
            protImportParticles.filesPath.set(self.ds.getFile('allMics'))
            project.saveProtocol(protImportParticles)
            
            protImportVol = project.getProtocolsByClass('ProtImportVolumes')[0]
            protImportVol.filesPath.set(self.ds.getFile('vol1'))
            project.saveProtocol(protImportVol)

        self.project = project
    
    
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
        
        projWindow = ProjectWindow(tutorial.project.getName())
        projWindow.show()
    else:
        pass #TODO: print list of all tutorials
