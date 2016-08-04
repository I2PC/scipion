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

import os
import sys
from collections import OrderedDict

import pyworkflow as pw
import pyworkflow.tests as tests
from pyworkflow.manager import Manager
import pyworkflow.utils as pwutils
from pyworkflow.gui.project import ProjectWindow


def getWorkflow(workflow):
    """ Return the full workflow path from
    the Scipion folder + config/workflows/
    """
    return pw.getConfigPath('workflows', workflow)
    

class Tutorial():
    """ Base class to implement some common functionalities. """
    def __init__(self):
        projName = self.__class__.__name__
        manager = Manager()
        if manager.hasProject(projName):
            self.project = manager.loadProject(projName)
        else:
            self.project = manager.createProject(projName)
            # Use graph view as default
            settings = self.project.getSettings()
            settings.setRunsView(1) # graph view
            settings.write()
            self.loadWorkflow()
    

class TutorialIntro(Tutorial):
    
    def loadWorkflow(self):            
        # Create a new project
        self.ds = tests.DataSet.getDataSet('xmipp_tutorial')
        self.project.loadProtocols(getWorkflow('workflow_tutorial_intro.json'))
        
        # Update the path of imports
        protImportMics = self.project.getProtocolsByClass('ProtImportMicrographs')[0]
        protImportMics.filesPath.set(self.ds.getFile('allMics'))
        self.project.saveProtocol(protImportMics)
        
        protImportVol = self.project.getProtocolsByClass('ProtImportVolumes')[0]
        protImportVol.filesPath.set(self.ds.getFile('vol110'))
        self.project.saveProtocol(protImportVol)


class TutorialBetagal(Tutorial):
    
    def loadWorkflow(self):            
        # Update the path of imports
        self.project.loadProtocols(getWorkflow('workflow_betagal1.json'))


ALL_TUTORIALS = OrderedDict([('intro', TutorialIntro),
                             ('betagal', TutorialBetagal)])

if __name__ == '__main__':

    def printUsage(msg):
        if msg:
            print "ERROR: ", msg
            
        print "\nUSAGE: scipion tutorial [TUTORIAL_NAME]"
        print "\nwhere TUTORIAL_NAME can be:"
        print "\n".join([' %s' % k for k in ALL_TUTORIALS.keys()])
        
    if pwutils.envVarOn('SCIPION_DEBUG'):
        # Add callback for remote debugging if available.
        try:
            from rpdb2 import start_embedded_debugger
            from signal import signal, SIGUSR2
            signal(SIGUSR2, lambda sig, frame: start_embedded_debugger('a'))
        except ImportError:
            pass

    if len(sys.argv) == 2:
        manager = Manager()
        tutorialName = sys.argv[1]
        
        if not tutorialName in ALL_TUTORIALS:
            printUsage("Invalid tutorial '%s'." % tutorialName)
        else:
            # Instanciate the proper tutorial class
            tutorial = ALL_TUTORIALS[tutorialName]()
        
            projWindow = ProjectWindow(tutorial.project.getName())
            projWindow.show()
    else:
        msg = 'Too many arguments.' if len(sys.argv) > 2 else ''
        printUsage(msg)
