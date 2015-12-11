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

import os

from tests import BaseTest, DataSet
from pyworkflow.manager import Manager
from pyworkflow.protocol import MODE_RESTART

import pyworkflow.em as em

    
class TestProject(BaseTest):

    def test_create(self):
        """Test the list with several Complex"""
        cwd = os.getcwd()
        manager = Manager()
        projectName = "TestProjectCreate"
        project = manager.createProject(projectName)
        self.assertNotEqual(cwd, os.getcwd())
        self.assertEqual(project.getPath(), os.getcwd())
        
        os.chdir(cwd)
        manager.deleteProject(projectName)
            
    def test_createNoChdir(self):
        """Test the list with several Complex"""
        cwd = os.getcwd()
        manager = Manager()
        projectName = "TestProjectCreate"
        project = manager.createProject(projectName, chdir=False)
        dataset = DataSet.getDataSet('emx')

        prot1 = project.newProtocol(em.ProtImportParticles,
                               objLabel='from emx (coordinatesT1)',
                               importFrom=em.ProtImportParticles.IMPORT_FROM_EMX,
                               emxFile=dataset.getFile('coordinatesT1'),
                               alignType=3,
                               voltage=100,
                               magnification=10000,
                               samplingRate=2.46)
        project.launchProtocol(prot1, wait=True)
        
        # Execute a copy of prot1
        prot2 = project.copyProtocol(prot1)
        project.launchProtocol(prot2, wait=True)
        
        # Now restart prot2
        prot2.runMode.set(MODE_RESTART)
        project.launchProtocol(prot2, wait=True)       
        
        #self.assertEqual(cwd, project.getPath())
        self.assertEqual(cwd, os.getcwd())

        #manager.deleteProject(projectName)
