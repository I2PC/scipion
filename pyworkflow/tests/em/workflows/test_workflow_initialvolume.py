# **************************************************************************
# *
# * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
# *  e-mail address 'jgomez@cnb.csic.es'
# *
# **************************************************************************

"""
Test the execution of the initial volume service with different datasets.
Execute:
1) Import averages
2) Run initial volume algorithms (ransac, eman, significant, prime)
3) Run align volumes

Datasets:
1) Groel
2) Ribosome
3) BPV
"""

import os

import pyworkflow.tests as tests
import pyworkflow.em as em
import pyworkflow.em.packages.xmipp3 as xmipp3
import pyworkflow.em.packages.eman2 as eman2
  
   
class TestGroel(tests.BaseTest):
    """ Check the images are converted properly to spider format. """
    
    @classmethod
    def setUpClass(cls):    
        # Create a new project
        tests.setupTestProject(cls)
        cls.ds = tests.DataSet.getDataSet('initial_volume')
    
    def test_groel(self):
        """ Run an Import particles protocol. """
        cpus = os.environ.get('SCIPION_TEST_CPU', 4)
        # 1. Run import of averages
        groelAvg = self.ds.getFile('groel')
        sym = 'd7'
        protImport = self.newProtocol(em.ProtImportAverages, 
                                      filesPath=groelAvg, 
                                      samplingRate=1)
        self.launchProtocol(protImport)
        
        # 2. Run initial models
        # 2a. Ransac 
        protRansac = self.newProtocol(xmipp3.XmippProtRansac,
                                      objLabel='xmipp - ransac',
                                      symmetryGroup=sym,
                                      nRansac=200,
                                      maxFreq=20,
                                      numberOfThreads=cpus
                                      )
        protRansac.inputSet.set(protImport.outputAverages)
        self.launchProtocol(protRansac)
        
        # 2b. Eman 
        protEmanInitVol = self.newProtocol(eman2.EmanProtInitModel,
                                           objLabel='eman - initial vol',
                                           numberOfThreads=cpus)
        protEmanInitVol.inputSet.set(protImport.outputAverages)
        self.launchProtocol(protEmanInitVol)
      
        
        
        
