# **************************************************************************
# *
# * Authors:    Erney Ramírez Aportela (eramirez@cnb.csic.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pyworkflow.em.packages.phenix.protocol_automated_sharpening import PhenixProtAutomatedSharpening
from pyworkflow.em.protocol import ProtImportVolumes, ProtImportPdb



class TestPhenixAutSharpBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='resmap'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.map = cls.dataset.getFile('betagal')
        cls.half1 = cls.dataset.getFile('betagal_half1')
        cls.half2 = cls.dataset.getFile('betagal_half2')

    @classmethod
    def runImportVolumes(cls, pattern, samplingRate):
        """ Run an Import volumes protocol. """
        cls.protImport = cls.newProtocol(ProtImportVolumes,
                                         filesPath=pattern,
                                         samplingRate=samplingRate
                                         )
        cls.launchProtocol(cls.protImport)
        return cls.protImport


class TestPhenixAutSharp(TestPhenixAutSharpBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestPhenixAutSharpBase.setData()
        cls.protImportMap = cls.runImportVolumes(cls.map, 3.54)        
        cls.protImportHalf1 = cls.runImportVolumes(cls.half1, 3.54)
        cls.protImportHalf2 = cls.runImportVolumes(cls.half2, 3.54)

    def testMapSharp(self):
        #print "Run autosharpen for input map"        
        autosharpen1 = self.newProtocol(PhenixProtAutomatedSharpening,
                                   inputMap=self.protImportMap.outputVolume,
                                   resolution=10)
        self.launchProtocol(autosharpen1)
        self.assertIsNotNone(autosharpen1.outputVol,
                        "autosharpen1 has failed")
         
    def testMapSharpHalfMap(self):
        #print "Run autosharpen for input map with half map"        
        autosharpen2 = self.newProtocol(PhenixProtAutomatedSharpening,
                                   inputMap=self.protImportMap.outputVolume,
                                   useSplitVolumes=True,
                                   volumeHalf1=self.protImportHalf1.outputVolume,
                                   volumeHalf2=self.protImportHalf2.outputVolume,
                                   resolution=10,
                                   sharpening=4) 
        self.launchProtocol(autosharpen2)
        self.assertIsNotNone(autosharpen2.outputVol,
                        "autosharpen2 has failed")    
   
    def testMapSharpPDB(self):
        #print "Run autosharpen for input map with atomic model'pdb'"        
        autosharpen3 = self.newProtocol(PhenixProtAutomatedSharpening,
                                   inputMap=self.protImportMap.outputVolume,
                                   usePDB=True,
                                   pdbId='5a1a',
                                   resolution=10,
                                   sharpening=5)
        self.launchProtocol(autosharpen3)
        self.assertIsNotNone(autosharpen3.outputVol,
                        "autosharpen3 has failed")          
             
