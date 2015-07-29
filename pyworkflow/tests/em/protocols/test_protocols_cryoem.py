# **************************************************************************
# *
# * Authors:    Laura del Cano (ldelcano@cnb.csic.es)
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

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.em.protocol import ProtImportAverages

from pyworkflow.em.packages.cryoem_scipion import ProtCryoem


class TestCryoEmBase(BaseTest):
    def test_initialmodel(self):
        #Import a set of averages
        print "Import Set of averages"
        protImportAvg = self.newProtocol(ProtImportAverages, 
                                         filesPath=self.averages, 
                                         checkStack=True, 
                                         samplingRate=2.1)
        self.launchProtocol(protImportAvg)
        self.assertIsNotNone(protImportAvg.getFiles(), "There was a problem with the import")

        return 
        print "Run CryoEM"
        protIniModel = self.newProtocol(ProtCryoem,
                                      symmetry=self.symmetry, 
                                      numberOfIterations=self.numberOfIterations, 
                                      numberOfModels=self.numberOfModels, 
                                      numberOfThreads=4)
        protIniModel.inputSet.set(protImportAvg.outputAverages)
        self.launchProtocol(protIniModel)
        self.assertIsNotNone(protIniModel.outputVolumes, "There was a problem with eman initial model protocol")


class TestCryoEmMDA(TestCryoEmBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('mda')
        cls.averages = cls.dataset.getFile('averages')
        cls.samplingRate = 3.5
        cls.symmetry = 'd6'
        cls.numberOfIterations = 5
        cls.numberOfModels = 2


class TestCryoEmGroel(TestCryoEmBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('groel')
        cls.averages = cls.dataset.getFile('averages')
        cls.samplingRate = 2.1
        cls.symmetry = 'd7'
        cls.numberOfIterations = 10
        cls.numberOfModels = 10


