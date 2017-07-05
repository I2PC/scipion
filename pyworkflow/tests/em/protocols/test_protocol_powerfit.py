# **************************************************************************
# *
# * Authors:    Carlos Oscar Sorzano (coss@cnb.csic.es)
# *             
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

from pyworkflow.tests import *
from pyworkflow.em.protocol.protocol_import import ProtImportPdb, ProtImportVolumes
from pyworkflow.em.packages.powerfit.protocol_powerfit import PowerfitProtRigidFit

class TestPowerFit(BaseTest):     
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('powerfit')
        cls.pdb = cls.dataset.getFile('pdb')
        cls.vol = cls.dataset.getFile('vol')
   
    def testPowerFit(self):
        print "Run Powerfit"
        protImportPDB = self.newProtocol(ProtImportPdb, 
                                      inputPdbData=ProtImportPdb.IMPORT_FROM_FILES, 
                                      pdbFile=self.pdb)
        self.launchProtocol(protImportPDB)
        self.assertIsNotNone(protImportPDB.outputPdb.getFileName(), 
                             "There was a problem with the import") 
               
        protImportVol = self.newProtocol(ProtImportVolumes,
                                      filesPath=self.vol, 
                                      samplingRate=4)
        self.launchProtocol(protImportVol)
        self.assertIsNotNone(protImportVol.outputVolume.getFileName(), 
                             "There was a problem with the import") 

        protPower = self.newProtocol(PowerfitProtRigidFit, angleStep=20)                                       
        protPower.inputVol.set(protImportVol.outputVolume)
        protPower.inputPDB.set(protImportPDB.outputPdb)
        self.launchProtocol(protPower)
        self.assertIsNotNone(protPower.outputPDBs.getFileName(), 
                             "There was a problem with the alignment") 
         
