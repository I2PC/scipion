# ***************************************************************************
# * Authors:    Marta Martinez (mmmtnez@cnb.csic.es)
# *             Roberto Marabini (roberto@cnb.csic.es)
# *
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
# ***************************************************************************/

# general protocol to test the model building workflow starting from a whole
# volume or a unit cell extracted from a map volume. Here we are going to test
# two protocols of rigid fitting (powerfit and chimera) and two protocols of
# flexible fitting (coot and refmac), as well as validation programs such as
# emringer and molprobity

from pyworkflow.tests import *
from pyworkflow.em.protocol.protocol_import import ProtImportPdb, \
    ProtImportVolumes
from pyworkflow.em.packages.powerfit.protocol_powerfit import \
    PowerfitProtRigidFit


class TestImportBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsModBuild = DataSet.getDataSet('model_building_tutorial')


class TestPowerFit(TestImportBase):
    """ Test the rigid fit of power fit protocol
    """
    def testPowerFitFromVolAndPDB(self):
        print "Run powerfit from imported volume and pdb"
        # This test check that prowerfit runs with a volume provided
        # directly as inputVol

        args = {'filesPath': self.dsModBuild.getFile('volumes/1ake_4-5A.mrc'),
                'samplingRate': 1.5,
                'setDefaultOrigin': True
                }
        protImportVol = self.newProtocol(ProtImportVolumes, **args)
        protImportVol.setObjLabel('volume 1ake_4-5A\ntest 1')
        self.launchProtocol(protImportVol)
        volume = protImportVol.outputVolume
        self.assertTrue(volume.getFileName())

        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile('PDBs/1ake_start.pdb'),
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('pdb 1ake_start')
        self.launchProtocol(protImportPDB)
        pdb = protImportPDB.outputPdb
        self.assertTrue(pdb.getFileName())

        args = {'resolution': 2.,
                'angleStep': 10.,
                'nModels': 3.
                }
        protPower = self.newProtocol(PowerfitProtRigidFit, **args)
        protPower.inputVol.set(volume)
        protPower.inputPDB.set(pdb)
        self.launchProtocol(protPower)
        self.assertIsNotNone(protPower.outputPDBs.getFileName(),
                             "There was a problem with the alignment")

    def testPowerFitFromVolAssocToPDB(self):
        print "Run powerfit from imported pdb and volume associated"
        # This test check that prowerfit runs when a volume is provided
        # associated to the imputPDB and not directly as inputVol

        args = {'filesPath': self.dsModBuild.getFile('volumes/1ake_4-5A.mrc'),
                'samplingRate': 1.5,
                'setDefaultOrigin': True
                }
        protImportVol = self.newProtocol(ProtImportVolumes, **args)
        protImportVol.setObjLabel('volume 1ake_4-5A\ntest 2')
        self.launchProtocol(protImportVol)
        volume = protImportVol.outputVolume
        self.assertTrue(volume.getFileName())

        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile('PDBs/1ake_start.pdb'),
                'inputVolume': volume
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('pdb 1ake_start')
        self.launchProtocol(protImportPDB)
        pdb = protImportPDB.outputPdb
        self.assertTrue(pdb.getFileName())

        args = {'resolution': 2.,
                'angleStep': 10.,
                'nModels': 3.
                }
        protPower = self.newProtocol(PowerfitProtRigidFit, **args)
        protPower.inputPDB.set(pdb)
        self.launchProtocol(protPower)
        self.assertIsNotNone(protPower.outputPDBs.getFileName(),
                             "There was a problem with the alignment")

    def testPowerFitFromPDBWithoutVol(self):
        print "Run powerfit from imported pdb without imported or " \
              "pdb-associated volume"
        # This test corroborates that prowerfit does not run unless a volume
        # is provided (directly as inputVol or associated to the imputPDB)

        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile('PDBs/1ake_start.pdb')
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('pdb 1ake_start\ntest 3\n')
        self.launchProtocol(protImportPDB)
        pdb = protImportPDB.outputPdb
        self.assertTrue(pdb.getFileName())
        self.assertFalse(pdb.getVolume())

        args = {'resolution': 2.,
                'angleStep': 10.,
                'nModels': 3.
                }
        protPower = self.newProtocol(PowerfitProtRigidFit, **args)
        protPower.inputPDB.set(pdb)
        self.launchProtocol(protPower)
