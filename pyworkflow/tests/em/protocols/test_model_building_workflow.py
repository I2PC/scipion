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

from pyworkflow.em.packages.chimera.protocol_fit import \
    ChimeraProtRigidFit
from pyworkflow.em.packages.powerfit.protocol_powerfit import \
    PowerfitProtRigidFit
from pyworkflow.em.protocol.protocol_import import ProtImportPdb, \
    ProtImportVolumes
from pyworkflow.tests import *



class TestImportBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsModBuild = DataSet.getDataSet('model_building_tutorial')


class TestImportData(TestImportBase):
    """ Import map volumes and atomic structures(PDBx/mmCIF files)
    """

    def _importVolume(self):
        args = {'filesPath': self.dsModBuild.getFile('volumes/1ake_4-5A.mrc'),
                'samplingRate': 1.5,
                'setDefaultOrigin': True
                }
        protImportVol = self.newProtocol(ProtImportVolumes, **args)
        protImportVol.setObjLabel('volume 1ake_4-5A\n')
        self.launchProtocol(protImportVol)
        volume = protImportVol.outputVolume
        self.assertTrue(volume.getFileName())
        return volume

    def _importStructurePDBWoVol(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile('PDBx_mmCIF/1ake_start.pdb'),
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('pdb 1ake_start')
        self.launchProtocol(protImportPDB)
        structure1_PDB = protImportPDB.outputPdb
        self.assertTrue(structure1_PDB.getFileName())
        return structure1_PDB

    def _importStructuremmCIFWoVol(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/1ake_start.pdb.cif'),
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('mmCIF 1ake_start')
        self.launchProtocol(protImportPDB)
        structure1_mmCIF = protImportPDB.outputPdb
        self.assertTrue(structure1_mmCIF.getFileName())
        return structure1_mmCIF

    def _importStructurePDBWithVol(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile('PDBx_mmCIF/1ake_start.pdb'),
                'inputVolume': self._importVolume()
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('pdb 1ake_start')
        self.launchProtocol(protImportPDB)
        structure2_PDB = protImportPDB.outputPdb
        self.assertTrue(structure2_PDB.getFileName())
        return structure2_PDB

    def _importStructuremmCIFWithVol(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile('PDBx_mmCIF/'
                                                   '1ake_start.pdb.cif'),
                'inputVolume': self._importVolume()
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('mmCIF 1ake_start')
        self.launchProtocol(protImportPDB)
        structure2_mmCIF = protImportPDB.outputPdb
        self.assertTrue(structure2_mmCIF.getFileName())
        return structure2_mmCIF

class TestPowerFit(TestImportData):
    """ Test the rigid fit of power fit protocol
    """

    def testPowerFitFromVolAndPDB(self):
        print "Run powerfit from imported volume and pdb file"
        # This test check that powerfit runs with a volume provided
        # directly as inputVol

        args = {'resolution': 2.,
                'angleStep': 10.,
                'nModels': 3.
                }
        protPower = self.newProtocol(PowerfitProtRigidFit, **args)
        volume = self._importVolume()
        protPower.inputVol.set(volume)
        structure1_PDB = self._importStructurePDBWoVol()
        protPower.inputPDB.set(structure1_PDB)
        self.launchProtocol(protPower)
        self.assertIsNotNone(protPower.outputPDBs.getFileName(),
                             "There was a problem with the alignment")

    def testPowerFitFromVolAndmmCIF(self):
        print "Run powerfit from imported volume and mmCIF file"
        # This test check that powerfit runs with a volume provided
        # directly as inputVol

        args = {'resolution': 2.,
                'angleStep': 10.,
                'nModels': 3.
                }
        protPower = self.newProtocol(PowerfitProtRigidFit, **args)
        volume = self._importVolume()
        protPower.inputVol.set(volume)
        structure1_mmCIF = self._importStructuremmCIFWoVol()
        protPower.inputPDB.set(structure1_mmCIF)
        self.launchProtocol(protPower)
        self.assertIsNotNone(protPower.outputPDBs.getFileName(),
                             "There was a problem with the alignment")

    def testPowerFitFromVolAssocToPDB(self):
        print "Run powerfit from imported pdb file and volume associated"
        # This test check that powerfit runs when a volume is provided
        # associated to the imputPDB and not directly as inputVol

        args = {'resolution': 2.,
                'angleStep': 10.,
                'nModels': 3.
                }
        protPower = self.newProtocol(PowerfitProtRigidFit, **args)
        structure2_PDB = self._importStructurePDBWithVol()
        protPower.inputPDB.set(structure2_PDB)
        self.launchProtocol(protPower)
        self.assertIsNotNone(protPower.outputPDBs.getFileName(),
                             "There was a problem with the alignment")

    def testPowerFitFromVolAssocTommCIF(self):
        print "Run powerfit from imported mmCIF file and volume associated"
        # This test check that powerfit runs when a volume is provided
        # associated to the imputPDB and not directly as inputVol

        args = {'resolution': 2.,
                'angleStep': 10.,
                'nModels': 3.
                }
        protPower = self.newProtocol(PowerfitProtRigidFit, **args)
        structure2_mmCIF = self._importStructuremmCIFWithVol()
        protPower.inputPDB.set(structure2_mmCIF)
        self.launchProtocol(protPower)
        self.assertIsNotNone(protPower.outputPDBs.getFileName(),
                             "There was a problem with the alignment")

    def testPowerFitFromPDBWithoutVol(self):
        print "Run powerfit from imported pdb file without imported or " \
              "pdb-associated volume"
        # This test corroborates that powerfit does not run unless a volume
        # is provided (directly as inputVol or associated to the imputPDB)

        structure1_PDB = self._importStructurePDBWoVol()
        self.assertTrue(structure1_PDB.getFileName())
        self.assertFalse(structure1_PDB.getVolume())

        args = {'resolution': 2.,
                'angleStep': 10.,
                'nModels': 3.
                }
        protPower = self.newProtocol(PowerfitProtRigidFit, **args)
        protPower.inputPDB.set(structure1_PDB)

        try:
            self.launchProtocol(protPower)
        except Exception as e:
            self.assertTrue(True)
            print "This test should return a error message as '" \
                  " ERROR running protocol scipion - rigid fit"

            return
        self.assertTrue(False)

    def testPowerFitFrommmCIFWithoutVol(self):
        print "Run powerfit from imported mmCIF file without imported or " \
              "mmCIF-associated volume"
        # This test corroborates that powerfit does not run unless a volume
        # is provided (directly as inputVol or associated to the imputPDB)

        structure1_mmCIF = self._importStructuremmCIFWoVol()
        self.assertTrue(structure1_mmCIF.getFileName())
        self.assertFalse(structure1_mmCIF.getVolume())

        args = {'resolution': 2.,
                'angleStep': 10.,
                'nModels': 3.
                }
        protPower = self.newProtocol(PowerfitProtRigidFit, **args)
        protPower.inputPDB.set(structure1_mmCIF)

        try:
            self.launchProtocol(protPower)
        except Exception as e:
            self.assertTrue(True)
            print "This test should return a error message as '" \
                  " ERROR running protocol scipion - rigid fit"

            return
        self.assertTrue(False)

# class TestChimeraFit(TestPowerFit):
#     """ Test the rigid fit of chimera fit protocol
#     """
#     def testChimeraFitFromVolAndPDBwoSavingVol(self):
#         print "Run Chimera fit from imported volume and pdb file"
#         # This test check that chimera runs with a volume provided
#         # directly as inputVol
#
#         volume = self._importVolume()
#         structure1_PDB = self._importStructurePDBWoVol()
#
#         args = {'inputVolume': volume,
#                 'pdbFileToBeRefined': structure1_PDB
#                 }
#         protChimera = self.newProtocol(ChimeraProtRigidFit, **args)
#         protChimera.inputVolume.set(volume)
#         protChimera.pdbFileToBeRefined.set(structure1_PDB)
#         self.launchProtocol(protChimera)
#         self.assertIsNotNone(protChimera.outputPdb.getFileName(),
#                              "There was a problem with the alignment")
#
#         import os
#         print "Current working directory: ", os.getcwd()
#         print "Listdir: ", os.listdir(os.getcwd())
#         print os.path.abspath("chimera_input.cmd")
#
#
#         tmpFileNameCMD = os.path.abspath( "chimera_input.cmd")/\
#                          "chimera_test.cmd"
#         print "tmpFileNameCMD: ", tmpFileNameCMD
#         f = open(tmpFileNameCMD, "w")
#         f.write("sleep 5")
#         f.write("move x -24.11")
#         f.write("move y -45.76")
#         f.write("move z -24.60")
#         f.close()



