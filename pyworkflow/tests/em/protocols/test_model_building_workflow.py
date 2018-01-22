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

from pyworkflow.em.packages.chimera.protocol_fit import ChimeraProtRigidFit 
from pyworkflow.em.packages.powerfit.protocol_powerfit import PowerfitProtRigidFit 
from pyworkflow.em.protocol.protocol_import import ProtImportPdb
from pyworkflow.em.protocol.protocol_import import ProtImportVolumes
from pyworkflow.tests import *
import os.path


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
        return volume

    def _importStructurePDBWoVol(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/1ake_start.pdb'),
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('pdb 1ake_start')
        self.launchProtocol(protImportPDB)
        structure1_PDB = protImportPDB.outputPdb
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
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/1ake_start.pdb'),
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

    def _importMut1StructurePDBWoVol(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/1ake_mut1.pdb'),
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('pdb 1ake_mut1')
        self.launchProtocol(protImportPDB)
        structure3_PDB = protImportPDB.outputPdb
        return structure3_PDB

    def _importMut2StructurePDBWoVol(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/1ake_mut2.pdb'),
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('pdb 1ake_mut2')
        self.launchProtocol(protImportPDB)
        structure4_PDB = protImportPDB.outputPdb
        return structure4_PDB


class TestPowerFit(TestImportData):
    """ Test the rigid fit of power fit protocol
    """

    def testPowerFitFromVolAndPDB(self):
        print "Run powerfit from imported volume and pdb file"
        # This test checks that powerfit runs with a volume provided
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
        # This test checks that powerfit runs with a volume provided
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
        # This test checks that powerfit runs when a volume is provided
        # associated to the imput PDB and not directly as inputVol

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
        # This test checks that powerfit runs when a volume is provided
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


class TestChimeraFit(TestImportData):
    """ Test the rigid fit of chimera fit protocol
    """

    def testChimeraFitFromVolAndPDBWithSavingVol(self):
        """ This test checks that chimera runs with a volume provided
        directly as inputVol, input PDB """
        print "Run Chimera fit from imported volume and pdb file"

        # Import Volume
        volume = self._importVolume()

        # import PDB
        structure1_PDB = self._importStructurePDBWoVol()

        # create auxiliary CMD file for chimera fit
        extraCommands = ""
        extraCommands += "runCommand('move -24.11,-45.76,-24.60 model #2 " \
                         "coord #1')\n"
        extraCommands += "runCommand('fitmap #2 #1')\n"
        extraCommands += "runCommand('scipionwrite model #2 refmodel #1 " \
                         "saverefmodel 1')\n"
        extraCommands += "runCommand('stop')\n"

        args = {'extraCommands': extraCommands,
                'inputVolume': volume,
                'pdbFileToBeRefined': structure1_PDB
                }
        protChimera = self.newProtocol(ChimeraProtRigidFit,
                                       **args)
        self.launchProtocol(protChimera)
        self.assertIsNotNone(protChimera.outputPdb.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protChimera.outputPdb.getFileName()))
        self.assertTrue(os.path.exists(protChimera.outputPdb.
                                       getVolume().getFileName()))

    def testChimeraFitFromVolAndmmCIFWithSavingVol(self):
        """ This test checks that chimera runs with a volume provided
        directly as inputVol, input CIF file """
        print "Run Chimera fit from imported volume and cif file"

        volume = self._importVolume()
        structure1_mmCIF = self._importStructuremmCIFWoVol()
        extraCommands = ""
        extraCommands += "runCommand('move -24.11,-45.76,-24.60 model #2 " \
                         "coord #1')\n"
        extraCommands += "runCommand('fitmap #2 #1')\n"
        extraCommands += "runCommand('scipionwrite model #2 refmodel #1 " \
                         "saverefmodel 1')\n"
        extraCommands += "runCommand('stop')\n"

        args = {'extraCommands': extraCommands,
                'inputVolume': volume,
                'pdbFileToBeRefined': structure1_mmCIF
                }
        protChimera = self.newProtocol(ChimeraProtRigidFit, **args)
        self.launchProtocol(protChimera)
        self.assertIsNotNone(protChimera.outputPdb.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protChimera.outputPdb.
                                       getVolume().getFileName()))

    def testChimeraFitFromVolAssocToPDBWithSavingVol(self):
        # This test checks that chimera runs when a volume is provided
        # associated to the input PDB and not directly as inputVol
        print "Run Chimera fit from imported pdb file and volume associated"

        structure2_PDB = self._importStructurePDBWithVol()
        extraCommands = ""
        extraCommands += "runCommand('move -24.11,-45.76,-24.60 model #2 " \
                         "coord #1')\n"
        extraCommands += "runCommand('fitmap #2 #1')\n"
        extraCommands += "runCommand('scipionwrite model #2 refmodel #1 " \
                         "saverefmodel 1')\n"
        extraCommands += "runCommand('stop')\n"

        args = {'extraCommands': extraCommands,
                'pdbFileToBeRefined': structure2_PDB
                }
        protChimera = self.newProtocol(ChimeraProtRigidFit, **args)
        self.launchProtocol(protChimera)
        self.assertIsNotNone(protChimera.outputPdb.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protChimera.outputPdb.
                                       getVolume().getFileName()))

    def testChimeraFitFromVolAssocTommCIFWithSavingVol(self):
        # This test checks that chimera runs when a volume is provided
        # associated to the imput mmCIF file and not directly as inputVol
        print "Run Chimera fit from imported mmCIF file and volume associated"

        structure2_mmCIF = self._importStructuremmCIFWithVol()
        extraCommands = ""
        extraCommands += "runCommand('move -24.11,-45.76,-24.60 model #2 " \
                         "coord #1')\n"
        extraCommands += "runCommand('fitmap #2 #1')\n"
        extraCommands += "runCommand('scipionwrite model #2 refmodel #1 " \
                         "saverefmodel 1')\n"
        extraCommands += "runCommand('stop')\n"

        args = {'extraCommands': extraCommands,
                'pdbFileToBeRefined': structure2_mmCIF
                }
        protChimera = self.newProtocol(ChimeraProtRigidFit, **args)
        self.launchProtocol(protChimera)
        self.assertIsNotNone(protChimera.outputPdb.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protChimera.outputPdb.
                                       getVolume().getFileName()))

    def testChimeraFitFromPDBWithoutVol(self):
        # This test corroborates that chimera does not run unless a volume
        # is provided (directly as inputVol or associated to the imputPDB)
        # protocol should raise an exception
        print "Run Chimera from imported pdb file without imported or " \
              "pdb-associated volume"

        structure1_PDB = self._importStructurePDBWoVol()
        self.assertTrue(structure1_PDB.getFileName())
        self.assertFalse(structure1_PDB.getVolume())

        protChimera = self.newProtocol(ChimeraProtRigidFit)
        protChimera.pdbFileToBeRefined.set(structure1_PDB)

        try:
            self.launchProtocol(protChimera)
        except Exception as e:
            self.assertTrue(True)
            print "This test should return a error message as '" \
                  " ERROR running protocol scipion - chimera rigid fit"

            return
        self.assertTrue(False)

    def testChimeraFitFrommmCIFWithoutVol(self):
        # This test corroborates that chimera does not run unless a volume
        # is provided (directly as inputVol or associated to the imputPDB)
        # protocol should raise an exception
        print "Run chimera from imported mmCIF file without imported or " \
              "mmCIF-associated volume"

        structure1_mmCIF = self._importStructuremmCIFWoVol()
        self.assertTrue(structure1_mmCIF.getFileName())
        self.assertFalse(structure1_mmCIF.getVolume())
        protChimera = self.newProtocol(ChimeraProtRigidFit)
        protChimera.pdbFileToBeRefined.set(structure1_mmCIF)

        try:
            self.launchProtocol(protChimera)
        except Exception as e:
            self.assertTrue(True)
            print "This test should return a error message as '" \
                  " ERROR running protocol scipion - chimera rigid fit"

            return
        self.assertTrue(False)

    def testChimeraFitFromPDBFromPowerFitWithSavingVol(self):
        # This test checks that chimera runs when a volume is provided
        # associated to the input PDB file and not directly as inputVol
        # PDB file created by powerfit as a set of PDBs
        print "Run Chimera fit from PDB file saved from PowerFit"

        # Powerfit
        args = {'resolution': 2.,
                'angleStep': 45.,
                'nModels': 3.
                }
        protPower = self.newProtocol(PowerfitProtRigidFit, **args)
        volume = self._importVolume()
        protPower.inputVol.set(volume)
        structure1_PDB = self._importStructurePDBWoVol()
        protPower.inputPDB.set(structure1_PDB)
        self.launchProtocol(protPower)

        # chimera fit
        structure2_PDB = protPower.outputPDBs.getFirstItem()

        # trick to save a single PDB from the set of PDBs
        protPower._defineOutputs(outputPDB2s=structure2_PDB)

        extraCommands = ""
        extraCommands += "runCommand('fitmap #2 #1')\n"
        extraCommands += "runCommand('scipionwrite model #2 refmodel #1 " \
                         "saverefmodel 1')\n"
        extraCommands += "runCommand('stop')\n"

        args = {'extraCommands': extraCommands,
                'pdbFileToBeRefined': structure2_PDB
                }

        protChimera = self.newProtocol(ChimeraProtRigidFit, **args)
        self.launchProtocol(protChimera)
        self.assertIsNotNone(protChimera.outputPdb.getFileName(),
                             "There was a problem with the alignment")

    def testChimeraFitFromVolAssocToPDBPlusPDBsWithSavingVol(self):
        # This test checks that chimera runs when a volume is provided
        # associated to the input PDB and several PDB files are added

        print "Run Chimera fit from imported pdb file and volume associated " \
              "and addition of two other pdb files"

        structure2_PDB = self._importStructurePDBWithVol()
        structure3_PDB = self._importMut1StructurePDBWoVol()
        structure4_PDB = self._importMut2StructurePDBWoVol()

        _pdbFiles = []
        _pdbFiles.append(structure3_PDB)
        _pdbFiles.append(structure4_PDB)

        extraCommands = ""
        extraCommands += "runCommand('move -24.11,-45.76,-24.60 model #2 " \
                         "coord #1')\n"
        extraCommands += "runCommand('fitmap #2 #1')\n"
        extraCommands += "runCommand('scipionwrite model #2 refmodel #1 " \
                         "saverefmodel 1')\n"
        extraCommands += "runCommand('stop')\n"

        args = {'extraCommands': extraCommands,
                'pdbFileToBeRefined': structure2_PDB,
                'inputPdbFiles': _pdbFiles
                }
        protChimera = self.newProtocol(ChimeraProtRigidFit, **args)
        self.launchProtocol(protChimera)
        self.assertIsNotNone(protChimera.outputPdb.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protChimera.outputPdb.
                                       getVolume().getFileName()))

    def testChimeraFitFromChimeraPDB(self):
        # This test checks that chimera runs with objects not imported
        # but generated in other programs
        print "Run Chimera fit using the pdb and its volume associated " \
              "generated in a previous protocol of Chimera rigid fit"

        volume = self._importVolume()
        structure1_PDB = self._importStructurePDBWoVol()
        extraCommands = ""
        extraCommands += "runCommand('move -24.11,-45.76,-24.60 model #2 " \
                         "coord #1')\n"
        extraCommands += "runCommand('fitmap #2 #1')\n"
        extraCommands += "runCommand('scipionwrite model #2 refmodel #1 " \
                         "saverefmodel 1')\n"
        extraCommands += "runCommand('stop')\n"

        args = {'extraCommands': extraCommands,
                'inputVolume': volume,
                'pdbFileToBeRefined': structure1_PDB
                }
        protChimera = self.newProtocol(ChimeraProtRigidFit, **args)
        self.launchProtocol(protChimera)
        structure2_PDB = protChimera.outputPdb

        extraCommands = ""
        extraCommands += "runCommand('move 24.11,45.76,24.60 model #2 " \
                         "coord #1')\n"
        extraCommands += "runCommand('scipionwrite model #2 refmodel #1 " \
                         "saverefmodel 1')\n"
        extraCommands += "runCommand('stop')\n"
        args = {'extraCommands': extraCommands,
                'pdbFileToBeRefined': structure2_PDB,
                }
        protChimera = self.newProtocol(ChimeraProtRigidFit, **args)
        self.launchProtocol(protChimera)
        structure3_PDB = protChimera.outputPdb
        extraCommands = ""
        extraCommands += "runCommand('move -24.11,-45.76,-24.60 model #2 " \
                         "coord #1')\n"
        extraCommands += "runCommand('fitmap #2 #1')\n"
        extraCommands += "runCommand('scipionwrite model #2 refmodel #1 " \
                         "saverefmodel 1')\n"
        extraCommands += "runCommand('stop')\n"
        args = {'extraCommands': extraCommands,
                'inputVolume': volume,
                'pdbFileToBeRefined': structure3_PDB,
                }
        protChimera = self.newProtocol(ChimeraProtRigidFit, **args)
        self.launchProtocol(protChimera)
        self.assertIsNotNone(protChimera.outputPdb.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protChimera.outputPdb.
                                       getVolume().getFileName()))

    def testChimeraFitromVolAssocToPDBPlusPDBsFromPowerfit(self):
        # This test checks that chimera runs when a volume is provided
        # associated to the input PDB and several PDB files are added,
        # one of them generated from powerfit rigid fit protocol
        print "Run Chimera fit from imported pdb file and volume associated " \
              "and addition of other three pdb files (one of them generated " \
              "by Powerfit)"

        structure2_PDB = self._importStructurePDBWithVol()
        structure3_PDB = self._importMut1StructurePDBWoVol()
        structure4_PDB = self._importMut2StructurePDBWoVol()

        # Powerfit
        args = {'resolution': 2.,
                'angleStep': 45.,
                'nModels': 3.
                }
        protPower = self.newProtocol(PowerfitProtRigidFit, **args)
        volume = self._importVolume()
        protPower.inputVol.set(volume)
        structure1_PDB = self._importStructurePDBWoVol()
        protPower.inputPDB.set(structure1_PDB)
        self.launchProtocol(protPower)

        # chimera fit
        structure5_PDB = protPower.outputPDBs.getFirstItem()
        protPower._defineOutputs(outputPDB2s=structure5_PDB)

        _pdbFiles = []
        _pdbFiles.append(structure3_PDB)
        _pdbFiles.append(structure4_PDB)
        _pdbFiles.append(structure5_PDB)

        extraCommands = ""
        extraCommands += "runCommand('move -24.11,-45.76,-24.60 model #2 " \
                         "coord #1')\n"
        extraCommands += "runCommand('fitmap #2 #1')\n"
        extraCommands += "runCommand('scipionwrite model #2 refmodel #1 " \
                         "saverefmodel 1')\n"
        extraCommands += "runCommand('stop')\n"

        args = {'extraCommands': extraCommands,
                'pdbFileToBeRefined': structure2_PDB,
                'inputPdbFiles': _pdbFiles
                }
        protChimera = self.newProtocol(ChimeraProtRigidFit, **args)
        self.launchProtocol(protChimera)
        self.assertIsNotNone(protChimera.outputPdb.getFileName(),
                             "There was a problem with the alignment")
