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


# protocol to test the chimera operate protocol starting from a
# volume, a pdb or both. Here we are going to test the suitability of chimera
# operate to save pdbs and, optionally, volumes after carrying out different
# manipulations with chimera

from pyworkflow.em.packages.chimera.protocol_operate import ChimeraProtOperate
from pyworkflow.em.protocol.protocol_import import ProtImportPdb, \
    ProtImportVolumes
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


class TestChimeraOperate(TestImportData):
    """ Test the chimera operate protocol
    """

    def testChimeraOperateFromVolAndPDB(self):
        """ This test checks that chimera runs with a volume provided
        directly as inputVol, input PDB """
        print "Run Chimera operate from imported volume and pdb file"

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
        protChimera = self.newProtocol(ChimeraProtOperate,
                                       **args)
        self.launchProtocol(protChimera)
        self.assertIsNotNone(protChimera.outputPdb_01.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protChimera.outputPdb_01.getFileName()))
        self.assertTrue(os.path.exists(protChimera.outputPdb_01.
                                       getVolume().getFileName()))

    def testChimeraOperateFromVolAndmmCIF(self):
        """ This test checks that chimera runs with a volume provided
        directly as inputVol, input CIF file """
        print "Run Chimera operate from imported volume and cif file"

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
        protChimera = self.newProtocol(ChimeraProtOperate, **args)
        self.launchProtocol(protChimera)
        self.assertIsNotNone(protChimera.outputPdb_01.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protChimera.outputPdb_01.
                                       getVolume().getFileName()))

    def testChimeraOperateFromVolAssocToPDB(self):
        # This test checks that chimera runs when a volume is provided
        # associated to the input PDB and not directly as inputVol
        print "Run Chimera operate from imported pdb file and volume associated"

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
        protChimera = self.newProtocol(ChimeraProtOperate, **args)
        self.launchProtocol(protChimera)
        self.assertIsNotNone(protChimera.outputPdb_01.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protChimera.outputPdb_01.
                                       getVolume().getFileName()))

    def testChimeraOperateFromVolAssocTommCIF(self):
        # This test checks that chimera runs when a volume is provided
        # associated to the imput mmCIF file and not directly as inputVol
        print "Run Chimera operate from imported mmCIF file and volume " \
              "associated"

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
        protChimera = self.newProtocol(ChimeraProtOperate, **args)
        self.launchProtocol(protChimera)
        self.assertIsNotNone(protChimera.outputPdb_01.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protChimera.outputPdb_01.
                                       getVolume().getFileName()))

    def testChimeraOperateFromPDBWithoutVol(self):
        # This test corroborates that chimera runs with a pdb and without
        # providing a volume

        print "Run Chimera operate from imported pdb file without imported or " \
              "pdb-associated volume"

        structure1_PDB = self._importStructurePDBWoVol()
        self.assertTrue(structure1_PDB.getFileName())
        self.assertFalse(structure1_PDB.getVolume())

        extraCommands = ""
        extraCommands += "runCommand('move -24.11,-45.76,-24.60 model #1 " \
                         "coord #1')\n"
        extraCommands += "runCommand('scipionwrite model #1 refmodel #0')\n"
        extraCommands += "runCommand('stop')\n"

        args = {'extraCommands': extraCommands,
                'pdbFileToBeRefined': structure1_PDB
                }
        protChimera = self.newProtocol(ChimeraProtOperate, **args)
        self.launchProtocol(protChimera)
        self.assertIsNotNone(protChimera.outputPdb_01.getFileName(),
                             "There was a problem with the alignment")

    def testChimeraOperateFrommmCIFWithoutVol(self):
        # This test corroborates that chimera runs with a mmCIF file and
        # without providing a volume
        print "Run chimera operate from imported mmCIF file without imported " \
              "or mmCIF-associated volume"

        structure1_mmCIF = self._importStructuremmCIFWoVol()
        self.assertTrue(structure1_mmCIF.getFileName())
        self.assertFalse(structure1_mmCIF.getVolume())

        extraCommands = ""
        extraCommands += "runCommand('move -24.11,-45.76,-24.60 model #1 " \
                         "coord #1')\n"
        extraCommands += "runCommand('scipionwrite model #1 refmodel #0')\n"
        extraCommands += "runCommand('stop')\n"

        args = {'extraCommands': extraCommands,
                'pdbFileToBeRefined': structure1_mmCIF
                }
        protChimera = self.newProtocol(ChimeraProtOperate, **args)
        self.launchProtocol(protChimera)
        self.assertIsNotNone(protChimera.outputPdb_01.getFileName(),
                             "There was a problem with the alignment")


    def testChimeraOperateFromChimeraPDB(self):
        # This test checks that chimera runs with objects not imported
        # but generated in other chimera programs
        print "Run Chimera operate using the pdb and its volume associated " \
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
        protChimera = self.newProtocol(ChimeraProtOperate, **args)
        self.launchProtocol(protChimera)
        structure2_PDB = protChimera.outputPdb_01

        extraCommands = ""
        extraCommands += "runCommand('move 24.11,45.76,24.60 model #2 " \
                         "coord #1')\n"
        extraCommands += "runCommand('scipionwrite model #2 refmodel #1 " \
                         "saverefmodel 1')\n"
        extraCommands += "runCommand('stop')\n"
        args = {'extraCommands': extraCommands,
                'pdbFileToBeRefined': structure2_PDB,
                }
        protChimera = self.newProtocol(ChimeraProtOperate, **args)
        self.launchProtocol(protChimera)
        structure3_PDB = protChimera.outputPdb_01
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
        protChimera = self.newProtocol(ChimeraProtOperate, **args)
        self.launchProtocol(protChimera)
        self.assertIsNotNone(protChimera.outputPdb_01.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protChimera.outputPdb_01.
                                       getVolume().getFileName()))

    def testChimeraOperateFromVolAssocToPDBPlusPDBs(self):
        # This test checks that chimera runs when a volume is provided
        # associated to the input PDB and several PDB files are added,
        # one of them generated from powerfit rigid fit protocol
        print "Run Chimera fit from imported pdb file and volume associated " \
              "and addition of other three pdb files"

        structure2_PDB = self._importStructurePDBWithVol()
        structure3_PDB = self._importMut1StructurePDBWoVol()
        structure4_PDB = self._importMut2StructurePDBWoVol()

        # chimera operate
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
        protChimera = self.newProtocol(ChimeraProtOperate, **args)
        self.launchProtocol(protChimera)
        self.assertIsNotNone(protChimera.outputPdb_01.getFileName(),
                             "There was a problem with the alignment")
