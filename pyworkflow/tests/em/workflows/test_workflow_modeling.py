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
# a protocol of rigid fitting (chimera) and two protocols of
# flexible fitting (coot and refmac), as well as validation programs such as
# emringer and molprobity

from pyworkflow.em.packages.chimera.protocol_fit import ChimeraProtRigidFit
from pyworkflow.em.protocol.protocol_import import ProtImportPdb, \
    ProtImportVolumes
from pyworkflow.em.packages.ccp4.protocol_coot import CootRefine
from pyworkflow.em.packages.ccp4.protocol_refmac import CCP4ProtRunRefmac
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
                'setOrigCoord': False
                }
        protImportVol = self.newProtocol(ProtImportVolumes, **args)
        protImportVol.setObjLabel('import volume 1ake_4-5A\n with default '
                                  'origin\n')
        self.launchProtocol(protImportVol)
        volume = protImportVol.outputVolume
        return volume

    def _importVolume2(self):
        args = {'filesPath': self.dsModBuild.getFile('volumes/1ake_4-5A.mrc'),
                'samplingRate': 1.5,
                'setOrigCoord': True,
                'x': 11.994,
                'y': -7.881,
                'z': 10.91
                }
        protImportVol = self.newProtocol(ProtImportVolumes, **args)
        protImportVol.setObjLabel('import volume 1ake_4-5A\n set origin in 11 '
                                  '-7 10\n')
        self.launchProtocol(protImportVol)
        volume2 = protImportVol.outputVolume
        return volume2

    def _importStructurePDBWoVol(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/1ake_start.pdb'),
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import pdb\n 1ake_start')
        self.launchProtocol(protImportPDB)
        structure1_PDB = protImportPDB.outputPdb
        return structure1_PDB

    def _importStructuremmCIFWoVol(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/1ake_start.pdb.cif'),
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import mmCIF\n 1ake_start')
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
        protImportPDB.setObjLabel('import pdb\n volume associated\n 1ake_start')
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
        protImportPDB.setObjLabel('import mmCIF\n volume associated\n '
                                  '1ake_start')
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
        protImportPDB.setObjLabel('import pdb\n 1ake_mut1')
        self.launchProtocol(protImportPDB)
        structure3_PDB = protImportPDB.outputPdb
        return structure3_PDB

    def _importMut2StructurePDBWoVol(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/1ake_mut2.pdb'),
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import pdb\n 1ake_mut2')
        self.launchProtocol(protImportPDB)
        structure4_PDB = protImportPDB.outputPdb
        return structure4_PDB

    def _importCootStructureWoVol(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                        'PDBx_mmCIF/scipionOut0001.pdb')
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import pdb\n coot')
        self.launchProtocol(protImportPDB)
        structureCoot_PDB = protImportPDB.outputPdb
        self.assertTrue(structureCoot_PDB.getFileName())
        return structureCoot_PDB

    def _createExtraCommandLine(self, x, y, z):
        if (x != 0. or y != 0. or z != 0.):
            return """translate_molecule_by(0, %f, %f, %f)
fit_molecule_to_map_by_random_jiggle(0,7000,2)
scipion_write()
coot_real_exit(0)
""" % (x, y, z)
        else:
            return """scipion_write()
coot_real_exit(0)
"""


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
        protChimera.setObjLabel('chimera fit\n volume and pdb\n save volume '
                                'and model')
        self.launchProtocol(protChimera)
        self.assertIsNotNone(protChimera.outputPdb_01.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protChimera.outputPdb_01.getFileName()))
        self.assertTrue(os.path.exists(protChimera.outputPdb_01.
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
        protChimera.setObjLabel('chimera fit\n volume and mmCIF\n save volume '
                                'and model')
        self.launchProtocol(protChimera)
        self.assertIsNotNone(protChimera.outputPdb_01.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protChimera.outputPdb_01.
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
        protChimera.setObjLabel('chimera fit\n pdb and associated volume\n '
                                'save volume and model')
        self.launchProtocol(protChimera)
        self.assertIsNotNone(protChimera.outputPdb_01.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protChimera.outputPdb_01.
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
        protChimera.setObjLabel('chimera fit\n mmCIF and associated volume\n '
                                'save volume and model')
        self.launchProtocol(protChimera)
        self.assertIsNotNone(protChimera.outputPdb_01.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protChimera.outputPdb_01.
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
        protChimera.setObjLabel('chimera fit\n no volume\n associated to pdb')

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
        protChimera.setObjLabel('chimera fit\n no volume\n associated to mmCIF')

        try:
            self.launchProtocol(protChimera)
        except Exception as e:
            self.assertTrue(True)
            print "This test should return a error message as '" \
                  " ERROR running protocol scipion - chimera rigid fit"

            return
        self.assertTrue(False)

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
        protChimera.setObjLabel('chimera fit\n pdb and associated volume\n '
                                'plus other pdbs\n save volume and model')
        self.launchProtocol(protChimera)
        self.assertIsNotNone(protChimera.outputPdb_01.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protChimera.outputPdb_01.
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
        protChimera.setObjLabel("chimera fit \npdb and volume\n save volume "
                                "and model")
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
        protChimera = self.newProtocol(ChimeraProtRigidFit, **args)
        protChimera.setObjLabel("chimera fit\n pdb and associated volume\n ("
                                "just moves the structure\nto "
                                "the starting position) "
                                "\nsave "
                                "volume and model")
        self.launchProtocol(protChimera)
        structure3_PDB = protChimera.outputPdb_01
        extraCommands = ""
        extraCommands += "runCommand('move -24.11,-45.76,-24.60 model #2 " \
                         "coord #1')\n"
        extraCommands += "runCommand('fitmap #2 #1')\n"
        extraCommands += "runCommand('scipionwrite model #2 refmodel #1 " \
                         "saverefmodel 0')\n"
        extraCommands += "runCommand('stop')\n"
        args = {'extraCommands': extraCommands,
                'pdbFileToBeRefined': structure3_PDB,
                }
        protChimera = self.newProtocol(ChimeraProtRigidFit, **args)
        protChimera.setObjLabel("chimera fit\n pdb and associated volume\n "
                                "save only model")
        self.launchProtocol(protChimera)

        structure4_PDB = protChimera.outputPdb_01
        extraCommands = ""
        extraCommands += "runCommand('volume #1 voxelSize 1.55')\n"
        extraCommands += "runCommand('scipionwrite model #2 refmodel #1 " \
                         "saverefmodel 1')\n"
        extraCommands += "runCommand('stop')\n"
        args = {'extraCommands': extraCommands,
                'pdbFileToBeRefined': structure4_PDB,
                }
        protChimera = self.newProtocol(ChimeraProtRigidFit, **args)
        protChimera.setObjLabel("chimera fit \n pdb and associated "
                                "volume\n volume resampled\n save volume "
                                "and model\n")
        self.launchProtocol(protChimera)
        self.assertIsNotNone(protChimera.outputPdb_01.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protChimera.outputPdb_01.
                                       getVolume().getFileName()))


class TestCootRefinement(TestImportData):
    """ Test the flexible fitting of coot refinement protocol
    """

    def testCootFlexibleFitFromPDB(self):
        """ This test checks that coot runs with an atomic structure;
         No Volume was provided and an error message is expected"""
        print "Run Coot fit from imported pdb file without imported or " \
              "pdb-associated volume"

        # import PDB
        structure_PDB = self._importStructurePDBWoVol()
        self.assertTrue(structure_PDB.getFileName())
        self.assertFalse(structure_PDB.getVolume())

        args = {'extraCommands': self._createExtraCommandLine(0., 0., 0.),
                'pdbFileToBeRefined': structure_PDB,
                'doInteractive': False
                }
        protCoot = self.newProtocol(CootRefine, **args)
        protCoot.setObjLabel('coot refinement\n no volume\n associated to pdb')

        try:
            self.launchProtocol(protCoot)
        except Exception as e:
            self.assertTrue(True)
            print "This test should return a error message as '" \
                  " ERROR running protocol scipion - coot refinement"

            return
        self.assertTrue(False)

    def testCootFlexibleFitFromUnfittedVolAndPDB(self):
        """ This test checks that coot runs with a volume provided
        directly as inputVol, input PDB (not previously fitted with Chimera)
         """
        print "Run Coot fit from imported volume and pdb file not fitted"

        # Import Volume
        volume = self._importVolume()

        # import PDB
        structure_PDB = self._importStructurePDBWoVol()

        listVolCoot = [volume]
        args = {'extraCommands': self._createExtraCommandLine(-24.11, -45.76,
                                                              -24.60),
                'inputVolumes': listVolCoot,
                'pdbFileToBeRefined': structure_PDB,
                'doInteractive': False
                }
        protCoot = self.newProtocol(CootRefine, **args)
        protCoot.setObjLabel('coot refinement\n volume and unfitted pdb\n '
                             'save model')
        self.launchProtocol(protCoot)
        self.assertIsNotNone(protCoot.outputPdb_0001.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protCoot.outputPdb_0001.getFileName()))
        self.assertTrue(
            os.path.exists(protCoot.output3DMap_0001.getFileName()))

    def testCootFlexibleFitFromVolAndPDB(self):
        """ This test checks that coot runs with a volume provided
        directly as inputVol, input PDB """
        print "Run Coot fit from imported volume and pdb file"

        # Import Volume
        volume = self._importVolume()

        # import PDB
        structure_PDB = self._importStructurePDBWoVol()

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
                'pdbFileToBeRefined': structure_PDB
                }
        protChimera = self.newProtocol(ChimeraProtRigidFit, **args)
        protChimera.setObjLabel('chimera fit\n volume and pdb\n '
                             'save volume and model')
        self.launchProtocol(protChimera)

        structure2_PDB = protChimera.outputPdb_01
        volume2 = protChimera.output3Dmap

        listVolCoot = [volume2]
        args = {'extraCommands': self._createExtraCommandLine(0., 0., 0.),
                'inputVolumes': listVolCoot,
                'pdbFileToBeRefined': structure2_PDB,
                'doInteractive': False
                }
        protCoot = self.newProtocol(CootRefine, **args)
        protCoot.setObjLabel('coot refinement\n volume and fitted pdb\n '
                             'save model')

        self.launchProtocol(protCoot)
        self.assertIsNotNone(protCoot.outputPdb_0001.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protCoot.outputPdb_0001.getFileName()))
        self.assertTrue(
            os.path.exists(protCoot.output3DMap_0001.getFileName()))

    def testCootFlexibleFitFromVolAssocToPDB(self):

        # This test checks that coot runs when a volume is provided
        # associated to the input PDB
        print "Run Coot fit from imported pdb file and volume associated "

        # Import Volume
        volume = self._importVolume()

        # import PDB
        structure_PDB = self._importStructurePDBWoVol()

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
                'pdbFileToBeRefined': structure_PDB
                }
        protChimera = self.newProtocol(ChimeraProtRigidFit, **args)
        protChimera.setObjLabel('chimera fit\n volume and pdb\n '
                             'save volume and model')
        self.launchProtocol(protChimera)

        structure2_PDB = protChimera.outputPdb_01

        args = {'extraCommands': self._createExtraCommandLine(0., 0., 0.),
                'pdbFileToBeRefined': structure2_PDB,
                'doInteractive': False
                }
        protCoot = self.newProtocol(CootRefine, **args)
        protCoot.setObjLabel('coot refinement\n pdb and associated volume\n '
                             'save model')

        self.launchProtocol(protCoot)
        self.assertIsNotNone(protCoot.outputPdb_0001.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protCoot.outputPdb_0001.getFileName()))
        self.assertTrue(
            os.path.exists(protCoot.output3DMap_0001.getFileName()))

    def testCootFlexibleFitFromtwoVolAndPDB(self):
        """ This test checks that coot runs with two volumes provided
        directly as inputVol, input PDB """
        print "Run Coot fit from imported volume and pdb file"

        # Import Volume
        volume = self._importVolume()

        # import PDB
        structure_PDB = self._importStructurePDBWoVol()

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
                'pdbFileToBeRefined': structure_PDB
                }
        protChimera = self.newProtocol(ChimeraProtRigidFit, **args)
        protChimera.setObjLabel('chimera fit\n volume and pdb\n save '
                                'volume and model')
        self.launchProtocol(protChimera)

        volume2 = protChimera.output3Dmap

        structureCoot_PDB = self._importCootStructureWoVol()

        listVolCoot = [volume, volume2]
        args = {'extraCommands': self._createExtraCommandLine(0., 0., 0.),
                'inputVolumes': listVolCoot,
                'pdbFileToBeRefined': structureCoot_PDB,
                'doInteractive': False
                }
        protCoot = self.newProtocol(CootRefine, **args)
        protCoot.setObjLabel('coot refinement\n two volumes and pdb\n '
                             'save model')

        self.launchProtocol(protCoot)
        self.assertIsNotNone(protCoot.outputPdb_0001.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protCoot.outputPdb_0001.getFileName()))
        self.assertTrue(
            os.path.exists(protCoot.output3DMap_0001.getFileName()))

    def testCootFitFromPDBFromChimera_2(self):
        # This test checks that coot runs once when a volume is provided
        # associated to the input PDB file after Chimera
        # workflow, and not directly as inputVol
        # starting volume with a different coordinate origin
        print "Run Coot fit from PDB file saved from Chimera_2"

        volume2 = self._importVolume2()
        structure1_PDB = self._importStructurePDBWoVol()

        # chimera fit

        extraCommands = ""
        extraCommands += "runCommand('fitmap #2 #1')\n"
        extraCommands += "runCommand('scipionwrite model #2 refmodel #1 " \
                         "saverefmodel 1')\n"
        extraCommands += "runCommand('stop')\n"

        args = {'extraCommands': extraCommands,
                'inputVolume': volume2,
                'pdbFileToBeRefined': structure1_PDB
                }

        protChimera = self.newProtocol(ChimeraProtRigidFit, **args)
        protChimera.setObjLabel('chimera fit\n pdb and volume associated\n '
                             'save volume and model')
        self.launchProtocol(protChimera)

        volume = protChimera.output3Dmap
        structure3_PDB = protChimera.outputPdb_01

        listVolCoot = [volume]
        args = {'extraCommands': self._createExtraCommandLine(0., 0., 0.),
                'inputVolumes': listVolCoot,
                'pdbFileToBeRefined': structure3_PDB,
                'doInteractive': False
                }
        protCoot = self.newProtocol(CootRefine, **args)
        protCoot.setObjLabel('coot refinement\n volume and pdb\n '
                             'save model')

        self.launchProtocol(protCoot)
        self.assertIsNotNone(protCoot.outputPdb_0001.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protCoot.outputPdb_0001.getFileName()))
        self.assertTrue(
            os.path.exists(protCoot.output3DMap_0001.getFileName()))

    def testMultipleCootFit(self):
        # This test checks that coot runs twice when a volume is provided
        # associated to the input PDB file after Chimera
        # workflow, and not directly as inputVol
        # starting volume with a different coordinate origin
        print "Run Coot fit from PDB file saved from Chimera_2"

        volume2 = self._importVolume2()
        structure1_PDB = self._importStructurePDBWoVol()

        # chimera fit

        extraCommands = ""
        extraCommands += "runCommand('fitmap #2 #1')\n"
        extraCommands += "runCommand('scipionwrite model #2 refmodel #1 " \
                         "saverefmodel 1')\n"
        extraCommands += "runCommand('stop')\n"

        args = {'extraCommands': extraCommands,
                'inputVolume': volume2,
                'pdbFileToBeRefined': structure1_PDB
                }

        protChimera = self.newProtocol(ChimeraProtRigidFit, **args)
        protChimera.setObjLabel('chimera fit\n pdb and associated volume\n '
                             'save volume and model')
        self.launchProtocol(protChimera)

        volume = protChimera.output3Dmap
        structure3_PDB = protChimera.outputPdb_01

        listVolCoot = [volume]
        args = {'extraCommands': self._createExtraCommandLine(0., 0., 0.),
                'inputVolumes': listVolCoot,
                'pdbFileToBeRefined': structure3_PDB,
                'doInteractive': True
                }
        protCoot = self.newProtocol(CootRefine, **args)
        protCoot.setObjLabel('coot refinement\n volume pdb\n '
                             'save model twice')
        try:
            self.launchProtocol(protCoot)
        except:
            print "first call to coot ended"
        self.assertIsNotNone(protCoot.outputPdb_0001.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protCoot.outputPdb_0001.getFileName()))
        self.assertTrue(
            os.path.exists(protCoot.output3DMap_0001.getFileName()))

        protCoot.doInteractive.set(False)
        self.launchProtocol(protCoot)
        self.assertIsNotNone(protCoot.outputPdb_0002.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protCoot.outputPdb_0002.getFileName()))
        self.assertTrue(
            os.path.exists(protCoot.output3DMap_0001.getFileName()))


class TestRefmacRefinement(TestImportData):
    """ Test the flexible fitting of refmac refinement protocol
    """

    def testRefmacFlexibleFitFromPDB(self):
        """ This test checks that refmac runs with an atomic structure;
         No Volume was provided and an error message is expected"""
        print "Run Refmac refinement from imported pdb file without imported " \
              "or pdb-associated volume"

        # import PDB
        structure_PDB = self._importStructurePDBWoVol()
        self.assertTrue(structure_PDB.getFileName())
        self.assertFalse(structure_PDB.getVolume())
        args = {'inputStructure': structure_PDB
                }

        protRefmac = self.newProtocol(CCP4ProtRunRefmac, **args)
        protRefmac.setObjLabel('refmac refinement\n no volume associated '
                               'to pdb\n save model')

        try:
            self.launchProtocol(protRefmac)
        except Exception as e:
            self.assertTrue(True)
            print "This test should return a error message as '" \
                  " ERROR running protocol scipion - coot refinement"

            return
        self.assertTrue(False)

    def testRefmacFlexibleFitAfterCoot(self):
        """ This test checks that refmac runs with a volume provided
        directly as inputVol, the input PDB was fitted to the volume and
        refined previously by coot
         """
        print "Run Refmac refinement from imported volume and pdb file " \
              "fitted and refined by Coot"

        # Import Volume
        volume = self._importVolume()

        # import PDB
        structure_PDB = self._importStructurePDBWoVol()

        # coot
        listVolCoot = [volume]
        args = {'extraCommands': self._createExtraCommandLine(-24.11, -45.76,
                                                              -24.60),
                'inputVolumes': listVolCoot,
                'pdbFileToBeRefined': structure_PDB,
                'doInteractive': False
                }
        protCoot = self.newProtocol(CootRefine, **args)
        protCoot.setObjLabel('coot refinement\n volume and pdb\n save model')
        self.launchProtocol(protCoot)

        # refmac
        coot_PDB = protCoot.outputPdb_0001
        args = {'inputVolume': volume,
                'inputStructure': coot_PDB
                }
        protRefmac = self.newProtocol(CCP4ProtRunRefmac, **args)
        protRefmac.setObjLabel('refmac refinement\n volume and pdb\n save '
                               'model')
        self.launchProtocol(protRefmac)
        self.assertIsNotNone(protRefmac.outputPdb.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protRefmac.outputPdb.getFileName()))

    def testRefmacFlexibleFitAfterChimeraAndCoot(self):
        """ This test checks that refmac runs with a volume provided
        by Chimera, the input PDB is provided by Coot """
        print "Run Refmac refinement from volume provided by Chimera " \
              "and pdb file provided by Coot"

        # Import Volume
        volume = self._importVolume()

        # import PDB
        structure_PDB = self._importStructurePDBWoVol()

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
                'pdbFileToBeRefined': structure_PDB
                }
        protChimera = self.newProtocol(ChimeraProtRigidFit, **args)
        protChimera.setObjLabel('chimera fit\n volume and pdb\n save volume '
                                'and model')
        self.launchProtocol(protChimera)

        structure2_PDB = protChimera.outputPdb_01
        volume2 = protChimera.output3Dmap

        listVolCoot = [volume2]
        args = {'extraCommands': self._createExtraCommandLine(0., 0., 0.),
                'inputVolumes': listVolCoot,
                'pdbFileToBeRefined': structure2_PDB,
                'doInteractive': False
                }
        protCoot = self.newProtocol(CootRefine, **args)
        protCoot.setObjLabel('coot refinement\n volume and pdb\n save model')
        self.launchProtocol(protCoot)

        # refmac
        coot_PDB = protCoot.outputPdb_0001
        args = {'inputVolume': volume2,
                'inputStructure': coot_PDB
                }

        protRefmac = self.newProtocol(CCP4ProtRunRefmac, **args)
        protRefmac.setObjLabel('refmac refinement\n volume and pdb\n save '
                               'model')
        self.launchProtocol(protRefmac)
        self.assertIsNotNone(protRefmac.outputPdb.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protRefmac.outputPdb.getFileName()))

    def testRefmacRefinementAfterMultipleCootFit(self):
        # This test checks that refmac runs when a volume provided
        # by Chimera workflow
        # the PDB is provided by Coot
        # starting volume with a different coordinate origin
        print "Run Refmac refinement from PDB file saved from " \
              "Chimera_2/Coot"

        volume2 = self._importVolume2()
        structure1_PDB = self._importStructurePDBWoVol()

        # chimera fit

        extraCommands = ""
        extraCommands += "runCommand('fitmap #2 #1')\n"
        extraCommands += "runCommand('scipionwrite model #2 refmodel #1 " \
                         "saverefmodel 1')\n"
        extraCommands += "runCommand('stop')\n"

        args = {'extraCommands': extraCommands,
                'inputVolume': volume2,
                'pdbFileToBeRefined': structure1_PDB
                }

        protChimera = self.newProtocol(ChimeraProtRigidFit, **args)
        protChimera.setObjLabel('chimera fit\n volume associated '
                               'to pdb\n save volume and model')
        self.launchProtocol(protChimera)

        volume = protChimera.output3Dmap
        structure3_PDB = protChimera.outputPdb_01

        listVolCoot = [volume]
        args = {'extraCommands': self._createExtraCommandLine(0., 0., 0.),
                'inputVolumes': listVolCoot,
                'pdbFileToBeRefined': structure3_PDB,
                'doInteractive': True
                }
        protCoot = self.newProtocol(CootRefine, **args)
        protCoot.setObjLabel('coot refinement\n volume and pdb\n 2 runs\n '
                             'save model')
        try:
            self.launchProtocol(protCoot)
        except:
            print "first call to coot ended"
        self.assertIsNotNone(protCoot.outputPdb_0001.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protCoot.outputPdb_0001.getFileName()))
        self.assertTrue(
            os.path.exists(protCoot.output3DMap_0001.getFileName()))

        protCoot.doInteractive.set(False)
        self.launchProtocol(protCoot)

        # refmac
        coot_PDB = protCoot.outputPdb_0002
        args = {'inputVolume': volume,
                'inputStructure': coot_PDB
                }

        protRefmac = self.newProtocol(CCP4ProtRunRefmac, **args)
        protRefmac.setObjLabel('refmac refinement\n volume and '
                               'pdb\n save model')
        self.launchProtocol(protRefmac)
        self.assertIsNotNone(protRefmac.outputPdb.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protRefmac.outputPdb.getFileName()))




