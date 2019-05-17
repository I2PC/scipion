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



import os.path
from pyworkflow.tests import *
import json
from pyworkflow.tests import *
from pyworkflow.utils import importFromPlugin
from pyworkflow.em.protocol.protocol_import import ProtImportPdb, \
    ProtImportVolumes
from pyworkflow.tests import *
import os.path
import json


ChimeraProtRigidFit = importFromPlugin('chimera.protocols',
                                       'ChimeraProtRigidFit', doRaise=True)
CootRefine = importFromPlugin('ccp4.protocols', 'CootRefine', doRaise=True)
CCP4ProtRunRefmac = importFromPlugin('ccp4.protocols', 'CCP4ProtRunRefmac')
PhenixProtRunEMRinger = importFromPlugin('phenix.protocols',
                                         'PhenixProtRunEMRinger', doRaise=True)
PhenixProtRunMolprobity = importFromPlugin('phenix.protocols',
                                           'PhenixProtRunMolprobity')

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

    def _importVolume3(self):
        args = {'filesPath': self.dsModBuild.getFile(
            'volumes/emd_4116.map'),
                'samplingRate': 0.637,
                'setOrigCoord': False
                }
        protImportVol = self.newProtocol(ProtImportVolumes, **args)
        protImportVol.setObjLabel('import volume emd_4116\nwith default '
                                  'origin\n')
        self.launchProtocol(protImportVol)
        volume3 = protImportVol.outputVolume
        return volume3

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

    def _importStructurePDBWoVol2(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/3i3e_fitted.pdb'),
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import pdb\n 3i3e_fitted')
        self.launchProtocol(protImportPDB)
        structure5_PDB = protImportPDB.outputPdb
        return structure5_PDB

    def _importStructureMolProbity1(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/jlv_chimeraOut0001.pdb'),
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import pdb\n jlv_chimeraOut0001')
        self.launchProtocol(protImportPDB)
        structure6_PDB = protImportPDB.outputPdb
        return structure6_PDB

    def _importStructureMolProbity2(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/jlv_cootOut0016.pdb'),
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import pdb\n jlv_cootOut0016')
        self.launchProtocol(protImportPDB)
        structure7_PDB = protImportPDB.outputPdb
        return structure7_PDB

    def _createExtraCommandLine(self, x, y, z, label=None):
        if label is not None:
            writeLine = 'scipion_write(0, "%s")' % label
        else:
            writeLine = 'scipion_write()'

        if (x != 0. or y != 0. or z != 0.):
            return """translate_molecule_by(0, %f, %f, %f)
fit_molecule_to_map_by_random_jiggle(0,7000,2)
%s
coot_real_exit(0)
""" % (x, y, z, writeLine)
        else:
            return """%s
coot_real_exit(0)
""" % writeLine


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
        label = 'testLabel'
        args = {'extraCommands': self._createExtraCommandLine(-24.11, -45.76,
                                                              -24.60, label),
                'inputVolumes': listVolCoot,
                'pdbFileToBeRefined': structure_PDB,
                'doInteractive': False
                }
        protCoot = self.newProtocol(CootRefine, **args)
        protCoot.setObjLabel('coot refinement\n volume and unfitted pdb\n '
                             'save model')

        self.launchProtocol(protCoot)
        self.assertIsNotNone(protCoot.testLabel.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protCoot.testLabel.getFileName()))
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
        label = 'testLabel2'
        args = {'extraCommands': self._createExtraCommandLine(0., 0., 0.,
                                                              label),
                'inputVolumes': listVolCoot,
                'pdbFileToBeRefined': structure2_PDB,
                'doInteractive': False
                }
        protCoot = self.newProtocol(CootRefine, **args)
        protCoot.setObjLabel('coot refinement\n volume and fitted pdb\n '
                             'save model')

        self.launchProtocol(protCoot)
        self.assertIsNotNone(protCoot.testLabel2.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protCoot.testLabel2.getFileName()))
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

        label = 'testLabel3'
        args = {'extraCommands': self._createExtraCommandLine(0., 0., 0.,
                                                              label),
                'pdbFileToBeRefined': structure2_PDB,
                'doInteractive': False
                }
        protCoot = self.newProtocol(CootRefine, **args)
        protCoot.setObjLabel('coot refinement\n pdb and associated volume\n '
                             'save model')

        self.launchProtocol(protCoot)
        self.assertIsNotNone(protCoot.testLabel3.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protCoot.testLabel3.getFileName()))
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

        label = 'testLabel4'
        listVolCoot = [volume, volume2]
        args = {'extraCommands': self._createExtraCommandLine(0., 0., 0.,
                                                              label),
                'inputVolumes': listVolCoot,
                'pdbFileToBeRefined': structureCoot_PDB,
                'doInteractive': False
                }
        protCoot = self.newProtocol(CootRefine, **args)
        protCoot.setObjLabel('coot refinement\n two volumes and pdb\n '
                             'save model')

        self.launchProtocol(protCoot)
        self.assertIsNotNone(protCoot.testLabel4.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protCoot.testLabel4.getFileName()))
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

        label = 'testLabel5'
        listVolCoot = [volume]
        args = {'extraCommands': self._createExtraCommandLine(0., 0., 0.,
                                                              label),
                'inputVolumes': listVolCoot,
                'pdbFileToBeRefined': structure3_PDB,
                'doInteractive': False
                }
        protCoot = self.newProtocol(CootRefine, **args)
        protCoot.setObjLabel('coot refinement\n volume and pdb\n '
                             'save model')

        self.launchProtocol(protCoot)
        self.assertIsNotNone(protCoot.testLabel5.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protCoot.testLabel5.getFileName()))
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

        #first coot
        label = 'testLabel6'
        listVolCoot = [volume]
        args = {'extraCommands': self._createExtraCommandLine(0., 0., 0.,
                                                              label),
                'inputVolumes': listVolCoot,
                'pdbFileToBeRefined': structure3_PDB,
                'doInteractive': True
                }
        protCoot = self.newProtocol(CootRefine, **args)
        protCoot.setObjLabel('coot refinement\n volume pdb\n '
                             'save model three times')

        try:
            self.launchProtocol(protCoot)
        except:
            print "first call to coot ended"
        self.assertIsNotNone(protCoot.testLabel6.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protCoot.testLabel6.getFileName()))
        self.assertTrue(
            os.path.exists(protCoot.output3DMap_0001.getFileName()))

        #second coot (label=None)
        newExtraCommands = self._createExtraCommandLine(0., 0., 0.)
        protCoot.extraCommands.set(newExtraCommands)

        try:
            self.launchProtocol(protCoot)
        except:
            print "second call to coot ended"
        self.assertIsNotNone(protCoot.cootOut0001.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protCoot.cootOut0001.getFileName()))
        self.assertTrue(
            os.path.exists(protCoot.output3DMap_0001.getFileName()))

        #third coot
        protCoot.doInteractive.set(False)
        label = 'lastTestLabel'
        lastExtraCommands = self._createExtraCommandLine(0., 0., 0., label)
        protCoot.extraCommands.set(lastExtraCommands)
        try:
            self.launchProtocol(protCoot)
        except:
            print "third call to coot ended"
        self.assertIsNotNone(protCoot.lastTestLabel.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protCoot.lastTestLabel.getFileName()))
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
                  " ERROR running protocol scipion - refmac refinement"

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
        label = 'testLabel'
        args = {'extraCommands': self._createExtraCommandLine(-24.11, -45.76,
                                                              -24.60, label),
                'inputVolumes': listVolCoot,
                'pdbFileToBeRefined': structure_PDB,
                'doInteractive': False
                }
        protCoot = self.newProtocol(CootRefine, **args)
        protCoot.setObjLabel('coot refinement\n volume and unfitted pdb\n '
                             'save model')
        self.launchProtocol(protCoot)

        # refmac
        coot_PDB = protCoot.testLabel
        args = {'inputVolume': volume,
                'inputStructure': coot_PDB,
                'generateMaskedVolume': False
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
        label = 'testLabel2'
        args = {'extraCommands': self._createExtraCommandLine(0., 0., 0.,
                                                              label),
                'inputVolumes': listVolCoot,
                'pdbFileToBeRefined': structure2_PDB,
                'doInteractive': False
                }
        protCoot = self.newProtocol(CootRefine, **args)
        protCoot.setObjLabel('coot refinement\n volume and pdb\n save model')
        self.launchProtocol(protCoot)

        # refmac
        coot_PDB = protCoot.testLabel2
        args = {'inputVolume': volume2,
                'inputStructure': coot_PDB,
                'generateMaskedVolume': False
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
        label = 'testLabel3'
        args = {'extraCommands': self._createExtraCommandLine(0., 0., 0.,
                                                              label),
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
        self.assertIsNotNone(protCoot.testLabel3.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protCoot.testLabel3.getFileName()))
        self.assertTrue(
            os.path.exists(protCoot.output3DMap_0001.getFileName()))

        protCoot.doInteractive.set(False)
        label = 'testLabel4'
        newExtraCommands = self._createExtraCommandLine(0., 0., 0., label)
        protCoot.extraCommands.set(newExtraCommands)
        self.launchProtocol(protCoot)

        # refmac
        coot_PDB = protCoot.testLabel4
        args = {'inputVolume': volume,
                'inputStructure': coot_PDB,
                'generateMaskedVolume': False
                }

        protRefmac = self.newProtocol(CCP4ProtRunRefmac, **args)
        protRefmac.setObjLabel('refmac refinement\n volume and '
                               'pdb\n save model')
        self.launchProtocol(protRefmac)
        self.assertIsNotNone(protRefmac.outputPdb.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protRefmac.outputPdb.getFileName()))


class TestEMRingerValidation(TestImportData):
    """ Test the protocol of EMRinger validation
    """
    def checkResults(self, optThresh, rotRatio, maxZscore, modLength,
                     EMScore, protEMRinger, places=2):
        # method to check EMRinger statistic results of the Final Results Table
        textFileName = protEMRinger._getExtraPath(
            protEMRinger.EMRINGERTRANSFERFILENAME.replace('py', 'txt'))
        with open(textFileName, "r") as f:
            self.resultsDict = json.loads(str(f.read()))
            self.assertAlmostEqual(self.resultsDict[
                                'Optimal Threshold'], optThresh, delta=0.5)
            self.assertAlmostEqual(self.resultsDict[
                                'Rotamer-Ratio'], rotRatio, delta=0.5)
            self.assertAlmostEqual(self.resultsDict[
                                'Max Zscore'], maxZscore, delta=0.5)
            self.assertAlmostEqual(self.resultsDict[
                                'Model Length'], modLength, places)
            self.assertAlmostEqual(self.resultsDict[
                                'EMRinger Score'], EMScore, delta=0.5)

    def testEMRingerValidationFromPDB(self):
        """ This test checks that EMRinger validation protocol runs with an
        atomic structure; No Volume was provided and an error message is
        expected"""
        print "Run EMRinger validation protocol from imported pdb file " \
              "without imported or pdb-associated volume"

        # import PDB
        structure_PDB = self._importStructurePDBWoVol()
        self.assertTrue(structure_PDB.getFileName())
        self.assertFalse(structure_PDB.getVolume())
        args = {'inputStructure': structure_PDB
                }

        protEMRinger = self.newProtocol(PhenixProtRunEMRinger, **args)
        protEMRinger.setObjLabel('EMRinger validation\n no volume associated '
                                 'to pdb\n')

        try:
            self.launchProtocol(protEMRinger)
        except Exception as e:
            self.assertTrue(True)
            print "This test should return a error message as '" \
                  " Error: You should provide a volume.\n"

            return
        self.assertTrue(False)

    def testEMRingerValidationAfterRefmacNoMask(self):
        """ This test checks that EMRinger validation protocol runs with a
        volume provided directly as inputVol, the input PDB was fitted to
        the volume and refined previously by coot and refmac(without mask)
         """
        print "Run EMRinger validation from imported volume and pdb file " \
              "fitted and refined by Coot/Refmac without mask"

        # Import Volume
        volume = self._importVolume()

        # import PDB
        structure_PDB = self._importStructurePDBWoVol()

        # coot
        label = 'testLabel1'
        listVolCoot = [volume]
        args = {'extraCommands': self._createExtraCommandLine(-24.11, -45.76,
                                                              -24.60, label),
                'inputVolumes': listVolCoot,
                'pdbFileToBeRefined': structure_PDB,
                'doInteractive': False
                }
        protCoot = self.newProtocol(CootRefine, **args)
        protCoot.setObjLabel('coot refinement\n volume and pdb\n save model')
        self.launchProtocol(protCoot)

        # refmac without mask
        coot_PDB = protCoot.testLabel1
        args = {'inputVolume': volume,
                'inputStructure': coot_PDB,
                'generateMaskedVolume': False
                }
        protRefmac = self.newProtocol(CCP4ProtRunRefmac, **args)
        protRefmac.setObjLabel('refmac refinement\n volume and pdb\n save '
                               'model')
        self.launchProtocol(protRefmac)
        self.assertIsNotNone(protRefmac.outputPdb.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protRefmac.outputPdb.getFileName()))

        # EMRinger
        refmac_PDB = protRefmac.outputPdb
        args = {'inputVolume': volume,
                'inputStructure': refmac_PDB,
                'doTest': True
                }
        protEMRinger = self.newProtocol(PhenixProtRunEMRinger, **args)
        protEMRinger.setObjLabel('EMRinger validation\n volume and pdb\n')
        self.launchProtocol(protEMRinger)

        # check EMRinger results
        self.checkResults(optThresh = 0.606299973426968,
                          rotRatio = 0.8181818181818182,
                          maxZscore = 5.521788316969326,
                          modLength = 121,
                          EMScore = 5.019807560881206,
                          protEMRinger = protEMRinger)

    def testEMRingerValidationAfterRefmacWithMask(self):
        """ This test checks that EMRinger validation protocol runs with a
        volume provided directly as inputVol, the input PDB was fitted to
        the volume and refined previously by coot and refmac(with mask)
        """
        print "Run EMRinger validation from imported volume and pdb file " \
              "fitted and refined by Coot/Refmac with mask"

        # Import Volume
        volume = self._importVolume()

        # import PDB
        structure_PDB = self._importStructurePDBWoVol()

        # coot
        label = 'testLabel2'
        listVolCoot = [volume]
        args = {
                'extraCommands': self._createExtraCommandLine(-24.11, -45.76,
                                                              -24.60, label),
                'inputVolumes': listVolCoot,
                'pdbFileToBeRefined': structure_PDB,
                'doInteractive': False
                }
        protCoot = self.newProtocol(CootRefine, **args)
        protCoot.setObjLabel(
                'coot refinement\n volume and pdb\n save model')
        self.launchProtocol(protCoot)

        # refmac with mask
        coot_PDB = protCoot.testLabel2
        args = {'inputVolume': volume,
                'inputStructure': coot_PDB,
                'generateMaskedVolume': True
                }
        protRefmac = self.newProtocol(CCP4ProtRunRefmac, **args)
        protRefmac.setObjLabel('MASK refmac refinement\n volume and pdb\n save '
                               'model')
        self.launchProtocol(protRefmac)
        self.assertIsNotNone(protRefmac.outputPdb.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protRefmac.outputPdb.getFileName()))

        # EMRinger
        refmac_PDB = protRefmac.outputPdb
        args = {'inputVolume': volume,
                'inputStructure': refmac_PDB,
                'doTest': True
                }
        protEMRinger = self.newProtocol(PhenixProtRunEMRinger, **args)
        protEMRinger.setObjLabel('EMRinger validation\n volume and pdb\n')
        self.launchProtocol(protEMRinger)

        # check EMRinger results
        self.checkResults(optThresh=0.6363458779254832,
                          rotRatio=0.8260869565217391,
                          maxZscore=5.475171299613324,
                          modLength=121,
                          EMScore=4.97742845419393,
                          protEMRinger=protEMRinger)

    def testEMRingerValidationAfterChimeraAndCootAndRefmacNoMask(self):
        """ This test checks that EMRinger validation protocol runs with a
        volume provided by Chimera, the input PDB is provided by Coot """
        print "Run EMRinger validation from volume provided by Chimera " \
              "and pdb file provided by Coot/Refmac without mask"

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

        # coot
        label = 'testLabel3'
        listVolCoot = [volume2]
        args = {'extraCommands': self._createExtraCommandLine(0., 0., 0.,
                                                              label),
                'inputVolumes': listVolCoot,
                'pdbFileToBeRefined': structure2_PDB,
                'doInteractive': False
                }
        protCoot = self.newProtocol(CootRefine, **args)
        protCoot.setObjLabel('coot refinement\n volume and pdb\n save model')
        self.launchProtocol(protCoot)
        self.assertIsNotNone(protCoot.testLabel3.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protCoot.testLabel3.getFileName()))

        # EMRinger
        coot_PDB = protCoot.testLabel3
        args = {'inputVolume': volume2,
                'inputStructure': coot_PDB,
                'doTest': True
                }

        protEMRinger = self.newProtocol(PhenixProtRunEMRinger, **args)
        protEMRinger.setObjLabel('EMRinger validation\n volume and pdb\n')
        self.launchProtocol(protEMRinger)

        # check EMRinger results
        self.checkResults(optThresh=0.8904416297331501,
                          rotRatio=0.9166666666666666,
                          maxZscore=2.607144567760445,
                          modLength=121,
                          EMScore=2.370131425236768,
                          protEMRinger=protEMRinger)

        # refmac without mask
        coot_PDB = protCoot.testLabel3
        args = {'inputVolume': volume2,
                'inputStructure': coot_PDB,
                'generateMaskedVolume': False
                }

        protRefmac = self.newProtocol(CCP4ProtRunRefmac, **args)
        protRefmac.setObjLabel('refmac refinement\n volume and pdb\n save '
                               'model')
        self.launchProtocol(protRefmac)
        self.assertIsNotNone(protRefmac.outputPdb.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protRefmac.outputPdb.getFileName()))

        # EMRinger
        refmac_PDB = protRefmac.outputPdb
        args = {'inputVolume': volume,
                'inputStructure': refmac_PDB,
                'doTest': True
                }
        protEMRinger = self.newProtocol(PhenixProtRunEMRinger, **args)
        protEMRinger.setObjLabel('EMRinger validation\n volume and pdb\n')
        self.launchProtocol(protEMRinger)

        # check EMRinger results
        self.checkResults(optThresh=0.6058890610462733,
                          rotRatio=0.8181818181818182,
                          maxZscore=5.521788316969326,
                          modLength=121,
                          EMScore=5.019807560881206,
                          protEMRinger=protEMRinger)

    def testEMRingerValidationAfterChimeraAndCootAndRefmacWithMask(self):
        """ This test checks that EMRinger validation protocol runs with a
        volume provided by Chimera, the input PDB is provided by Coot """
        print "Run EMRinger validation from volume provided by Chimera " \
              "and pdb file provided by Coot/Refmac with mask"

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
        protChimera.setObjLabel(
            'chimera fit\n volume and pdb\n save volume and model')
        self.launchProtocol(protChimera)

        structure2_PDB = protChimera.outputPdb_01
        volume2 = protChimera.output3Dmap

        # coot
        label = 'testLabel4'
        listVolCoot = [volume2]
        args = {'extraCommands': self._createExtraCommandLine(0., 0., 0.,
                                                              label),
                'inputVolumes': listVolCoot,
                'pdbFileToBeRefined': structure2_PDB,
                'doInteractive': False
                }
        protCoot = self.newProtocol(CootRefine, **args)
        protCoot.setObjLabel(
            'coot refinement\n volume and pdb\n save model')
        self.launchProtocol(protCoot)
        self.assertIsNotNone(protCoot.testLabel4.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(
            os.path.exists(protCoot.testLabel4.getFileName()))

        # EMRinger
        coot_PDB = protCoot.testLabel4
        args = {'inputVolume': volume2,
                'inputStructure': coot_PDB,
                'doTest': True
                }

        protEMRinger = self.newProtocol(PhenixProtRunEMRinger, **args)
        protEMRinger.setObjLabel('EMRinger validation\n volume and pdb\n')
        self.launchProtocol(protEMRinger)

        # check EMRinger results
        self.checkResults(optThresh=0.8904416297331501,
                          rotRatio=0.9166666666666666,
                          maxZscore=2.607144567760445,
                          modLength=121,
                          EMScore=2.370131425236768,
                          protEMRinger=protEMRinger)

        # refmac with mask
        coot_PDB = protCoot.testLabel4
        args = {'inputVolume': volume,
                'inputStructure': coot_PDB,
                'generateMaskedVolume': True
                }
        protRefmac = self.newProtocol(CCP4ProtRunRefmac, **args)
        protRefmac.setObjLabel('MASK refmac refinement\n volume and pdb\n save '
                               'model')
        self.launchProtocol(protRefmac)
        self.assertIsNotNone(protRefmac.outputPdb.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protRefmac.outputPdb.getFileName()))

        # EMRinger
        refmac_PDB = protRefmac.outputPdb
        args = {'inputVolume': volume,
                'inputStructure': refmac_PDB,
                'doTest': True
                }
        protEMRinger = self.newProtocol(PhenixProtRunEMRinger, **args)
        protEMRinger.setObjLabel('EMRinger validation\n volume and pdb\n')
        self.launchProtocol(protEMRinger)

        # check EMRinger results
        self.checkResults(optThresh=0.6357197891869828,
                          rotRatio=0.8260869565217391,
                          maxZscore=5.475171299613324,
                          modLength=121,
                          EMScore=4.97742845419393,
                          protEMRinger=protEMRinger)

    def testEMRingerValidationAfterMultipleCootAndRefmacFitNoMask(self):
        # This test checks that EMRinger runs when a volume provided
        # by Chimera workflow; the PDB is provided by Coot/Refmac without mask
        # starting volume with a different coordinate origin
        print "Run EMRinger validation from PDB file saved from " \
              "Chimera_2/Coot/Refmac without mask"

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

        # coot
        label = 'testLabel5'
        listVolCoot = [volume]
        args = {'extraCommands': self._createExtraCommandLine(0., 0., 0.,
                                                              label),
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
        self.assertIsNotNone(protCoot.testLabel5.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protCoot.testLabel5.getFileName()))
        self.assertTrue(
            os.path.exists(protCoot.output3DMap_0001.getFileName()))

        protCoot.doInteractive.set(False)
        newExtraCommands = self._createExtraCommandLine(0., 0., 0.)
        protCoot.extraCommands.set(newExtraCommands)
        self.launchProtocol(protCoot)

        # EMRinger
        coot_PDB = protCoot.cootOut0001
        args = {'inputVolume': volume,
                'inputStructure': coot_PDB,
                'doTest': True
                }
        protEMRinger = self.newProtocol(PhenixProtRunEMRinger, **args)
        protEMRinger.setObjLabel('EMRinger validation\n volume and pdb\n')
        self.launchProtocol(protEMRinger)

        # check EMRinger results
        self.checkResults(optThresh=0.671007165163965,
                          rotRatio=0.6756756756756757,
                          maxZscore=2.313625586144688,
                          modLength=121,
                          EMScore=2.103295987404262,
                          protEMRinger=protEMRinger)

        # refmac withouth mask
        coot_PDB = protCoot.cootOut0001
        args = {'inputVolume': volume,
                'inputStructure': coot_PDB,
                'generateMaskedVolume': False
                }
        protRefmac = self.newProtocol(CCP4ProtRunRefmac, **args)
        protRefmac.setObjLabel('refmac refinement\n volume and pdb\n save '
                               'model')
        self.launchProtocol(protRefmac)
        self.assertIsNotNone(protRefmac.outputPdb.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protRefmac.outputPdb.getFileName()))

        # EMRinger
        refmac_PDB = protRefmac.outputPdb
        args = {'inputVolume': volume,
                'inputStructure': refmac_PDB,
                'doTest': True
                }
        protEMRinger = self.newProtocol(PhenixProtRunEMRinger, **args)
        protEMRinger.setObjLabel('EMRinger validation\n volume and pdb\n')
        self.launchProtocol(protEMRinger)

        # check EMRinger results
        self.checkResults(optThresh=0.6118165038763741,
                          rotRatio=0.8181818181818182,
                          maxZscore=5.521788316969326,
                          modLength=121,
                          EMScore=5.019807560881206,
                          protEMRinger=protEMRinger)

    def testEMRingerValidationAfterMultipleCootAndRefmacFitWithMask(self):
        # This test checks that EMRinger runs when a volume provided
        # by Chimera workflow; the PDB is provided by Coot/Refmac with mask
        # starting volume with a different coordinate origin
        print "Run EMRinger validation from PDB file saved from " \
              "Chimera_2/Coot/Refmac with mask"

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

        # coot
        label = 'testLabel6'
        listVolCoot = [volume]
        args = {'extraCommands': self._createExtraCommandLine(0., 0., 0.,
                                                              label),
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
        self.assertIsNotNone(protCoot.testLabel6.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(
            os.path.exists(protCoot.testLabel6.getFileName()))
        self.assertTrue(
            os.path.exists(protCoot.output3DMap_0001.getFileName()))

        protCoot.doInteractive.set(False)
        newExtraCommands = self._createExtraCommandLine(0., 0., 0.)
        protCoot.extraCommands.set(newExtraCommands)
        self.launchProtocol(protCoot)

        # EMRinger
        coot_PDB = protCoot.cootOut0001
        args = {'inputVolume': volume,
                'inputStructure': coot_PDB,
                'doTest': True
                }
        protEMRinger = self.newProtocol(PhenixProtRunEMRinger, **args)
        protEMRinger.setObjLabel('EMRinger validation\n volume and pdb\n')
        self.launchProtocol(protEMRinger)

        # check EMRinger results
        self.checkResults(optThresh=0.671007165163965,
                          rotRatio=0.6756756756756757,
                          maxZscore=2.313625586144688,
                          modLength=121,
                          EMScore=2.103295987404262,
                          protEMRinger=protEMRinger)

        # refmac with mask
        coot_PDB = protCoot.cootOut0001
        args = {'inputVolume': volume,
                'inputStructure': coot_PDB,
                'generateMaskedVolume': True
                }
        protRefmac = self.newProtocol(CCP4ProtRunRefmac, **args)
        protRefmac.setObjLabel('MASK refmac refinement\n volume and pdb\n save '
                               'model')
        self.launchProtocol(protRefmac)
        self.assertIsNotNone(protRefmac.outputPdb.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protRefmac.outputPdb.getFileName()))

        # EMRinger
        refmac_PDB = protRefmac.outputPdb
        args = {'inputVolume': volume,
                'inputStructure': refmac_PDB,
                'doTest': True
                }
        protEMRinger = self.newProtocol(PhenixProtRunEMRinger, **args)
        protEMRinger.setObjLabel('EMRinger validation\n volume and pdb\n')
        self.launchProtocol(protEMRinger)

        # check EMRinger results
        self.checkResults(optThresh=0.5640620255535769,
                          rotRatio=0.780952380952381,
                          maxZscore=4.921014490568343,
                          modLength=121,
                          EMScore=4.473649536880312,
                          protEMRinger=protEMRinger)

    def testEMRingerValidationSeveralChains(self):
        """ This test checks that EMRinger validation protocol runs with a
        volume provided directly as inputVol, the input PDB was fitted to
        the volume and refined previously by coot and refmac in another project
         """
        print "Run EMRinger validation from imported volume and pdb file " \
              "already refined by Coot and Refmac in another project"

        # Import Volume
        volume3 = self._importVolume3()

        # import PDB
        structure5_PDB = self._importStructurePDBWoVol2()

        # EMRinger
        args = {'inputVolume': volume3,
                'inputStructure': structure5_PDB,
                'doTest': True
                }
        protEMRinger = self.newProtocol(PhenixProtRunEMRinger, **args)
        protEMRinger.setObjLabel('EMRinger validation\n volume and pdb\n')
        self.launchProtocol(protEMRinger)

        # check EMRinger results
        self.checkResults(optThresh=0.016682716149338486,
                          rotRatio=0.8158347676419966,
                          maxZscore=26.526393559321498,
                          modLength=2587,
                          EMScore=5.21530839370391,
                          protEMRinger=protEMRinger)

class TestMolprobityValidation(TestImportData):
    """ Test the protocol of MolProbity validation
    """
    def checkResults(self, ramOutliers, ramFavored, rotOutliers, cbetaOutliers,
                     clashScore, overallScore, protMolProbity, places=0):
        # method to check MolProbity statistic results of the Final Results
        # Table
        self.assertAlmostEqual(protMolProbity.ramachandranOutliers.get(),
                               ramOutliers, places)
        self.assertAlmostEqual(protMolProbity.ramachandranFavored.get(),
                               ramFavored, places=0)
        self.assertAlmostEqual(protMolProbity.rotamerOutliers.get(),
                               rotOutliers, places)
        self.assertAlmostEqual(protMolProbity.cbetaOutliers.get(),
                               cbetaOutliers, delta=2)
        self.assertAlmostEqual(protMolProbity.clashscore.get(),
                               clashScore, places)
        self.assertAlmostEqual(protMolProbity.overallScore.get(),
                               overallScore, places)

    def testMolProbityValidationFromPDB(self):
        """ This test checks that EMRinger validation protocol runs with an
        atomic structure; No Volume was provided and no error message is
        expected"""
        print "Run MolProbity validation protocol from imported pdb file " \
              "without imported or pdb-associated volume"

        # import PDB
        structure_PDB = self._importStructurePDBWoVol()
        self.assertTrue(structure_PDB.getFileName())
        self.assertFalse(structure_PDB.getVolume())
        args = {'inputStructure': structure_PDB
               }

        protMolProbity = self.newProtocol(PhenixProtRunMolprobity, **args)
        protMolProbity.setObjLabel('PhenixProtRunMolprobity validation\n '
                                   'no volume associated to pdb\n')
        self.launchProtocol(protMolProbity)

        # check MolProbity results
        self.checkResults(ramOutliers=0.94,
                          ramFavored=81.60,
                          rotOutliers=3.98,
                          cbetaOutliers=1,  # 0
                          clashScore=4.77,
                          overallScore=2.42,
                          protMolProbity=protMolProbity)

    def testMolProbityValidationFromVolume(self):
        """ This test checks that MolProbity validation protocol runs with a
        density volume; No atomic structure was provided and a error message is
        expected"""

        print "Run MolProbity validation protocol from imported volume file " \
          "without imported pdb"

        # import volume
        volume = self._importVolume()
        self.assertTrue(volume.getFileName())

        args = {'inputVolume': volume,
                'resolution': 3.5
               }

        protMolProbity = self.newProtocol(PhenixProtRunMolprobity, **args)
        protMolProbity.setObjLabel(
        'PhenixProtRunMolprobity validation\n volume and no pdb\n')

        try:
            self.launchProtocol(protMolProbity)
        except Exception as e:
            self.assertTrue(True)
            print "This test should return a error message as '" \
              " Input atomic structure cannot be EMPTY.\n"

            return
        self.assertTrue(False)

    def testMolprobityValidationAfterRefmacNoMask(self):
        """ This test checks that MolProbity validation protocol runs with a
        volume provided directly as inputVol, the input PDB was fitted to
        the volume and refined previously by coot and refmac without mask
         """
        print "Run MolProbity validation from imported volume and pdb file " \
              "fitted and refined by Coot/Refmac"

        # Import Volume
        volume = self._importVolume()

        # import PDB
        structure_PDB = self._importStructurePDBWoVol()

        # coot
        label = 'testLabel1'
        listVolCoot = [volume]
        args = {'extraCommands': self._createExtraCommandLine(-24.11, -45.76,
                                                              -24.60, label),
                'inputVolumes': listVolCoot,
                'pdbFileToBeRefined': structure_PDB,
                'doInteractive': False
                }
        protCoot = self.newProtocol(CootRefine, **args)
        protCoot.setObjLabel('coot refinement\n volume and pdb\n save model')
        self.launchProtocol(protCoot)

        # refmac withouth mask
        coot_PDB = protCoot.testLabel1
        args = {'inputVolume': volume,
                'inputStructure': coot_PDB,
                'generateMaskedVolume': False
                }
        protRefmac = self.newProtocol(CCP4ProtRunRefmac, **args)
        protRefmac.setObjLabel('refmac refinement\n volume and pdb\n save '
                               'model')
        self.launchProtocol(protRefmac)
        self.assertIsNotNone(protRefmac.outputPdb.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protRefmac.outputPdb.getFileName()))

        # MolProbity
        refmac_PDB = protRefmac.outputPdb
        args = {'inputVolume': volume,
                'resolution': 3.5,
                'inputStructure': refmac_PDB
                }
        protMolProbity = self.newProtocol(PhenixProtRunMolprobity, **args)
        protMolProbity.setObjLabel('MolProbity validation\n volume and pdb\n')
        self.launchProtocol(protMolProbity)

        # check MolProbity results
        self.checkResults(ramOutliers=0.47,
                          ramFavored=83.49,
                          rotOutliers=6.25,
                          cbetaOutliers=2,
                          clashScore=4.77,
                          overallScore=2.54,
                          protMolProbity=protMolProbity)

    def testMolprobityValidationAfterRefmacWithMask(self):
        """ This test checks that MolProbity validation protocol runs with a
        volume provided directly as inputVol, the input PDB was fitted to
        the volume and refined previously by coot and refmac with mask
         """
        print "Run MolProbity validation from imported volume and pdb file " \
              "fitted and refined by Coot/Refmac with mask"

        # Import Volume
        volume = self._importVolume()

        # import PDB
        structure_PDB = self._importStructurePDBWoVol()

        # coot
        label = 'testLabel2'
        listVolCoot = [volume]
        args = {'extraCommands': self._createExtraCommandLine(-24.11, -45.76,
                                                              -24.60, label),
                'inputVolumes': listVolCoot,
                'pdbFileToBeRefined': structure_PDB,
                'doInteractive': False
                }
        protCoot = self.newProtocol(CootRefine, **args)
        protCoot.setObjLabel('coot refinement\n volume and pdb\n save model')
        self.launchProtocol(protCoot)

        # refmac with mask
        coot_PDB = protCoot.testLabel2
        args = {'inputVolume': volume,
                'inputStructure': coot_PDB,
                'generateMaskedVolume': True
                }
        protRefmac = self.newProtocol(CCP4ProtRunRefmac, **args)
        protRefmac.setObjLabel('MASK refmac refinement\n volume and pdb\n save '
                               'model')
        self.launchProtocol(protRefmac)
        self.assertIsNotNone(protRefmac.outputPdb.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protRefmac.outputPdb.getFileName()))

        # MolProbity
        refmac_PDB = protRefmac.outputPdb
        args = {'inputVolume': volume,
                'resolution': 3.5,
                'inputStructure': refmac_PDB
                }
        protMolProbity = self.newProtocol(PhenixProtRunMolprobity, **args)
        protMolProbity.setObjLabel('MolProbity validation\n volume and pdb\n')
        self.launchProtocol(protMolProbity)

        # check MolProbity results
        self.checkResults(ramOutliers=0.47,
                          ramFavored=83.02,  # 83.49 TestMolprobityValidation.testMolProbityValidationAfterMultipleCootAndRefmacFitWithMask
                          rotOutliers=5.68,
                          cbetaOutliers=0,
                          clashScore=4.47,
                          overallScore=2.48,
                          protMolProbity=protMolProbity)

    def testMolProbityValidationAfterChimeraAndCootAndRefmacNoMask(self):
        """ This test checks that MolProbity validation protocol runs with a
        volume provided by Chimera, the input PDB is provided by Coot """
        print "Run MolProbity validation from volume provided by Chimera " \
              "and pdb file provided by Coot/Refmac without mask"

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

        # coot
        label = 'testLabel3'
        listVolCoot = [volume2]
        args = {'extraCommands': self._createExtraCommandLine(0., 0., 0.,
                                                              label),
                'inputVolumes': listVolCoot,
                'pdbFileToBeRefined': structure2_PDB,
                'doInteractive': False
                }
        protCoot = self.newProtocol(CootRefine, **args)
        protCoot.setObjLabel('coot refinement\n volume and pdb\n save model')
        self.launchProtocol(protCoot)
        self.assertIsNotNone(protCoot.testLabel3.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protCoot.testLabel3.getFileName()))

        # MolProbity
        coot_PDB = protCoot.testLabel3
        args = {'inputVolume': volume2,
                'resolution': 3.5,
                'inputStructure': coot_PDB
                }

        protMolProbity = self.newProtocol(PhenixProtRunMolprobity, **args)
        protMolProbity.setObjLabel('MolProbity validation\n volume and pdb\n')
        self.launchProtocol(protMolProbity)

        # check MolProbity results
        self.checkResults(ramOutliers=0.94,
                          ramFavored=81.60,
                          rotOutliers=3.98,
                          cbetaOutliers=0,
                          clashScore=4.77,
                          overallScore=2.42,
                          protMolProbity=protMolProbity)

        # refmac withouth mask
        coot_PDB = protCoot.testLabel3
        args = {'inputVolume': volume2,
                'inputStructure': coot_PDB,
                'generateMaskedVolume': False
                }

        protRefmac = self.newProtocol(CCP4ProtRunRefmac, **args)
        protRefmac.setObjLabel('refmac refinement\n volume and pdb\n save '
                               'model')
        self.launchProtocol(protRefmac)
        self.assertIsNotNone(protRefmac.outputPdb.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protRefmac.outputPdb.getFileName()))

        # MolProbity
        refmac_PDB = protRefmac.outputPdb
        args = {'inputVolume': volume2,
                'resolution': 3.5,
                'inputStructure': refmac_PDB
                }

        protMolProbity = self.newProtocol(PhenixProtRunMolprobity, **args)
        protMolProbity.setObjLabel('MolProbity validation\n volume and pdb\n')
        self.launchProtocol(protMolProbity)

        # check MolProbity results
        self.checkResults(ramOutliers=0.47,
                          ramFavored=83.49,
                          rotOutliers=5.68,
                          cbetaOutliers=2,
                          clashScore=4.77,
                          overallScore=2.51,
                          protMolProbity=protMolProbity)

    def testMolProbityValidationAfterChimeraAndCootAndRefmacWithMask(self):
        """ This test checks that MolProbity validation protocol runs with a
        volume provided by Chimera, the input PDB is provided by Coot """
        print "Run MolProbity validation from volume provided by Chimera " \
              "and pdb file provided by Coot/Refmac with mask"

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

        # coot
        label = 'testLabel4'
        listVolCoot = [volume2]
        args = {'extraCommands': self._createExtraCommandLine(0., 0., 0.,
                                                              label),
                'inputVolumes': listVolCoot,
                'pdbFileToBeRefined': structure2_PDB,
                'doInteractive': False
                }
        protCoot = self.newProtocol(CootRefine, **args)
        protCoot.setObjLabel('coot refinement\n volume and pdb\n save model')
        self.launchProtocol(protCoot)
        self.assertIsNotNone(protCoot.testLabel4.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protCoot.testLabel4.getFileName()))

        # MolProbity
        coot_PDB = protCoot.testLabel4
        args = {'inputVolume': volume2,
                'resolution': 3.5,
                'inputStructure': coot_PDB
                }

        protMolProbity = self.newProtocol(PhenixProtRunMolprobity, **args)
        protMolProbity.setObjLabel('MolProbity validation\n volume and pdb\n')
        self.launchProtocol(protMolProbity)

        # check MolProbity results
        self.checkResults(ramOutliers=0.94,
                          ramFavored=81.60,
                          rotOutliers=3.98,
                          cbetaOutliers=0,
                          clashScore=4.77,
                          overallScore=2.42,
                          protMolProbity=protMolProbity)

        # refmac with mask
        coot_PDB = protCoot.testLabel4
        args = {'inputVolume': volume2,
                'inputStructure': coot_PDB,
                'generateMaskedVolume': True
                }

        protRefmac = self.newProtocol(CCP4ProtRunRefmac, **args)
        protRefmac.setObjLabel('MASK refmac refinement\n volume and pdb\n save '
                               'model')
        self.launchProtocol(protRefmac)
        self.assertIsNotNone(protRefmac.outputPdb.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protRefmac.outputPdb.getFileName()))

        # MolProbity
        refmac_PDB = protRefmac.outputPdb
        args = {'inputVolume': volume2,
                'resolution': 3.5,
                'inputStructure': refmac_PDB
                }

        protMolProbity = self.newProtocol(PhenixProtRunMolprobity, **args)
        protMolProbity.setObjLabel('MolProbity validation\n volume and pdb\n')
        self.launchProtocol(protMolProbity)

        # check MolProbity results
        self.checkResults(ramOutliers=0.47,
                          ramFavored=83.96,
                          rotOutliers=5.68,
                          cbetaOutliers=0,
                          clashScore=4.47,
                          overallScore=2.47,
                          protMolProbity=protMolProbity)

    def testMolProbityValidationAfterMultipleCootAndRefmacFitNoMask(self):
        # This test checks that MolProbity runs when a volume provided
        # by Chimera workflow; the PDB is provided by Coot/Refmac
        # starting volume with a different coordinate origin
        print "Run MolProbity validation from PDB file saved from " \
              "Chimera_2/Coot/Refmac withouth mask"

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

        # coot
        label = 'testLabel5'
        listVolCoot = [volume]
        args = {'extraCommands': self._createExtraCommandLine(0., 0., 0.,
                                                              label),
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
        self.assertIsNotNone(protCoot.testLabel5.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protCoot.testLabel5.getFileName()))
        self.assertTrue(
            os.path.exists(protCoot.output3DMap_0001.getFileName()))

        protCoot.doInteractive.set(False)
        newExtraCommands = self._createExtraCommandLine(0., 0., 0.)
        protCoot.extraCommands.set(newExtraCommands)
        self.launchProtocol(protCoot)

        # MolProbity
        coot_PDB = protCoot.cootOut0001
        args = {'inputVolume': volume,
                'resolution': 3.5,
                'inputStructure': coot_PDB
                }
        protMolProbity = self.newProtocol(PhenixProtRunMolprobity, **args)
        protMolProbity.setObjLabel('MolProbity validation\n volume and pdb\n')
        self.launchProtocol(protMolProbity)

        # check MolProbity results
        self.checkResults(ramOutliers=0.94,
                          ramFavored=81.60,
                          rotOutliers=3.98,
                          cbetaOutliers=0,
                          clashScore=4.77,
                          overallScore=2.42,
                          protMolProbity=protMolProbity)

        # refmac withouth mask
        coot_PDB = protCoot.cootOut0001
        args = {'inputVolume': volume,
                'inputStructure': coot_PDB,
                'generateMaskedVolume': False
                }

        protRefmac = self.newProtocol(CCP4ProtRunRefmac, **args)
        protRefmac.setObjLabel('refmac refinement\n volume and '
                               'pdb\n save model')
        self.launchProtocol(protRefmac)
        self.assertIsNotNone(protRefmac.outputPdb.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protRefmac.outputPdb.getFileName()))

        # MolProbity
        refmac_PDB = protRefmac.outputPdb
        args = {'inputVolume': volume,
                'resolution': 3.5,
                'inputStructure': refmac_PDB
                }
        protMolProbity = self.newProtocol(PhenixProtRunMolprobity, **args)
        protMolProbity.setObjLabel('MolProbity validation\n volume and pdb\n')
        self.launchProtocol(protMolProbity)

        # check MolProbity results
        self.checkResults(ramOutliers=0.47,
                          ramFavored=83.96,
                          rotOutliers=5.68,
                          cbetaOutliers=1,
                          clashScore=5.07,
                          overallScore=2.52,
                          protMolProbity=protMolProbity)

    def testMolProbityValidationAfterMultipleCootAndRefmacFitWithMask(self):
        # This test checks that MolProbity runs when a volume provided
        # by Chimera workflow; the PDB is provided by Coot/Refmac
        # starting volume with a different coordinate origin
        print "Run MolProbity validation from PDB file saved from " \
              "Chimera_2/Coot/Refmac with mask"

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

        # coot
        label = 'testLabel6'
        listVolCoot = [volume]
        args = {'extraCommands': self._createExtraCommandLine(0., 0., 0.,
                                                              label),
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
        self.assertIsNotNone(protCoot.testLabel6.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protCoot.testLabel6.getFileName()))
        self.assertTrue(
            os.path.exists(protCoot.output3DMap_0001.getFileName()))

        protCoot.doInteractive.set(False)
        newExtraCommands = self._createExtraCommandLine(0., 0., 0.)
        protCoot.extraCommands.set(newExtraCommands)
        self.launchProtocol(protCoot)

        # MolProbity
        coot_PDB = protCoot.cootOut0001
        args = {'inputVolume': volume,
                'resolution': 3.5,
                'inputStructure': coot_PDB
                }
        protMolProbity = self.newProtocol(PhenixProtRunMolprobity, **args)
        protMolProbity.setObjLabel('MolProbity validation\n volume and pdb\n')
        self.launchProtocol(protMolProbity)

        # check MolProbity results
        self.checkResults(ramOutliers=0.94,
                          ramFavored=81.60,
                          rotOutliers=3.98,
                          cbetaOutliers=0,
                          clashScore=4.77,
                          overallScore=2.42,
                          protMolProbity=protMolProbity)

        # refmac with mask
        coot_PDB = protCoot.cootOut0001
        args = {'inputVolume': volume,
                'inputStructure': coot_PDB,
                'generateMaskedVolume': True
                }

        protRefmac = self.newProtocol(CCP4ProtRunRefmac, **args)
        protRefmac.setObjLabel('MASK refmac refinement\n volume and '
                               'pdb\n save model')
        self.launchProtocol(protRefmac)
        self.assertIsNotNone(protRefmac.outputPdb.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protRefmac.outputPdb.getFileName()))

        # MolProbity
        refmac_PDB = protRefmac.outputPdb
        args = {'inputVolume': volume,
                'resolution': 3.5,
                'inputStructure': refmac_PDB
                }
        protMolProbity = self.newProtocol(PhenixProtRunMolprobity, **args)
        protMolProbity.setObjLabel('MolProbity validation\n volume and pdb\n')
        self.launchProtocol(protMolProbity)

        # check MolProbity results
        self.checkResults(ramOutliers=0.47,
                          ramFavored=83.96,
                          rotOutliers=5.68,
                          cbetaOutliers=1,
                          clashScore=4.47,
                          overallScore=2.47,
                          protMolProbity=protMolProbity)

    def testMolProbityValidationSeveralChains(self):
        """ This test checks that MolProbity validation protocol runs with a
        volume provided directly as inputVol, the input PDB was fitted to
        the volume and refined previously by coot and refmac in another project
         """
        print "Run MolProbity validation from imported volume and pdb file " \
              "already refined by Coot and Refmac in another project"

        # Import Volume
        volume3 = self._importVolume3()

        # import PDB
        structure5_PDB = self._importStructurePDBWoVol2()

        # MolProbity
        args = {'inputVolume': volume3,
                'resolution': 2.2,
                'inputStructure': structure5_PDB
                }
        protMolProbity = self.newProtocol(PhenixProtRunMolprobity, **args)
        protMolProbity.setObjLabel('MolProbity validation\n volume and pdb\n')
        self.launchProtocol(protMolProbity)

        # check MolProbity results
        self.checkResults(ramOutliers=0.12,
                          ramFavored=95.86,
                          rotOutliers=0.52,
                          cbetaOutliers=0,
                          clashScore=9.74,
                          overallScore=1.80,
                          protMolProbity=protMolProbity)

    def testMolProbityValidationManyOutliers1(self):
        """ This test checks that MolProbity validation protocol runs with
        two independent pdbs provided directly without inputVol associated and
        allows the comparison of results
         """
        print "Run MolProbity validation to compare an imported pdb " \
              "files obtained in another project"

        # import first PDB
        structure6_PDB = self._importStructureMolProbity1()

        # MolProbity
        args = {
                'inputStructure': structure6_PDB
                }
        protMolProbity = self.newProtocol(PhenixProtRunMolprobity, **args)
        protMolProbity.setObjLabel('MolProbity validation\n pdb\n')
        self.launchProtocol(protMolProbity)

        # check MolProbity results
        self.checkResults(ramOutliers=0.20,
                          ramFavored=97.35,
                          rotOutliers=12.24,
                          cbetaOutliers=0,
                          clashScore=130.72,
                          overallScore=3.53,
                          protMolProbity=protMolProbity)

    def testMolProbityValidationManyOutliers2(self):
        """ This test checks that MolProbity validation protocol runs with
        two independent pdbs provided directly without inputVol associated and
        allows the comparison of results
            """
        print "Run MolProbity validation to compare an imported pdb " \
                  "files obtained in another project"

        # import second PDB (with higher number of Outliers)
        structure7_PDB = self._importStructureMolProbity2()

        # MolProbity
        args = {
            'inputStructure': structure7_PDB
        }
        protMolProbity = self.newProtocol(PhenixProtRunMolprobity, **args)
        protMolProbity.setObjLabel('MolProbity validation\n pdb\n')
        self.launchProtocol(protMolProbity)

        # check MolProbity results
        self.checkResults(ramOutliers=3.82,
                          ramFavored=89.09,
                          rotOutliers=31.35,
                          cbetaOutliers=746,
                          clashScore=276.52,
                          overallScore=4.61,
                          protMolProbity=protMolProbity)

