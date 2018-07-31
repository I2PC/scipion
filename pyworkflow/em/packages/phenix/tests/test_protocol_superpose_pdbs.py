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

# protocol to test the phenix protocol superpose_pdbs

from pyworkflow.em.protocol.protocol_import import ProtImportPdb,\
    ProtImportVolumes
from pyworkflow.em.packages.phenix.protocol_superpose_pdbs import \
    PhenixProtRunSuperposePDBs
from pyworkflow.tests import *


class TestImportBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsModBuild = DataSet.getDataSet('model_building_tutorial')


class TestImportData(TestImportBase):
    """ Import atomic structures(PDBx/mmCIF files)
    """

    def _importVolume2(self):
        args = {'filesPath': self.dsModBuild.getFile('volumes/1ake_4-5A.mrc'),
                'samplingRate': 1.5,
                'setOrigCoord': True,
                'x': 11.994,
                'y': -7.881,
                'z': 10.91
                }
        protImportVol = self.newProtocol(ProtImportVolumes, **args)
        protImportVol.setObjLabel('import volume 1ake_4-5A\n'
                                  'set origin in 11 -7 10\n')
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

    def _importStructurePDBWithVol2(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/1ake_start.pdb'),
                'inputVolume': self._importVolume2()
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import pdb\nvolume associated\n1ake_start')
        self.launchProtocol(protImportPDB)
        structure3_PDB = protImportPDB.outputPdb
        self.assertTrue(structure3_PDB.getFileName())
        return structure3_PDB

    def _importStructuremmCIFWithVol2(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile('PDBx_mmCIF/'
                                                   '1ake_start.pdb.cif'),
                'inputVolume': self._importVolume2()
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import mmCIF\n volume associated\n '
                                  '1ake_start')
        self.launchProtocol(protImportPDB)
        structure3_mmCIF = protImportPDB.outputPdb
        self.assertTrue(structure3_mmCIF.getFileName())
        return structure3_mmCIF

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

class TestProtSuperposePdbs(TestImportData):
    """ Test the protocol of MolProbity validation
    """
    def checkResults(self, startRMSD, finalRMSD, protSuperposePdbs, places=2):
        # method to check start and final RMSD values showed in Summary
        logFile = os.path.abspath(protSuperposePdbs._getLogsPath()) +\
                  "/run.stdout"
        protSuperposePdbs._parseLogFile(logFile)
        self.assertAlmostEqual(protSuperposePdbs.startRMSD,
                               startRMSD, places)
        self.assertAlmostEqual(protSuperposePdbs.finalRMSD,
                               finalRMSD, places)

    def testSuperposePdbsFromPDBAndCIF(self):
        """ This test checks that phenix superpose_pdbs protocol runs with
        two atomic structures (pdb and cif)"""
        print "Run phenix superpose_pdbs protocol from an imported pdb file " \
              "and an imported cif file without volumes associated"

        # import PDBs
        structure1_PDB = self._importStructurePDBWoVol()
        self.assertTrue(structure1_PDB.getFileName())
        self.assertFalse(structure1_PDB.getVolume())

        structure1_mmCIF = self._importStructuremmCIFWoVol()
        self.assertTrue(structure1_mmCIF.getFileName())
        self.assertFalse(structure1_mmCIF.getVolume())

        args = {
                'inputStructureFixed': structure1_PDB,
                'inputStructureMoving': structure1_mmCIF
               }

        protSuperposePdbs = self.newProtocol(PhenixProtRunSuperposePDBs, **args)
        protSuperposePdbs.setObjLabel('SuperposePDBs\n'
                                   'no volumes associated to pdb/cif\n')
        self.launchProtocol(protSuperposePdbs)
        self.assertTrue(os.path.exists(
            protSuperposePdbs.outputPdb.getFileName()))

        # check RMSD results
        self.checkResults(startRMSD=0.000,
                          finalRMSD=0.000,
                          protSuperposePdbs=protSuperposePdbs)

    def testSuperposePdbsFromPDBAndPDB(self):
        """ This test checks that phenix superpose_pdbs protocol runs with
        two atomic structures (pdbs)"""
        print "Run phenix superpose_pdbs protocol from two imported pdb" \
              "files without volumes associated"

        # import PDBs
        structure6_PDB = self._importStructureMolProbity1()
        self.assertTrue(structure6_PDB.getFileName())
        self.assertFalse(structure6_PDB.getVolume())

        structure7_PDB = self._importStructureMolProbity2()
        self.assertTrue(structure7_PDB.getFileName())
        self.assertFalse(structure7_PDB.getVolume())

        args = {
                'inputStructureFixed': structure6_PDB,
                'inputStructureMoving': structure7_PDB
                }

        protSuperposePdbs = self.newProtocol(PhenixProtRunSuperposePDBs,
                                                 **args)
        protSuperposePdbs.setObjLabel('SuperposePDBs\n'
                                          'no volumes associated to pdbs\n')
        self.launchProtocol(protSuperposePdbs)
        self.assertTrue(os.path.exists(
            protSuperposePdbs.outputPdb.getFileName()))

        # check RMSD results
        self.checkResults(startRMSD=0.464,
                          finalRMSD=0.462,
                          protSuperposePdbs=protSuperposePdbs)

    def testSuperposePdbsFromPDBAndCIFWithVol(self):
        """ This test checks that phenix superpose_pdbs protocol runs with
        two atomic structures (pdb and cif) with associated volumes"""
        print "Run phenix superpose_pdbs protocol from an imported pdb file " \
              "and a imported cif file with volumes associated"

        # import PDBs
        structure2_PDB = self._importStructurePDBWithVol2()
        self.assertTrue(structure2_PDB.getFileName())
        self.assertTrue(structure2_PDB.getVolume())

        structure2_mmCIF = self._importStructuremmCIFWithVol2()
        self.assertTrue(structure2_mmCIF.getFileName())
        self.assertTrue(structure2_mmCIF.getVolume())

        args = {
                'inputStructureFixed': structure2_PDB,
                'inputStructureMoving': structure2_mmCIF
                }

        protSuperposePdbs = self.newProtocol(PhenixProtRunSuperposePDBs,
                                                 **args)
        protSuperposePdbs.setObjLabel('SuperposePDBs\n'
                                      'volumes associated to pdb/cif\n')
        self.launchProtocol(protSuperposePdbs)
        self.assertTrue(os.path.exists(
            protSuperposePdbs.outputPdb.getFileName()))
        self.assertTrue(os.path.exists(
            protSuperposePdbs.outputPdb.getVolume().getFileName()))

        # check RMSD results
        self.checkResults(startRMSD=0.000,
                          finalRMSD=0.000,
                          protSuperposePdbs=protSuperposePdbs)

    def testSuperposePdbsFromPDBWithVolAndPDBWoVol1(self):
        """ This test checks that phenix superpose_pdbs protocol runs with
        two atomic structures (pdb with associated volumes and pdb without
        associated volume"""
        print "Run phenix superpose_pdbs protocol from an imported pdb file " \
              "with associated volume and an imported pdb without volume " \
              "associated"

        # import PDBs
        structure2_PDB = self._importStructurePDBWithVol2()
        self.assertTrue(structure2_PDB.getFileName())
        self.assertTrue(structure2_PDB.getVolume())

        structure3_PDB = self._importMut1StructurePDBWoVol()
        self.assertTrue(structure3_PDB.getFileName())
        self.assertFalse(structure3_PDB.getVolume())

        args = {
                'inputStructureFixed': structure2_PDB,
                'inputStructureMoving': structure3_PDB
                }

        protSuperposePdbs = self.newProtocol(PhenixProtRunSuperposePDBs,
                                                 **args)
        protSuperposePdbs.setObjLabel('SuperposePDBs\n'
                                          'volume associated to one pdb\n')
        self.launchProtocol(protSuperposePdbs)
        self.assertTrue(os.path.exists(
            protSuperposePdbs.outputPdb.getFileName()))
        self.assertTrue(os.path.exists(
            protSuperposePdbs.outputPdb.getVolume().getFileName()))

        # check RMSD results
        self.checkResults(startRMSD=0.000,
                          finalRMSD=0.000,
                          protSuperposePdbs=protSuperposePdbs)

    def testSuperposePdbsFromPDBWithVolAndPDBWoVol2(self):
        """ This test checks that phenix superpose_pdbs protocol runs with
        two atomic structures (pdb with associated volumes and pdb without
        associated volume"""
        print "Run phenix superpose_pdbs protocol from an imported pdb file " \
              "with associated volume and an imported pdb without volume " \
              "associated"

        # import PDBs
        structure2_PDB = self._importStructurePDBWithVol2()
        self.assertTrue(structure2_PDB.getFileName())
        self.assertTrue(structure2_PDB.getVolume())

        structure3_PDB = self._importMut2StructurePDBWoVol()
        self.assertTrue(structure3_PDB.getFileName())
        self.assertFalse(structure3_PDB.getVolume())

        args = {
                'inputStructureFixed': structure2_PDB,
                'inputStructureMoving': structure3_PDB
                }

        protSuperposePdbs = self.newProtocol(PhenixProtRunSuperposePDBs,
                                                 **args)
        protSuperposePdbs.setObjLabel('SuperposePDBs\n'
                                          'volume associated to one pdb\n')
        self.launchProtocol(protSuperposePdbs)
        self.assertTrue(os.path.exists(
            protSuperposePdbs.outputPdb.getFileName()))
        self.assertTrue(os.path.exists(
            protSuperposePdbs.outputPdb.getVolume().getFileName()))

        # check RMSD results
        self.checkResults(startRMSD=0.000,
                          finalRMSD=0.000,
                          protSuperposePdbs=protSuperposePdbs)

    def testSuperposePdbsFromPDBWithVolAndPDBWoVol3(self):
        """ This test checks that phenix superpose_pdbs protocol runs with
        two atomic structures (pdb with associated volumes and pdb without
        associated volume"""
        print "Run phenix superpose_pdbs protocol from an imported pdb file " \
              "with associated volume and an imported pdb without volume " \
              "associated"

        # import PDBs
        structure2_PDB = self._importStructurePDBWithVol2()
        self.assertTrue(structure2_PDB.getFileName())
        self.assertTrue(structure2_PDB.getVolume())

        structure6_PDB = self._importStructureMolProbity1()
        self.assertTrue(structure6_PDB.getFileName())
        self.assertFalse(structure6_PDB.getVolume())

        args = {
                'inputStructureFixed': structure2_PDB,
                'inputStructureMoving': structure6_PDB
                }

        protSuperposePdbs = self.newProtocol(PhenixProtRunSuperposePDBs,
                                                 **args)
        protSuperposePdbs.setObjLabel('SuperposePDBs\n'
                                      'volume associated to one pdb\n'
                                      'unrelated pdbs')
        self.launchProtocol(protSuperposePdbs)
        self.assertTrue(os.path.exists(
            protSuperposePdbs.outputPdb.getFileName()))
        self.assertTrue(os.path.exists(
            protSuperposePdbs.outputPdb.getVolume().getFileName()))

        # check RMSD results
        self.checkResults(startRMSD=60.703,
                          finalRMSD=19.500,
                          protSuperposePdbs=protSuperposePdbs)
