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

# protocol to test the validation method EMRinger

from pyworkflow.em.protocol.protocol_import import ProtImportPdb, \
    ProtImportVolumes
from pyworkflow.em.packages.phenix.protocol_emringer import \
    PhenixProtRunEMRinger
from pyworkflow.tests import *
import json


class TestImportBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsModBuild = DataSet.getDataSet('model_building_tutorial')


class TestImportData(TestImportBase):
    """ Import map volumes and atomic structures(PDBx/mmCIF files)
    """

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

    def _importVolCoot1(self):
        args = {'filesPath': self.dsModBuild.getFile(
            'volumes/coot1.mrc'),
            'samplingRate': 1.5,
            'setOrigCoord': False
            }
        protImportVol = self.newProtocol(ProtImportVolumes, **args)
        protImportVol.setObjLabel('import volume coot1.mrc\set'
                                  'origin default\n')
        self.launchProtocol(protImportVol)
        volume_coot1 = protImportVol.outputVolume
        return volume_coot1

    def _importVolCoot2(self):
        args = {'filesPath': self.dsModBuild.getFile(
            'volumes/coot2.mrc'),
            'samplingRate': 1.5,
            'setOrigCoord': True,
            'x': 12.0,
            'y': -7.5,
            'z': 10.5
            }
        protImportVol = self.newProtocol(ProtImportVolumes, **args)
        protImportVol.setObjLabel('import volume coot2.mrc\nset '
                                  'origin 12.0 -7.5 10.5\n')
        self.launchProtocol(protImportVol)
        volume_coot2 = protImportVol.outputVolume
        return volume_coot2

    def _importVolRefmac1(self):
        args = {'filesPath': self.dsModBuild.getFile(
            'volumes/refmac1.mrc'),
            'samplingRate': 1.5,
            'setOrigCoord': False
            }
        protImportVol = self.newProtocol(ProtImportVolumes, **args)
        protImportVol.setObjLabel('import volume refmac1.mrc\nset '
                                  'origin default\n')
        self.launchProtocol(protImportVol)
        volume_refmac1 = protImportVol.outputVolume
        return volume_refmac1

    def _importVolRefmac2(self):
        args = {'filesPath': self.dsModBuild.getFile(
            'volumes/refmac2.mrc'),
            'samplingRate': 1.5,
            'setOrigCoord': True,
            'x': 12.0,
            'y': -7.5,
            'z': 10.5
        }
        protImportVol = self.newProtocol(ProtImportVolumes, **args)
        protImportVol.setObjLabel('import volume refmac2.mrc\nset '
                                  'origin 12.0 -7.5 10.5\n')
        self.launchProtocol(protImportVol)
        volume_refmac2 = protImportVol.outputVolume
        return volume_refmac2

    def _importVolRefmac3(self):
        args = {'filesPath': self.dsModBuild.getFile(
            'volumes/refmac3.mrc'),
            'samplingRate': 1.5,
            'setOrigCoord': True,
            'x': 37.5,
            'y': 37.5,
            'z': 37.5
        }
        protImportVol = self.newProtocol(ProtImportVolumes, **args)
        protImportVol.setObjLabel('import volume refmac3.mrc\nset '
                                  'origin 37.5 37.5 37.5\n')
        self.launchProtocol(protImportVol)
        volume_refmac3 = protImportVol.outputVolume
        return volume_refmac3

    def _importVolEMRinger1(self):
        args = {'filesPath': self.dsModBuild.getFile(
                'volumes/emd_5995.map'),
                'samplingRate': 0.637,
                'setOrigCoord': False
                }
        protImportVol = self.newProtocol(ProtImportVolumes, **args)
        protImportVol.setObjLabel('import volume emd_5995.map\nset '
                                  'origin default\n')
        self.launchProtocol(protImportVol)
        volume_emringer1 = protImportVol.outputVolume
        return volume_emringer1

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

    def _importStructurePDBWithVol(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/coot1.pdb'),
                'inputVolume': self._importVolCoot1()
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
                'inputVolume': self._importVolCoot2()
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import mmCIF\n volume associated\n '
                                  '1ake_start')
        self.launchProtocol(protImportPDB)
        structure2_mmCIF = protImportPDB.outputPdb
        self.assertTrue(structure2_mmCIF.getFileName())
        return structure2_mmCIF

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

    def _importStructCoot1(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/coot1.pdb'),
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import pdb\n coot1.pdb')
        self.launchProtocol(protImportPDB)
        structure_coot1 = protImportPDB.outputPdb
        return structure_coot1

    def _importStructCoot2(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/coot2.pdb'),
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import pdb\n coot2.pdb')
        self.launchProtocol(protImportPDB)
        structure_coot2 = protImportPDB.outputPdb
        return structure_coot2

    def _importStructRefmac1(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/refmac1.pdb'),
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import pdb\n refmac1.pdb')
        self.launchProtocol(protImportPDB)
        structure_refmac1 = protImportPDB.outputPdb
        return structure_refmac1

    def _importStructRefmac2(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/refmac2.pdb'),
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import pdb\n refmac2.pdb')
        self.launchProtocol(protImportPDB)
        structure_refmac2 = protImportPDB.outputPdb
        return structure_refmac2

    def _importStructRefmac3(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/refmac3.pdb'),
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import pdb\n refmac3.pdb')
        self.launchProtocol(protImportPDB)
        structure_refmac3 = protImportPDB.outputPdb
        return structure_refmac3

class TestEMRingerValidation2(TestImportData):
    """ Test the protocol of EMRinger validation
    """
    def checkResults(self, optThresh, rotRatio, maxZscore, modLength,
                     EMScore, protEMRinger, places=4):
        # method to check EMRinger statistic results of the Final Results Table
        textFileName = protEMRinger._getExtraPath(
            protEMRinger.EMRINGERTRANSFERFILENAME.replace('py', 'txt'))
        with open(textFileName, "r") as f:
            self.resultsDict = json.loads(str(f.read()))
            self.assertAlmostEqual(self.resultsDict[
                                'Optimal Threshold'], optThresh, places)
            self.assertAlmostEqual(self.resultsDict[
                                'Rotamer-Ratio'], rotRatio, places)
            self.assertAlmostEqual(self.resultsDict[
                                'Max Zscore'], maxZscore, places)
            self.assertAlmostEqual(self.resultsDict[
                                'Model Length'], modLength, places)
            self.assertAlmostEqual(self.resultsDict[
                                'EMRinger Score'], EMScore, places)

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
        args = {'inputStructure': structure_PDB,
                'doTest': True
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

    def testEMRingerValidationFromVolAssociatedToPDB(self):
        """ This test checks that EMRinger validation protocol runs with an
        atomic structure; No Volume was provided and an error message is
        expected"""
        print "Run EMRinger validation protocol from imported pdb file " \
              "with pdb-associated volume"

        # import PDB
        structure2_PDB = self._importStructurePDBWithVol()
        self.assertTrue(structure2_PDB.getFileName())
        self.assertTrue(structure2_PDB.getVolume())

        # EMRinger
        args = {'inputStructure': structure2_PDB,
                'doTest': True
                }
        protEMRinger = self.newProtocol(PhenixProtRunEMRinger, **args)
        protEMRinger.setObjLabel('EMRinger validation\n volume associated to '
                                 'pdb\n')
        self.launchProtocol(protEMRinger)

        # check EMRinger results
        self.checkResults(optThresh=0.5988849483678802,
                          rotRatio=0.9166666666666666,
                          maxZscore=2.607144567760445,
                          modLength=121,
                          EMScore=2.370131425236768,
                          protEMRinger=protEMRinger)

    def testEMRingerValidationFromVolAssociatedToCIF(self):
        """ This test checks that EMRinger validation protocol runs with an
        atomic structureand associated
        """
        print "Run EMRinger validation protocol from imported cif file " \
              "with cif-associated volume"

        # import PDB
        structure2_mmCIF = self._importStructuremmCIFWithVol()
        self.assertTrue(structure2_mmCIF.getFileName())
        self.assertTrue(structure2_mmCIF.getVolume())

        # EMRinger
        args = {'inputStructure': structure2_mmCIF,
                'doTest': True
                }
        protEMRinger = self.newProtocol(PhenixProtRunEMRinger, **args)
        protEMRinger.setObjLabel('EMRinger validation\n volume associated to '
                                 'pdb\n')
        self.launchProtocol(protEMRinger)

        # check EMRinger results
        self.checkResults(optThresh=0.46013039428861446,
                          rotRatio=0.6944444444444444,
                          maxZscore=2.6017745423519636,
                          modLength=121,
                          EMScore=2.3652495839563303,
                          protEMRinger=protEMRinger)

    def testEMRingerValidationFromVolumeAndPDB1(self):
        """ This test checks that EMRinger validation protocol runs with a
        volume provided directly as inputVol, the input PDB was fitted to
        the volume and refined previously by coot
         """
        print "Run EMRinger validation from imported volume and pdb file " \
              "previously fitted and refined by Coot"

        # Import Volume
        volume_coot1 = self._importVolCoot1()

        # import PDB
        structure_coot1 = self._importStructCoot1()

        # EMRinger
        args = {'inputVolume': volume_coot1,
                'inputStructure': structure_coot1,
                'doTest': True
                }
        protEMRinger = self.newProtocol(PhenixProtRunEMRinger, **args)
        protEMRinger.setObjLabel('EMRinger validation\n coot1.mrc coot1.pdb\n')
        self.launchProtocol(protEMRinger)

        # check EMRinger results
        self.checkResults(optThresh=0.5988849483678802,
                          rotRatio=0.9166666666666666,
                          maxZscore=2.607144567760445,
                          modLength=121,
                          EMScore=2.370131425236768,
                          protEMRinger=protEMRinger)

    def testEMRingerValidationFromVolumeAndPDB2(self):
        """ This test checks that EMRinger validation protocol runs with a
        volume provided directly as inputVol, the input PDB was fitted to
        the volume and refined previously by coot
         """
        print "Run EMRinger validation from imported volume and pdb file " \
              "previously fitted and refined by Coot"

        # Import Volume
        volume_coot2 = self._importVolCoot2()

        # import PDB
        structure_coot2 = self._importStructCoot2()

        # EMRinger
        args = {'inputVolume': volume_coot2,
                'inputStructure': structure_coot2,
                'doTest': True
                }
        protEMRinger = self.newProtocol(PhenixProtRunEMRinger, **args)
        protEMRinger.setObjLabel('EMRinger validation\n coot2.mrc coot2.pdb\n')
        self.launchProtocol(protEMRinger)

        # check EMRinger results
        self.checkResults(optThresh=0.451299768696973,
                          rotRatio=0.6756756756756757,
                          maxZscore=2.313625586144688,
                          modLength=121,
                          EMScore=2.103295987404262,
                          protEMRinger=protEMRinger)

    def testEMRingerValidationFromVolumeAndPDB3(self):
        """ This test checks that EMRinger validation protocol runs with a
        volume provided directly as inputVol, the input PDB was fitted to
        the volume and refined previously by coot and refmac without mask
         """
        print "Run EMRinger validation from imported volume and pdb file " \
              "previously fitted and refined by Coot and Refmac without mask"

        # Import Volume
        volume_refmac1 = self._importVolRefmac1()

        # import PDB
        structure_refmac1 = self._importStructRefmac1()

        # EMRinger
        args = {'inputVolume': volume_refmac1,
                'inputStructure': structure_refmac1,
                'doTest': True
                }
        protEMRinger = self.newProtocol(PhenixProtRunEMRinger, **args)
        protEMRinger.setObjLabel('EMRinger validation\n refmac1.mrc '
                                 'refmac1.pdb\n')
        self.launchProtocol(protEMRinger)

        # check EMRinger results
        self.checkResults(optThresh=0.6095088692601145,
                          rotRatio=0.826530612244898,
                          maxZscore=5.659704361609791,
                          modLength=121,
                          EMScore=5.145185783281628,
                          protEMRinger=protEMRinger)

    def testEMRingerValidationFromVolumeAndPDB4(self):
        """ This test checks that EMRinger validation protocol runs with a
        volume provided directly as inputVol, the input PDB was fitted to
        the volume and refined previously by coot and refmac without mask
         """
        print "Run EMRinger validation from imported volume and pdb file " \
              "previously fitted and refined by Coot and Refmac without mask"

        # Import Volume
        volume_refmac2 = self._importVolRefmac2()

        # import PDB
        structure_refmac2 = self._importStructRefmac2()

        # EMRinger
        args = {'inputVolume': volume_refmac2,
                'inputStructure': structure_refmac2,
                'doTest': True
                }
        protEMRinger = self.newProtocol(PhenixProtRunEMRinger, **args)
        protEMRinger.setObjLabel('EMRinger validation\n refmac2.mrc '
                                 'refmac2.pdb\n')
        self.launchProtocol(protEMRinger)

        # check EMRinger results
        self.checkResults(optThresh=0.6134629219492418,
                          rotRatio=0.8247422680412371,
                          maxZscore=5.59540502751673,
                          modLength=121,
                          EMScore=5.086731843197027,
                          protEMRinger=protEMRinger)

    def testEMRingerValidationFromVolumeAndPDB5(self):
        """ This test checks that EMRinger validation protocol runs with a
        volume provided directly as inputVol, the input PDB was fitted to
        the volume and refined previously by coot and refmac without mask
         """
        print "Run EMRinger validation from imported volume and pdb file " \
              "previously fitted and refined by Coot and Refmac without mask"

        # Import Volume
        volume_refmac3 = self._importVolRefmac3()

        # import PDB
        structure_refmac3 = self._importStructRefmac3()

        # EMRinger
        args = {'inputVolume': volume_refmac3,
                'inputStructure': structure_refmac3,
                'doTest': True
                }
        protEMRinger = self.newProtocol(PhenixProtRunEMRinger, **args)
        protEMRinger.setObjLabel('EMRinger validation\n refmac3.mrc '
                                 'refmac3.pdb\n')
        self.launchProtocol(protEMRinger)

        # check EMRinger results
        self.checkResults(optThresh=0.609911277932157,
                          rotRatio=0.8350515463917526,
                          maxZscore=5.799183055832969,
                          modLength=121,
                          EMScore=5.27198459621179,
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

    def testEMRingerValidationManyOutliers1(self):
        """ This test checks that EMRinger validation protocol runs with
        a pdb and inputVol associated
         """
        print "Run EMRinger validation to compare an imported pdb " \
              "file obtained in another project and volume"

        # Import Volume
        volume_emringer1 = self._importVolEMRinger1()

        # import first PDB
        structure6_PDB = self._importStructureMolProbity1()

        # EMRinger
        args = {'inputVolume': volume_emringer1,
                'inputStructure': structure6_PDB,
                'doTest': True
                }
        protEMRinger = self.newProtocol(PhenixProtRunEMRinger, **args)
        protEMRinger.setObjLabel('EMRinger validation\n emd_5995.map and '
                                 'pdb\n')
        self.launchProtocol(protEMRinger)

        # check EMRinger results
        self.checkResults(optThresh=0.037343173353483194,
                          rotRatio=0.5714285714285714,
                          maxZscore=0.1580348853102536,
                          modLength=2616,
                          EMScore=0.03089826515461954,
                          protEMRinger=protEMRinger)

    def testEMRingerValidationManyOutliers2(self):
        """ This test checks that EMRinger validation protocol runs with
        a pdb and inputVol associated
         """
        print "Run EMRinger validation to compare an imported pdb " \
              "file obtained in another project and volume"

        # Import Volume
        volume_emringer1 = self._importVolEMRinger1()

        # import first PDB
        structure7_PDB = self._importStructureMolProbity2()

        # EMRinger
        args = {'inputVolume': volume_emringer1,
                'inputStructure': structure7_PDB,
                'doTest': True
                }
        protEMRinger = self.newProtocol(PhenixProtRunEMRinger, **args)
        protEMRinger.setObjLabel('EMRinger validation\n emd_5995.map and '
                                 'pdb\n')
        self.launchProtocol(protEMRinger)

        # check EMRinger results
        self.checkResults(optThresh=0.027276289739258344,
                          rotRatio=0.5809199318568995,
                          maxZscore=1.9087017071890435,
                          modLength=2616,
                          EMScore=0.37318071471385195,
                          protEMRinger=protEMRinger)

