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

# protocol to test the phenix method real_spacerefine

from pyworkflow.em.protocol.protocol_import import ProtImportPdb, \
    ProtImportVolumes
from pyworkflow.em.packages.phenix.protocol_real_space_refine import \
    PhenixProtRunRSRefine
from pyworkflow.tests import *


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

    def _importVolHemoOrg(self):
        args = {'filesPath': self.dsModBuild.getFile(
            'volumes/emd_3488.map'),
            'samplingRate': 1.05,
            'setOrigCoord': True,
            'x': 0.0,
            'y': 0.0,
            'z': 0.0
        }
        protImportVol = self.newProtocol(ProtImportVolumes, **args)
        protImportVol.setObjLabel('import volume emd_3488.map\nset '
                                  'origin 0.0 0.0 0.0\n')
        self.launchProtocol(protImportVol)
        volume_hemo_orig = protImportVol.outputVolume
        return volume_hemo_orig

    def _importStructurePDBWoVol2(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/3i3e_fitted.pdb')
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import pdb\n 3i3e_fitted')
        self.launchProtocol(protImportPDB)
        structure5_PDB = protImportPDB.outputPdb
        return structure5_PDB

    def _importStructurePDBWithVol2(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/3i3e_fitted.pdb'),
                'inputVolume': self._importVolume3()
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import pdb\n 3i3e_fitted with\n associated '
                                  'volume')
        self.launchProtocol(protImportPDB)
        structure6_PDB = protImportPDB.outputPdb
        return structure6_PDB

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

    def _importStructHemoPDB(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/5ni1.pdb'),
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import pdb\n 5ni1.pdb')
        self.launchProtocol(protImportPDB)
        structure_hemo_pdb = protImportPDB.outputPdb
        return structure_hemo_pdb

    def _importStructHemoCIF(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/5ni1.cif'),
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import pdb\n 5ni1.cif')
        self.launchProtocol(protImportPDB)
        structure_hemo_cif = protImportPDB.outputPdb
        return structure_hemo_cif

class TestPhenixRSRefine(TestImportData):
    """ Test the protocol of MolProbity validation
    """
    def checkResults(self, ramOutliers, ramFavored, rotOutliers, cbetaOutliers,
                     clashScore, overallScore, protRSRefine, places=2):
        # method to check MolProbity statistic results of the Final Results
        # Table
        self.assertAlmostEqual(protRSRefine.ramachandranOutliers,
                               ramOutliers, places)
        self.assertAlmostEqual(protRSRefine.ramachandranFavored,
                               ramFavored, places)
        self.assertAlmostEqual(protRSRefine.rotamerOutliers,
                               rotOutliers, places)
        self.assertAlmostEqual(protRSRefine.cbetaOutliers,
                               cbetaOutliers, places)
        self.assertAlmostEqual(protRSRefine.clashscore,
                               clashScore, places)
        protRSRefine.overallScore = round(protRSRefine.overallScore, 2)
        self.assertAlmostEqual(protRSRefine.overallScore,
                               overallScore, places)
    #
    def testPhenixRSRefineFromPDB(self):
        """ This test checks that phenix real space refine protocol runs
        with an atomic structure; No Volume was provided and an error message
        is expected
        """
        print "Run phenix real space refine protocol from imported pdb file " \
              "without imported or pdb-associated volume"

        # import PDB
        structure_refmac3 = self._importStructRefmac3()
        self.assertTrue(structure_refmac3.getFileName())
        self.assertFalse(structure_refmac3.getVolume())

        # real_space_refine
        args = {
                'resolution': 3.5,
                'inputStructure': structure_refmac3
                }
        protRSRefine = self.newProtocol(PhenixProtRunRSRefine, **args)
        protRSRefine.setObjLabel('RSRefine\nrefmac3.pdb\n')
        try:
            self.launchProtocol(protRSRefine)
        except Exception as e:
            self.assertTrue(True)
            print "This test should return a error message as '" \
                  " Input volume cannot be EMPTY.\n"

            return
        self.assertTrue(False)

    def testPhenixRSRefineFromVolumeAndPDB1(self):
        """ This test checks that phenix real_space_refine protocol runs
        with a volume provided directly as inputVol, the input PDB was fitted
        to the volume and refined previously by coot and refmac withouth mask
        in another project
        """
        print "Run phenix real_space_refine from imported volume and pdb file " \
              "previously fitted and refined by Coot and Refmac without mask"

        # Import Volume
        volume_refmac3 = self._importVolRefmac3()

        # import PDB
        structure_refmac3 = self._importStructRefmac3()

        # real_space_refine
        args = {'inputVolume': volume_refmac3,
                'resolution': 3.5,
                'inputStructure': structure_refmac3
                }
        protRSRefine = self.newProtocol(PhenixProtRunRSRefine, **args)
        protRSRefine.setObjLabel('RSRefine\n refmac3.mrc and '
                                   'refmac3.pdb\n')
        self.launchProtocol(protRSRefine)

        # check real_space_refine results
        self.checkResults(ramOutliers=0.00,
                          ramFavored=96.23,
                          rotOutliers=0.00,
                          cbetaOutliers=0,
                          clashScore=1.79,
                          overallScore=1.19,
                          protRSRefine=protRSRefine)

    def testPhenixRSRefineFromVolumeAndPDB4(self):
        """ This test checks that phenix real_space_refine protocol runs
        with a volume provided directly as inputVol and the input PDB from
        data banks
        """
        print "Run phenix real_space_refine from imported volume and pdb file " \
              "from data banks (vol origin 0.0, 0.0, 0.0)"

        # Import Volume
        volume_hemo_org = self._importVolHemoOrg()

        # import PDB
        structure_hemo_pdb = self._importStructHemoPDB()

        # real_space_refine
        args = {'inputVolume': volume_hemo_org,
                'resolution': 3.2,
                'inputStructure': structure_hemo_pdb
                }
        protRSRefine = self.newProtocol(PhenixProtRunRSRefine, **args)
        protRSRefine.setObjLabel('RSRefine hemo\n emd_3488.map and '
                                 '5ni1.pdb\n')
        self.launchProtocol(protRSRefine)

        # check real_space_refine results
        self.checkResults(ramOutliers=0.00,
                          ramFavored=97.17,
                          rotOutliers=2.17,
                          cbetaOutliers=0,
                          clashScore=2.32,
                          overallScore=1.42,
                          protRSRefine=protRSRefine)

    def testPhenixRSRefineFromVolumeAndPDB5(self):
        """ This test checks that phenix real_space_refine protocol runs
        with a volume provided directly as inputVol and the input PDB from
        data banks
        """
        print "Run phenix real_space_refine from imported volume and cif file " \
              "from data banks (vol origin 0.0, 0.0, 0.0)"

        # Import Volume
        volume_hemo_org = self._importVolHemoOrg()

        # import PDB
        structure_hemo_cif = self._importStructHemoCIF()

        # real_space_refine
        args = {'inputVolume': volume_hemo_org,
                'resolution': 3.2,
                'inputStructure': structure_hemo_cif
                }
        protRSRefine = self.newProtocol(PhenixProtRunRSRefine, **args)
        protRSRefine.setObjLabel('RSRefine hemo\n emd_3488.map and '
                                 '5ni1.cif\n')
        self.launchProtocol(protRSRefine)

        # check real_space_refine results
        self.checkResults(ramOutliers=0.00,
                          ramFavored=97.17,
                          rotOutliers=2.17,
                          cbetaOutliers=0,
                          clashScore=2.32,
                          overallScore=1.42,
                          protRSRefine=protRSRefine)

    # def testPhenixRSRefineSeveralChains(self):
    #     """ This test checks that the phenix real_space_refine protocol runs
    #     with a volume provided directly as inputVol, the input PDB was fitted
    #     to the volume and refined previously by coot and refmac in another
    #     project
    #     """
    #     print "Run phenix real_space_refine protocol from imported volume " \
    #           "and pdb file already refined by Coot and Refmac in another " \
    #           "project"
    #
    #     # Import Volume
    #     volume3 = self._importVolume3()
    #
    #     # import PDB
    #     structure5_PDB = self._importStructurePDBWoVol2()
    #
    #     # real_space_refine
    #     args = {'inputVolume': volume3,
    #             'resolution': 2.2,
    #             'inputStructure': structure5_PDB
    #             }
    #     protRSRefine = self.newProtocol(PhenixProtRunRSRefine, **args)
    #     protRSRefine.setObjLabel('RSRefine\n volume and pdb\n')
    #     self.launchProtocol(protRSRefine)
    #
    #     # check real_space_refine results
    #     self.checkResults(ramOutliers=0.12,
    #                       ramFavored=95.86,
    #                       rotOutliers=0.52,
    #                       cbetaOutliers=0,
    #                       clashScore=9.74,
    #                       overallScore=1.80,
    #                       protRSRefine=protRSRefine)
    #
    # def testPhenixRSRefineSeveralChains2(self):
    #     """ This test checks that the phenix real_space_refine protocol runs
    #     with an input PDB that has a volume associated (the pdb was fitted
    #     to the volume and refined previously by coot and refmac in another
    #     project
    #     """
    #     print "Run phenix real_space_refine protocol from imported pdb " \
    #           "and and volume associated and already refined by Coot and " \
    #           "Refmac in another project"
    #
    #     # import PDB
    #     structure6_PDB = self._importStructurePDBWithVol2()
    #
    #     # real_space_refine
    #     args = {
    #             'resolution': 2.2,
    #             'inputStructure': structure6_PDB
    #             }
    #     protRSRefine = self.newProtocol(PhenixProtRunRSRefine, **args)
    #     protRSRefine.setObjLabel('RSRefine\n volume and pdb\n')
    #     self.launchProtocol(protRSRefine)
    #
    #     # check real_space_refine results
    #     self.checkResults(ramOutliers=0.12,
    #                       ramFavored=95.86,
    #                       rotOutliers=0.52,
    #                       cbetaOutliers=0,
    #                       clashScore=9.74,
    #                       overallScore=1.80,
    #                       protRSRefine=protRSRefine)

