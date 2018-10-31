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

# protocol to test the two ccp4 protocols of flexible fitting (coot and
# refmac)

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
        protImportPDB.setObjLabel('import pdb\nvolume associated\n1ake_start')
        self.launchProtocol(protImportPDB)
        structure2_PDB = protImportPDB.outputPdb
        self.assertTrue(structure2_PDB.getFileName())
        return structure2_PDB

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


class TestCootRefinement2(TestImportData):
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

    def testCootFlexibleFitFromVolAssocToPDB(self):

        # This test checks that coot runs when a volume is provided
        # associated to the input PDB
        print "Run Coot fit from imported pdb file and volume associated "

        # import PDB
        structure2_PDB = self._importStructurePDBWithVol()

        label = 'testLabel2'
        args = {'extraCommands': self._createExtraCommandLine(-24.11, -45.76,
                                                              -24.60, label),
                'pdbFileToBeRefined': structure2_PDB,
                'doInteractive': False
                }
        protCoot = self.newProtocol(CootRefine, **args)
        protCoot.setObjLabel('coot refinement\n pdb and associated volume\n '
                             'save model')

        self.launchProtocol(protCoot)
        self.assertIsNotNone(protCoot.testLabel2.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protCoot.testLabel2.getFileName()))
        self.assertTrue(
            os.path.exists(protCoot.output3DMap_0001.getFileName()))

    def testCootFlexibleFitFromVolAssocToCIF(self):

        # This test checks that coot runs when a volume is provided
        # associated to the input CIF
        print "Run Coot fit from imported cif file and volume associated "

        # import PDB
        structure_mmCIF = self._importStructuremmCIFWithVol()

        label = 'testLabel3'
        args = {'extraCommands': self._createExtraCommandLine(-24.11, -45.76,
                                                              -24.60, label),
                'pdbFileToBeRefined': structure_mmCIF,
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
        volume2 = self._importVolume2()

        # import PDB
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

    def testMultipleCootFit(self):
        # This test checks that coot runs three times when a volume is provided
        # associated to the input PDB file
        # starting volume with a different coordinate origin
        print "Run Coot fit from PDB file multiple times"

        volume2 = self._importVolume2()
        structure1_PDB = self._importStructurePDBWoVol()

        # first coot
        label = 'testLabel5'
        listVolCoot = [volume2]
        args = {'extraCommands': self._createExtraCommandLine(0., 0., 0.,
                                                              label),
                'inputVolumes': listVolCoot,
                'pdbFileToBeRefined': structure1_PDB,
                'doInteractive': True
                }
        protCoot = self.newProtocol(CootRefine, **args)
        protCoot.setObjLabel('coot refinement\n volume pdb\n '
                             'save model three times')

        try:
            self.launchProtocol(protCoot)
        except:
            print "first call to coot ended"
        self.assertIsNotNone(protCoot.testLabel5.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protCoot.testLabel5.getFileName()))
        self.assertTrue(
            os.path.exists(protCoot.output3DMap_0001.getFileName()))

        # second coot (label=None)
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

        # third coot
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


class TestRefmacRefinement2(TestImportData):
    """ Test the flexible fitting of refmac refinement protocol
    """

    def testRefmacFlexibleFitFromPDB(self):
        """ This test checks that refmac runs with an atomic structure;
         No Volume was provided and an error message is expected"""
        print "Run Refmac refinement from imported pdb file without" \
              "imported or pdb-associated volume"

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

    def testRefmacFlexibleFit1(self):
        """ This test checks that refmac runs with a volume provided
        directly as inputVol and an input PDB (refmac without mask)
         """
        print "Run Refmac refinement withouth mask from imported volume and " \
              "pdb file"

        # Import Volume
        volume = self._importVolume2()

        # import PDB
        structure_PDB = self._importStructurePDBWoVol()

        self.assertTrue(structure_PDB.getFileName())
        self.assertFalse(structure_PDB.getVolume())
        args = {'inputStructure': structure_PDB,
                'inputVolume': volume,
                'generateMaskedVolume': False
                }

        protRefmac = self.newProtocol(CCP4ProtRunRefmac, **args)
        protRefmac.setObjLabel('refmac refinement\n'
                               'pdb and volume\n save model')
        self.launchProtocol(protRefmac)
        self.assertIsNotNone(protRefmac.outputPdb.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protRefmac.outputPdb.getFileName()))

    def testRefmacFlexibleFit2(self):
        """ This test checks that refmac runs with a volume associated to an
        input CIF (refmac with mask)
         """
        print "Run MASK Refmac refinement withouth mask from associated " \
              "volume to a cif file"

        # import PDB
        structure_PDB = self._importStructuremmCIFWithVol2()

        self.assertTrue(structure_PDB.getFileName())
        self.assertTrue(structure_PDB.getVolume())
        args = {'inputStructure': structure_PDB,
                'generateMaskedVolume': True
                }

        protRefmac = self.newProtocol(CCP4ProtRunRefmac, **args)
        protRefmac.setObjLabel('MASK refmac refinement\n'
                               'pdb and volume\n save model')
        self.launchProtocol(protRefmac)
        self.assertIsNotNone(protRefmac.outputPdb.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protRefmac.outputPdb.getFileName()))

    def testRefmacFlexibleFitAfterCoot(self):
        """ This test checks that refmac runs with a volume provided
        directly as inputVol, the input PDB was fitted to the volume and
        refined previously by coot (refmac without mask)
         """
        print "Run Refmac refinement withouth mask from imported volume and " \
              "pdb file fitted and refined by Coot"

        # Import Volume
        volume = self._importVolume()

        # import PDB
        structure_PDB = self._importStructurePDBWoVol()

        # coot
        listVolCoot = [volume]
        label = 'testLabel1'
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

        # refmac without MASK
        coot_PDB = protCoot.testLabel1
        args = {'inputVolume': volume,
                'inputStructure': coot_PDB,
                'generateMaskedVolume': False
                }
        protRefmac = self.newProtocol(CCP4ProtRunRefmac, **args)
        protRefmac.setObjLabel('refmac refinement\nvolume and pdb\n'
                               'save model')
        self.launchProtocol(protRefmac)
        self.assertIsNotNone(protRefmac.outputPdb.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protRefmac.outputPdb.getFileName()))

    def testMASKRefmacFlexibleFitAfterCoot(self):
        """ This test checks that refmac runs with a volume provided
        directly as inputVol, the input PDB was fitted to the volume and
        refined previously by coot (refmac with mask)
         """
        print "Run MASK Refmac refinement withouth mask from imported " \
              "volume and pdb file fitted and refined by Coot"

        # Import Volume
        volume = self._importVolume()

        # import PDB
        structure_PDB = self._importStructurePDBWoVol()

        # coot
        listVolCoot = [volume]
        label = 'testLabel2'
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

        # refmac with MASK
        coot_PDB = protCoot.testLabel2
        args = {'inputVolume': volume,
                'inputStructure': coot_PDB,
                'generateMaskedVolume': True
                }
        protRefmac = self.newProtocol(CCP4ProtRunRefmac, **args)
        protRefmac.setObjLabel('MASK refmac refinement\nvolume and pdb\n'
                               'save model')
        self.launchProtocol(protRefmac)
        self.assertIsNotNone(protRefmac.outputPdb.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protRefmac.outputPdb.getFileName()))

    def testRefmacFlexibleFitAfterCoot2(self):
        """ This test checks that refmac runs with a volume provided
        directly as inputVol, the input PDB was fitted to the
        associated volume and
        refined previously by coot (refmac without mask)
         """
        print "Run Refmac refinement withouth mask from imported " \
              "pdb file (volume associated) fitted and refined by Coot"

        # import PDB
        structure_PDB = self._importStructurePDBWithVol2()

        # coot
        label = 'testLabel3'
        args = {'extraCommands': self._createExtraCommandLine(0., 0., 0.,
                                                              label),
                'pdbFileToBeRefined': structure_PDB,
                'doInteractive': False
                }
        protCoot = self.newProtocol(CootRefine, **args)
        protCoot.setObjLabel('coot refinement\n volume and unfitted pdb\n '
                             'save model')
        self.launchProtocol(protCoot)

        # refmac without MASK
        coot_PDB = protCoot.testLabel3
        volume = protCoot.pdbFileToBeRefined.get().getVolume()
        args = {'inputVolume': volume,
                'inputStructure': coot_PDB,
                'generateMaskedVolume': False
                }
        protRefmac = self.newProtocol(CCP4ProtRunRefmac, **args)
        protRefmac.setObjLabel('refmac refinement\nvolume and pdb\n'
                               'save model')
        self.launchProtocol(protRefmac)
        self.assertIsNotNone(protRefmac.outputPdb.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protRefmac.outputPdb.getFileName()))

    def testMASKRefmacFlexibleFitAfterCoot2(self):
        """ This test checks that refmac runs with a volume provided
        directly as inputVol, the input PDB was fitted to the volume and
        refined previously by coot (refmac with mask)
         """
        print "Run MASK Refmac refinement withouth mask from imported " \
              "volume and pdb file fitted and refined by Coot"

        # import PDB
        structure_PDB = self._importStructurePDBWithVol2()

        # coot
        label = 'testLabel4'
        args = {'extraCommands': self._createExtraCommandLine(0., 0., 0.,
                                                              label),
                'pdbFileToBeRefined': structure_PDB,
                'doInteractive': False
                }
        protCoot = self.newProtocol(CootRefine, **args)
        protCoot.setObjLabel('coot refinement\n volume and unfitted pdb\n '
                             'save model')
        self.launchProtocol(protCoot)

        # refmac with MASK
        coot_PDB = protCoot.testLabel4
        volume = protCoot.pdbFileToBeRefined.get().getVolume()
        args = {'inputVolume': volume,
                'inputStructure': coot_PDB,
                'generateMaskedVolume': False
                }
        protRefmac = self.newProtocol(CCP4ProtRunRefmac, **args)
        protRefmac.setObjLabel('MASK refmac refinement\nvolume and pdb\n'
                               'save model')
        self.launchProtocol(protRefmac)
        self.assertIsNotNone(protRefmac.outputPdb.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protRefmac.outputPdb.getFileName()))

    def testRefmacRefinementAfterMultipleCootFit(self):
        # This test checks that refmac runs when a volume and a pdb provided
        # by Coot (three runs) with or without mask
        print "Run Refmac refinement from PDB file saved from " \
              "Coot (3 runs) with or without mask"

        volume2 = self._importVolume2()
        structure1_PDB = self._importStructurePDBWoVol()

        # first coot
        label = 'testLabel5'
        listVolCoot = [volume2]
        args = {'extraCommands': self._createExtraCommandLine(0., 0., 0.,
                                                              label),
                'inputVolumes': listVolCoot,
                'pdbFileToBeRefined': structure1_PDB,
                'doInteractive': True
                }
        protCoot = self.newProtocol(CootRefine, **args)
        protCoot.setObjLabel('coot refinement\n volume pdb\n '
                             'save model three times')

        try:
            self.launchProtocol(protCoot)
        except:
            print "first call to coot ended"
        self.assertIsNotNone(protCoot.testLabel5.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(protCoot.testLabel5.getFileName()))
        self.assertTrue(
            os.path.exists(protCoot.output3DMap_0001.getFileName()))

        # second coot (label=None)
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

        # third coot
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

        # refmac without mask
        coot_PDB = protCoot.lastTestLabel
        args = {'inputVolume': volume2,
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

        # refmac with mask
        coot_PDB = protCoot.lastTestLabel
        args = {'inputVolume': volume2,
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
