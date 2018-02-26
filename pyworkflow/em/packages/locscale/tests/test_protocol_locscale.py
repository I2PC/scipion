# ***************************************************************************
# * Authors:     David Maluenda (dmaluenda@cnb.csic.es) (2018)
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

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.em.protocol import ProtImportVolumes, ProtImportPdb, \
                                                                  ProtImportMask
from pyworkflow.em.packages.locscale import ProtLocScale
from pyworkflow.em.packages.xmipp3.protocol_preprocess import \
                                                      XmippProtCropResizeVolumes
from pyworkflow.utils import redStr, greenStr, magentaStr

class TestProtLocscale(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataSet = DataSet.getDataSet('nma')

        #
        # Imports
        #
        print magentaStr("\n==> Importing data - Input data")
        new = cls.proj.newProtocol  # short notation
        launch = cls.proj.launchProtocol

        # Volume
        print magentaStr("\nImporting Volume:")
        pImpVolume = new(ProtImportVolumes, samplingRate=1,
                         filesPath=cls.dataSet.getFile('vol'))
        launch(pImpVolume, wait=True)
        pResizeVol = new(XmippProtCropResizeVolumes,
                         inputVolumes=pImpVolume.outputVolume,
                         doResize=True,
                         resizeOption=1,
                         resizeDim=32)
        launch(pResizeVol, wait=True)
        cls.inputVol = pResizeVol.outputVol


        # PDB
        print magentaStr("\nImporting PDB:")
        pImpPdb = new(ProtImportPdb,
                      inputPdbData=ProtImportPdb.IMPORT_FROM_FILES,
                      pdbFile=cls.dataSet.getFile('pdb'))
        launch(pImpPdb, wait=True)
        cls.inputPdb = pImpPdb.outputPdb

        # # Mask  ---   FIXME : vol is not a mask...    ---
        print magentaStr("\nImporting Mask:")
        pImpMask = new(ProtImportMask,
                       maskPath=cls.dataSet.getFile('vol'),
                       samplingRate=1)
        launch(pImpMask, wait=True)
        pResizeMask = new(XmippProtCropResizeVolumes,
                          inputVolumes=pImpMask.outputMask,
                          doResize=True,
                          resizeOption=1,
                          resizeDim=32)
        launch(pResizeMask, wait=True)
        cls.mask = pResizeMask.outputVol

        
    def testLocscaleSimple(self):
        """ Check that an output was generated and the condition is valid.
            In addition, returns the size of the set.
        """
        def launchTest(pdbFrom, label, pdbId=None, refObj=None, inputPath=None):
            print magentaStr("\n%s test:" % label)
            pLocScale = self.proj.newProtocol(ProtLocScale,
                                              objLabel='locscale - ' + label,
                                              inputVolume=self.inputVol,
                                              inputPdbData=pdbFrom,
                                              pdbId=pdbId,
                                              refObj=refObj,
                                              inputPath=inputPath,
                                              patchSize=16)
            self.proj.launchProtocol(pLocScale, wait=True)
            self.assertIsNotNone(pLocScale.outputVolume,
                                 "outputVolume is None for %s test." % label)
            # return pLocScale.outputVolume

        # pdb IMPORT_FROM_ID   -  FIXME: this pdbID is not for this volume!!!
        launchTest(ProtLocScale.IMPORT_FROM_ID, 'pdbFromID', pdbId='3j5p')

        # pdb IMPORT_FROM_OBJ
        launchTest(ProtLocScale.IMPORT_FROM_OBJ, 'pdbFromObject',
                                     refObj=self.inputPdb)

        # ref IMPORT_FROM_OBJ
        launchTest(ProtLocScale.IMPORT_FROM_OBJ, 'refFromObject',
                                     refObj=self.inputVol)

        # pdb IMPORT_FROM_FILE
        launchTest(ProtLocScale.IMPORT_FROM_FILES, 'pdbFromFile',
                                   inputPath=self.dataSet.getFile('pdb'))


    def testLocscaleMask(self):
        """ Check that an output was generated and the condition is valid.
            In addition, returns the size of the set.
        """
        def launchTest(pdbFrom, label, pdbId=None, refObj=None, inputPath=None):
            print magentaStr("\n%s test with mask:" % label)
            pLocScale = self.proj.newProtocol(ProtLocScale,
                                              objLabel='locscale Mask - '+label,
                                              inputVolume=self.inputVol,
                                              inputPdbData=pdbFrom,
                                              pdbId=pdbId,
                                              refObj=refObj,
                                              inputPath=inputPath,
                                              patchSize=16,
                                              binaryMask=self.mask)
            self.proj.launchProtocol(pLocScale, wait=True)
            self.assertIsNotNone(pLocScale.outputVolume,
                                 "outputVolume is None for %s test." % label)
            # return pLocScale.outputVolume

        # pdb IMPORT_FROM_ID   -  FIXME: this pdbID is not for this volume!!!
        launchTest(ProtLocScale.IMPORT_FROM_ID, 'pdbFromID', pdbId='3j5p')

        # pdb IMPORT_FROM_OBJ
        launchTest(ProtLocScale.IMPORT_FROM_OBJ, 'pdbFromObject',
                   refObj=self.inputPdb)

        # ref IMPORT_FROM_OBJ
        launchTest(ProtLocScale.IMPORT_FROM_OBJ, 'refFromObject',
                   refObj=self.inputVol)

        # pdb IMPORT_FROM_FILE
        launchTest(ProtLocScale.IMPORT_FROM_FILES, 'pdbFromFile',
                   inputPath=self.dataSet.getFile('pdb'))


    def testLocscaleMPI(self):
        """ Check that an output was generated and the condition is valid.
            In addition, returns the size of the set.
        """
        def launchTest(pdbFrom, label, pdbId=None, refObj=None, inputPath=None):
            print magentaStr("\n%s test in MPI:" % label)
            pLocScale = self.proj.newProtocol(ProtLocScale,
                                              objLabel='locscale MPI - '+label,
                                              inputVolume=self.inputVol,
                                              inputPdbData=pdbFrom,
                                              pdbId=pdbId,
                                              refObj=refObj,
                                              inputPath=inputPath,
                                              patchSize=16,
                                              doParalelize=True)
            self.proj.launchProtocol(pLocScale, wait=True)
            self.assertIsNotNone(pLocScale.outputVolume,
                                 "outputVolume is None for %s test." % label)
            # return pLocScale.outputVolume

        # pdb IMPORT_FROM_ID   -  FIXME: this pdbID is not for this volume!!!
        launchTest(ProtLocScale.IMPORT_FROM_ID, 'pdbFromID', pdbId='3j5p')

        # pdb IMPORT_FROM_OBJ
        launchTest(ProtLocScale.IMPORT_FROM_OBJ, 'pdbFromObject', 
                   refObj=self.inputPdb)

        # ref IMPORT_FROM_OBJ
        launchTest(ProtLocScale.IMPORT_FROM_OBJ, 'refFromObject',
                   refObj=self.inputVol)

        # pdb IMPORT_FROM_FILE
        launchTest(ProtLocScale.IMPORT_FROM_FILES, 'pdbFromFile',
                   inputPath=self.dataSet.getFile('pdb'))
