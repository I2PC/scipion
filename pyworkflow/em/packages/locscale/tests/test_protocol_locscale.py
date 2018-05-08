# ***************************************************************************
# *
# * Authors:     David Maluenda (dmaluenda@cnb.csic.es) (2018)
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
from pyworkflow.em.protocol import ProtImportVolumes, ProtImportMask
from pyworkflow.em.packages.locscale import ProtLocScale
from pyworkflow.utils import magentaStr

class TestProtLocscale(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataSet = DataSet.getDataSet('xmipp_tutorial')

        #
        # Imports
        #
        print magentaStr("\n==> Importing data - Input data")
        new = cls.proj.newProtocol  # short notation
        launch = cls.proj.launchProtocol

        # Volume
        print magentaStr("\nImporting Volume:")
        pImpVolume = new(ProtImportVolumes, samplingRate=1,
                         filesPath=cls.dataSet.getFile('vol2'))
        launch(pImpVolume, wait=True)
        cls.inputVol = pImpVolume.outputVolume

        # Reference
        print magentaStr("\nImporting Reference:")
        pImpRef = new(ProtImportVolumes, samplingRate=1,
                      filesPath=cls.dataSet.getFile('vol3'))
        launch(pImpRef, wait=True)
        cls.inputRef = pImpRef.outputVolume

        # Mask
        print magentaStr("\nImporting Mask:")
        pImpMask = new(ProtImportMask,
                       maskPath=cls.dataSet.getFile('mask3d'),
                       samplingRate=1)
        launch(pImpMask, wait=True)
        cls.mask = pImpMask.outputMask

        
    def testLocscale(self):
        """ Check that an output was generated and the condition is valid.
            In addition, returns the size of the set.
        """
        print magentaStr("\n==> Testing locscale:")
        def launchTest(label, mpi=4, mask=None):
            print magentaStr("\nTest %s:" % label)
            pLocScale = self.proj.newProtocol(ProtLocScale,
                                              objLabel='locscale - ' + label,
                                              inputVolume=self.inputVol,
                                              refObj=self.inputRef,
                                              patchSize=16,
                                              binaryMask=mask,
                                              numberOfMpi=mpi)
            self.proj.launchProtocol(pLocScale, wait=True)

            self.assertIsNotNone(pLocScale.outputVolume,
                                 "outputVolume is None for %s test." % label)

            self.assertEqual(self.inputVol.getDim(),
                             pLocScale.outputVolume.getDim(),
                             "outputVolume has diferent size than inputVol "
                             "for %s test" % label)

            self.assertEqual(self.inputVol.getSamplingRate(),
                             pLocScale.outputVolume.getSamplingRate(),
                             "outputVolume has diferent sampling rate than "
                             "inputVol for %s test" % label)

        # default test
        launchTest('with MPI + noMask')

        # with mask test
        launchTest('with MPI + Mask', mask=self.mask)

        # with mask test
        launchTest('with Mask + noMPI', mask=self.mask, mpi=1)

        # without MPI
        launchTest('noMask + noMPI', mpi=1)
