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
from pyworkflow.em.protocol import ProtImportVolumes, ProtImportPdb
from pyworkflow.em.packages.locscale import ProtLocScale
from pyworkflow.utils import redStr, greenStr, magentaStr

class TestProtLocscale(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsXmipp = DataSet.getDataSet('xmipp_tutorial')

        #
        # Imports
        #
        new = cls.proj.newProtocol  # short notation
        launch = cls.proj.launchProtocol

        # Volumes
        print magentaStr("\n==> Importing data - Input data")
        p_imp_volume = new(ProtImportVolumes,
                           filesPath=cls.dsXmipp.getFile('volumes/'),
                           filesPattern='BPV_scale_filtered_windowed_64.vol',
                           samplingRate=12)
        launch(p_imp_volume, wait=True)
        cls.inputVol = p_imp_volume.outputVolume

        
    def testLocscale(self):
        """ Check that an ouput was generated and the condition is valid.
            In addition, returns the size of the set.
        """
        # IMPORT_FROM_ID = 0
        # IMPORT_OBJ = 1
        # IMPORT_FROM_FILES = 2 
        new = self.proj.newProtocol  # short notation
        launch = self.proj.launchProtocol

        p_locscale = new(ProtLocScale,
                         inputVolume=self.inputVol,
                         inputPdbData=ProtLocScale.IMPORT_FROM_ID,
                         pdbId='3j5p', #  FIXME: this pdbID is not for this volume!!!
                         patchSize=8)
        launch(p_locscale, wait=True)

        self.assertIsNotNone(p_locscale.outputVolume, "outputVolume is None")

