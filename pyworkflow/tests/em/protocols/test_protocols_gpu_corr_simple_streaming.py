# ***************************************************************************
# *
# * Authors:     Amaya Jimenez (ajimenez@cnb.csic.es)
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

import time
import os

from pyworkflow.tests import BaseTest, setupTestProject
from pyworkflow.em.protocol import ProtCreateStreamData
from pyworkflow.em.protocol.protocol_create_stream_data import \
    SET_OF_PARTICLES
from pyworkflow.em.packages.grigoriefflab import ProtCTFFind
from pyworkflow.protocol import getProtocolFromDb
from pyworkflow.em.packages.xmipp3 import XmippProtCTFSelection
from pyworkflow.em.data import SetOfCTF, SetOfMicrographs
from pyworkflow.em.protocol import ProtImportAverages


# Load the number of movies for the simulation, by default equal 5, but
# can be modified in the environment
NUM_PART = 100

CTF_SQLITE = "ctfs.sqlite"
MIC_SQLITE = "micrographs.sqlite"
MIC_DISCARDED_SQLITE = "micrographsDiscarded.sqlite"


class TestGpuCorrSimpleStreaming(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    def _updateProtocol(self, prot):
        prot2 = getProtocolFromDb(prot.getProject().path,
                                  prot.getDbPath(),
                                  prot.getObjId())
        # Close DB connections
        prot2.getProject().closeMapper()
        prot2.closeMappers()
        return prot2

    def test_pattern(self):
        """ Import several Particles from a given pattern.
        """

        args = {'importFrom': ProtImportAverages.IMPORT_FROM_FILES,
                'filesPath': self.dsXmipp.getFile('/home/ajimenez/Desktop/averages/'),
                'filesPattern': '*stk',
                'samplingRate': 1.5,
                }

        # Id's should be set increasing from 1 if ### is not in the
        # pattern
        protAvgImport = self.newProtocol(ProtImportAverages, **args)
        protAvgImport.setObjLabel('improt averages')
        self.launchProtocol(protAvgImport, wait=False)


        kwargs = {'xDim': 64,
                  'yDim': 64,
                  'nDim': NUM_PART,
                  'samplingRate': 1.25,
                  'creationInterval': 15,
                  'delay': 0,
                  'setof': SET_OF_PARTICLES  # SetOfParticles
                  }

        # create input particles
        protStream = self.newProtocol(ProtCreateStreamData, **kwargs)
        protStream.setObjLabel('create Stream Mic')
        self.proj.launchProtocol(protStream, wait=False)

        counter = 1
        while not protStream.hasAttribute('outputMicrographs'):
            time.sleep(2)
            protStream = self._updateProtocol(protStream)
            if counter > 100:
                break
            counter += 1

