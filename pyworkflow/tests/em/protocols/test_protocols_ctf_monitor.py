# ***************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
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
from pyworkflow.em.protocol import ProtCreateStreamData, ProtMonitorSystem
from pyworkflow.em.packages.grigoriefflab import ProtCTFFind
from pyworkflow.protocol import getProtocolFromDb
from pyworkflow.em.protocol import ProtMonitorCTF


# Load the number of movies for the simulation, by default equal 5, but
# can be modified in the environement
MICS = os.environ.get('SCIPION_TEST_MICS', 3)


class TestCtfStream(BaseTest):
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
        kwargs = {'xDim': 4096,
                  'yDim': 4096,
                  'nDim': MICS,
                  'samplingRate': 1.25,
                  'creationInterval': 5,
                  'delay':0,
                  'setof': 0 # SetOfMicrographs
                }
        
        #put some stress on the system
        protStream = self.newProtocol(ProtCreateStreamData, **kwargs)
        protStream.setObjLabel('create Stream Mic')
        self.proj.launchProtocol(protStream,wait=False)

        counter=1
        while not protStream.hasAttribute('outputMicrographs'):
            time.sleep(10)
            protStream = self._updateProtocol(protStream)
            if counter > 1000:
                break
            counter += 1

        #then introduce monitor, checking all the time ctf and sving to database
        protCTF = ProtCTFFind(useCftfind4=True)
        protCTF.inputMicrographs.set(protStream.outputMicrographs)
        protCTF.ctfDownFactor.set(2)
        protCTF.highRes.set(0.4)
        protCTF.lowRes.set(0.05)
        protCTF.numberOfThreads.set(4)
        self.proj.launchProtocol(protCTF, wait=False)

        counter=1
        print("while")
        while not protCTF.hasAttribute('outputCTF'):
            print("whilenot")
            time.sleep(10)
            protCTF = self._updateProtocol(protCTF)
            if counter > 1000:
                break
            counter += 1
        print("afterwhile")
        kwargs = {'samplingInterval':10,
                  'interval': 300,
                  'maxDefocus': 40000,
                  'minDefocus': 1000,
                  'astigmatism': 0.2,
                }

        protMonitor = self.newProtocol(ProtMonitorCTF, **kwargs)
        protMonitor.inputProtocol.set(protCTF)
        self.launchProtocol(protMonitor)

        baseFn = protMonitor._getPath(protMonitor.dataBase)
        self.assertTrue(os.path.isfile(baseFn))