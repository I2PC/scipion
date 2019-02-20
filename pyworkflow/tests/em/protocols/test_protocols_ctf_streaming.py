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
from pyworkflow.em.data import SetOfCTF
from pyworkflow.tests import BaseTest, setupTestProject
from pyworkflow.em.protocol import ProtCreateStreamData, ProtMonitorSystem
from pyworkflow.em.packages.grigoriefflab import ProtCTFFind
from pyworkflow.protocol import getProtocolFromDb
from pyworkflow.em.protocol.protocol_create_stream_data import \
    SET_OF_RANDOM_MICROGRAPHS
from pyworkflow.em.packages.xmipp3.protocol_ctf_micrographs import\
    XmippProtCTFMicrographs
from pyworkflow.em.packages.gctf import ProtGctf
from pyworkflow.em.protocol.monitors.pynvml import nvmlInit, NVMLError
# Load the number of movies for the simulation, by default equal 5, but
# can be modified in the environment
MICS = os.environ.get('SCIPION_TEST_MICS', 6)
CTF_SQLITE = "ctfs.sqlite"


class TestCtfStreaming(BaseTest):
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
        def checkOutputs(prot):

            t0 = time.time()

            while not (prot.isFinished() or prot.isFailed()):

                # Time out 6 minutes, just in case
                tdelta = time.time() -t0
                if tdelta > 6*60:
                    break

                time.sleep(10)

                prot = self._updateProtocol(prot)

                # Check if the protocol is still launched
                if prot.isLaunched():
                    break
                elif prot.isScheduled():
                    continue

                if prot.hasAttribute("outputCTF"):
                    ctfSet = SetOfCTF(filename=prot._getPath(CTF_SQLITE))
                    baseFn = prot._getPath(CTF_SQLITE)
                    self.assertTrue(os.path.isfile(baseFn))
                    counter = 0
                    if ctfSet.getSize() > counter:
                        counter += 1
                        for ctf in ctfSet:
                            self.assertNotEqual(ctf._resolution.get(), None)
                            self.assertNotEqual(ctf._fitQuality.get(), None)
                            self.assertNotEqual(ctf.isEnabled(), None)
                            self.assertNotEqual(ctf._defocusU.get(), None)
                            self.assertNotEqual(ctf._defocusV.get(), None)
                            self.assertNotEqual(ctf._defocusRatio.get(), None)
                            if ctf.getPhaseShift():
                                self.assertNotEqual(ctf.getPhaseShift(), None)
            self.assertIsNotNone(prot.outputCTF,
                                 "Error: outputCTF is not produced "
                                 "in %s." % prot.getClassName())
            self.assertEqual(prot.outputCTF.getSize(), MICS)

        # Simulate streaming with create stream data protocol
        kwargs = {'xDim': 4096,
                  'yDim': 4096,
                  'nDim': MICS,
                  'samplingRate': 1.25,
                  'creationInterval': 15,
                  'delay': 0,
                  'setof': SET_OF_RANDOM_MICROGRAPHS}  # SetOfMicrographs

        # create mic in streaming mode
        protStream = self.newProtocol(ProtCreateStreamData, **kwargs)
        protStream.setObjLabel('create Stream Mic')
        self.proj.launchProtocol(protStream, wait=False)

        # 1st ctf - ctffind4 in streaming
        protCTF = ProtCTFFind(useCftfind4=True)

        protCTF.inputMicrographs.set(protStream)
        protCTF.inputMicrographs.setExtended('outputMicrographs')
        protCTF.ctfDownFactor.set(2)
        protCTF.findPhaseShift.set(True)
        protCTF.slowSearch.set(False)
        protCTF.highRes.set(0.4)
        protCTF.lowRes.set(0.05)
        protCTF.numberOfThreads.set(4)
        self.proj.scheduleProtocol(protCTF)

        # 2nd ctf - Xmipp CTF
        kwargs = {'ctfDownFactor': 2,
                  'numberOfThreads': 4
                  }
        protCTF2 = self.newProtocol(XmippProtCTFMicrographs, **kwargs)
        protCTF2.inputMicrographs.set(protStream)
        protCTF2.inputMicrographs.setExtended('outputMicrographs')
        self.proj.scheduleProtocol(protCTF2)

        
        # 3rd ctf - Gctf
        protCTF3 = None
        try:
            # check if box has nvidia cuda libs.
            nvmlInit()  # fails if not GPU attached
            protCTF3 = ProtGctf()
            protCTF3.inputMicrographs.set(protStream)
            protCTF3.inputMicrographs.setExtended('outputMicrographs')
            protCTF3.ctfDownFactor.set(2)
            self.proj.scheduleProtocol(protCTF3)

        except NVMLError, err:
            print("Cannot find GPU."
                  "I assume that no GPU is connected to this machine")

        # Check first ctf output - Ctffind
        checkOutputs(protCTF)

        # Check 2nd ctf output - Xmipp
        checkOutputs(protCTF2)

        # If GPU and therefore GCTF protocl
        if protCTF3 is not None:
            checkOutputs(protCTF3)