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
# can be modified in the environement
MICS = os.environ.get('SCIPION_TEST_MICS', 3)
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
            while not (prot.isFinished() or prot.isFailed()):
                time.sleep(10)
                prot = self._updateProtocol(prot)
                if prot.hasAttribute("outputCTF"):
                    ctfSet = SetOfCTF(filename=prot._getPath(CTF_SQLITE))
                    baseFn = prot._getPath(CTF_SQLITE)
                    self.assertTrue(os.path.isfile(baseFn))
                    counter = 0
                    if ctfSet.getSize() > counter:
                        counter += 1
                        for ctf in ctfSet:
                            self.assertNotEqual(ctf._phaseShift.get(), None)
                            self.assertNotEqual(ctf._resolution.get(), None)
                            self.assertNotEqual(ctf._fitQuality.get(), None)
                            self.assertNotEqual(ctf.isEnabled(), None)
                            self.assertNotEqual(ctf._defocusU.get(), None)
                            self.assertNotEqual(ctf._defocusV.get(), None)
                            self.assertNotEqual(ctf._defocusRatio.get(), None)
            self.assertIsNotNone(prot.outputCTF,
                                 "Error: outputCTF is not produced "
                                 "in %s." % prot.getClassName())
            self.assertEqual(prot.outputCTF.getSize(), MICS)

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

        # wait until a micrograph has been created
        counter = 1
        while not protStream.hasAttribute('outputMicrographs'):
            time.sleep(10)
            protStream = self._updateProtocol(protStream)
            if counter > 10:
                break
            counter += 1

        # run ctffind4
        # then introduce monitor, checking all the time ctf and saving to
        # database
        protCTF = ProtCTFFind(useCftfind4=True)
        #time.sleep(10)
        protCTF.inputMicrographs.set(protStream.outputMicrographs)
        protCTF.ctfDownFactor.set(2)
        protCTF.highRes.set(0.4)
        protCTF.lowRes.set(0.05)
        protCTF.numberOfThreads.set(4)
        self.proj.launchProtocol(protCTF, wait=True)
        checkOutputs(protCTF)

        kwargs = {'ctfDownFactor': 2,
                  'numberOfThreads': 3
                  }
        protCTF2 = self.newProtocol(XmippProtCTFMicrographs, **kwargs)
        protCTF2.inputMicrographs.set(protStream.outputMicrographs)
        self.proj.launchProtocol(protCTF2)
        checkOutputs(protCTF2)
        
        # run gctf
        # check if box has nvidia cuda libs.
        try:
            nvmlInit()  # fails if not GPU attached
            protCTF3 = ProtGctf()
            protCTF3.inputMicrographs.set(protStream.outputMicrographs)
            protCTF3.ctfDownFactor.set(2)
            self.proj.launchProtocol(protCTF3, wait=False)
            checkOutputs(protCTF3)

        except NVMLError, err:
            print("Cannot find GPU."
                  "I assume that no GPU is connected to this machine")


        # run xmipp ctf. Since this is the slower method wait until finish
        # before running asserts
        kwargs = {
            'numberOfThreads': MICS + 1}
        protCTF2 = self.newProtocol(XmippProtCTFMicrographs, **kwargs)
        protCTF2.inputMicrographs.set(protStream.outputMicrographs)
        protCTF2.ctfDownFactor.set(2)

        self.proj.launchProtocol(protCTF2, wait=True)

        ctfSet = SetOfCTF(filename=protCTF._getPath(CTF_SQLITE))

        baseFn = protCTF._getPath(CTF_SQLITE)
        self.assertTrue(os.path.isfile(baseFn))

        self.assertSetSize(ctfSet, MICS, "Ctffind4 output size does not match")

        for ctf in ctfSet:
            self.assertNotEqual(ctf._phaseShift.get(), None)
            self.assertNotEqual(ctf._resolution.get(), None)
            self.assertNotEqual(ctf._fitQuality.get(), None)
            self.assertNotEqual(ctf.isEnabled(), None)
            self.assertNotEqual(ctf._defocusU.get(), None)
            self.assertNotEqual(ctf._defocusV.get(), None)
            self.assertNotEqual(ctf._defocusRatio.get(), None)

