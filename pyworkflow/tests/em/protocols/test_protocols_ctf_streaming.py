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
import sys

from pyworkflow.em.protocol.monitors.protocol_monitor_ctf import CTF_LOG_SQLITE
from pyworkflow.tests import BaseTest, setupTestProject
from pyworkflow.em.protocol import ProtCreateStreamData, ProtMonitorSystem
from pyworkflow.em.packages.grigoriefflab import ProtCTFFind
from pyworkflow.protocol import getProtocolFromDb
from pyworkflow.em.packages.xmipp3 import XmippProtCTFMicrographsStr, XmippProtCTFSelection, XmippProtCTFMicrographs
from pyworkflow.em.data import SetOfCTF, SetOfMicrographs


# Load the number of movies for the simulation, by default equal 5, but
# can be modified in the environement

MICS = os.environ.get('SCIPION_TEST_MICS', 10)

CTF_SQLITE = "ctfs.sqlite"
MIC_SQLITE = "micrographs.sqlite"


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
        kwargs = {'xDim': 1024,
                  'yDim': 1024,
                  'nDim': MICS,
                  'samplingRate': 3.0,
                  'creationInterval': 30,
                  'delay':0,
                  'setof': 0 # SetOfMicrographs
                  }

        #create input micrographs
        protStream = self.newProtocol(ProtCreateStreamData, **kwargs)
        protStream.setObjLabel('create Stream Mic')
        self.proj.launchProtocol(protStream,wait=False)

        counter=1
        while not protStream.hasAttribute('outputMicrographs'):
            time.sleep(10)
            protStream = self._updateProtocol(protStream)
            if counter > 100:
                break
            counter += 1




######################### AJ START NEW STREAMING PROTOCOL #########################################

        kwargs = {
            'numberOfThreads': 5}

        protCTF = self.newProtocol(XmippProtCTFMicrographsStr, **kwargs)
        protCTF.inputMicrographs.set(protStream.outputMicrographs)
        self.proj.launchProtocol(protCTF)

####################################################################################################

        while not protCTF.hasAttribute('outputCTF'):
            time.sleep(10)
            protCTF = self._updateProtocol(protCTF)

        ctfSet = SetOfCTF(filename=protCTF._getPath(CTF_SQLITE))
        while not ctfSet.getSize()>0:
            time.sleep(10)
            protCTF = self._updateProtocol(protCTF)
            ctfSet = SetOfCTF(filename=protCTF._getPath(CTF_SQLITE))

        #kwargs = {
        #    'maxDefocus': 28000,
        #    'minDefocus': 1000,
        #    'astigmatism': 1000,
        #    'resolution': 7
        #}

        #protCTFSel = self.newProtocol(XmippProtCTFSelection, **kwargs)
        #protCTFSel.inputCTFs.set(protCTF.outputCTF)
        #self.proj.launchProtocol(protCTFSel)


        kwargs = {
            'numberOfThreads': 5}
        protCTF_ref = self.newProtocol(XmippProtCTFMicrographs, **kwargs)
        protCTF_ref.inputMicrographs.set(protStream.outputMicrographs)
        self.proj.launchProtocol(protCTF_ref)



######################### AJ CHECKING OUTPUT OF NEW STREAMING PROTOCOL ##############################


        micSet = SetOfMicrographs(filename=protStream._getPath(MIC_SQLITE))
        ctfSet = SetOfCTF(filename=protCTF._getPath(CTF_SQLITE))
        ctfSet_ref = SetOfCTF(filename=protCTF_ref._getPath(CTF_SQLITE))

        while not (ctfSet.getSize() == micSet.getSize() and ctfSet_ref.getSize() == micSet.getSize()):
            time.sleep(10)
            protCTF = self._updateProtocol(protCTF)
            protCTF_ref = self._updateProtocol(protCTF_ref)
            micSet = SetOfMicrographs(filename=protStream._getPath(MIC_SQLITE))
            ctfSet = SetOfCTF(filename=protCTF._getPath(CTF_SQLITE))
            ctfSet_ref = SetOfCTF(filename=protCTF_ref._getPath(CTF_SQLITE))


        ctfSet = SetOfCTF(filename=protCTF._getPath(CTF_SQLITE))
        ctfSet_ref = SetOfCTF(filename=protCTF_ref._getPath(CTF_SQLITE))

        baseFn = protCTF._getPath(CTF_SQLITE)
        self.assertTrue(os.path.isfile(baseFn))

        self.assertEqual(ctfSet.getSize(), MICS)

        for ctf, ctfRef in zip(ctfSet, ctfSet_ref):
            self.assertEqual(ctf._xmipp_ctfCritMaxFreq.get(), ctfRef._xmipp_ctfCritMaxFreq.get())
            self.assertEqual(ctf.isEnabled(), ctfRef.isEnabled())
            self.assertEqual(ctf._defocusU.get(), ctfRef._defocusU.get())
            self.assertEqual(ctf._defocusV.get(), ctfRef._defocusV.get())
            self.assertEqual(ctf._defocusRatio.get(), ctfRef._defocusRatio.get())

####################################################################################################