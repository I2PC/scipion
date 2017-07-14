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
from pyworkflow.em.protocol import ProtCreateStreamData
from pyworkflow.em.packages.grigoriefflab import ProtCTFFind
from pyworkflow.protocol import getProtocolFromDb
from pyworkflow.em.packages.xmipp3 import XmippProtCTFSelection
from pyworkflow.em.data import SetOfCTF, SetOfMicrographs


# Load the number of movies for the simulation, by default equal 5, but
# can be modified in the environement
MICS = os.environ.get('SCIPION_TEST_MICS', 10)

CTF_SQLITE = "ctfs.sqlite"
MIC_SQLITE = "micrographs.sqlite"


class TestCtfSelection(BaseTest):
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
                  'samplingRate': 1.25,
                  'creationInterval': 5,
                  'delay':0,
                  'setof': 0 # SetOfMicrographs
                  }

        #create input micrographs
        protStream = self.newProtocol(ProtCreateStreamData, **kwargs)
        protStream.setObjLabel('create Stream Mic')
        self.proj.launchProtocol(protStream,wait=False)

        counter=1
        while not protStream.hasAttribute('outputMicrographs'):
            time.sleep(2)
            protStream = self._updateProtocol(protStream)
            if counter > 100:
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

        counter = 1

        while not protCTF.hasAttribute('outputCTF'):

            time.sleep(2)
            protCTF = self._updateProtocol(protCTF)
            if counter > 100:
                break
            counter += 1

        kwargs = {
            'maxDefocus': 28000,
            'minDefocus': 1000,
            'astigmatism': 1000,
            'resolution': 7
        }

        protCTFSel = self.newProtocol(XmippProtCTFSelection, **kwargs)
        protCTFSel.inputCTF.set(protCTF.outputCTF)
        self.proj.launchProtocol(protCTFSel,  wait=False)

        counter = 1

        while not protCTFSel.hasAttribute('outputCTF'):

            time.sleep(2)
            protCTFSel = self._updateProtocol(protCTFSel)
            if counter > 100:
                break
            counter += 1

        kwargs = {
            'maxDefocus': 40000,
            'minDefocus': 1000,
            'astigmatism': 1000,
            'resolution': 3.7
        }

        protCTFSel2 = self.newProtocol(XmippProtCTFSelection, **kwargs)
        protCTFSel2.inputCTF.set(protCTFSel.outputCTF)
        self.proj.launchProtocol(protCTFSel2)

        counter=1

        while not (protCTFSel2.hasAttribute('outputCTF') and
                   protCTFSel2.hasAttribute('outputMicrograph')):

            time.sleep(2)
            protCTFSel2 = self._updateProtocol(protCTFSel2)
            if counter > 100:
                self.assertTrue(False)
            counter += 1

        baseFn = protCTFSel2._getPath(CTF_SQLITE)
        self.assertTrue(os.path.isfile(baseFn))
        baseFn = protCTFSel2._getPath(MIC_SQLITE)
        self.assertTrue(os.path.isfile(baseFn))

        micSet = SetOfMicrographs(filename=protCTFSel2._getPath(MIC_SQLITE))
        ctfSet = SetOfCTF(filename=protCTFSel2._getPath(CTF_SQLITE))
        counter = 1
        while not (ctfSet.getSize()==10 and micSet.getSize()==10):
            time.sleep(2)
            micSet = SetOfMicrographs(filename=protCTFSel2._getPath(MIC_SQLITE))
            ctfSet = SetOfCTF(filename=protCTFSel2._getPath(CTF_SQLITE))
            if counter > 100:
                self.assertTrue(False)
            counter += 1

        ctfSetIn = SetOfCTF(filename=protCTF._getPath(CTF_SQLITE))
        for ctf, ctfOut in zip(ctfSetIn, ctfSet):
            defocusU = ctf.getDefocusU()
            defocusV = ctf.getDefocusV()
            astigm = defocusU - defocusV
            resol = ctf._ctffind4_ctfResolution.get() # TODO
            if defocusU > 1000 and defocusU < 28000 and \
            defocusV > 1000 and defocusV < 28000 and \
            astigm < 1000 and resol < 3.7:
                self.assertTrue(ctfOut.isEnabled(),
                                "ctf with id %d enabled is "
                                "False and should be True" % ctfOut.getObjId())
            else:
                self.assertFalse(ctfOut.isEnabled(),
                                 "ctf with id %d enabled is "
                                 "True and should be False" % ctfOut.getObjId())
