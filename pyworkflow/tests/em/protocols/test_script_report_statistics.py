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
import subprocess
from pyworkflow.tests import BaseTest, setupTestProject
from pyworkflow.em.protocol import ProtCreateStreamData
from pyworkflow.em.protocol.protocol_create_stream_data import \
    SET_OF_RANDOM_MICROGRAPHS
from pyworkflow.em.packages.grigoriefflab import ProtCTFFind
from pyworkflow.protocol import getProtocolFromDb
from os.path import dirname, realpath, join


# Load the number of movies for the simulation, by default equal 5, but
# can be modified in the environement
MICS = os.environ.get('SCIPION_TEST_MICS', 3)

CTF_SQLITE = "ctfs.sqlite"
MIC_SQLITE = "micrographs.sqlite"
MIC_DISCARDED_SQLITE = "micrographsDiscarded.sqlite"


class TestScriptReportStatistics(BaseTest):
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
        """ Create several micrographs
        """
        kwargs = {'xDim': 1024,
                  'yDim': 1024,
                  'nDim': MICS,
                  'samplingRate': 1.25,
                  'creationInterval': 5,
                  'delay': 0,
                  'setof': SET_OF_RANDOM_MICROGRAPHS  # SetOfMicrographs
                  }

        # create input micrographs
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

        # Compute CTF
        protCTF = ProtCTFFind(useCftfind4=True)
        protCTF.inputMicrographs.set(protStream.outputMicrographs)
        protCTF.ctfDownFactor.set(2)
        protCTF.highRes.set(0.4)
        protCTF.lowRes.set(0.05)
        protCTF.numberOfThreads.set(4)
        self.proj.launchProtocol(protCTF, wait=True)

        # execute script
        # cd tmp dir
        # popen
        # capture stderr
        # open setofCTF and compare
        self.assertTrue(False)
        SCIPION_HOME = dirname(dirname(realpath(__file__)))
        script = os.path.join(SCIPION_HOME,
                              'scripts/scipionbox_report_statistics.py')
        if not os.path.exists(script):
            print "HORROR script %d does not exist" % script
            exit(-1)
        args = ["python"]
        args += [script]
        args += ['-p', self.projName]
        scipion = os.path.join(SCIPION_HOME,'scipion')

        p = subprocess.Popen([scipion] +  args, stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        output, err = p.communicate()
        # if err contains something else (in addition to the dict)
        # remove it
        err = err[err.find("{"):err.find("}")+1]
        rc = p.returncode
        d = json.loads(err)
        if d:
            print "d", d
            # statistic.averageResolution =  d['averageResolution']
            # statistic.resolutionData = json.dumps(d['resolutionData'])
            # statistic.defocusData = json.dumps(d['defocusData'])
            # statistic.numberMovies = d['numberMovies']
            # statistic.save()
