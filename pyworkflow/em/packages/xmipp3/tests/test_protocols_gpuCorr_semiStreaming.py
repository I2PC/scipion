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

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.em.protocol import ProtImportAverages, ProtImportMicrographs, \
    ProtCreateStreamData
from pyworkflow.em.protocol.protocol_create_stream_data import \
    SET_OF_MICROGRAPHS
from pyworkflow.protocol import getProtocolFromDb
from pyworkflow.em.packages.grigoriefflab import ProtCTFFind
from pyworkflow.em.packages.eman2.protocol_autopick import *
from pyworkflow.em.packages.xmipp3.protocol_extract_particles import *
from pyworkflow.em.packages.xmipp3.protocol_classification_gpuCorr_semi \
    import *
import time


# Number of mics to be processed
NUM_MICS = 5

class TestGpuCorrSemiStreaming(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsRelion = DataSet.getDataSet('relion_tutorial')

    def importAverages(self):
        prot = self.newProtocol(ProtImportAverages,
                                filesPath=self.dsRelion.getFile(
                                    'import/averages.mrcs'),
                                samplingRate=1.0)
        self.launchProtocol(prot)

        return prot

    def importMicrographs(self):
        prot = self.newProtocol(ProtImportMicrographs,
                                filesPath=self.dsRelion.getFile('micrographs'),
                                filesPattern='*.mrc',
                                samplingRateMode=1,
                                magnification=79096,
                                scannedPixelSize=56, voltage=300,
                                sphericalAberration=2.0)
        self.launchProtocol(prot)

        return prot

    def importMicrographsStr(self, fnMics):
        kwargs = {'inputMics': fnMics,
                  'nDim': NUM_MICS,
                  'creationInterval': 30,
                  'delay': 10,
                  'setof': SET_OF_MICROGRAPHS  # SetOfMics
                  }
        protStream = self.newProtocol(ProtCreateStreamData, **kwargs)
        protStream.setObjLabel('create Stream Mic')
        self.proj.launchProtocol(protStream, wait=False)

        return protStream

    def calculateCtf(self, inputMics):
        protCTF = ProtCTFFind(useCftfind4=True)
        protCTF.inputMicrographs.set(inputMics)
        protCTF.ctfDownFactor.set(1.0)
        protCTF.lowRes.set(0.05)
        protCTF.highRes.set(0.5)
        self.proj.launchProtocol(protCTF, wait=False)

        return protCTF


    def runPicking(self, inputMicrographs):
        """ Run a particle picking. """
        protPicking = SparxGaussianProtPicking(boxSize=64)
        protPicking.inputMicrographs.set(inputMicrographs)
        self.proj.launchProtocol(protPicking, wait=False)

        return protPicking

    def runExtractParticles(self, inputCoord, setCtfs):
        protExtract = self.newProtocol(XmippProtExtractParticles,
                                       boxSize=64,
                                       doInvert = True,
                                       doFlip = False)

        protExtract.inputCoordinates.set(inputCoord)
        protExtract.ctfRelations.set(setCtfs)

        self.proj.launchProtocol(protExtract, wait=False)

        return protExtract

    def runClassify(self, inputParts, inputAvgs):
        protClassify = self.newProtocol(XmippProtStrGpuCrrSimple,
                                        useAsRef=REF_AVERAGES)

        protClassify.inputParticles.set(inputParts)
        protClassify.inputRefs.set(inputAvgs)
        self.proj.launchProtocol(protClassify, wait=False)

        return protClassify


    def _updateProtocol(self, prot):
            prot2 = getProtocolFromDb(prot.getProject().path,
                                      prot.getDbPath(),
                                      prot.getObjId())
            # Close DB connections
            prot2.getProject().closeMapper()
            prot2.closeMappers()
            return prot2



    def test_pattern(self):

        protImportAvgs = self.importAverages()
        if protImportAvgs.isFailed():
            self.assertTrue(False)
        protImportMics = self.importMicrographs()
        if protImportMics.isFailed():
            self.assertTrue(False)

        protImportMicsStr = self.importMicrographsStr\
            (protImportMics.outputMicrographs)
        counter = 1
        while not protImportMicsStr.hasAttribute('outputMicrographs'):
            time.sleep(2)
            protImportMicsStr = self._updateProtocol(protImportMicsStr)
            if counter > 100:
                break
            counter += 1
        if protImportMicsStr.isFailed():
            self.assertTrue(False)

        protCtf = self.calculateCtf(protImportMicsStr.outputMicrographs)
        counter = 1
        while not protCtf.hasAttribute('outputCTF'):
            time.sleep(2)
            protCtf = self._updateProtocol(protCtf)
            if counter > 100:
                break
            counter += 1
        if protCtf.isFailed():
            self.assertTrue(False)

        protPicking = self.runPicking(protImportMicsStr.outputMicrographs)
        counter = 1
        while not protPicking.hasAttribute('outputCoordinates'):
            time.sleep(2)
            protPicking = self._updateProtocol(protPicking)
            if counter > 100:
                break
            counter += 1
        if protPicking.isFailed():
            self.assertTrue(False)

        protExtract = self.runExtractParticles(protPicking.outputCoordinates,
                                               protCtf.outputCTF)
        counter = 1
        while not protExtract.hasAttribute('outputParticles'):
            time.sleep(2)
            protExtract = self._updateProtocol(protExtract)
            if counter > 100:
                break
            counter += 1
        if protExtract.isFailed():
            self.assertTrue(False)

        protClassify = self.runClassify(protExtract.outputParticles,
                                        protImportAvgs.outputAverages)
        while protClassify.getStatus()!=STATUS_FINISHED:
            protClassify = self._updateProtocol(protClassify)
            if protClassify.isFailed():
                self.assertTrue(False)
                break
        if not protClassify.hasAttribute('outputClasses'):
            self.assertTrue(False)
        if protClassify.outputClasses.getSize() != \
                protImportAvgs.outputAverages.getSize():
            self.assertTrue(False)
