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
from pyworkflow.em.protocol import ProtImportAverages, ProtImportMicrographs
from pyworkflow.em.protocol.protocol_sets import ProtSubSet
from pyworkflow.em.packages.grigoriefflab import ProtCTFFind
from pyworkflow.em.packages.eman2.protocol_autopick import *
from pyworkflow.em.packages.xmipp3.protocol_extract_particles import *
from pyworkflow.em.packages.xmipp3.protocol_classification_gpuCorr import *


# Number of mics to be processed
NUM_MICS = 5 #maximum the number of mics in the relion set

class TestGpuCorrClassifier(BaseTest):
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


    def subsetMics(self, inputMics):
        protSubset = ProtSubSet()
        protSubset.inputFullSet.set(inputMics)
        protSubset.chooseAtRandom.set(True)
        protSubset.nElements.set(NUM_MICS)
        self.launchProtocol(protSubset)

        return protSubset


    def calculateCtf(self, inputMics):
        protCTF = ProtCTFFind(useCftfind4=True)
        protCTF.inputMicrographs.set(inputMics)
        protCTF.ctfDownFactor.set(1.0)
        protCTF.lowRes.set(0.05)
        protCTF.highRes.set(0.5)
        self.launchProtocol(protCTF)

        return protCTF


    def runPicking(self, inputMicrographs):
        """ Run a particle picking. """
        protPicking = SparxGaussianProtPicking(boxSize=64,
                                               numberOfThreads=1,
                                               numberOfMpi=1)
        protPicking.inputMicrographs.set(inputMicrographs)
        self.launchProtocol(protPicking)

        return protPicking

    def runExtractParticles(self, inputCoord, setCtfs):
        protExtract = self.newProtocol(XmippProtExtractParticles,
                                       boxSize=64,
                                       doInvert = True,
                                       doFlip = False)

        protExtract.inputCoordinates.set(inputCoord)
        protExtract.ctfRelations.set(setCtfs)

        self.launchProtocol(protExtract)

        return protExtract

    def runClassify(self, inputParts):
        numClasses = int(inputParts.getSize()/1000)
        if numClasses<2:
            numClasses=4
        protClassify = self.newProtocol(XmippProtGpuCrrCL2D,
                                        useReferenceImages=False,
                                        numberOfClasses=numClasses)
        protClassify.inputParticles.set(inputParts)
        self.launchProtocol(protClassify)

        return protClassify, numClasses

    def runClassify2(self, inputParts, inputRefs):
        numClasses = int(inputParts.getSize()/1000)
        if numClasses<inputRefs.getSize():
            numClasses=inputRefs.getSize()
        protClassify = self.newProtocol(XmippProtGpuCrrCL2D,
                                        useReferenceImages=True,
                                        numberOfClasses=numClasses)
        protClassify.inputParticles.set(inputParts)
        protClassify.referenceImages.set(inputRefs)
        self.launchProtocol(protClassify)

        return protClassify, numClasses

    def runClassify3(self, inputParts, inputRefs):
        numClasses = inputRefs.getSize()+2
        protClassify = self.newProtocol(XmippProtGpuCrrCL2D,
                                        useReferenceImages=True,
                                        numberOfClasses=numClasses)
        protClassify.inputParticles.set(inputParts)
        protClassify.referenceImages.set(inputRefs)
        protClassify.useAttraction.set(False)
        self.launchProtocol(protClassify)

        return protClassify, numClasses


    def test_pattern(self):

        protImportMics = self.importMicrographs()
        if protImportMics.isFailed():
            self.assertTrue(False)

        protImportAvgs = self.importAverages()
        if protImportAvgs.isFailed():
            self.assertTrue(False)
        outMics = protImportMics.outputMicrographs

        if NUM_MICS<20:
            protSubsetMics = self.subsetMics(protImportMics.outputMicrographs)
            if protSubsetMics.isFailed():
                self.assertTrue(False)
            outMics = protSubsetMics.outputMicrographs


        protCtf = self.calculateCtf(outMics)
        if protCtf.isFailed():
            self.assertTrue(False)

        protPicking = self.runPicking(outMics)
        if protPicking.isFailed():
            self.assertTrue(False)

        protExtract = self.runExtractParticles(protPicking.outputCoordinates,
                                               protCtf.outputCTF)
        if protExtract.isFailed():
            self.assertTrue(False)

        protClassify, numClasses = self.runClassify(protExtract.outputParticles)
        if protClassify.isFailed():
            self.assertTrue(False)
        if not protClassify.hasAttribute('outputClasses'):
            self.assertTrue(False)
        if protClassify.outputClasses.getSize() != numClasses:
            self.assertTrue(False)

        protClassify2, numClasses2 = self.runClassify2(
            protExtract.outputParticles, protImportAvgs.outputAverages)
        if protClassify2.isFailed():
            self.assertTrue(False)
        if not protClassify2.hasAttribute('outputClasses'):
            self.assertTrue(False)
        if protClassify2.outputClasses.getSize() != numClasses2:
            self.assertTrue(False)

        protClassify3, numClasses3 = self.runClassify3(
            protExtract.outputParticles, protImportAvgs.outputAverages)
        if protClassify3.isFailed():
            self.assertTrue(False)
        if not protClassify3.hasAttribute('outputClasses'):
            self.assertTrue(False)
        if protClassify3.outputClasses.getSize() != numClasses3:
            self.assertTrue(False)

