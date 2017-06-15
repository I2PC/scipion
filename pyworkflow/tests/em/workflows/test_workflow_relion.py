# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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
# *
# **************************************************************************

import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3 import *
from pyworkflow.em.packages.xmipp3.constants import SAME_AS_PICKING
from pyworkflow.em.packages.grigoriefflab import *
from pyworkflow.em.packages.relion import *
from pyworkflow.em.packages.relion.protocol_autopick_v2 import *
from test_workflow import TestWorkflow



class TestWorkflowRelionPick(TestWorkflow):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('relion_tutorial')

    def _runPickWorkflow(self):
        #First, import a set of micrographs
        print "Importing a set of micrographs..."
        protImport = self.newProtocol(ProtImportMicrographs,
                                      filesPath=self.ds.getFile('micrographs'),
                                      filesPattern='*.mrc',
                                      samplingRateMode=1,
                                      magnification=79096,
                                      scannedPixelSize=56, voltage=300,
                                      sphericalAberration=2.0)
        protImport.setObjLabel('import 20 mics')
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputMicrographs,
                             "There was a problem with the import")
        
        print "Preprocessing the micrographs..."
        protCropMics = self.newProtocol(XmippProtPreprocessMicrographs,
                                          doCrop=True, cropPixels=25)
        protCropMics.inputMicrographs.set(protImport.outputMicrographs)
        protCropMics.setObjLabel('crop 50px')
        self.launchProtocol(protCropMics)
        self.assertIsNotNone(protCropMics.outputMicrographs,
                             "There was a problem with the downsampling")
        self.protCropMics = protCropMics

        # Now estimate CTF on the micrographs with ctffind
        print "Performing CTFfind..."
        protCTF = self.newProtocol(ProtCTFFind,
                                   useCtffind4=True,
                                   lowRes=0.02, highRes=0.45,
                                   minDefocus=1.2, maxDefocus=3,
                                   runMode=1,
                                   numberOfMpi=1, numberOfThreads=3)
        protCTF.inputMicrographs.set(protCropMics.outputMicrographs)
        protCTF.setObjLabel('CTF ctffind')
        self.protCTF = protCTF
        self.launchProtocol(protCTF)

        print "Importing 2D averages (subset of 4)"
        ih = ImageHandler()
        classesFn = self.ds.getFile('import/classify2d/extra/'
                                    'relion_it015_classes.mrcs')

        outputName = 'input_averages.mrcs'
        inputTmp = os.path.abspath(self.proj.getTmpPath())
        outputFn = self.proj.getTmpPath(outputName)

        for i, index in enumerate([5, 16, 17, 18, 24]):
            ih.convert((index, classesFn), (i+1, outputFn))

        protAvgs = self.newProtocol(ProtImportAverages,
                                    objLabel='avgs - 5',
                                    filesPath=inputTmp,
                                    filesPattern=outputName,
                                    samplingRate=7.08
                                    )

        self.launchProtocol(protAvgs)
        # Select some good averages from the iterations mrcs a

        protPick1 = self.newProtocol(ProtRelion2Autopick,
                                     objLabel='autopick refs (optimize)',
                                     runType=RUN_OPTIMIZE,
                                     micrographsList=5,
                                     referencesType=REF_AVERAGES,
                                     refsHaveInvertedContrast=True,
                                     particleDiameter=380
                                     )

        protPick1.inputMicrographs.set(protCropMics.outputMicrographs)
        protPick1.ctfRelations.set(protCTF.outputCTF)
        protPick1.inputReferences.set(protAvgs.outputAverages)

        self.launchProtocol(protPick1)

        return protPick1

    def test_ribo(self):
        protPick1 = self._runPickWorkflow()

        # Launch the same picking run but now for all micrographs
        protPick2 = self.proj.copyProtocol(protPick1)
        protPick2.setObjLabel('autopick refs (all)')
        protPick2.runType.set(RUN_COMPUTE)
        self.launchProtocol(protPick2)

        # Launch now using the Gaussian as references
        protPick3 = self.proj.copyProtocol(protPick1)
        protPick3.setObjLabel('autopick gauss (optimize)')
        protPick3.referencesType.set(REF_BLOBS)
        protPick3.inputReferences.set(None)
        self.launchProtocol(protPick3)

        # Launch the same picking run but now in 1 GPU.
        protPick4 = self.proj.copyProtocol(protPick1)
        protPick4.setObjLabel('autopick refs (optimize) 1 GPU')
        protPick4.gpusToUse.set('0:0:0:0')
        self.launchProtocol(protPick4)


class TestWorkflowRelionExtract(TestWorkflowRelionPick):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('relion_tutorial')

    def _checkOutput(self, prot, **kwargs):
        # Read expected parameters
        size = kwargs.get('size', 2618)
        sampling = kwargs.get('sampling', 7.08)
        dim = kwargs.get('dim', 64)

        outputParts = getattr(prot, 'outputParticles', None)

        self.assertIsNotNone(outputParts)
        self.assertEqual(outputParts.getSize(), size)

        first = outputParts.getFirstItem()
        ctfModel = first.getCTF()
        
        self.assertEqual(first.getDim(), (dim, dim, 1))
        self.assertAlmostEqual(first.getSamplingRate(), sampling, delta=0.001)
        self.assertAlmostEqual(ctfModel.getDefocusU(), 23467, delta=10)
        self.assertAlmostEqual(ctfModel.getDefocusV(), 23308, delta=10)

    def test_ribo(self):
        """ Reimplement this test to run several extract cases. """
        protPick1 = self._runPickWorkflow()
        proj = protPick1.getProject()
        size = protPick1.outputCoordinates.getSize()

        protExtract = self.newProtocol(ProtRelionExtractParticles,
                                       objLabel='extract - box=64',
                                       boxSize=64,
                                       doInvert=True
                                       )

        protExtract.inputCoordinates.set(protPick1.outputCoordinates)
        protExtract.ctfRelations.set(self.protCTF.outputCTF)

        self.launchProtocol(protExtract)
        self._checkOutput(protExtract, size=size)

        # Now test the re-scale option
        protExtract2 = self.proj.copyProtocol(protExtract)
        protExtract2.setObjLabel('extract - rescale 32')
        protExtract2.doRescale.set(True)
        protExtract2.rescaledSize.set(32)
        self.launchProtocol(protExtract2)

        self._checkOutput(protExtract2, size=size, dim=32, sampling=14.16)
        
        # Now test changing micrographs source option
        subsetProt = self.newProtocol(em.ProtSubSet,
                                      chooseAtRandom=True,
                                      nElements=10)
        subsetProt.inputFullSet.set(self.protCropMics.outputMicrographs)
        self.launchProtocol(subsetProt)

        protExtract3 = self.proj.copyProtocol(protExtract)
        protExtract3.setObjLabel('extract - Other')
        protExtract3.downsampleType.set(1)
        protExtract3.inputMicrographs.set(subsetProt.outputMicrographs)
        self.launchProtocol(protExtract3)
        
        # The size of the set is different every time is executed the test
        # due to the seleccion of the micrographs is random.
        
        partSize = protExtract3.outputParticles.getSize()
        self._checkOutput(protExtract3, size=partSize)