# **************************************************************************
# *
# * Authors:    Laura del Cano (ldelcano@cnb.csic.es)
# *             Josue Gomez Blanco (josue.gomez-blanco@mcgill.ca)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
from os.path import join, basename

from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3 import *
from pyworkflow.protocol import getProtocolFromDb
import pyworkflow.utils as pwutils

# Some utility functions to import micrographs that are used
# in several tests.
class TestXmippBase(BaseTest):
    @classmethod
    def setData(cls):
        cls.dataset = DataSet.getDataSet('xmipp_tutorial')
        cls.micFn = cls.dataset.getFile('mic1')
        cls.micsFn = cls.dataset.getFile('allMics')
        cls.coordsDir = cls.dataset.getFile('posSupervisedDir')
        cls.allCrdsDir = cls.dataset.getFile('posAllDir')

    @classmethod
    def runImportMicrograph(cls, pattern, samplingRate, voltage,
                            scannedPixelSize, magnification,
                            sphericalAberration):
        """ Run an Import micrograph protocol. """

        # We have two options:
        # 1) pass the SamplingRate or
        # 2) the ScannedPixelSize + microscope magnification
        if not samplingRate is None:
            cls.protImport = cls.newProtocol(ProtImportMicrographs,
                                             samplingRateMode=0,
                                             filesPath=pattern,
                                             samplingRate=samplingRate,
                                             magnification=magnification,
                                             voltage=voltage,
                                             sphericalAberration=sphericalAberration)
        else:
            cls.protImport = cls.newProtocol(ProtImportMicrographs,
                                             samplingRateMode=1,
                                             filesPath=pattern,
                                             scannedPixelSize=scannedPixelSize,
                                             voltage=voltage,
                                             magnification=magnification,
                                             sphericalAberration=sphericalAberration)

        cls.protImport.setObjLabel('import mics')
        cls.launchProtocol(cls.protImport)
        if cls.protImport.isFailed():
            raise Exception("Protocol has failed. Error: ", cls.protImport.getErrorMessage())
        # check that input micrographs have been imported (a better way to do this?)
        if cls.protImport.outputMicrographs is None:
            raise Exception('Import of micrograph: %s, failed. outputMicrographs is None.' % pattern)
        return cls.protImport

    @classmethod
    def runImportMicrographBPV(cls, pattern):
        """ Run an Import micrograph protocol. """
        return cls.runImportMicrograph(pattern, samplingRate=1.237,
                                       voltage=300, sphericalAberration=2,
                                       scannedPixelSize=None, magnification=56000)

    @classmethod
    def runDownsamplingMicrographs(cls, mics, downFactorValue, threads=1):
        # test downsampling a set of micrographs
        cls.protDown = XmippProtPreprocessMicrographs(doDownsample=True,
                                                      downFactor=downFactorValue,
                                                      numberOfThreads=threads)
        cls.protDown.inputMicrographs.set(mics)
        cls.proj.launchProtocol(cls.protDown, wait=True)
        return cls.protDown

    @classmethod
    def runFakedPicking(cls, mics, pattern):
        """ Run a faked particle picking. Coordinates already existing. """
        cls.protPP = XmippProtParticlePicking(importFolder=pattern, runMode=1)
        cls.protPP.inputMicrographs.set(mics)
        cls.proj.launchProtocol(cls.protPP, wait=True)
        # check that faked picking has run ok
        if cls.protPP.outputCoordinates is None:
            raise Exception('Faked particle picking failed. outputCoordinates is None.')
        return cls.protPP


class TestImportMicrographs(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData()

    """This class check if a set of micrographs is imported properly"""
    def testImport1(self):
        pattern = self.micsFn
        samplingRate = None
        scannedPixelSize = 7
        magnification = 56000
        voltage = 300
        sphericalAberration = 2

        protImport = self.runImportMicrograph(pattern, samplingRate=samplingRate,
                                              scannedPixelSize=scannedPixelSize,
                                              magnification=magnification, voltage=voltage,
                                              sphericalAberration=sphericalAberration)
        if protImport.isFailed():
            raise Exception(protImport.getError())

        m = protImport.outputMicrographs.getAcquisition()
        # Check that sampling rate on output micrographs is equal to
        self.assertEquals(protImport.outputMicrographs.getScannedPixelSize(), scannedPixelSize, "Incorrect ScannedPixelSize on output micrographs.")
        self.assertEquals(m.getMagnification(), magnification, "Incorrect Magnification on output micrographs.")
        self.assertEquals(m.getVoltage(), voltage, "Incorrect Voltage on output micrographs.")
        self.assertEquals(m.getSphericalAberration(), sphericalAberration, "Incorrect SphericalAberration on output micrographs.")

    def testImport2(self):
        pattern = self.micsFn
        samplingRate = 2.56
        scannedPixelSize = 7
        magnification = 56000
        voltage = 400
        sphericalAberration = 2.5

        protImport = self.runImportMicrograph(pattern, samplingRate=samplingRate,
                                              scannedPixelSize=scannedPixelSize,
                                              magnification=magnification, voltage=voltage,
                                              sphericalAberration=sphericalAberration)
        m = protImport.outputMicrographs.getAcquisition()
        # Check that sampling rate on output micrographs is equal to
        self.assertEquals(protImport.outputMicrographs.getSamplingRate(), samplingRate, "Incorrect SamplingRate on output micrographs.")
        self.assertEquals(m.getVoltage(), voltage, "Incorrect Voltage on output micrographs.")
        self.assertEquals(m.getSphericalAberration(), sphericalAberration, "Incorrect Spherical aberration on output micrographs.")

    def testImport3(self):
        pattern = self.dataset.getFile('micrographs/BPV_####.mrc')
        samplingRate = 2.56
        scannedPixelSize = 7
        magnification = 56000
        voltage = 400
        sphericalAberration = 2.5

        protImport = self.runImportMicrograph(pattern, samplingRate=samplingRate,
                                              scannedPixelSize=scannedPixelSize,
                                              magnification=magnification, voltage=voltage,
                                              sphericalAberration=sphericalAberration)
        m = protImport.outputMicrographs.getAcquisition()
        # Check that sampling rate on output micrographs is equal to
        self.assertEquals(protImport.outputMicrographs.getSamplingRate(), samplingRate, "Incorrect SamplingRate on output micrographs.")
        self.assertEquals(m.getVoltage(), voltage, "Incorrect Voltage on output micrographs.")
        self.assertEquals(m.getSphericalAberration(), sphericalAberration, "Incorrect Spherical aberration on output micrographs.")


class TestXmippPreprocessMicrographs(TestXmippBase):
    """This class check if the preprocessing micrographs protocol in Xmipp works properly."""
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData()
        cls.protImport = cls.runImportMicrographBPV(cls.micFn)

    def testDownsampling(self):
        # test downsampling a set of micrographs
        downFactorValue = 2
        protDown = XmippProtPreprocessMicrographs(doDownsample=True, downFactor=downFactorValue)
        protDown.inputMicrographs.set(self.protImport.outputMicrographs)
        self.proj.launchProtocol(protDown, wait=True)
        # check that output micrographs have double sampling rate than input micrographs
        self.assertEquals(protDown.outputMicrographs.getSamplingRate(), self.protImport.outputMicrographs.getSamplingRate()*downFactorValue, "Micrographs incorrectly downsampled")
        self.assertTrue(protDown.isFinished(), "Downsampling failed")

    def testPreprocessing(self):
        # test Crop, Take logarithm and Remove bad pixels on a set of micrographs
        cropPixels = 100
        protPreprocess = XmippProtPreprocessMicrographs(doCrop=True, doLog=True, doRemoveBadPix=True, cropPixels=cropPixels)
        protPreprocess.inputMicrographs.set(self.protImport.outputMicrographs)
        self.proj.launchProtocol(protPreprocess, wait=True)
        self.assertIsNotNone(protPreprocess.outputMicrographs, "SetOfMicrographs has not been preprocessed.")

    def testInvertNormalize(self):
        # test invert and normalize a set of micrographs
        protInvNorm = XmippProtPreprocessMicrographs(doInvert=True, doNormalize=True)
        protInvNorm.inputMicrographs.set(self.protImport.outputMicrographs)
        self.proj.launchProtocol(protInvNorm, wait=True)
        self.assertIsNotNone(protInvNorm.outputMicrographs, "SetOfMicrographs has not been preprocessed.")

    def testSmooth(self):
        # test smooth a set of micrographs
        protSmooth = XmippProtPreprocessMicrographs(doSmooth=True, sigmaConvolution=3)
        protSmooth.inputMicrographs.set(self.protImport.outputMicrographs)
        self.proj.launchProtocol(protSmooth, wait=True)
        self.assertIsNotNone(protSmooth.outputMicrographs, "SetOfMicrographs has not been preprocessed.")

    def _updateProtocol(self, prot):
        prot2 = getProtocolFromDb(prot.getProject().path,
                                  prot.getDbPath(),
                                  prot.getObjId())
        # Close DB connections
        prot2.getProject().closeMapper()
        prot2.closeMappers()
        return prot2

    def testStreaming(self):
        """ Import several Particles from a given pattern.
        """
        kwargs = {'xDim': 1024,
                  'yDim': 1024,
                  'nDim': 5,
                  'samplingRate': 3.0,
                  'creationInterval': 3,
                  'delay': 0,
                  'setof': 2   # RandomMicrographs
                  }

        # create input micrographs
        protStream = self.newProtocol(ProtCreateStreamData, **kwargs)
        self.proj.launchProtocol(protStream, wait=False)

        while not protStream.hasAttribute('outputMicrographs'):
            time.sleep(3)
            protStream = self._updateProtocol(protStream)

        protDenoise = self.newProtocol(XmippProtPreprocessMicrographs,
                                       doDenoise=True, maxIteration=50)
        protDenoise.inputMicrographs.set(protStream.outputMicrographs)
        self.proj.launchProtocol(protDenoise)


        micSet = SetOfMicrographs(
            filename=protStream._getPath("micrographs.sqlite"))
        denoiseSet = len(protDenoise._readDoneList())

        while not (denoiseSet == micSet.getSize()):
            time.sleep(5)
            print("Imported mics: %d, processed mics: %d" % (
            micSet.getSize(), len(protDenoise._readDoneList())))
            micSet = SetOfMicrographs(
                filename=protStream._getPath("micrographs.sqlite"))
            denoiseSet = len(protDenoise._readDoneList())


class TestXmippCTFEstimation(TestXmippBase):
    """This class check if the protocol to determine the CTF in Xmipp works properly."""
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData()
        cls.protImport = cls.runImportMicrographBPV(cls.micFn)
#         cls.protDown = cls.runDownsamplingMicrographs(cls.protImport.outputMicrographs, 3)

    def testCTF(self):
        # Estimate CTF on the downsampled micrographs
        print "Performing CTF..."
        protCTF = XmippProtCTFMicrographs()
        protCTF.inputMicrographs.set(self.protImport.outputMicrographs)
        protCTF.ctfDownFactor.set(2)
        self.proj.launchProtocol(protCTF, wait=True)
        self.assertIsNotNone(protCTF.outputCTF, "SetOfCTF has not been produced.")
        ctfModel = protCTF.outputCTF.getFirstItem()
        self.assertAlmostEquals(ctfModel.getDefocusU(),23928.4, delta=500)
        self.assertAlmostEquals(ctfModel.getDefocusV(),23535.2, delta=500)
        self.assertAlmostEquals(ctfModel.getDefocusAngle(), 63.669, delta=5)
        sampling = ctfModel.getMicrograph().getSamplingRate()
        self.assertAlmostEquals(sampling, 2.474, delta=0.001)


class TestXmippAutomaticPicking(TestXmippBase):
    """This class check if the protocol to pick the micrographs automatically in Xmipp works properly."""
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData()
        cls.protImport1 = cls.runImportMicrographBPV(cls.micsFn)
        cls.protImport2 = cls.runImportMicrographBPV(cls.micFn)
        cls.protDown1 = cls.runDownsamplingMicrographs(cls.protImport1.outputMicrographs, 5)
        cls.protDown2 = cls.runDownsamplingMicrographs(cls.protImport2.outputMicrographs, 5)
        cls.protPP = cls.runFakedPicking(cls.protDown1.outputMicrographs, cls.coordsDir)

    def testAutomaticPicking(self):
        print "Run automatic particle picking"
        protAutomaticPP = XmippParticlePickingAutomatic()
        protAutomaticPP.xmippParticlePicking.set(self.protPP)
        self.proj.launchProtocol(protAutomaticPP, wait=True)
        self.assertIsNotNone(protAutomaticPP.outputCoordinates,
                             "There was a problem with the automatic particle picking")

    def testAutomaticPickingOther(self):
        print "Run automatic particle picking"
        protAutomaticPP = XmippParticlePickingAutomatic()
        protAutomaticPP.xmippParticlePicking.set(self.protPP)
        protAutomaticPP.inputMicrographs.set(self.protDown2.outputMicrographs)
        protAutomaticPP.micsToPick.set(1)
        self.proj.launchProtocol(protAutomaticPP, wait=True)
        self.assertIsNotNone(protAutomaticPP.outputCoordinates,
                             "There was a problem with the automatic particle picking")


class TestXmippExtractParticles(TestXmippBase):
    """This class check if the protocol to extract particles
    in Xmipp works properly.
    """
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData()
        cls.DOWNSAMPLING = 5.0
        cls.protImport = cls.runImportMicrographBPV(cls.micsFn)
        cls.protDown = cls.runDownsamplingMicrographs(cls.protImport.outputMicrographs,
                                                      cls.DOWNSAMPLING)

        cls.protCTF = cls.newProtocol(ProtImportCTF,
                                      importFrom=ProtImportCTF.IMPORT_FROM_XMIPP3,
                                      filesPath=cls.dataset.getFile('ctfsDir'),
                                      filesPattern='*.ctfparam')
        cls.protCTF.inputMicrographs.set(cls.protImport.outputMicrographs)
        cls.proj.launchProtocol(cls.protCTF, wait=True)

        cls.protPP = cls.runFakedPicking(cls.protDown.outputMicrographs, cls.allCrdsDir)

    def _checkSamplingConsistency(self, outputSet):
        """ Check that the set sampling is the same as item sampling. """
        first = outputSet.getFirstItem()

        self.assertAlmostEqual(outputSet.getSamplingRate(),
                               first.getSamplingRate())

    def _checkVarianceAndGiniCoeff(self, particle, varianceScore, giniScore):
        """ Check the Variance and Gini coeff. added to a certain particle
        """
        self.assertTrue(particle.hasAttribute('_xmipp_scoreByVariance'),
                        'Particle has not scoreByVariance attribute.')
        self.assertAlmostEqual(particle._xmipp_scoreByVariance.get(), varianceScore,
                               3, "The was a problem with the varianceScore")

        self.assertTrue(particle.hasAttribute('_xmipp_scoreByGiniCoeff'),
                        'Particle has not scoreByGiniCoeff attribute.')
        self.assertAlmostEqual(particle._xmipp_scoreByGiniCoeff.get(), giniScore,
                               3, "The was a problem with the giniCoeffScore")

    def testExtractSameAsPicking(self):
        print "Run extract particles from same micrographs as picking"
        protExtract = self.newProtocol(XmippProtExtractParticles,
                                       boxSize=110,
                                       downsampleType=SAME_AS_PICKING,
                                       doInvert=False,
                                       doFlip=False)
        protExtract.setObjLabel("extract-same as picking")
        protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
        self.launchProtocol(protExtract)


        inputCoords = protExtract.inputCoordinates.get()
        outputParts = protExtract.outputParticles
        micSampling = protExtract.inputCoordinates.get().getMicrographs().getSamplingRate()

        self.assertIsNotNone(outputParts,
                             "There was a problem generating the output.")
        self.assertAlmostEqual(outputParts.getSamplingRate()/micSampling,
                               1, 1,
                               "There was a problem generating the output.")
        self._checkSamplingConsistency(outputParts)

        def compare(objId, delta=0.001):
            cx, cy = inputCoords[objId].getPosition()
            px, py = outputParts[objId].getCoordinate().getPosition()
            micNameCoord = inputCoords[objId].getMicName()
            micNamePart = outputParts[objId].getCoordinate().getMicName()
            self.assertAlmostEquals(cx, px, delta=delta)
            self.assertAlmostEquals(cy, py, delta=delta)
            self.assertEqual(micNameCoord, micNamePart,
                             "The micName should be %s and its %s"
                             %(micNameCoord, micNamePart))
        compare(83)
        compare(228)
        self._checkVarianceAndGiniCoeff(outputParts[170], 1.1640, 0.5190)

    def testExtractOriginal(self):
        print "Run extract particles from the original micrographs"
        protExtract = self.newProtocol(XmippProtExtractParticles,
                                       boxSize=550,
                                       downsampleType=OTHER,
                                       doInvert=False,
                                       doFlip=False)
        protExtract.setObjLabel("extract-original")
        protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
        protExtract.inputMicrographs.set(self.protImport.outputMicrographs)
        self.launchProtocol(protExtract)

        inputCoords = protExtract.inputCoordinates.get()
        outputParts = protExtract.outputParticles
        samplingCoords = self.protPP.outputCoordinates.getMicrographs().getSamplingRate()
        samplingFinal = self.protImport.outputMicrographs.getSamplingRate()
        samplingMics = protExtract.inputMicrographs.get().getSamplingRate()
        factor = samplingFinal / samplingCoords

        def compare(objId, delta=1.0):
            cx, cy = inputCoords[objId].getPosition()
            px, py = outputParts[objId].getCoordinate().getPosition()
            micNameCoord = inputCoords[objId].getMicName()
            micNamePart = outputParts[objId].getCoordinate().getMicName()
            self.assertAlmostEquals(cx/factor, px, delta=delta)
            self.assertAlmostEquals(cy/factor, py, delta=delta)
            self.assertEqual(micNameCoord, micNamePart,
                             "The micName should be %s and its %s"
                             %(micNameCoord, micNamePart))

        compare(111)
        compare(7)        
        
        self.assertIsNotNone(outputParts,
                             "There was a problem generating the output.")
        self.assertEqual(outputParts.getSamplingRate(), samplingMics,
                         "Output sampling rate should be equal to input "
                         "sampling rate.")
        self._checkSamplingConsistency(outputParts)
        self._checkVarianceAndGiniCoeff(outputParts[170], 1.2081, 0.5754)

    def testNoExtractBorders(self):
        print "Run extract particles avoiding extract in borders"
        protExtract = self.newProtocol(XmippProtExtractParticles,
                                       boxSize=750,
                                       downsampleType=OTHER,
                                       doInvert=False,
                                       doBorders=False,
                                       doFlip=False)

        protExtract.setObjLabel("extract-avoid borders")
        protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
        protExtract.inputMicrographs.set(self.protImport.outputMicrographs)
        self.launchProtocol(protExtract)

        inputCoords = protExtract.inputCoordinates.get()
        outputParts = protExtract.outputParticles
        samplingCoords = self.protPP.outputCoordinates.getMicrographs().getSamplingRate()
        samplingFinal = self.protImport.outputMicrographs.getSamplingRate()
        samplingMics = protExtract.inputMicrographs.get().getSamplingRate()
        factor = samplingFinal / samplingCoords

        def compare(objId, delta=1.0):
            cx, cy = inputCoords[objId].getPosition()
            px, py = outputParts[objId].getCoordinate().getPosition()
            micNameCoord = inputCoords[objId].getMicName()
            micNamePart = outputParts[objId].getCoordinate().getMicName()
            self.assertAlmostEquals(cx / factor, px, delta=delta)
            self.assertAlmostEquals(cy / factor, py, delta=delta)
            self.assertEqual(micNameCoord, micNamePart,
                             "The micName should be %s and its %s"
                             % (micNameCoord, micNamePart))

        compare(111)
        compare(7)

        self.assertIsNotNone(outputParts,
                             "There was a problem generating the output.")
        self.assertEqual(outputParts.getSamplingRate(), samplingMics,
                         "Output sampling rate should be equal to input "
                         "sampling rate.")
        self.assertAlmostEquals(outputParts.getSize(), 399, delta=1)
        self._checkSamplingConsistency(outputParts)
        self._checkVarianceAndGiniCoeff(outputParts[170], 1.2120, 0.5275)

    def testExtractOther(self):
        print "Run extract particles from original micrographs, with downsampling"
        downFactor = 3.0
        protExtract = self.newProtocol(XmippProtExtractParticles,
                                       boxSize=183, downsampleType=OTHER,
                                       doDownsample=True,
                                       downFactor=downFactor,
                                       doInvert=False,
                                       doFlip=False)
        # Get all the micrographs ids to validate that all particles
        # has the micId properly set
        micsId = [mic.getObjId() for mic in
                  self.protPP.outputCoordinates.getMicrographs()]
        
        protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
        protExtract.inputMicrographs.set(self.protImport.outputMicrographs)
        protExtract.setObjLabel("extract-other")
        self.launchProtocol(protExtract)

        inputCoords = protExtract.inputCoordinates.get()
        outputParts = protExtract.outputParticles
        samplingCoords = self.protPP.outputCoordinates.getMicrographs().getSamplingRate()
        samplingFinal = self.protImport.outputMicrographs.getSamplingRate() * downFactor
        samplingMics = protExtract.inputMicrographs.get().getSamplingRate()
        factor = samplingFinal / samplingCoords
        self.assertIsNotNone(outputParts,
                             "There was a problem generating the output.")
        
        def compare(objId, delta=1.0):
            cx, cy = inputCoords[objId].getPosition()
            px, py = outputParts[objId].getCoordinate().getPosition()
            micNameCoord = inputCoords[objId].getMicName()
            micNamePart = outputParts[objId].getCoordinate().getMicName()
            self.assertAlmostEquals(cx/factor, px, delta=delta)
            self.assertAlmostEquals(cy/factor, py, delta=delta)
            self.assertEqual(micNameCoord, micNamePart,
                             "The micName should be %s and its %s"
                             %(micNameCoord, micNamePart))

        compare(45)
        compare(229)

        outputSampling = outputParts.getSamplingRate()
        self.assertAlmostEqual(outputSampling/samplingMics,
                               downFactor, 1,
                               "There was a problem generating the output.")
        for particle in outputParts:
            self.assertTrue(particle.getCoordinate().getMicId() in micsId)
            self.assertAlmostEqual(outputSampling, particle.getSamplingRate())
        self._checkVarianceAndGiniCoeff(outputParts[170], 1.2472, 0.6052)

    def testExtractNoise(self):
        # here we will try a different patchSize than the default
        print "Run extract particles from original micrographs, with downsampling"
        downFactor = 5.0
        protExtract = self.newProtocol(XmippProtExtractParticles,
                                       boxSize=183, downsampleType=OTHER,
                                       doDownsample=True,
                                       downFactor=downFactor,
                                       doInvert=False,
                                       doFlip=False,
                                       extractNoise=True,
                                       patchSize=500)

        protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
        protExtract.inputMicrographs.set(self.protImport.outputMicrographs)
        protExtract.setObjLabel("extract-noise")
        self.launchProtocol(protExtract)

        outputParts = protExtract.outputParticles
        self.assertIsNotNone(outputParts, "There was a problem generating the output.")
        self.assertAlmostEquals(outputParts.getSize(), 403, delta=1)
        self._checkVarianceAndGiniCoeff(outputParts[170], 1.1594, 0.5702)

    def testExtractCTF(self):
        print "Run extract particles with CTF"
        protExtract = self.newProtocol(XmippProtExtractParticles,
                                       boxSize=110,
                                       downsampleType=SAME_AS_PICKING,
                                       doInvert=False,
                                       doFlip=True)
        protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
        protExtract.ctfRelations.set(self.protCTF.outputCTF)
        protExtract.setObjLabel("extract-ctf")
        self.launchProtocol(protExtract)

        inputCoords = protExtract.inputCoordinates.get()
        outputParts = protExtract.outputParticles

        def compare(objId, delta=0.001):
            cx, cy = inputCoords[objId].getPosition()
            px, py = outputParts[objId].getCoordinate().getPosition()
            micNameCoord = inputCoords[objId].getMicName()
            micNamePart = outputParts[objId].getCoordinate().getMicName()
            self.assertAlmostEquals(cx, px, delta=delta)
            self.assertAlmostEquals(cy, py, delta=delta)
            self.assertEqual(micNameCoord, micNamePart,
                             "The micName should be %s and its %s"
                             %(micNameCoord, micNamePart))

        compare(228)
        compare(83)

        def compareCTF(partId, ctfId):
            partDefU = outputParts[partId].getCTF().getDefocusU()
            defU = protExtract.ctfRelations.get()[ctfId].getDefocusU()
            self.assertAlmostEquals(partDefU, defU, delta=1)

        compareCTF(1, 1)
        compareCTF(150, 2)
        compareCTF(300, 3)
        
        self.assertIsNotNone(outputParts,
                             "There was a problem generating the output.")
        self.assertTrue(outputParts.hasCTF(), "Output does not have CTF.")
        self._checkSamplingConsistency(outputParts)
        self._checkVarianceAndGiniCoeff(outputParts[170], 1.1640, 0.5190)


class TestXmippEliminatingEmptyParticles(TestXmippBase):
    """This class check if the protocol for eliminating
     empty particles in Xmipp works properly."""
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData()
        cls.protImport1 = cls.runImportMicrographBPV(cls.micsFn)
        cls.protDown1 = cls.runDownsamplingMicrographs(
            cls.protImport1.outputMicrographs, 5)
        cls.protPP = cls.runFakedPicking(cls.protDown.outputMicrographs,
                                         cls.allCrdsDir)

    def _updateProtocol(self, prot):
        prot2 = getProtocolFromDb(prot.getProject().path,
                                  prot.getDbPath(),
                                  prot.getObjId())
        # Close DB connections
        prot2.getProject().closeMapper()
        prot2.closeMappers()
        return prot2

    def testStreamingAndNonStreaming(self):
        protExtract = self.newProtocol(XmippProtExtractParticles,
                                       boxSize=110,
                                       downsampleType=SAME_AS_PICKING,
                                       doInvert=False,
                                       doFlip=False)
        protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
        self.launchProtocol(protExtract)

        protElimination1 = self.newProtocol(XmippProtEliminateEmptyParticles)
        protElimination1.inputParticles.set(protExtract.outputParticles)
        self.launchProtocol(protElimination1)

        outSet = SetOfParticles(
            filename=protElimination1._getPath('outputParticles.sqlite'))
        elimSet = SetOfParticles(
            filename=protElimination1._getPath('eliminatedParticles.sqlite'))

        self.assertTrue(outSet.getSize() + elimSet.getSize() ==
                        protExtract.outputParticles.getSize(),
                        "Output sets size does not much the input set size.")

        kwargs = {'nDim': 20,  # 20 objects/particles
                  'creationInterval': 1,  # wait 1 sec. after creation
                  'setof': 3,  # create SetOfParticles
                  'inputParticles': protExtract.outputParticles
                  }

        protStream = self.newProtocol(ProtCreateStreamData, **kwargs)
        self.proj.launchProtocol(protStream, wait=False)

        while not protStream.hasAttribute('outputParticles'):
            time.sleep(1)
            protStream = self._updateProtocol(protStream)

        protElimination2 = self.newProtocol(XmippProtEliminateEmptyParticles)
        protElimination2.inputParticles.set(protStream.outputParticles)
        self.launchProtocol(protElimination2)

        partSet = SetOfParticles(
            filename=protStream._getPath("particles.sqlite"))
        outSet = SetOfParticles(
            filename=protElimination2._getPath('outputParticles.sqlite'))
        elimSet = SetOfParticles(
            filename=protElimination2._getPath('eliminatedParticles.sqlite'))
        self.assertTrue(outSet.getSize() + elimSet.getSize() ==
                        partSet.getSize(),
                        "Output sets size does not much the input set size.")


class TestXmippParticlesPickConsensus(TestXmippBase):
    """This class check if the protocol for particle
    picking consensus in Xmipp works properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData()
        cls.protImport = cls.runImportMicrographBPV(cls.micsFn)
        cls.protDown = cls.runDownsamplingMicrographs(
            cls.protImport.outputMicrographs, 5)
        cls.protFaPi = cls.runFakedPicking(cls.protDown.outputMicrographs,
                                           cls.allCrdsDir)

    def _updateProtocol(self, prot):
        prot2 = getProtocolFromDb(prot.getProject().path,
                                  prot.getDbPath(),
                                  prot.getObjId())
        # Close DB connections
        prot2.getProject().closeMapper()
        prot2.closeMappers()
        return prot2

    def testStreamingAndNonStreaming(self):
        protAutomaticPP = XmippParticlePickingAutomatic()
        protAutomaticPP.xmippParticlePicking.set(self.protFaPi)
        protAutomaticPP.inputMicrographs.set(self.protDown.outputMicrographs)
        protAutomaticPP.micsToPick.set(1)
        self.proj.launchProtocol(protAutomaticPP, wait=True)

        protCons1 = self.newProtocol(XmippProtConsensusPicking)
        protCons1.inputCoordinates.set([self.protFaPi.outputCoordinates,
                                        protAutomaticPP.outputCoordinates])
        self.launchProtocol(protCons1)

        self.assertTrue(protCons1.isFinished(), "Consensus failed")
        self.assertTrue(protCons1.consensusCoordinates.getSize() == 390,
                        "Output coordinates size does not is wrong.")

        kwargs = {'nDim': 3,  # 3 objects
                  'creationInterval': 20,  # wait 1 sec. after creation
                  'setof': 1,  # create SetOfMicrographs
                  'inputMics': self.protDown.outputMicrographs
                  }

        protStream = self.newProtocol(ProtCreateStreamData, **kwargs)
        self.proj.launchProtocol(protStream, wait=False)

        while not protStream.hasAttribute('outputMicrographs'):
            time.sleep(1)
            protStream = self._updateProtocol(protStream)

        protAutoPP = XmippParticlePickingAutomatic()
        protAutoPP.xmippParticlePicking.set(self.protFaPi)
        protAutoPP.micsToPick.set(1)
        protAutoPP.inputMicrographs.set(protStream.outputMicrographs)
        self.proj.launchProtocol(protAutoPP, wait=False)

        while not protAutoPP.hasAttribute('outputCoordinates'):
            time.sleep(1)
            protAutoPP = self._updateProtocol(protAutoPP)

        protCons2 = self.newProtocol(XmippProtConsensusPicking)
        protCons2.inputCoordinates.set([self.protFaPi.outputCoordinates,
                                        protAutoPP.outputCoordinates])
        self.launchProtocol(protCons2)
        self.assertTrue(protCons2.isFinished(), "Consensus failed")
        self.assertTrue(protCons2.consensusCoordinates.getSize() == 390,
                        "Output coordinates size does not is wrong.")


    # Sorting particles is not possible in streaming mode. Thus, all params
    # related with was removed from extract particle protocol. There exists
    # another protocol (screen particles) to do it.

    # def testExtractSort(self):
    #     print "Run extract particles with sort by statistics"
    #     protExtract = self.newProtocol(XmippProtExtractParticles,
    #                                    boxSize=110,
    #                                    downsampleType=SAME_AS_PICKING,
    #                                    doFlip=True, doSort=True,
    #                                    doInvert=False,
    #                                    rejectionMethod=1, maxZscore=2)
    #     protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
    #     protExtract.ctfRelations.set(self.protCTF.outputCTF)
    #     protExtract.setObjLabel("extract-sort")
    #     self.launchProtocol(protExtract)
    #
    #     inputCoords = protExtract.inputCoordinates.get()
    #     outputParts = protExtract.outputParticles
    #
    #     def compare(objId, delta=0.001):
    #         cx, cy = inputCoords[objId].getPosition()
    #         px, py = outputParts[objId].getCoordinate().getPosition()
    #         micNameCoord = inputCoords[objId].getMicName()
    #         micNamePart = outputParts[objId].getCoordinate().getMicName()
    #         self.assertAlmostEquals(cx, px, delta=delta)
    #         self.assertAlmostEquals(cy, py, delta=delta)
    #         self.assertEqual(micNameCoord, micNamePart,
    #                          "The micName should be %s and its %s"
    #                          %(micNameCoord, micNamePart))
    #
    #     compare(228)
    #     compare(83)
    #
    #     self.assertIsNotNone(outputParts, "There was a problem generating"
    #                                       " the output.")
    #     self.assertAlmostEquals(outputParts.getSize(), 267, delta=2)
    #     self._checkSamplingConsistency(outputParts)
    #
    # def testExtractSortSmall(self):
    #     print "Run extract small particles sort by statistics"
    #     downFactor = 3.0
    #     protExtract = self.newProtocol(XmippProtExtractParticles,
    #                                    boxSize=40,
    #                                    downsampleType=SAME_AS_PICKING,
    #                                    doDownsample=True,
    #                                    downFactor=downFactor,
    #                                    doFlip=True, doSort=True,
    #                                    doInvert=False,
    #                                    rejectionMethod=1, maxZscore=2)
    #     protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
    #     protExtract.ctfRelations.set(self.protCTF.outputCTF)
    #     protExtract.setObjLabel("extract-sort small")
    #     self.launchProtocol(protExtract)
    #
    #     inputCoords = protExtract.inputCoordinates.get()
    #     outputParts = protExtract.outputParticles
    #
    #     def compare(objId, delta=1):
    #         cx, cy = inputCoords[objId].getPosition()
    #         px, py = outputParts[objId].getCoordinate().getPosition()
    #         micNameCoord = inputCoords[objId].getMicName()
    #         micNamePart = outputParts[objId].getCoordinate().getMicName()
    #         self.assertAlmostEquals(cx/downFactor, px, delta=delta)
    #         self.assertAlmostEquals(cy/downFactor, py, delta=delta)
    #         self.assertEqual(micNameCoord, micNamePart,
    #                          "The micName should be %s and its %s"
    #                          %(micNameCoord, micNamePart))
    #
    #     compare(228)
    #     compare(82)
    #
    #     self.assertIsNotNone(outputParts, "There was a problem generating"
    #                                       " the output.")
    #     self.assertAlmostEquals(outputParts.getSize(), 280, delta=2)
    #     self._checkSamplingConsistency(outputParts)


    # We changed the behaviour for extract particles in streaming mode,
    # and, from now, extract particles from a SetOfMicrographs with a
    # different micName but same ids is not supported

    # def testAssignCTF(self):
    #     """ Test the particle extraction after importing another
    #     SetOfMicrographs with a different micName but same ids.
    #     We will use assign-ctf protocol and extract from the
    #     newly imported mics with the assigned CTF.
    #     For the other mics, we will just create symbolic links.
    #     """
    #     # Create the links with a different micrograph name
    #     micsPath = self.proj.getPath('otherMicrographs')
    #     pwutils.makePath(micsPath)
    #     for i in [6, 7, 8]:
    #         micPath = self.dataset.getFile('micrographs/BPV_138%d.mrc' % i)
    #         micLink = join(micsPath, basename(micPath).replace('.mrc', '_DW.mrc'))
    #         pwutils.createAbsLink(micPath, micLink)
    #
    #     protImportDW = self.proj.copyProtocol(self.protImport)
    #     protImportDW.setObjLabel('import -mics DW')
    #     protImportDW.filesPath.set(os.path.abspath(micsPath))
    #     protImportDW.filesPattern.set('*_DW.mrc')
    #     self.launchProtocol(protImportDW)
    #
    #     protAssignCTF = self.newProtocol(ProtCTFAssign)
    #     protAssignCTF.inputSet.set(protImportDW.outputMicrographs)
    #     protAssignCTF.inputCTF.set(self.protCTF.outputCTF)
    #     self.launchProtocol(protAssignCTF)
    #     downFactor = 3.0
    #
    #     protExtract = self.newProtocol(XmippProtExtractParticles,
    #                                    boxSize=183, downsampleType=OTHER,
    #                                    doDownsample=True,
    #                                    downFactor=downFactor,
    #                                    doInvert=False,
    #                                    doFlip=False)
    #     protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
    #     protExtract.inputMicrographs.set(protAssignCTF.outputMicrographs)
    #     protExtract.ctfRelations.set(self.protCTF.outputCTF)
    #     protExtract.setObjLabel("extract-other (DW mics)")
    #     self.launchProtocol(protExtract)
    #
    #     inputCoords = protExtract.inputCoordinates.get()
    #     outputParts = protExtract.outputParticles
    #     samplingCoords = inputCoords.getMicrographs().getSamplingRate()
    #     samplingFinal = protImportDW.outputMicrographs.getSamplingRate() * downFactor
    #     samplingMics = protExtract.inputMicrographs.get().getSamplingRate()
    #     factor = samplingFinal / samplingCoords
    #     self.assertIsNotNone(outputParts, "There was a problem generating the output.")
    #
    #     def compare(objId, delta=1.0):
    #         cx, cy = inputCoords[objId].getPosition()
    #         px, py = outputParts[objId].getCoordinate().getPosition()
    #         micNameCoord = inputCoords[objId].getMicName()
    #         micNamePart = outputParts[objId].getCoordinate().getMicName()
    #         self.assertAlmostEquals(cx / factor, px, delta=delta)
    #         self.assertAlmostEquals(cy / factor, py, delta=delta)
    #         self.assertEqual(micNameCoord, micNamePart,
    #                          "The micName should be %s and its %s"
    #                          %(micNameCoord, micNamePart))
    #
    #     compare(45)
    #     compare(229)
    #
    #     self.assertAlmostEqual(outputParts.getSamplingRate() / samplingMics,
    #                            downFactor, 1)
    #     self._checkSamplingConsistency(outputParts)