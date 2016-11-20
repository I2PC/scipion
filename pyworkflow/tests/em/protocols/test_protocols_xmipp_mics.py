# **************************************************************************
# *
# * Authors:    Laura del Cano (ldelcano@cnb.csic.es)
# *             Josue Gomez Blanco (jgomez@cnb.csic.es)
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
    def runImportMicrograph(cls, pattern, samplingRate, voltage, scannedPixelSize, magnification, sphericalAberration):
        """ Run an Import micrograph protocol. """
        
        # We have two options: pass the SamplingRate or the ScannedPixelSize + microscope magnification
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
        self.assertEquals(protDown.outputMicrographs.getSamplingRate(), self.protImport.outputMicrographs.getSamplingRate()*downFactorValue, "Micrographs uncorrectly downsampled")
    
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
        self.assertAlmostEquals(ctfModel.getDefocusU(),23825.9, delta=500)
        self.assertAlmostEquals(ctfModel.getDefocusV(),23520.3, delta=500)
        self.assertAlmostEquals(ctfModel.getDefocusAngle(), 49.2882, delta=5)
        sampling = ctfModel.getMicrograph().getSamplingRate() * protCTF.ctfDownFactor.get()
        self.assertAlmostEquals(sampling, 2.474, delta=0.001)


# class TestXmippCTFRestimation(TestXmippBase):
#     """This class check if the protocol to determine the CTF in Xmipp works properly."""
#     @classmethod
#     def setUpClass(cls):
#         setupTestProject(cls)
#         TestXmippBase.setData()
#         cls.protImport = cls.runImportMicrographBPV(cls.micsFn)
#         cls.protDown = cls.runDownsamplingMicrographs(cls.protImport.outputMicrographs, 3, 4)
#     
#     def testCTF(self):
#         # Estimate CTF on the downsampled micrographs
#         print "Performing CTF..."
#         protCTF = XmippProtCTFMicrographs(numberOfThreads=4)
#         protCTF.inputMicrographs.set(self.protDown.outputMicrographs)        
#         self.proj.launchProtocol(protCTF, wait=True)
#         self.assertIsNotNone(protCTF.outputCTF, "SetOfCTF has not been produced.")
#         
        print "Performing CTF Recalculation..."
#         str = "1,22000,24000,0,0.05,0.26; 3,21000,23000,0,0.04,0.3"
#         protReCTF = XmippProtRecalculateCTF(numberOfThreads=3, numberOfMpi=1)
#         protReCTF.inputCtf.set(protCTF.outputCTF)
#         protReCTF.inputValues.set(str)
#         self.proj.launchProtocol(protReCTF, wait=True)
#         self.assertIsNotNone(protReCTF.outputCTF, "SetOfCTF has not been produced in CTF Recalculation.")


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
    """This class check if the protocol to extract particles in Xmipp works properly."""
    
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
        def compare(objId, delta=0.001):
            cx, cy = inputCoords[objId].getPosition()
            px, py = outputParts[objId].getCoordinate().getPosition()
            self.assertAlmostEquals(cx, px, delta=delta)
            self.assertAlmostEquals(cy, py, delta=delta)

        compare(228)
        compare(83)
    
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
            self.assertAlmostEquals(cx/factor, px, delta=delta)
            self.assertAlmostEquals(cy/factor, py, delta=delta)
             
        compare(111)
        compare(7)        
        
        self.assertIsNotNone(outputParts, "There was a problem generating the output.")
        self.assertEqual(outputParts.getSamplingRate(), samplingMics, 
                         "Output sampling rate should be equal to input sampling rate.")


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
        micsId = [mic.getObjId() for mic in self.protPP.outputCoordinates.getMicrographs()]
        
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
        self.assertIsNotNone(outputParts, "There was a problem generating the output.")
        
        def compare(objId, delta=1.0):
            cx, cy = inputCoords[objId].getPosition()
            px, py = outputParts[objId].getCoordinate().getPosition()
            self.assertAlmostEquals(cx/factor, px, delta=delta)
            self.assertAlmostEquals(cy/factor, py, delta=delta)
            
        compare(45)
        compare(229)  

        self.assertAlmostEqual(outputParts.getSamplingRate()/samplingMics,
                               downFactor, 1, "There was a problem generating the output.")
        for particle in outputParts:
            self.assertTrue(particle.getCoordinate().getMicId() in micsId)
    
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
            self.assertAlmostEquals(cx, px, delta=delta)
            self.assertAlmostEquals(cy, py, delta=delta)
            
        compare(228)
        compare(83) 
                        
        def compareCTF(partId, ctfId):
            partDefU = outputParts[partId].getCTF().getDefocusU()
            defU = protExtract.ctfRelations.get()[ctfId].getDefocusU()
            self.assertAlmostEquals(partDefU, defU, delta=1)
            
        compareCTF(1, 1)
        compareCTF(150, 2)
        compareCTF(300, 3)
        
        self.assertIsNotNone(outputParts, "There was a problem generating the output.")
        self.assertTrue(outputParts.hasCTF(), "Output does not have CTF.")
    
    def testExtractSort(self):
        print "Run extract particles with sort by statistics"
        protExtract = self.newProtocol(XmippProtExtractParticles, 
                                       boxSize=110, 
                                       downsampleType=SAME_AS_PICKING,
                                       doFlip=True, doSort=True,
                                       doInvert=False,
                                       rejectionMethod=1, maxZscore=2)
        protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
        protExtract.ctfRelations.set(self.protCTF.outputCTF)
        protExtract.setObjLabel("extract-sort")
        self.launchProtocol(protExtract)
 
        inputCoords = protExtract.inputCoordinates.get()
        outputParts = protExtract.outputParticles 
        
        def compare(objId, delta=0.001):
            cx, cy = inputCoords[objId].getPosition()
            px, py = outputParts[objId].getCoordinate().getPosition()
            self.assertAlmostEquals(cx, px, delta=delta)
            self.assertAlmostEquals(cy, py, delta=delta)
            
        compare(228)
        compare(83) 
               
        self.assertIsNotNone(outputParts, "There was a problem generating the output.")
        self.assertAlmostEquals(outputParts.getSize(), 249, delta=10)
        
    def testAssignCTF(self):
        """ Test the particle extraction after importing another
        SetOfMicrographs with a different micName but same ids.
        We will use assign-ctf protocol and extract from the 
        newly imported mics with the assigned CTF.
        For the other mics, we will just create symbolic links.
        """
        # Create the links with a different micrograph name
        micsPath = self.proj.getPath('otherMicrographs')
        pwutils.makePath(micsPath)
        for i in [6, 7, 8]:
            micPath = self.dataset.getFile('micrographs/BPV_138%d.mrc' % i)
            micLink = join(micsPath, basename(micPath).replace('.mrc', '_DW.mrc'))
            pwutils.createAbsLink(micPath, micLink)
            
        protImportDW = self.proj.copyProtocol(self.protImport)
        protImportDW.setObjLabel('import -mics DW')
        protImportDW.filesPath.set(os.path.abspath(micsPath))
        protImportDW.filesPattern.set('*_DW.mrc')
        self.launchProtocol(protImportDW)
        
        protAssignCTF = self.newProtocol(ProtCTFAssign)
        protAssignCTF.inputSet.set(protImportDW.outputMicrographs)
        protAssignCTF.inputCTF.set(self.protCTF.outputCTF)
        self.launchProtocol(protAssignCTF)
        downFactor = 3.0
        
        protExtract = self.newProtocol(XmippProtExtractParticles, 
                                       boxSize=183, downsampleType=OTHER,
                                       doDownsample=True,
                                       downFactor=downFactor,
                                       doInvert=False,
                                       doFlip=False)
        protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
        protExtract.inputMicrographs.set(protAssignCTF.outputMicrographs)
        protExtract.ctfRelations.set(self.protCTF.outputCTF)
        protExtract.setObjLabel("extract-other (DW mics)")
        self.launchProtocol(protExtract)

        inputCoords = protExtract.inputCoordinates.get()
        outputParts = protExtract.outputParticles
        samplingCoords = inputCoords.getMicrographs().getSamplingRate()
        samplingFinal = protImportDW.outputMicrographs.getSamplingRate() * downFactor
        samplingMics = protExtract.inputMicrographs.get().getSamplingRate()
        factor = samplingFinal / samplingCoords
        self.assertIsNotNone(outputParts, "There was a problem generating the output.")

        def compare(objId, delta=1.0):
            cx, cy = inputCoords[objId].getPosition()
            px, py = outputParts[objId].getCoordinate().getPosition()
            self.assertAlmostEquals(cx / factor, px, delta=delta)
            self.assertAlmostEquals(cy / factor, py, delta=delta)

        compare(45)
        compare(229)

        self.assertAlmostEqual(outputParts.getSamplingRate() / samplingMics,
                               downFactor, 1)
