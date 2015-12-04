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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3 import *


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
        # We have two options: passe the SamplingRate or the ScannedPixelSize + microscope magnification
        if not samplingRate is None:
            cls.protImport = ProtImportMicrographs(samplingRateMode=0, filesPath=pattern, samplingRate=samplingRate, magnification=magnification, 
                                                   voltage=voltage, sphericalAberration=sphericalAberration)
        else:
            cls.protImport = ProtImportMicrographs(samplingRateMode=1, filesPath=pattern, scannedPixelSize=scannedPixelSize, 
                                                   voltage=voltage, magnification=magnification, sphericalAberration=sphericalAberration)
            
        cls.proj.launchProtocol(cls.protImport, wait=True)
        if cls.protImport.isFailed():
            raise Exception("Protocol has failed. Error: ", cls.protImport.getErrorMessage())
        # check that input micrographs have been imported (a better way to do this?)
        if cls.protImport.outputMicrographs is None:
            raise Exception('Import of micrograph: %s, failed. outputMicrographs is None.' % pattern)
        return cls.protImport
    
    @classmethod
    def runImportMicrographBPV(cls, pattern):
        """ Run an Import micrograph protocol. """
        return cls.runImportMicrograph(pattern, samplingRate=1.237, voltage=300, sphericalAberration=2, scannedPixelSize=None, magnification=56000)
    
    @classmethod
    def runDownsamplingMicrographs(cls, mics, downFactorValue, threads=1):
        # test downsampling a set of micrographs
        cls.protDown = XmippProtPreprocessMicrographs(doDownsample=True, downFactor=downFactorValue, numberOfThreads=threads)
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
        self.assertIsNotNone(protAutomaticPP.outputCoordinates, "There was a problem with the automatic particle picking")
    
    def testAutomaticPickingOther(self):
        print "Run automatic particle picking"
        protAutomaticPP = XmippParticlePickingAutomatic()
        protAutomaticPP.xmippParticlePicking.set(self.protPP)
        protAutomaticPP.inputMicrographs.set(self.protDown2.outputMicrographs)
        protAutomaticPP.micsToPick.set(1)
        self.proj.launchProtocol(protAutomaticPP, wait=True)
        self.assertIsNotNone(protAutomaticPP.outputCoordinates, "There was a problem with the automatic particle picking")


class TestXmippExtractParticles(TestXmippBase):
    """This class check if the protocol to extract particles in Xmipp works properly."""
    
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData()
        cls.DOWNSAMPLING = 5.0
        cls.protImport = cls.runImportMicrographBPV(cls.micsFn)
        cls.protDown = cls.runDownsamplingMicrographs(cls.protImport.outputMicrographs, cls.DOWNSAMPLING)
        
        cls.protCTF = cls.newProtocol(ProtImportCTF,
                                 importFrom=ProtImportCTF.IMPORT_FROM_XMIPP3,
                                 filesPath=cls.dataset.getFile('ctfsDir'),
                                 filesPattern='*.ctfparam')
        cls.protCTF.inputMicrographs.set(cls.protImport.outputMicrographs)
        cls.proj.launchProtocol(cls.protCTF, wait=True)
         
        cls.protPP = cls.runFakedPicking(cls.protDown.outputMicrographs, cls.allCrdsDir)
    
    def testExtractSameAsPicking(self):
        print "Run extract particles with downsampling factor equal to the one at picking"
        protExtract = XmippProtExtractParticles(boxSize=110, downsampleType=SAME_AS_PICKING, doFlip=False)
        protExtract.inputMicrographs.set(self.protImport.outputMicrographs)
        protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
        protExtract.setObjLabel("extract-same as picking")
        self.proj.launchProtocol(protExtract, wait=True)
        partCrdX1 = protExtract.outputParticles[228].getCoordinate().getX()
        crdX1 = protExtract.inputCoordinates.get()[228].getX()
        partCrdY1 = protExtract.outputParticles[228].getCoordinate().getY()
        crdY1 = protExtract.inputCoordinates.get()[228].getY()
        partCrdX2 = protExtract.outputParticles[83].getCoordinate().getX()
        crdX2 = protExtract.inputCoordinates.get()[83].getX()
        partCrdY2 = protExtract.outputParticles[83].getCoordinate().getY()
        crdY2 = protExtract.inputCoordinates.get()[83].getY()
        self.assertAlmostEquals(partCrdX1, crdX1, delta=0.001)
        self.assertAlmostEquals(partCrdY1, crdY1, delta=0.001)
        self.assertAlmostEquals(partCrdX2, crdX2, delta=0.001)
        self.assertAlmostEquals(partCrdY2, crdY2, delta=0.001)
        self.assertIsNotNone(protExtract.outputParticles, "There was a problem generating the output.")
        self.assertAlmostEqual(protExtract.outputParticles.getSamplingRate()/protExtract.inputMicrographs.get().getSamplingRate(), self.DOWNSAMPLING, 1, "There was a problem generating the output.")
        
    
    def testExtractOriginal(self):
        print "Run extract particles with downsampling factor equal to the original micrographs"
        protExtract = XmippProtExtractParticles(boxSize=550, downsampleType=ORIGINAL, doFlip=False)
        protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
        protExtract.inputMicrographs.set(self.protImport.outputMicrographs)
        protExtract.setObjLabel("extract-original")
        self.proj.launchProtocol(protExtract, wait=True)
        
        samplingCoords = self.protPP.outputCoordinates.getMicrographs().getSamplingRate()
        samplingFinal = self.protImport.outputMicrographs.getSamplingRate()
        factor = samplingFinal / samplingCoords
        
        partCrdX1 = protExtract.outputParticles[111].getCoordinate().getX()
        crdX1 = int(protExtract.inputCoordinates.get()[111].getX() / factor)
        partCrdY1 = protExtract.outputParticles[111].getCoordinate().getY()
        crdY1 = int(protExtract.inputCoordinates.get()[111].getY() / factor)
        partCrdX2 = protExtract.outputParticles[7].getCoordinate().getX()
        crdX2 = int(protExtract.inputCoordinates.get()[7].getX() / factor)
        partCrdY2 = protExtract.outputParticles[7].getCoordinate().getY()
        crdY2 = int(protExtract.inputCoordinates.get()[7].getY() / factor)
        self.assertAlmostEquals(partCrdX1, crdX1, delta=0.001)
        self.assertAlmostEquals(partCrdY1, crdY1, delta=0.001)
        self.assertAlmostEquals(partCrdX2, crdX2, delta=0.001)
        self.assertAlmostEquals(partCrdY2, crdY2, delta=0.001)

        self.assertIsNotNone(protExtract.outputParticles, "There was a problem generating the output.")
        self.assertEqual(protExtract.outputParticles.getSamplingRate(), protExtract.inputMicrographs.get().getSamplingRate(), "Output sampling rate should be equal to input sampling rate.")

    def testExtractOther(self):
        print "Run extract particles with downsampling factor equal to other"
        downFactor=3.0
        protExtract = XmippProtExtractParticles(boxSize=183, downsampleType=OTHER, downFactor=downFactor,doFlip=False)
        # Get all the micrographs ids to validate that all particles
        # has the micId properly set
        micsId = [mic.getObjId() for mic in self.protPP.outputCoordinates.getMicrographs()]
        
        protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
        protExtract.inputMicrographs.set(self.protImport.outputMicrographs)
        protExtract.setObjLabel("extract-other")
        self.proj.launchProtocol(protExtract, wait=True)
        self.assertIsNotNone(protExtract.outputParticles, "There was a problem generating the output.")
        
        samplingCoords = self.protPP.outputCoordinates.getMicrographs().getSamplingRate()
        samplingFinal = self.protImport.outputMicrographs.getSamplingRate() * downFactor
        factor = samplingFinal / samplingCoords
        partCrdX1 = protExtract.outputParticles[45].getCoordinate().getX()
        crdX1 = int(protExtract.inputCoordinates.get()[45].getX() / factor)
        partCrdY1 = protExtract.outputParticles[45].getCoordinate().getY()
        crdY1 = int(protExtract.inputCoordinates.get()[45].getY() / factor)
        partCrdX2 = protExtract.outputParticles[229].getCoordinate().getX()
        crdX2 = int(protExtract.inputCoordinates.get()[229].getX() / factor)
        partCrdY2 = protExtract.outputParticles[229].getCoordinate().getY()
        crdY2 = int(protExtract.inputCoordinates.get()[229].getY() / factor)
        self.assertAlmostEquals(partCrdX1, crdX1, delta=0.001)
        self.assertAlmostEquals(partCrdY1, crdY1, delta=0.001)
        self.assertAlmostEquals(partCrdX2, crdX2, delta=0.001)
        self.assertAlmostEquals(partCrdY2, crdY2, delta=0.001)

        self.assertAlmostEqual(protExtract.outputParticles.getSamplingRate()/ protExtract.inputMicrographs.get().getSamplingRate(), downFactor, 1, "There was a problem generating the output.")
        for particle in protExtract.outputParticles:
            self.assertTrue(particle.getCoordinate().getMicId() in micsId)
    
    def testExtractCTF(self):
        print "Run extract particles with CTF"#        
        protExtract = XmippProtExtractParticles(boxSize=110, downsampleType=SAME_AS_PICKING,doFlip=True)
        protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
        protExtract.inputMicrographs.set(self.protCTF.inputMicrographs.get())
        protExtract.ctfRelations.set(self.protCTF.outputCTF)
        protExtract.setObjLabel("extract-ctf")
        self.proj.launchProtocol(protExtract, wait=True)
        partDef1 = protExtract.outputParticles[1].getCTF().getDefocusU()
        defU1 = protExtract.ctfRelations.get()[1].getDefocusU()
        partDef2 = protExtract.outputParticles[150].getCTF().getDefocusU()
        defU2 = protExtract.ctfRelations.get()[2].getDefocusU()
        partDef3 = protExtract.outputParticles[300].getCTF().getDefocusU()
        defU3 = protExtract.ctfRelations.get()[3].getDefocusU()
        partCrdX1 = protExtract.outputParticles[228].getCoordinate().getX()
        crdX1 = protExtract.inputCoordinates.get()[228].getX()
        partCrdY1 = protExtract.outputParticles[228].getCoordinate().getY()
        crdY1 = protExtract.inputCoordinates.get()[228].getY()
        partCrdX2 = protExtract.outputParticles[83].getCoordinate().getX()
        crdX2 = protExtract.inputCoordinates.get()[83].getX()
        partCrdY2 = protExtract.outputParticles[83].getCoordinate().getY()
        crdY2 = protExtract.inputCoordinates.get()[83].getY()
        self.assertAlmostEquals(partDef1, defU1, delta=1)
        self.assertAlmostEquals(partDef2, defU2, delta=1)
        self.assertAlmostEquals(partDef3, defU3, delta=1)
        self.assertAlmostEquals(partCrdX1, crdX1, delta=0.001)
        self.assertAlmostEquals(partCrdY1, crdY1, delta=0.001)
        self.assertAlmostEquals(partCrdX2, crdX2, delta=0.001)
        self.assertAlmostEquals(partCrdY2, crdY2, delta=0.001)
        self.assertIsNotNone(protExtract.outputParticles, "There was a problem generating the output.")
        self.assertTrue(protExtract.outputParticles.hasCTF(), "Output does not have CTF.")
    
    def testExtractSort(self):
        print "Run extract particles with sort by statistics"#
        protExtract = XmippProtExtractParticles(boxSize=110, downsampleType=SAME_AS_PICKING,doFlip=True,
                                                doSort=True,rejectionMethod=1, maxZscore=2)
        protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
        protExtract.inputMicrographs.set(self.protCTF.inputMicrographs.get())
        protExtract.ctfRelations.set(self.protCTF.outputCTF)
        protExtract.setObjLabel("extract-sort")
        self.proj.launchProtocol(protExtract, wait=True)
        partCrdX1 = protExtract.outputParticles[228].getCoordinate().getX()
        crdX1 = protExtract.inputCoordinates.get()[228].getX()
        partCrdY1 = protExtract.outputParticles[228].getCoordinate().getY()
        crdY1 = protExtract.inputCoordinates.get()[228].getY()
        partCrdX2 = protExtract.outputParticles[83].getCoordinate().getX()
        crdX2 = protExtract.inputCoordinates.get()[83].getX()
        partCrdY2 = protExtract.outputParticles[83].getCoordinate().getY()
        crdY2 = protExtract.inputCoordinates.get()[83].getY()
        self.assertAlmostEquals(partCrdX1, crdX1, delta=0.001)
        self.assertAlmostEquals(partCrdY1, crdY1, delta=0.001)
        self.assertAlmostEquals(partCrdX2, crdX2, delta=0.001)
        self.assertAlmostEquals(partCrdY2, crdY2, delta=0.001)
        self.assertIsNotNone(protExtract.outputParticles, "There was a problem generating the output.")
        self.assertAlmostEquals(protExtract.outputParticles.getSize(), 249, delta=10)