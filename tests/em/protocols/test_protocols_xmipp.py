'''
Created on May 8, 2013

@author: laura
'''
import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3 import *


# Some utility functions to import micrographs that are used
# in several tests.
class TestXmippBase(unittest.TestCase):
    
    @classmethod
    def runImportMicrograph(cls, pattern, samplingRate, voltage):
        """ Run an Import micrograph protocol. """
        cls.protImport = ProtImportMicrographs(pattern=pattern, samplingRate=samplingRate, voltage=voltage)
        cls.proj.launchProtocol(cls.protImport, wait=True)
        # check that input micrographs have been imported (a better way to do this?)
        if cls.protImport.outputMicrographs is None:
            raise Exception('Import of micrograph: %s, failed. outputMicrographs is None.' % pattern)
        return cls.protImport
        
    @classmethod
    def runImportMicrographBPV1(cls):
        """ Run an Import micrograph protocol. """
        pattern = getInputPath('Micrographs_BPV1', '*.mrc')
        return cls.runImportMicrograph(pattern, samplingRate=1.237, voltage=300)
    
    @classmethod
    def runFakedPicking(cls, mics, pattern):
        """ Run a faked particle picking. Coordinates already existing. """
        coordsFolder = getInputPath(pattern)
        cls.protPP = XmippProtParticlePicking(importFolder=coordsFolder)                
        cls.protPP.inputMicrographs.set(mics)        
        cls.proj.launchProtocol(cls.protPP, wait=True)
        # check that faked picking has run ok
        if cls.protPP.outputCoordinates is None:
            raise Exception('Faked particle picking: %s, failed. outputCoordinates is None.' % coordsFolder)
        return cls.protPP       

    @classmethod
    def runImportParticles(cls, pattern, samplingRate):
        """ Run an Import particles protocol. """
        cls.protImport = ProtImportParticles(pattern=pattern, samplingRate=samplingRate)
        cls.proj.launchProtocol(cls.protImport, wait=True)
        # check that input images have been imported (a better way to do this?)
        if cls.protImport.outputImages is None:
            raise Exception('Import of images: %s, failed. outputImages is None.' % pattern)
        return cls.protImport        
 
class TestXmippPreprocessMicrographs(TestXmippBase):

    @classmethod
    def setUpClass(cls):
        setupProject(cls)
        cls.protImport = cls.runImportMicrographBPV1()

    def testDownsampling(self):
        # test downsampling a set of micrographs
        downFactor = 2
        protDown = XmippProtPreprocessMicrographs(doDownsample=True, downFactor=downFactor)
        protDown.inputMicrographs.set(self.protImport.outputMicrographs)
        self.proj.launchProtocol(protDown, wait=True)
        
        # check that output micrographs have double sampling rate than input micrographs
        self.assertTrue(protDown.outputMicrographs.samplingRate.get() == self.protImport.outputMicrographs.samplingRate.get()*downFactor, "Micrographs uncorrectly downsampled")
        
    def testCrop(self):
        # test crop on a set of micrographs
        cropPixels = 100
        protCrop = XmippProtPreprocessMicrographs(doCrop=True, cropPixels=cropPixels)
        protCrop.inputMicrographs.set(self.protImport.outputMicrographs)
        self.proj.launchProtocol(protCrop, wait=True)
        
        
class TestXmippCTFEstimation(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupProject(cls)
            
    def doCTF(self, pattern):
        #First, import a set of micrographs
        protImport = self.runImportMicrograph(pattern, samplingRate=3.711, voltage=300)
        
        # Now estimate CTF on the downsampled micrographs
        print "Performing CTF..."   
        protCTF = XmippProtCTFMicrographs()                
        protCTF.inputMicrographs.set(protImport.outputMicrographs)        
        self.proj.launchProtocol(protCTF, wait=True)
        
        self.assertTrue(protCTF.outputMicrographs.hasCTF(), "CTF estimation has not been performed.")
        
    def test_Micrographs_BPV1_Down3(self):
        self.doCTF(pattern = getInputPath('Micrographs_BPV1_Down3', '*.mrc'))
        
    def atest_Micrographs_BPV3_Down3(self):
        self.doCTF(pattern = getInputPath('Micrographs_BPV3_Down3', '*.mrc')) 
    
    
class TestXmippExtractParticles(TestXmippBase):
    
    SAME_AS_PICKING = 1
    ORIGINAL = 0
    OTHER = 2
    
    @classmethod
    def setUpClass(cls):
        setupProject(cls)    
        pattern = getInputPath('Micrographs_BPV1_Down3', '*.mrc')
        protImport = cls.runImportMicrograph(pattern, samplingRate=3.711, voltage=300)
        pattern = getInputPath('Micrographs_BPV1', '*.mrc')
        cls.protImport_ori = cls.runImportMicrograph(pattern, samplingRate=1.237, voltage=300)        
        cls.protPP = cls.runFakedPicking(protImport.outputMicrographs, 'Picking_XmippBPV1_Down3')
    
    def doExtract(self, boxSize, downsampleType):
        print "Run extract particles"
        protExtract = XmippProtExtractParticles(boxSize=boxSize, downsampleType=downsampleType)
        protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
        protExtract.inputMicrographs.set(self.protImport_ori.outputMicrographs)
        self.proj.launchProtocol(protExtract, wait=True)
        
        self.assertIsNotNone(protExtract.outputImages, "There was a problem with the extract particles")
        
    def testExtractSameAsPicking(self):
        print "Run extract particles with downsampling factor equal to the one at picking"
        protExtract = XmippProtExtractParticles(boxSize=171, downsampleType=self.SAME_AS_PICKING)
        protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
        self.proj.launchProtocol(protExtract, wait=True)
        
        self.assertIsNotNone(protExtract.outputImages, "There was a problem with the extract particles")
        
    def testExtractOriginal(self):
        print "Run extract particles with downsampling factor equal to the original micrographs"
        protExtract = XmippProtExtractParticles(boxSize=512, downsampleType=self.ORIGINAL)
        protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
        protExtract.inputMicrographs.set(self.protImport_ori.outputMicrographs)
        self.proj.launchProtocol(protExtract, wait=True)
        
        self.assertIsNotNone(protExtract.outputImages, "There was a problem with the extract particles")

    def testExtractOther(self):
        print "Run extract particles with downsampling factor equal to other"
        protExtract = XmippProtExtractParticles(boxSize=256, downsampleType=self.OTHER, downFactor=2)
        protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
        protExtract.inputMicrographs.set(self.protImport_ori.outputMicrographs)
        self.proj.launchProtocol(protExtract, wait=True)
        
        self.assertIsNotNone(protExtract.outputImages, "There was a problem with the extract particles")


     
def setupClassification(cls):
    """ Method to setup classification Test Cases. """
    setupProject(cls)
    #TODO: Find a set of images to make this work, with this it does not
    pattern = getInputPath('images_LTA', '*.xmp')
    cls.protImport = cls.runImportParticles(pattern=pattern, samplingRate=5.6)
    
       
class TestXmippML2D(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupClassification(cls)
        
    def testML2D(self):
        print "Run ML2D"
        protML2D = XmippProtML2D(numberOfReferences=2, maxIters=3, 
                                 numberOfMpi=2, numberOfThreads=2)
        protML2D.inputImages.set(self.protImport.outputImages)
        self.proj.launchProtocol(protML2D, wait=True)        
        
        self.assertIsNotNone(protML2D.outputClassification, "There was a problem with ML2D")  
        
class TestXmippCL2D(TestXmippBase):

    @classmethod
    def setUpClass(cls):
        setupClassification(cls)
        
    def testCL2D(self):
        print "Run CL2D"
        protCL2D = XmippProtCL2D(numberOfReferences=2, numberOfInitialReferences=1, 
                                 numberOfIterations=3, numberOfMpi=4)
        protCL2D.inputImages.set(self.protImport.outputImages)
        self.proj.launchProtocol(protCL2D, wait=True)        
        
        self.assertIsNotNone(protCL2D.outputClassification, "There was a problem with CL2D")  

        

if __name__ == "__main__":
    #suite = unittest.TestLoader().loadTestsFromTestCase(TestXmippExtractParticles)
    #unittest.TextTestRunner(verbosity=2).run(suite)
    unittest.main()
