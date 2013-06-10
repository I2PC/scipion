import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3 import *
    
    
class TestXmippWorkflow(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):    
        # Create a new project
        setupProject(cls)
        cls.pattern = getInputPath('Micrographs_BPV3', '*.mrc')        
        cls.importFolder = getInputPath('Picking_XmippBPV3')

    def validateFiles(self, prot, filesSet):
        """ Validate if the produced files are the expected ones.
        Params:
            prot: the protocol to validate. 
            filesSet: the known files that should be produced (set)
        """
        self.assertEqual(prot.getFiles(), filesSet)
        
    def testXmippWorkflow(self):
        #First, import a set of micrographs
        protImport = ProtImportMicrographs(pattern=self.pattern, samplingRate=1.237, voltage=300)
        self.proj.launchProtocol(protImport, wait=True)
        
        self.assertIsNotNone(protImport.outputMicrographs, "There was a problem with the import")
        
        # Perform a downsampling on the micrographs

        print "Downsampling..."
        protDownsampling = XmippProtPreprocessMicrographs(doDownsample=True, downFactor=3, doCrop=False)
        protDownsampling.inputMicrographs.set(protImport.outputMicrographs)
        self.proj.launchProtocol(protDownsampling, wait=True)
          
        self.assertIsNotNone(protDownsampling.outputMicrographs, "There was a problem with the downsampling")
          
        # Now estimate CTF on the downsampled micrographs 
        print "Performing CTF..."   
        protCTF = XmippProtCTFMicrographs(numberOfThreads=3)                
        protCTF.inputMicrographs.set(protDownsampling.outputMicrographs)        
        self.proj.launchProtocol(protCTF, wait=True)
        
        print "Running fake particle picking..."   
        protPP = XmippProtParticlePicking(importFolder=self.importFolder)                
        protPP.inputMicrographs.set(protCTF.outputMicrographs)        
        self.proj.launchProtocol(protPP, wait=True)
            
        self.assertIsNotNone(protPP.outputCoordinates, "There was a problem with the faked picking")
            
        print "Run extract particles with Same as picking"
        protExtract = XmippProtExtractParticles(boxSize=64, downsampleType=1)
        protExtract.inputCoordinates.set(protPP.outputCoordinates)
        #protExtract.inputMicrographs.set(protDownsampling.outputMicrographs)
        self.proj.launchProtocol(protExtract, wait=True)
        
        self.assertIsNotNone(protExtract.outputImages, "There was a problem with the extract particles")
        
        print "Run ML2D"
        protML2D = XmippProtML2D(numberOfReferences=1, maxIters=4, 
                                 numberOfMpi=2, numberOfThreads=1)
        protML2D.inputImages.set(protExtract.outputImages)
        self.proj.launchProtocol(protML2D, wait=True)        
        
        self.assertIsNotNone(protML2D.outputClassification, "There was a problem with ML2D")  
        
        print "Run CL2D"
        protCL2D = XmippProtCL2D(numberOfReferences=2, numberOfInitialReferences=1, 
                                 numberOfIterations=4, numberOfMpi=2)
        protCL2D.inputImages.set(protExtract.outputImages)
        self.proj.launchProtocol(protCL2D, wait=True)        
        
        self.assertIsNotNone(protCL2D.outputClassification, "There was a problem with CL2D")         
        
if __name__ == "__main__":
    unittest.main()