
import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3 import *


class TestXmippTiltedWorkflow(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):    
        # Create a new project
        setupProject(cls)
        cls.pattern = getInputPath('Micrographs_TiltedPhantom', '*.mrc')        
        cls.importFolder = getInputPath('Picking_TiltedPhantom')

    def testXmippTiltedWorkflow(self):
        #First, import a set of micrographs
        protImport = ProtImportMicrographs(pattern=self.pattern, samplingRate=1.237, 
                                           voltage=300, tiltPairs=True)
        self.proj.launchProtocol(protImport, wait=True)
        
        self.assertTrue(protImport.outputMicrographs is not None, "There was a problem with the import")
        
        # Perform a downsampling on the micrographs

        print "Downsampling..."
        protDownsampling = XmippProtPreprocessMicrographs(doDownsample=True, downFactor=2)
        protDownsampling.inputMicrographs.set(protImport.outputMicrographs)
        self.proj.launchProtocol(protDownsampling, wait=True)
          
        self.assertTrue(protDownsampling.outputMicrographs is not None, "There was a problem with the downsampling")
        
        print "Running fake particle picking..."   
        protPP = XmippProtParticlePicking(importFolder=self.importFolder)
                
        protPP.inputMicrographs.set(protDownsampling.outputMicrographs)
        
        self.proj.launchProtocol(protPP, wait=True)
            
        self.assertTrue(protPP.outputCoordinates is not None, "There was a problem with the faked picking")
            
        print "Run extract particles with Same as picking"
        protExtract = XmippProtExtractParticles(boxSize=20, downsampleType=1)
        protExtract.inputCoordinates.set(protPP.outputCoordinates)

        self.proj.launchProtocol(protExtract, wait=True)
        
        self.assertIsNotNone(protExtract.outputImages, "There was a problem with the extract particles")
        
        print "Run ML2D"
        protML2D = XmippProtML2D(numberOfReferences=1, maxIters=4, numberOfMpi=1)
        protML2D.inputImages.set(protExtract.outputImages)
        self.proj.launchProtocol(protML2D, wait=True)        
        
        self.assertIsNotNone(protML2D.outputClassification, "There was a problem with ML2D")  
        
#        print "Run CL2D"
#        protCL2D = XmippProtCL2D(numberOfReferences=2, numberOfInitialReferences=1, 
#                                 numberOfIterations=4, numberOfMpi=1)
#        protCL2D.inputImages.set(protExtract.outputImages)
#        self.proj.launchProtocol(protCL2D, wait=True)        
#        
#        self.assertIsNotNone(protCL2D.outputClassification, "There was a problem with CL2D")         
        
        
if __name__ == "__main__":
#    suite = unittest.TestLoader().loadTestsFromTestCase(TestXmippTiltedWorkflow)    
#    unittest.TextTestRunner(verbosity=2).run(suite)
    unittest.main()