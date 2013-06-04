
import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3 import *


class TestXmippTiltedWorkflow(unittest.TestCase):
    
    def setUp(self):    
        setupProject(self, 'TestProject_XmippTilted')
        
        self.pattern = getInputPath('Micrographs_TiltedPhantom', '*.mrc')        
        self.importFolder = getInputPath('Picking_XmippBPV3')

    def testWorkflow(self):
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
        #protExtract.inputMicrographs.set(protDownsampling.outputMicrographs)
        self.proj.launchProtocol(protExtract, wait=True)
        
        
if __name__ == "__main__":
#    suite = unittest.TestLoader().loadTestsFromTestCase(TestXmippTiltedWorkflow)    
#    unittest.TextTestRunner(verbosity=2).run(suite)
    unittest.main()