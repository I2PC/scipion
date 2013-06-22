import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3 import *
from pyworkflow.em.packages.brandeis import *
from pyworkflow.em.packages.eman2 import *
from test_workflow import TestWorkflow
    
    
class TestMixedWorkflow_1(TestWorkflow):
    
    @classmethod
    def setUpClass(cls):    
        # Create a new project
        setupProject(cls)
        cls.pattern = getInputPath('Micrographs_BPV1', '*.mrc')        
        cls.importFolder = getInputPath('Picking_XmippBPV1')
        
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
        print "Performing CTFfind..."   
        protCTF = ProtCTFFind(runMode=1)         
        protCTF.inputMicrographs.set(protDownsampling.outputMicrographs)        
        self.proj.launchProtocol(protCTF, wait=True)
        
        print "Running fake particle picking..."   
        protPP = XmippProtParticlePicking(importFolder=self.importFolder, runMode=1)                
        protPP.inputMicrographs.set(protCTF.outputMicrographs)        
        self.proj.launchProtocol(protPP, wait=True)
            
        self.assertIsNotNone(protPP.outputCoordinates, "There was a problem with the faked picking")
            
        print "Run extract particles with Same as picking"
        protExtract = XmippProtExtractParticles(boxSize=258, downsampleType=1, runMode=1)
        protExtract.inputCoordinates.set(protPP.outputCoordinates)
        #protExtract.inputMicrographs.set(protDownsampling.outputMicrographs)
        self.proj.launchProtocol(protExtract, wait=True)
        
        self.assertIsNotNone(protExtract.outputImages, "There was a problem with the extract particles")
               
        
if __name__ == "__main__":
    unittest.main()