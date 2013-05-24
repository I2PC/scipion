import sys

import unittest, sys
from pyworkflow.manager import Manager
from pyworkflow.project import *
from pyworkflow.em import *
from pyworkflow.em.packages.xmipp3 import *

class TestNormalWorkflow(unittest.TestCase):
    
    def setUp(self):    
        # Create or load project
        projName = 'tests'
        self.pattern = '/home/laura/Scipion_Projects/InputData/*.mrc'
        self.importFolder ='/home/laura/Scipion_Projects/PosData' 
        
        manager = Manager()
        self.proj = manager.createProject(projName) # Now it will be loaded if exists

    def testXmippWorkflow(self):
        #First, import a set of micrographs
        protImport = ProtImportMicrographs(pattern=self.pattern, samplingRate=1.237, voltage=300)
        self.proj.launchProtocol(protImport, wait=True)
        
        self.assertTrue(protImport.outputMicrographs is not None, "There was a problem with the import")
        
        # Perform a downsampling on the micrographs

        print "Downsampling..."
        protDownsampling = XmippProtPreprocessMicrographs(doDownsample=True, downFactor=2, doCrop=False)
        protDownsampling.inputMicrographs.set(protImport.outputMicrographs)
        self.proj.launchProtocol(protDownsampling, wait=True)
          
        self.assertTrue(protDownsampling.outputMicrographs is not None, "There was a problem with the downsampling")
          
        # Now estimate CTF on the downsampled micrographs    

        print "Performing CTF..."   
        protCTF = XmippProtCTFMicrographs()
                
        protCTF.inputMicrographs.set(protDownsampling.outputMicrographs)
        
        self.proj.launchProtocol(protCTF, wait=True)
        
        print "Running fake particle picking..."   
        protPP = XmippProtParticlePicking(importFolder=self.importFolder)
                
        protPP.inputMicrographs.set(protCTF.outputMicrographs)
        
        self.proj.launchProtocol(protPP, wait=True)
            
        self.assertTrue(protPP.outputCoordinates is not None, "There was a problem with the faked picking")
            
        print "Run extract particles with Same as picking"
        protExtract = XmippProtExtractParticles(boxSize=258, downsampleType=1)
        protExtract.inputCoordinates.set(protPP.outputCoordinates)
        #protExtract.inputMicrographs.set(protDownsampling.outputMicrographs)
        self.proj.launchProtocol(protExtract, wait=True)
        
class TestXmippTiltedWorkflow(unittest.TestCase):
    
    def setUp(self):    
        # Create or load project
        projName = 'tests'
        self.pattern = '/home/laura/Scipion_Projects/TiltedData/*.mrc'
        self.importFolder ='/home/laura/Scipion_Projects/PosTiltedData' 
        
        manager = Manager()
        self.proj = manager.createProject(projName) # Now it will be loaded if exists

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
    suite = unittest.TestLoader().loadTestsFromTestCase(TestXmippTiltedWorkflow)    
    unittest.TextTestRunner(verbosity=2).run(suite)