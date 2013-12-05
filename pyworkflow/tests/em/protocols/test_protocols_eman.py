'''
Created on Jun 6, 2013

@author: laura
'''
import unittest, sys
from pyworkflow.tests import *
from pyworkflow.em import *
from pyworkflow.em.packages.xmipp3 import *
from pyworkflow.em.packages.eman2 import *

class TestEmanBoxing(unittest.TestCase):

    @classmethod
    def setUpClass(cls):    
        # Create a new project
        setupProject(cls)
        cls.pattern = getInputPath('Micrographs_BPV3', '*.mrc')        
        cls.importFolder = getInputPath('EmanTestProject2')
        
    def testCreateOutput(self):    
        #First, import a set of micrographs
        protImport = ProtImportMicrographs(pattern=self.pattern, samplingRate=1.237, voltage=300)
        self.proj.launchProtocol(protImport, wait=True)
        
        self.assertIsNotNone(protImport.outputMicrographs, "There was a problem with the import")
        
        print "Running Eman fake particle picking..."   
        protPP = EmanProtBoxing(importFolder=self.importFolder, runMode=1)                
#        protPP.inputMicrographs.set(protCTF.outputMicrographs)        
        protPP.inputMicrographs.set(protImport.outputMicrographs)
        protPP.boxSize.set(110)
        self.proj.launchProtocol(protPP, wait=True)
            
        self.assertIsNotNone(protPP.outputCoordinates, "There was a problem with the faked picking")

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestEmanBoxing)
    unittest.TextTestRunner(verbosity=2).run(suite)