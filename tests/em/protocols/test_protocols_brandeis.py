'''
Created on Jun 6, 2013

@author: laura
'''
import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.brandeis import *


class TestBrandeisCtffind(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):    
        # Create a new project for the tests
        setupProject(cls)
        
    @classmethod
    def runImportMicrograph(cls, pattern, samplingRate, voltage):
        """ Run an Import micrograph protocol. """
        cls.protImport = ProtImportMicrographs(pattern=pattern, samplingRate=samplingRate, voltage=voltage)
        cls.proj.launchProtocol(cls.protImport, wait=True)
        # check that input micrographs have been imported (a better way to do this?)
        if cls.protImport.outputMicrographs is None:
            raise Exception('Import of micrograph: %s, failed. outputMicrographs is None.' % pattern)
        return cls.protImport
    
    def testCtffind(self):
        #First, import a set of micrographs
        pattern = getInputPath('Micrographs_BPV1', '*.mrc')
        protImport = self.runImportMicrograph(pattern, samplingRate=1.273, voltage=300)
        
        protCTF = ProtCTFFind()
        protCTF.inputMicrographs.set(protImport.outputMicrographs)
        self.proj.launchProtocol(protCTF, wait=True)
        
        self.assertTrue(protCTF.outputMicrographs.hasCTF(), "CTF estimation has not been performed.")
        
if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestBrandeisCtffind)
    #suite = unittest.TestLoader().loadTestsFromName('test_protocols_brandeis.TestBrandeisCtffind.testCtffind')
    unittest.TextTestRunner(verbosity=2).run(suite)