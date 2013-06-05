"""
Created on Apr 9, 2013

@author: antonio
"""
import unittest
import argparse
import sys
from pyworkflow.manager import Manager
from pyworkflow.em.packages.xmipp3 import *
from pyworkflow.utils.utils import executeRemote, executeLongRemote


projName = None

class TestRemoteProtocol(unittest.TestCase):
    
    def setUp(self):
        manager = Manager()
        self.proj = manager.createProject(projName)
      
    def testRemoteCTF(self):
        # Recover import micrograph protocol to get output micrographs
        result = self.proj.mapper.selectByClass('ProtImportMicrographs')
        if len(result):    
            # We only need one for this test
            protImportMicrographs = result[0]
            protCTF = XmippProtCTFMicrographs(numberOfThreads=1)
            protCTF.inputMicrographs.set(protImportMicrographs.outputMicrographs)
            # We create the protocol
            protCTF.setHostName("crunchy")
            self.proj._storeProtocol(protCTF)
            # We send the protocol
            self.proj.launchRemoteProtocol(protCTF)
        else:
            print "Not ProtImportMicrographs found"        
        self.assertTrue(True)  
        
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('projName')
    parser.add_argument('unittest_args', nargs='*')
    args = parser.parse_args()
    projName = args.projName
    sys.argv[1:] = args.unittest_args
    unittest.main()