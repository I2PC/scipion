"""
Created on Apr 9, 2013

@author: antonio
"""
import sys
import unittest
from pyworkflow.manager import Manager

projName = "TestProject"

class TestTransfer(unittest.TestCase):
    
    def setUp(self):
        manager = Manager()
        self.proj = manager.createProject(projName)  # Now it will be loaded if exists
    
    def testTransferProtocols(self):
        result = self.proj.mapper.selectByClass('EmanProtBoxing')
        if len(result):    
            for emanProtBoxing in result:
                emanProtBoxing.setHostName('glassfishdev')
                self.proj.sendProtocol(emanProtBoxing)
        self.assertTrue(True)  

if __name__ == "__main__":
    #import sys;sys.argv = ["", "Test.testName"]
    print "HOLAAAA"
    unittest.main()