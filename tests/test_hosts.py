"""
Created on Apr 9, 2013

@author: antonio
"""
import unittest
import sys
from pyworkflow.manager import Manager
import argparse

projName = None

class TestHosts(unittest.TestCase):
    
    def setUp(self):
        manager = Manager()
        self.proj = manager.createProject(projName)
      
    def testAllProjectHosts(self):
        for executionHostConfig in self.proj.getHosts():
            print (executionHostConfig.getLabel())     
        self.assertTrue(True)  
        
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('projName')
    parser.add_argument('unittest_args', nargs='*')
    args = parser.parse_args()
    projName = args.projName
    sys.argv[1:] = args.unittest_args
    unittest.main()