'''
Created on May 8, 2013

@author: antonio
'''
import sys

from pyworkflow.manager import Manager
from pyworkflow.project import *
from pyworkflow.em import *
import unittest
import argparse

projName = None

class TestIterEman(unittest.TestCase):
    
    def setUp(self):
        manager = Manager()
        self.proj = manager.createProject(projName)
        
    def testIterEman(self):   
        print ('*************************************************************************************************************************')
        result = self.proj.mapper.selectByClass('EmanSetOfCoordinates')
        if len(result):    
            for emanSetOfCoordinates in result:
                for emanCoordinate in emanSetOfCoordinates.iterCoordinates():
                    print ("Coordinate: " + str(emanCoordinate.getPosition()))
        else:
            print "Not EmanSetOfCoordinates found"
        print ('*************************************************************************************************************************')
        self.assertTrue(True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('projName')
    parser.add_argument('unittest_args', nargs='*')
    args = parser.parse_args()
    projName = args.projName
    sys.argv[1:] = args.unittest_args
    unittest.main()