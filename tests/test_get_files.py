'''
Created on May 8, 2013

@author: antonio
'''
import sys

from pyworkflow.manager import Manager
from pyworkflow.project import *
from pyworkflow.em import *
from em.packages.xmipp3.data import XmippSetOfMicrographs
import unittest
import argparse

projName = None

class TestGetFiles(unittest.TestCase):
    
    def setUp(self):
        manager = Manager()
        self.proj = manager.createProject(projName)
        
    def testXmippSetOfMicrographs(self): 
        print ('*****************************************   XmippSetOfMicrographs   ***************************************************')
        result = self.proj.mapper.selectByClass('XmippSetOfMicrographs')
        if len(result):    
            for xmippSetOfMicrographs in result:
                print ("XmippSetOfMicrographs files: " + str(xmippSetOfMicrographs.getFiles()))
        else:
            print "Not XmippSetOfMicrographs found"
        print ('*************************************************************************************************************************')    
    
    def testXmippSetOfCoordinates(self):
        print ('****************************************  XmippSetOfCoordinates   ****************************************************')
        result = self.proj.mapper.selectByClass('XmippSetOfCoordinates')
        if len(result):    
            for xmippSetOfCoordinates in result:
                print ("XmippSetOfCoordinates files: " + str(xmippSetOfCoordinates.getFiles()))
        else:
            print "Not XmippSetOfCoordinates found"
        print ('*************************************************************************************************************************')
        self.assertTrue(True)    
      
    def testEmanSetOfCoordinates(self): 
        print ('****************************************  EmanSetOfCoordinates  ******************************************************')
        result = self.proj.mapper.selectByClass('EmanSetOfCoordinates')
        if len(result):    
            for emanSetOfCoordinates in result:
                print ("EmanSetOfCoordinates files: " + str(emanSetOfCoordinates.getFiles()))
        else:
            print "Not EmanSetOfCoordinates found"  
        print ('*************************************************************************************************************************')
        self.assertTrue(True)
    
    def testProtImportMicrographs(self):
        print ('*****************************************  ProtImportMicrographs  *********************************************************')   
        result = self.proj.mapper.selectByClass('ProtImportMicrographs')
        if len(result):    
            for protImportMicrographs in result:
                print ("ProtImportMicrographs files: " + str(protImportMicrographs.getFiles()))
        else:
            print "Not ProtImportMicrographs found"
        print ('*************************************************************************************************************************')
        self.assertTrue(True)
        
    def testXmippProtPreprocessMicrographs(self):
        print ('****************************************  XmippProtPreprocessMicrographs  ******************************************************')
        result = self.proj.mapper.selectByClass('XmippProtPreprocessMicrographs')
        if len(result):    
            for xmippProtPreprocessMicrographs in result:
                print ("XmippProtPreprocessMicrographs files: " + str(xmippProtPreprocessMicrographs.getFiles()))
        else:
            print "Not XmippProtPreprocessMicrographs found"
        print ('*************************************************************************************************************************')  
        self.assertTrue(True) 
    
    def testXmippProtParticlePicking(self):
        print ('****************************************  XmippProtParticlePicking  ******************************************************')
        result = self.proj.mapper.selectByClass('XmippProtParticlePicking')
        if len(result):    
            for xmippProtParticlePicking in result:
                print ("XmippProtParticlePicking files: " + str(xmippProtParticlePicking.getFiles()))
        else:
            print "Not XmippProtParticlePicking found"
        print ('*************************************************************************************************************************')  
        self.assertTrue(True) 
    
    def testXmippProtExtractParticles(self):
        print ('****************************************  XmippProtExtractParticles  ******************************************************')
        result = self.proj.mapper.selectByClass('XmippProtExtractParticles')
        if len(result):    
            for xmippProtExtractParticles in result:
                print ("XmippProtExtractParticles files: " + str(xmippProtExtractParticles.getFiles()))
        else:
            print "Not XmippProtExtractParticles found"
        print ('*************************************************************************************************************************')  
        self.assertTrue(True)        
        
    def testXmippCTFModel(self):
        print ('****************************************  XmippCTFModel  ******************************************************')
        result = self.proj.mapper.selectByClass('XmippCTFModel')
        if len(result):    
            for xmippCTFModel in result:
                print ("XmippCTFModel files: " + str(xmippCTFModel.getFiles()))
        else:
            print "Not XmippCTFModel found"
        print ('*************************************************************************************************************************')  
        self.assertTrue(True) 
    
    def testEmanProtBoxing(self): 
        print ('***************************************  EmanProtBoxing  ********************************************************')       
        result = self.proj.mapper.selectByClass('EmanProtBoxing')
        if len(result):    
            for emanProtBoxing in result:
                print ("EmanProtBoxing files: " + str(emanProtBoxing.getFiles()))
        else:
            print "Not EmanProtBoxing found"
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
