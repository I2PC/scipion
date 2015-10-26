'''
Created on May 20, 2013

@author: laura
'''

import unittest
from pyworkflow.em.packages.eman2 import *
import pyworkflow.em.packages.eman2.convert as e2convert
from pyworkflow.tests import *


class TestEmanSetOfMicrographs(BaseTest):
        
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        cls.dataset = DataSet.getDataSet('xmipp_tutorial')  
        cls.dbGold = cls.dataset.getFile( 'micsGoldSqlite')
        cls.mdGold = cls.dataset.getFile('micsGoldXmd')
        cls.relionDataset = DataSet.getDataSet('relion_tutorial')  
                

    #Fixme: cache output is written in test input files        
    def testConvert(self):
        micSet = SetOfMicrographs(filename=self.dbGold)
        
        
        # Test writing micrgraphs to an hdf        
        filename = self.getOutputPath('micrographs.hdf')
        e2convert.writeSetOfParticles(micSet, filename, self.dataset.getPath())
        
        # Test writing a set of particles
#        partSet = SetOfParticles()
#        e2convert.readSetOfParticles(fnImages, imgSet)
        

class TestEmanSetOfCoordinates(BaseTest):
    
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        cls.dataset = DataSet.getDataSet('xmipp_tutorial')  
        cls.micBd = cls.dataset.getFile( 'micsSqlite')
        cls.emanCoord = cls.dataset.getFile('boxingFile')

        
    def testIterate(self):
        """ Test reading an EmanSetOfCoordinates from an existing  directory """
#         eman2.getEnviron()
        emanSetCoords = EmanSetOfCoordinates(self.emanCoord)     
#        #Set micrographs associated to coordinates
        emanSetMics = SetOfMicrographs(filename=self.micBd)
        emanSetCoords.setMicrographs(emanSetMics)   
        for coord in emanSetCoords.iterCoordinates():
            (x, y) = coord.getPosition()
            print ("Coordinate: x=%s y=%s" %(x,y))     
               
                    
if __name__ == '__main__':
#    suite = unittest.TestLoader().loadTestsFromName('test_data_enab.TestEmanSetOfMicrographs.testConvert')
#    unittest.TextTestRunner(verbosity=2).run(suite)
    
    unittest.main()
