'''
Created on May 20, 2013

@author: laura
'''

import unittest
from glob import glob
from pyworkflow.em.packages.eman2 import *
from pyworkflow.tests import *
   
class TestEmanSetOfCoordinates(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.outputPath = getOutputPath('test_data_eman')
        cls.emanCoord = getInputPath('EmanTestProject2', 'scipion_micrographs_coordinates.json')
        cls.micBd = getInputPath('Micrographs_BPV3', 'micrographs.sqlite')
        
    def setUp(self):
        cleanPath(self.outputPath)
        makePath(self.outputPath)
        
    def testIterate(self):
        """ Test reading an EmanSetOfCoordinates from an existing  directory """
        #eman2.loadEnvironment()
        emanSetCoords = EmanSetOfCoordinates(self.emanCoord)     
        #Set micrographs associated to coordinates
        emanSetMics = SetOfMicrographs(self.micBd)
        emanSetCoords.setMicrographs(emanSetMics)   
        for coord in emanSetCoords.iterCoordinates():
            (x, y) = coord.getPosition()
            print ("Coordinate: x=%s y=%s" %(x,y))       
            
            
if __name__ == '__main__':
#    suite = unittest.TestLoader().loadTestsFromName('test_eman_data.TestEmanSetOfImages.testReadParams')
#    unittest.TextTestRunner(verbosity=2).run(suite)
    
    unittest.main()