'''
Created on May 20, 2013

@author: laura
'''

import unittest
from glob import glob
from pyworkflow.em.packages.eman2 import *
from pyworkflow.tests import *

class TestEmanUtils(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):        
        cls.outputPath = getOutputPath('test_data_eman')
        
        cls.emanDir = getInputPath('EmanTestProject')
        
    def setUp(self):
        cleanPath(self.outputPath)
        makePath(self.outputPath)
                
    def testReadParams(self):
        """ Test reading parameters from eman bdb """
        # Enter eman directory
        os.chdir(self.emanDir)
        value = self.__getEmanParamValue('box_size')
        self.assertTrue(int(value) == 110, "Incorrect particle size")
            
        self.assertTrue(self.__getEmanParamValue('write_particles'), "Not particles written")

        
    def __getEmanParamValue(self, paramName):
        """ Recover a parameter value from EMAN Berkeley data base. """        
        command = "e2bdb.py -D bdb:emboxerbase"
        pipe = os.popen(command)
        stOutput = pipe.readlines()
        pipe.close()
        auxValue = None
        for line in stOutput:
            if (paramName in line):
                auxValue = line.split(" : ")[1]
        if auxValue is None:
            return None
    
        return auxValue
    
class TestEmanSetOfCoordinates(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.outputPath = getOutputPath('test_data_eman')
        cls.emanDir = getInputPath('EmanTestProject')
        cls.micBd = getInputPath('Micrographs_BPV3', 'micrographs.sqlite')
        
    def setUp(self):
        cleanPath(self.outputPath)
        makePath(self.outputPath)
        
    def testIterate(self):
        """ Test reading an EmanSetOfCoordinates from an existing  directory """
        emanSetCoords = EmanSetOfCoordinates(self.emanDir)     
        #Set micrographs associated to coordinates
        emanSetMics = SetOfMicrographs(self.micBd)
        emanSetCoords.setMicrographs(emanSetMics)   
        for coord in emanSetCoords.iterCoordinates():
            (x, y) = coord.getPosition()
            #print ("Coordinate: x=%s y=%s" %(x,y))       
            
            
if __name__ == '__main__':
#    suite = unittest.TestLoader().loadTestsFromName('test_eman_data.TestEmanSetOfImages.testReadParams')
#    unittest.TextTestRunner(verbosity=2).run(suite)
    
    unittest.main()