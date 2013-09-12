'''
Created on May 20, 2013

@author: laura
'''

import unittest

from pyworkflow.em.packages.xmipp3 import *
from pyworkflow.tests import *


class TestXmippSetOfMicrographs(unittest.TestCase):
        
    @classmethod
    def setUpClass(cls):
        cls.outputPath = getOutputPath('test_data_xmipp')
        
        cls.dbGold = getGoldPath('Micrographs_BPV3', 'micrographs_gold.sqlite')
        cls.mdGold = getGoldPath('Micrographs_BPV3', 'micrographs_gold.xmd')
        
                
    def setUp(self):
#        cleanPath(self.outputPath)
        makePath(self.outputPath)
            
    def testConvert(self):
        """ Test the convertion of a SetOfMicrographs to Xmipp"""
        micSet = SetOfMicrographs()
        dbMic = getInputPath('Micrographs_BPV3', 'micrographs.sqlite')
        micSet.setFileName(dbMic)
        mdFn = getOutputPath(self.outputPath, 'micrographs.xmd')
        
        writeSetOfMicrographs(micSet, mdFn)
                
        self.assertEqual(xmipp.MetaData(mdFn), xmipp.MetaData(self.mdGold), "xmipp metadata wrong")
                  
    
class TestXmippSetOfCoordinates(unittest.TestCase):  
    
    @classmethod
    def setUpClass(cls):
        cls.outputPath = getOutputPath('test_data_xmipp')   
        
        cls.dbGold = getGoldPath('Micrographs_BPV3', 'micrographs_gold.sqlite')
        cls.mdGold = getGoldPath('Picking_XmippBPV3_Down3', 'coordinates_gold.sqlite')

        
    def setUp(self):
        cleanPath(self.outputPath)
        makePath(self.outputPath)
        
    def testConvert(self):
        """ Test converting a SetOfCoordinates to Xmipp """
        micSet = SetOfMicrographs()
        micSet.setFileName(self.dbGold)
        
        # Test reading a set of coordinates
        posDir = getInputPath('Picking_XmippBPV3_Down3')
        print "reading pos from :", posDir
        coordSet = SetOfCoordinates()
        fn = getOutputPath(self.outputPath, 'coordinates.sqlite')
        coordSet.setFileName(fn)
        
        readSetOfCoordinates(posDir, micSet, coordSet)
        coordSet.write()
        
        # Test writing a set of coordinates to xmipp format
        coordSet = SetOfCoordinates()
        dbCoord = getInputPath('Picking_XmippBPV3_Down3', 'coordinates.sqlite')
        coordSet.setFileName(dbCoord)
        coordSet.setMicrographs(micSet)
        
        posDir = getOutputPath(self.outputPath)
        writePosCoordinates(posDir, coordSet)
        
    
class TestXmippCTFModel(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.outputPath = getOutputPath('test_data_xmipp')   
        
    def setUp(self):
        cleanPath(self.outputPath)
        makePath(self.outputPath)
                    
    def testConvert(self):
        """ Test converting a CTFModel to Xmipp """
        pass

        
if __name__ == '__main__':
#    suite = unittest.TestLoader().loadTestsFromName('test_data_xmipp.TestXmippSetOfMicrographs.testCopy')
#    unittest.TextTestRunner(verbosity=2).run(suite)
    
    unittest.main()