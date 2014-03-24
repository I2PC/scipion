'''
Created on May 20, 2013

@author: laura
'''

import unittest

from pyworkflow.em.packages.xmipp3 import *
from pyworkflow.tests import *


class TestXmippSetOfMicrographs(BaseTest):
        
    @classmethod
    def setUpClass(cls):
        super(TestXmippSetOfMicrographs, cls).setUpClass()
        cls.dataset = DataSet.getDataSet('xmipp_tutorial')  
        cls.dbGold = cls.dataset.getFile( 'micsGoldSqlite')
        cls.mdGold = cls.dataset.getFile('micsGoldXmd')
        
            
    def testConvert(self):
        """ Test the convertion of a SetOfMicrographs to Xmipp metadata"""
        micSet = SetOfMicrographs(filename=self.dbGold)
        mdFn = self.getOutputPath('micrographs.xmd')
        
        writeSetOfMicrographs(micSet, mdFn)
                
        self.assertEqual(xmipp.MetaData(mdFn), xmipp.MetaData(self.mdGold), "xmipp metadata wrong")
                  
    
class TestXmippSetOfCoordinates(BaseTest):  
    
    @classmethod
    def setUpClass(cls):
        super(TestXmippSetOfCoordinates, cls).setUpClass()
        cls.dataset = DataSet.getDataSet('xmipp_tutorial')  
        cls.dbGold = cls.dataset.getFile( 'micsGoldSqlite')
        cls.posDir = cls.dataset.getFile('posDir')
        cls.dbCoord = cls.dataset.getFile('coordsGoldSqlite')
        
        
    def testConvert(self):
        """ Test converting a SetOfCoordinates to Xmipp """
        micSet = SetOfMicrographs(filename=self.dbGold)
        
        # Test reading a set of coordinates
        fn = self.getOutputPath('coordinates.sqlite')
        coordSet = SetOfCoordinates(filename=fn)
        readSetOfCoordinates(self.posDir, micSet, coordSet)
        coordSet.write()
        
        # Test writing a set of coordinates to xmipp format
        coordSet = SetOfCoordinates(filename=self.dbCoord)
        coordSet.setBoxSize(512)
        coordSet.setMicrographs(micSet)
        writeSetOfCoordinates(self.outputPath, coordSet)
        
    
class TestXmippCTFModel(BaseTest):
    
    @classmethod
    def setUpClass(cls):
        super(TestXmippCTFModel, cls).setUpClass()
        cls.dataset = DataSet.getDataSet('xmipp_tutorial')  

                    
    def testConvert(self):
        """ Test converting a CTFModel to Xmipp """
        pass

        
if __name__ == '__main__':
#    suite = unittest.TestLoader().loadTestsFromName('test_data_xmipp.TestXmippSetOfMicrographs.testCopy')
#    unittest.TextTestRunner(verbosity=2).run(suite)
    
    unittest.main()