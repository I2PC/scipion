'''
Created on May 20, 2013

@author: laura
'''

import unittest
from pyworkflow.em.packages.eman2 import *
import pyworkflow.em.packages.eman2.convert as e2convert
from pyworkflow.tests import *


class TestEmanSetOfMicrographs(unittest.TestCase):
        
    @classmethod
    def setUpClass(cls):
        cls.outputPath = getOutputPath('test_data_eman')
        
        cls.dbGold = getGoldPath('Micrographs_BPV3', 'micrographs_gold.sqlite')
        cls.mdGold = getGoldPath('Micrographs_BPV3', 'micrographs_gold.xmd')
        
                
    def setUp(self):
#        cleanPath(self.outputPath)
        makePath(self.outputPath)
            
    def testConvert(self):
        micSet = SetOfMicrographs()
        dbMic = getInputPath('Micrographs_BPV3', 'micrographs.sqlite')
        micSet.setFileName(dbMic)
        
        cwd = os.getcwd()
        # Change to test path
        os.chdir(getPath())
        
        # Test writing micrgraphs to an hdf        
        filename = getOutputPath('test_data_eman', 'micrographs.hdf')
        e2convert.writeSetOfParticles(micSet, filename)
        
        # Test writing a set of particles
        #partSet = SetOfParticles()
        #readSetOfParticles(fnImages, imgSet)
        
        os.chdir(cwd)

        
if __name__ == '__main__':
#    suite = unittest.TestLoader().loadTestsFromName('test_data_enab.TestEmanSetOfMicrographs.testConvert')
#    unittest.TextTestRunner(verbosity=2).run(suite)
    
    unittest.main()