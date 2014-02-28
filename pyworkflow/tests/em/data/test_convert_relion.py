'''
Created on 5th Feb, 2014

@author: Roberto Marabini
         J.M. de la Rosa Trevin
'''

import os
import unittest

import xmipp
from pyworkflow.em.packages.xmipp3 import *
from pyworkflow.tests import *
import pyworkflow.em.packages.relion as relion


class TestConversions(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.outputPath = getOutputPath('test_convert_relion')   

    def setUp(self):
        cleanPath(self.outputPath)
        makePath(self.outputPath)
        
    def getPath(self, fn):
        return os.path.join(self.outputPath, fn)
        
    def test_particlesToRelion(self):
        """ Test the convertion of a SetOfParticles to Xmipp metadata. """
        imgSet = SetOfParticles(filename=self.getPath("particles.sqlite"))
        n = 10
        fn = "particles.mrc"
        ctfs = [CTFModel(defocusU=10000, defocusV=15000, defocusAngle=15),
                CTFModel(defocusU=20000, defocusV=25000, defocusAngle=25)
               ]
        acquisition = Acquisition(magnification=60000, voltage=300,
                                  sphericalAberration=2., amplitudeContrast=0.07)
        mdXmipp = xmipp.MetaData()
        imgSet.setAcquisition(acquisition)

        for i in range(n):
            p = Particle()
            p.setLocation(i+1, fn)
            ctf = ctfs[i%2]
            p.setCTF(ctf)
            p.setAcquisition(acquisition)
            imgSet.append(p)
            
        fnStar = self.getPath('particles.star')
        fnStk = self.getPath('particles.stk')
        
        relion.writeSetOfParticles(imgSet, fnStar, fnStk)
        

       
if __name__ == '__main__':
#    suite = unittest.TestLoader().loadTestsFromName('test_data_xmipp.TestXmippSetOfMicrographs.testCopy')
#    unittest.TextTestRunner(verbosity=2).run(suite)
    
    unittest.main()