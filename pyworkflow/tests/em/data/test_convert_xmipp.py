'''
Created on 5th Feb, 2014

@author: Roberto Marabini
         J.M. de la Rosa Trevin
'''

import unittest

import xmipp
from pyworkflow.em.packages.xmipp3 import *
from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3.convert import *

class TestConversions(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.outputPath = getOutputPath('test_convert_xmipp')   

    def setUp(self):
        cleanPath(self.outputPath)
        makePath(self.outputPath)
        
    def getOutputFn(self, fn):
        return os.path.join(self.outputPath, fn)
        
    def test_particlesToMd(self):
        """ Test the convertion of a SetOfParticles to Xmipp metadata. """
        imgFn = self.getOutputFn("particles.sqlite")
        imgSet = SetOfParticles(filename=imgFn)
        n = 10
        fn = "images.stk"
        ctfs = [CTFModel(defocusU=10000, defocusV=15000, defocusAngle=15),
                CTFModel(defocusU=20000, defocusV=25000, defocusAngle=25)
               ]
        acq = Acquisition(magnification=60000, voltage=300,
                          sphericalAberration=2.,amplitudeContrast=0.07)
        mdXmipp = xmipp.MetaData()

        for i in range(n):
            p = Particle()
            p.setLocation(i+1, fn)
            ctf = ctfs[i%2]
            p.setCTF(ctf)
            p.setAcquisition(acq)
            imgSet.append(p)
            id = mdXmipp.addObject()
            mdXmipp.setValue(xmipp.MDL_ITEM_ID, long(i+1), id)
            mdXmipp.setValue(xmipp.MDL_IMAGE,locationToXmipp(i+1, fn), id)
            mdXmipp.setValue(xmipp.MDL_CTF_DEFOCUSU,ctf.getDefocusU(), id)
            mdXmipp.setValue(xmipp.MDL_CTF_DEFOCUSV,ctf.getDefocusV(), id)
            mdXmipp.setValue(xmipp.MDL_CTF_DEFOCUS_ANGLE,ctf.getDefocusAngle(), id)
            mdXmipp.setValue(xmipp.MDL_CTF_DEFOCUSU,acq.getA, id)

        mdScipion = xmipp.MetaData()
        setOfParticlesToMd(imgSet, mdScipion)
        self.assertEqual(mdScipion, mdXmipp, "metadata are not the same")
        
if __name__ == '__main__':
#    suite = unittest.TestLoader().loadTestsFromName('test_data_xmipp.TestXmippSetOfMicrographs.testCopy')
#    unittest.TextTestRunner(verbosity=2).run(suite)
    
    unittest.main()