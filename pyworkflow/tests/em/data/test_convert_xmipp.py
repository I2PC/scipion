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
from pyworkflow.em.packages.xmipp3.convert import *

class TestConversions(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.outputPath = getOutputPath('test_convert_xmipp')   

    def setUp(self):
        cleanPath(self.outputPath)
        makePath(self.outputPath)
        
    def getPath(self, fn):
        return os.path.join(self.outputPath, fn)
        
    def test_micrographsToMd(self):
        """ Test the convertion of a SetOfMicrographs to Xmipp metadata. """
        micSet = SetOfMicrographs(filename=self.getPath("micrographs.sqlite"))
        n = 10
        fn = "images.stk"
        ctfs = [CTFModel(defocusU=10000, defocusV=15000, defocusAngle=15),
                CTFModel(defocusU=20000, defocusV=25000, defocusAngle=25)
               ]
        acquisition = Acquisition(magnification=60000, 
                                  voltage=300,
                                  sphericalAberration=2., 
                                  amplitudeContrast=0.07)
        micSet.setAcquisition(acquisition)
        micSet.setSamplingRate(1.)
        mdXmipp = xmipp.MetaData()
        

        for i in range(n):
            p = Micrograph()
            p.setLocation(i+1, fn)
            ctf = ctfs[i%2]
            p.setCTF(ctf)
            micSet.append(p)
            id = mdXmipp.addObject()
            mdXmipp.setValue(xmipp.MDL_ITEM_ID, long(i+1), id)
            mdXmipp.setValue(xmipp.MDL_MICROGRAPH, locationToXmipp(i+1, fn), id)
            # set CTFModel params
            mdXmipp.setValue(xmipp.MDL_CTF_DEFOCUSU, ctf.getDefocusU(), id)
            mdXmipp.setValue(xmipp.MDL_CTF_DEFOCUSV, ctf.getDefocusV(), id)
            mdXmipp.setValue(xmipp.MDL_CTF_DEFOCUS_ANGLE, ctf.getDefocusAngle(), id)
            # set Acquisition params
            mdXmipp.setValue(xmipp.MDL_CTF_Q0, acquisition.getAmplitudeContrast(), id)
            mdXmipp.setValue(xmipp.MDL_CTF_CS, acquisition.getSphericalAberration(), id)
            mdXmipp.setValue(xmipp.MDL_CTF_VOLTAGE, acquisition.getVoltage(), id)
            
        mdScipion = xmipp.MetaData()
        setOfMicrographsToMd(micSet, mdScipion)
        writeSetOfMicrographs(micSet, self.getPath("micrographs.xmd"))
        self.assertEqual(mdScipion, mdXmipp, "metadata are not the same")
        
    def test_particlesToMd(self):
        """ Test the convertion of a SetOfParticles to Xmipp metadata. """
        imgSet = SetOfParticles(filename=self.getPath("particles.sqlite"))
        n = 10
        fn = "images.stk"
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
            id = mdXmipp.addObject()
            mdXmipp.setValue(xmipp.MDL_ITEM_ID, long(i+1), id)
            mdXmipp.setValue(xmipp.MDL_IMAGE, locationToXmipp(i+1, fn), id)
            # set CTFModel params
            mdXmipp.setValue(xmipp.MDL_CTF_DEFOCUSU, ctf.getDefocusU(), id)
            mdXmipp.setValue(xmipp.MDL_CTF_DEFOCUSV, ctf.getDefocusV(), id)
            mdXmipp.setValue(xmipp.MDL_CTF_DEFOCUS_ANGLE, ctf.getDefocusAngle(), id)
            # set Acquisition params
            mdXmipp.setValue(xmipp.MDL_CTF_Q0, acquisition.getAmplitudeContrast(), id)
            mdXmipp.setValue(xmipp.MDL_CTF_CS, acquisition.getSphericalAberration(), id)
            mdXmipp.setValue(xmipp.MDL_CTF_VOLTAGE, acquisition.getVoltage(), id)
            
        mdScipion = xmipp.MetaData()
        setOfParticlesToMd(imgSet, mdScipion)
        self.assertEqual(mdScipion, mdXmipp, "metadata are not the same")

    def test_writeSetOfDefocusGroups(self):
        #reference metadata
        md = xmipp.MetaData()
        objId = md.addObject()
        defocusGroupRow = XmippMdRow()

        defocusGroupRow.setValue(xmipp.MDL_CTF_GROUP, 1)
        defocusGroupRow.setValue(xmipp.MDL_MIN, 2000.)
        defocusGroupRow.setValue(xmipp.MDL_MAX, 2500.)
        defocusGroupRow.setValue(xmipp.MDL_AVG, 2100.)
        defocusGroupRow.writeToMd(md, objId)

        objId = md.addObject()
        defocusGroupRow.setValue(xmipp.MDL_CTF_GROUP, 2)
        defocusGroupRow.setValue(xmipp.MDL_MIN, 3000.)
        defocusGroupRow.setValue(xmipp.MDL_MAX, 5500.)
        defocusGroupRow.setValue(xmipp.MDL_AVG, 5000.)
        defocusGroupRow.writeToMd(md, objId)
        #
        fnScipion=self.getPath("writeSetOfDefocusGroups.sqlite")
        setOfDefocus = SetOfDefocusGroup(filename=fnScipion)

        df = DefocusGroup()
        df.setDefocusMin(2000.)
        df.setDefocusMax(2500.)
        df.setDefocusAvg(2100.)
        setOfDefocus.append(df)
        
        df.cleanObjId()
        df.setDefocusMin(3000)
        df.setDefocusMax(5500)
        df.setDefocusAvg(5000)
        setOfDefocus.append(df)
        
        fnXmipp=self.getPath("writeSetOfDefocusGroups.xmd")
        writeSetOfDefocusGroups(setOfDefocus, fnXmipp)
        mdAux = xmipp.MetaData(fnXmipp)
        self.assertEqual(md,mdAux, "test writeSetOfDefocusGroups fails")
       
if __name__ == '__main__':
#    suite = unittest.TestLoader().loadTestsFromName('test_data_xmipp.TestXmippSetOfMicrographs.testCopy')
#    unittest.TextTestRunner(verbosity=2).run(suite)
    
    unittest.main()