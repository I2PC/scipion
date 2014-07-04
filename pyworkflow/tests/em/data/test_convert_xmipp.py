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

class TestConversions(BaseTest):
    
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        cls.dataset = DataSet.getDataSet('xmipp_tutorial')  
        cls.dbGold = cls.dataset.getFile( 'micsGoldSqlite')
        cls.particles = cls.dataset.getFile( 'particles1')

        
    def test_micrographsToMd(self):
        """ Test the convertion of a SetOfMicrographs to Xmipp metadata. """
        micSet = SetOfMicrographs(filename=self.getOutputPath("micrographs.sqlite"))
        n = 3
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
            file = self.dataset.getFile("mic%s"%(i + 1))
            p.setLocation(file)
            ctf = ctfs[i%2]
            p.setCTF(ctf)
            micSet.append(p)
            id = mdXmipp.addObject()
            mdXmipp.setValue(xmipp.MDL_ITEM_ID, long(i+1), id)
            mdXmipp.setValue(xmipp.MDL_MICROGRAPH, file, id)
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
        writeSetOfMicrographs(micSet, self.getOutputPath("micrographs.xmd"))
        self.assertEqual(mdScipion, mdXmipp, "metadata are not the same")
        
    def test_particlesToMd(self):
        """ Test the convertion of a SetOfParticles to Xmipp metadata. """
        imgSet = SetOfParticles(filename=self.getOutputPath("particles.sqlite"))
        n = 10
        fn = self.particles
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
        
    def test_CTF(self):
        """ Test the convertion of a SetOfParticles to Xmipp metadata. """
        mdCtf = xmipp.MetaData(self.dataset.getFile('ctfGold'))
        objId = mdCtf.firstObject()
        ctf = rowToCtfModel(mdCtf, objId)        
        
        ALL_CTF_LABELS = [   
            xmipp.MDL_CTF_CA,
            xmipp.MDL_CTF_ENERGY_LOSS,
            xmipp.MDL_CTF_LENS_STABILITY,
            xmipp.MDL_CTF_CONVERGENCE_CONE,
            xmipp.MDL_CTF_LONGITUDINAL_DISPLACEMENT,
            xmipp.MDL_CTF_TRANSVERSAL_DISPLACEMENT,
            xmipp.MDL_CTF_K,
            xmipp.MDL_CTF_BG_GAUSSIAN_K,
            xmipp.MDL_CTF_BG_GAUSSIAN_SIGMAU,
            xmipp.MDL_CTF_BG_GAUSSIAN_SIGMAV,
            xmipp.MDL_CTF_BG_GAUSSIAN_CU,
            xmipp.MDL_CTF_BG_GAUSSIAN_CV,
            xmipp.MDL_CTF_BG_SQRT_K,
            xmipp.MDL_CTF_BG_SQRT_U,
            xmipp.MDL_CTF_BG_SQRT_V,
            xmipp.MDL_CTF_BG_SQRT_ANGLE,
            xmipp.MDL_CTF_BG_BASELINE,
            xmipp.MDL_CTF_BG_GAUSSIAN2_K,
            xmipp.MDL_CTF_BG_GAUSSIAN2_SIGMAU,
            xmipp.MDL_CTF_BG_GAUSSIAN2_SIGMAV,
            xmipp.MDL_CTF_BG_GAUSSIAN2_CU,
            xmipp.MDL_CTF_BG_GAUSSIAN2_CV,
            xmipp.MDL_CTF_BG_GAUSSIAN2_ANGLE,
            xmipp.MDL_CTF_CRIT_FITTINGSCORE,
            xmipp.MDL_CTF_CRIT_FITTINGCORR13,
            xmipp.MDL_CTF_DOWNSAMPLE_PERFORMED,
            xmipp.MDL_CTF_CRIT_PSDVARIANCE,
            xmipp.MDL_CTF_CRIT_PSDPCA1VARIANCE,
            xmipp.MDL_CTF_CRIT_PSDPCARUNSTEST,
            xmipp.MDL_CTF_CRIT_FIRSTZEROAVG,
            xmipp.MDL_CTF_CRIT_DAMPING,
            xmipp.MDL_CTF_CRIT_FIRSTZERORATIO,
            xmipp.MDL_CTF_CRIT_PSDCORRELATION90,
            xmipp.MDL_CTF_CRIT_PSDRADIALINTEGRAL,
            xmipp.MDL_CTF_CRIT_NORMALITY,  
        ]      

        for label in ALL_CTF_LABELS:
            attrName = '_xmipp_%s' % xmipp.label2Str(label)
            self.assertAlmostEquals(mdCtf.getValue(label, objId), ctf.getAttributeValue(attrName))
            
        
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
        fnScipion=self.getOutputPath("writeSetOfDefocusGroups.sqlite")
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
        
        fnXmipp=self.getOutputPath("writeSetOfDefocusGroups.xmd")
        writeSetOfDefocusGroups(setOfDefocus, fnXmipp)
        mdAux = xmipp.MetaData(fnXmipp)
        self.assertEqual(md,mdAux, "test writeSetOfDefocusGroups fails")
       
if __name__ == '__main__':
#    suite = unittest.TestLoader().loadTestsFromName('test_data_xmipp.TestXmippSetOfMicrographs.testCopy')
#    unittest.TextTestRunner(verbosity=2).run(suite)
    
    unittest.main()