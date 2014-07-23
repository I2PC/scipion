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


class TestConversions(BaseTest):
    
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        cls.dataset = DataSet.getDataSet('xmipp_tutorial')  
        cls.particles = cls.dataset.getFile( 'particles1')
        
    def test_particlesToRelion(self):
        """ Test the convertion of a SetOfParticles to Xmipp metadata. """
        imgSet = SetOfParticles(filename=self.getOutputPath("particles.sqlite"))
        n = 10
        fn = self.particles
        ctfs = [CTFModel(defocusU=10000, defocusV=15000, defocusAngle=15),
                CTFModel(defocusU=20000, defocusV=25000, defocusAngle=25)
               ]
        acquisition = Acquisition(magnification=60000, voltage=300,
                                  sphericalAberration=2., amplitudeContrast=0.07)
        imgSet.setAcquisition(acquisition)
        coord = Coordinate()
        coord.setMicId(1)

        for i in range(n):
            p = Particle()
            p.setLocation(i+1, fn)
            ctf = ctfs[i%2]
            p.setCTF(ctf)
            p.setAcquisition(acquisition)
            p._xmipp_zScore = Float(i)
            coord.setX(i*10)
            coord.setY(i*10)
            p.setCoordinate(coord)
            imgSet.append(p)
            
        fnStar = self.getOutputPath('particles.star')
        fnStk = self.getOutputPath('particles.stk')
        
        print ">>> Writing to file: %s" % fnStar
        relion.writeSetOfParticles(imgSet, fnStar, fnStk)
        
        relion.addRelionLabels()
        md = xmipp.MetaData(fnStar)
        self.assertTrue(md.containsLabel(xmipp.MDL_XCOOR))
        self.assertTrue(md.containsLabel(xmipp.MDL_YCOOR))
        self.assertFalse(md.containsLabel(xmipp.MDL_ZSCORE))
