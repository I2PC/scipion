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
        cls.dataset = DataSet.getDataSet('relion_tutorial')  
        cls.getFile = cls.dataset.getFile
        
        cls.ds = DataSet.getDataSet('relion_tutorial')
        
        
    def aaatest_particlesToStar(self):
        """ Write a SetOfParticles to Relion star input file. """
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
        
    def test_particlesFromStar(self):
        """ Read a set of particles from an .star file.  """
        fnStar = self.getFile('relion_it020_data')
        
        relion.addRelionLabels()
        print ">>> Reading star file: ", fnStar
        md = xmipp.MetaData(fnStar)
        print "labels: ", [xmipp.label2Str(l) for l in md.getActiveLabels()]
        print "size: ", md.size()
        
        
    def test_particlesToStar(self):
        fnSqlite = self.getFile('particles')
        partSet = SetOfParticles(filename=fnSqlite)
        
        for i, img in enumerate(partSet):
            img.printAll()
            if i > 10: 
                break
        partSet.printAll()
        #print md
        
