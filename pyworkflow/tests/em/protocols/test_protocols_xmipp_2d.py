# **************************************************************************
# *
# * Authors:    Laura del Cano (ldelcano@cnb.csic.es)
# *             Josue Gomez Blanco (jgomez@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3 import *


# Some utility functions to import particles that are used
# in several tests.
class TestXmippBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='xmipp_tutorial'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.particlesFn = cls.dataset.getFile('particles')
    
    @classmethod
    def runImportParticles(cls, pattern, samplingRate, checkStack=False):
        """ Run an Import particles protocol. """
        cls.protImport = cls.newProtocol(ProtImportParticles, 
                                         pattern=pattern, samplingRate=samplingRate, 
                                         checkStack=checkStack)
        cls.launchProtocol(cls.protImport)
        # check that input images have been imported (a better way to do this?)
        if cls.protImport.outputParticles is None:
            raise Exception('Import of images: %s, failed. outputParticles is None.' % pattern)
        return cls.protImport
    
    @classmethod
    def runCL2DAlign(cls, particles):
        cls.CL2DAlign = cls.newProtocol(XmippProtCL2DAlign, 
                                        maximumShift=2, numberOfIterations=2,
                                        numberOfMpi=4, numberOfThreads=1, useReferenceImage=False)
        cls.CL2DAlign.inputParticles.set(particles)
        cls.launchProtocol(cls.CL2DAlign)
        return cls.CL2DAlign
    
    @classmethod
    def runClassify(cls, particles):
        cls.ProtClassify = cls.newProtocol(XmippProtML2D, 
                                           numberOfReferences=8, maxIters=4, doMlf=False,
                                           numberOfMpi=2, numberOfThreads=2)
        cls.ProtClassify.inputParticles.set(particles)
        cls.launchProtocol(cls.ProtClassify)
        return cls.ProtClassify


class TestXmippPreprocessParticles(TestXmippBase):
    """This class check if the protocol to preprocess particles in Xmipp works properly."""
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
    
    def test_preprocessPart(self):
        print "Run Preprocess particles"
        protPreproc = self.newProtocol(XmippProtPreprocessParticles, 
                                      doRemoveDust=True, doNormalize=True, 
                                      backRadius=48, doInvert=True,
                                      doThreshold=True, thresholdType=1)
        
        protPreproc.inputParticles.set(self.protImport.outputParticles)
        self.launchProtocol(protPreproc)
        
        self.assertIsNotNone(protPreproc.outputParticles, "There was a problem with preprocess particles")


class TestXmippCropResizeParticles(TestXmippBase):
    """This class check if the protocol to crop/resize particles in Xmipp works properly."""
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData('xmipp_tutorial')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 1.237, True)
        cls.acquisition = cls.protImport.outputParticles.getAcquisition()
        
    def validateAcquisition(self, particles):
        acquisition = particles.getAcquisition()
        self.assertAlmostEqual(acquisition.getVoltage(), self.acquisition.getVoltage())
    
    def test_cropResizePart(self):
        print "Run crop/resize particles"
        protCropResize = self.newProtocol(XmippProtCropResizeParticles, 
                                         doResize=True, resizeOption=1, resizeDim=128, 
                                         doWindow=True, windowOperation=1, windowSize=256)
        protCropResize.inputParticles.set(self.protImport.outputParticles)
        self.launchProtocol(protCropResize)
        
        self.assertIsNotNone(protCropResize.outputParticles, "There was a problem with resize/crop the particles")
        self.validateAcquisition(protCropResize.outputParticles)
    
    def test_cropResizePart2(self):
        print "Run crop/resize particles v2"
        protCropResize = self.newProtocol(XmippProtCropResizeParticles,  
                                         doResize=True, resizeOption=3, factor=0.5, 
                                         doWindow=True, windowOperation=1, windowSize=256)
        protCropResize.inputParticles.set(self.protImport.outputParticles)
        self.launchProtocol(protCropResize)
        
        self.assertIsNotNone(protCropResize.outputParticles, "There was a problem with resize/crop v2 the particles")
        self.validateAcquisition(protCropResize.outputParticles)


class TestXmippML2D(TestXmippBase):
    """This class check if the protocol to classify with ML2D in Xmipp works properly."""
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
    
    def test_ml2d(self):
        print "Run ML2D"
        protML2D = self.newProtocol(XmippProtML2D, 
                                   numberOfReferences=2, maxIters=3, 
                                   numberOfMpi=2, numberOfThreads=2)
        protML2D.inputParticles.set(self.protImport.outputParticles)
        self.launchProtocol(protML2D)        
        
        self.assertIsNotNone(protML2D.outputClasses, "There was a problem with ML2D")  


class TestXmippCL2D(TestXmippBase):
    """This class check if the protocol to classify with CL2D in Xmipp works properly."""
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
    
    def test_cl2d(self):
        print "Run CL2D"
        protCL2D = self.newProtocol(XmippProtCL2D, 
                                   numberOfReferences=2, numberOfInitialReferences=1, 
                                   numberOfIterations=4, numberOfMpi=2)
        protCL2D.inputParticles.set(self.protImport.outputParticles)
        self.launchProtocol(protCL2D)      
        self.assertIsNotNone(protCL2D.outputClasses, "There was a problem with CL2D")


class TestXmippProtCL2DAlign(TestXmippBase):
    """This class check if the protocol to align particles in Xmipp works properly."""
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
    
    def test_xmippProtCL2DAlign(self):
        print "Run Only Align"
        CL2DAlign = self.newProtocol(XmippProtCL2DAlign, 
                                    maximumShift=5, numberOfIterations=5,
                                    numberOfMpi=4, numberOfThreads=1, useReferenceImage=False)
        CL2DAlign.inputParticles.set(self.protImport.outputParticles)
        self.launchProtocol(CL2DAlign)
        self.assertIsNotNone(CL2DAlign.outputParticles, "There was a problem with Only align2d")    


class TestXmippRotSpectra(TestXmippBase):
    """This class check if the protocol to calculate the rotational spectra from particles in Xmipp works properly."""
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
        cls.align2D = cls.runCL2DAlign(cls.protImport.outputParticles)
         
    def test_rotSpectra(self):
        print "Run Rotational Spectra"
        xmippProtRotSpectra = self.newProtocol(XmippProtRotSpectra, SomXdim=2, SomYdim=2)
        xmippProtRotSpectra.inputImages.set(self.align2D.outputParticles)
        self.launchProtocol(xmippProtRotSpectra)        
        self.assertIsNotNone(xmippProtRotSpectra.outputClasses, "There was a problem with Rotational Spectra")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        className = sys.argv[1]
        cls = globals().get(className, None)
        if cls:
            suite = unittest.TestLoader().loadTestsFromTestCase(cls)
            unittest.TextTestRunner(verbosity=2).run(suite)
        else:
            print "Test: '%s' not found." % className
    else:
        unittest.main()
