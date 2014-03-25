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
        cls.protImport = ProtImportParticles(pattern=pattern, samplingRate=samplingRate, checkStack=checkStack)
        cls.proj.launchProtocol(cls.protImport, wait=True)
        # check that input images have been imported (a better way to do this?)
        if cls.protImport.outputParticles is None:
            raise Exception('Import of images: %s, failed. outputParticles is None.' % pattern)
        return cls.protImport
    
#     def setupClassification(cls):
#         """ Method to setup classification Test Cases. """
#         setupProject(cls)
#         #TODO: Find a set of images to make this work, with this it does not
#         #pattern = getInputPath('images_LTA', '*.xmp')
#         #cls.protImport = cls.runImportParticles(pattern=pattern, samplingRate=5.6, checkStack=False)
#         
#         images = getInputPath('Images_Vol_ML3D/phantom_images', '*.xmp')
#         cls.protImport = cls.runImportParticles(pattern=images, samplingRate=1, checkStack=False)
    
    def launchXmippProtCL2DAlign(cls):
        xmippProtCL2DAlign = XmippProtCL2DAlign(maximumShift=5, numberOfIterations=5, 
                                 numberOfMpi=2, numberOfThreads=1, useReferenceImage=False)
        xmippProtCL2DAlign.inputParticles.set(test.protImport.outputParticles)
        cls.proj.launchProtocol(xmippProtCL2DAlign, wait=True)
        return xmippProtCL2DAlign


class TestXmippPreprocessParticles(TestXmippBase):
    """This class check if the protocol to preprocess particles in Xmipp works properly."""
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
    
    def testPreprocessPart(self):
        print "Run Preprocess particles"
        protPreproc = XmippProtPreprocessParticles(doRemoveDust=True, doNormalize=True, backRadius=48, doInvert=True,
                                              doThreshold=True, thresholdType=1)
        protPreproc.inputParticles.set(self.protImport.outputParticles)
        self.proj.launchProtocol(protPreproc, wait=True)
        
        self.assertIsNotNone(protPreproc.outputParticles, "There was a problem with preprocess particles")


class TestXmippCropResizeParticles(TestXmippBase):
    """This class check if the protocol to crop/resize particles in Xmipp works properly."""
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData('xmipp_tutorial')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 1.237, True)
    
    def testPreprocessPart(self):
        print "Run Preprocess particles"
        protCropResize = XmippProtCropResizeParticles(doResize=True, resizeOption=1, resizeDim=128, doWindow=True,
                                                      windowOperation=1, windowSize=256)
        protCropResize.inputParticles.set(self.protImport.outputParticles)
        self.proj.launchProtocol(protCropResize, wait=True)
        
        self.assertIsNotNone(protCropResize.outputParticles, "There was a problem with resize/crop the particles")


class TestXmippML2D(TestXmippBase):
    """This class check if the protocol to classify with ML2D in Xmipp works properly."""
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
    
    def testML2D(self):
        print "Run ML2D"
        protML2D = XmippProtML2D(numberOfReferences=2, maxIters=3, 
                                 numberOfMpi=2, numberOfThreads=2)
        protML2D.inputParticles.set(self.protImport.outputParticles)
        self.proj.launchProtocol(protML2D, wait=True)        
        
        self.assertIsNotNone(protML2D.outputClasses, "There was a problem with ML2D")  


class TestXmippCL2D(TestXmippBase):
    """This class check if the protocol to classify with CL2D in Xmipp works properly."""
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
    
    def testCL2D(self):
        print "Run CL2D"
        protCL2D = XmippProtCL2D(numberOfReferences=2, numberOfInitialReferences=1, 
                                 numberOfIterations=4, numberOfMpi=2)
        protCL2D.inputImages.set(self.protImport.outputParticles)
        self.proj.launchProtocol(protCL2D, wait=True)        
        self.assertIsNotNone(protCL2D.outputClasses, "There was a problem with CL2D")


class TestXmippProtCL2DAlign(TestXmippBase):
    """This class check if the protocol to align particles in Xmipp works properly."""
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
    
    def testXmippProtCL2DAlign(self):
        print "Run Only Align"
        xmippProtCL2DAlign = self.launchXmippProtCL2DAlign()
        self.assertIsNotNone(xmippProtCL2DAlign.outputParticles, "There was a problem with Only align2d")    


# class TestXmippRotSpectra(TestXmippBase):
#     """This class check if the protocol to calculate the rotational spectra from particles in Xmipp works properly."""
#     @classmethod
#     def setUpClass(cls):
#         setupClassification(cls)
#         
#     def testRotSpectra(self):
#         print "Run Rotational Spectra"
#         xmippProtCL2DAlign = launchXmippProtCL2DAlign(self)
#         
#         xmippProtRotSpectra = XmippProtRotSpectra(SomXdim=2, SomYdim=2)
#         xmippProtRotSpectra.inputImages.set(xmippProtCL2DAlign.outputParticles)
#         self.proj.launchProtocol(xmippProtRotSpectra, wait=True)        
#         self.assertIsNotNone(xmippProtRotSpectra.outputClasses, "There was a problem with Rotational Spectra")


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
