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
        print '_label: ', cls.protImport._label
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


class TestXmippCreateMask2D(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData()
        cls.protImport = cls.runImportParticles(cls.particlesFn, 1.237, True)
        cls.samplingRate = cls.protImport.getSamplingRate()
        cls.size = 20
        
    def testCreateCircularMask(self):
        print "Run create circular mask for particles"
        protMask1 = self.newProtocol(XmippProtCreateMask2D,
                                     samplingRate = self.samplingRate, 
                                     size= self.size, 
                                     geo=0, radius=-1 )
        protMask1.setObjLabel('circular mask')
        self.launchProtocol(protMask1)
        self.assertIsNotNone(protMask1.outputMask, "There was a problem with create circular mask for particles")
    
    def testCreateBoxMask(self):
        print "Run create box mask for particles"
        protMask2 = self.newProtocol(XmippProtCreateMask2D,
                                     samplingRate = self.samplingRate, 
                                     size= self.size, 
                                     geo=1, boxSize=-1 )
        protMask2.setObjLabel('box mask')
        self.launchProtocol(protMask2)
        self.assertIsNotNone(protMask2.outputMask, "There was a problem with create boxed mask for particles")
    
    def testCreateCrownMask(self):
        print "Run create crown mask for particles"
        protMask3 = self.newProtocol(XmippProtCreateMask2D,
                                     samplingRate = self.samplingRate, 
                                     size= self.size, 
                                     geo=2, innerRadius=2, outerRadius=12 )
        protMask3.setObjLabel('crown mask')
        self.launchProtocol(protMask3)
        self.assertIsNotNone(protMask3.outputMask, "There was a problem with create crown mask for particles")
    
    def testCreateGaussianMask(self):
        print "Run create gaussian mask for particles"
        protMask4 = self.newProtocol(XmippProtCreateMask2D,
                                     samplingRate = self.samplingRate, 
                                     size= self.size, 
                                     geo=3, sigma=-1 )
        protMask4.setObjLabel('gaussian mask')
        self.launchProtocol(protMask4)
        self.assertIsNotNone(protMask4.outputMask, "There was a problem with create gaussian mask for particles")
    
    def testCreateRaisedCosineMask(self):
        print "Run create raised cosine mask for particles"
        protMask5 = self.newProtocol(XmippProtCreateMask2D,
                                     samplingRate = self.samplingRate, 
                                     size= self.size,
                                     geo=4, innerRadius=2, outerRadius=12 )
        protMask5.setObjLabel('raised cosine mask')
        self.launchProtocol(protMask5)
        self.assertIsNotNone(protMask5.outputMask, "There was a problem with create raised cosine mask for particles")
    
    def testCreateRaisedCrownMask(self):
        print "Run create raised crown mask for particles"
        protMask6 = self.newProtocol(XmippProtCreateMask2D,
                                     samplingRate = self.samplingRate, 
                                     size= self.size, 
                                     geo=5, innerRadius=2, outerRadius=12, decay=2 )
        protMask6.setObjLabel('raised crown mask')
        self.launchProtocol(protMask6)
        self.assertIsNotNone(protMask6.outputMask, "There was a problem with create raised crown mask for particles")
    

class TestXmippApplyMask2D(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData()
        cls.protImport = cls.runImportParticles(cls.particlesFn, 1.237, True)

    def testApplyCircularMask(self):
        print "Run apply circular mask for particles"
        protMask1 = self.newProtocol(XmippProtMaskParticles,
                                     source=0, geo=0, radius=-1,
                                     fillType=0, fillValue=5 )
        protMask1.inputParticles.set(self.protImport.outputParticles)
        protMask1.setObjLabel('circular mask')
        self.launchProtocol(protMask1)
        self.assertIsNotNone(protMask1.outputParticles, "There was a problem with apply circular mask for particles")
    
    def testApplyBoxMask(self):
        print "Run apply box mask for particles"
        protMask2 = self.newProtocol(XmippProtMaskParticles,
                                     source=0, geo=1, boxSize=-1,
                                     fillType=1 )
        protMask2.inputParticles.set(self.protImport.outputParticles)
        protMask2.setObjLabel('box mask')
        self.launchProtocol(protMask2)
        self.assertIsNotNone(protMask2.outputParticles, "There was a problem with apply boxed mask for particles")
    
    def testapplyCrownMask(self):
        print "Run apply crown mask for particles"
        protMask3 = self.newProtocol(XmippProtMaskParticles,
                                     source=0, geo=2, innerRadius=2, outerRadius=12,
                                     fillType=2 )
        protMask3.inputParticles.set(self.protImport.outputParticles)
        protMask3.setObjLabel('crown mask')
        self.launchProtocol(protMask3)
        self.assertIsNotNone(protMask3.outputParticles, "There was a problem with apply crown mask for particles")
    
    def testApplyGaussianMask(self):
        print "Run apply gaussian mask for particles"
        protMask4 = self.newProtocol(XmippProtMaskParticles,
                                     source=0, geo=3, sigma=-1,
                                     fillType=3 )
        protMask4.inputParticles.set(self.protImport.outputParticles)
        protMask4.setObjLabel('gaussian mask')
        self.launchProtocol(protMask4)
        self.assertIsNotNone(protMask4.outputParticles, "There was a problem with apply gaussian mask for particles")
    
    def testApplyRaisedCosineMask(self):
        print "Run apply raised cosine mask for particles"
        protMask5 = self.newProtocol(XmippProtMaskParticles,
                                     source=0, geo=4, innerRadius=2, outerRadius=12,
                                     fillType=0, fillValue=5 )
        protMask5.inputParticles.set(self.protImport.outputParticles)
        protMask5.setObjLabel('raised cosine mask')
        self.launchProtocol(protMask5)
        self.assertIsNotNone(protMask5.outputParticles, "There was a problem with apply raised cosine mask for particles")
    
    def testApplyRaisedCrownMask(self):
        print "Run apply raised crown mask for particles"
        protMask6 = self.newProtocol(XmippProtMaskParticles,
                                     source=0, geo=5, innerRadius=2, outerRadius=12, decay=2,
                                     fillType=1 )
        protMask6.inputParticles.set(self.protImport.outputParticles)
        protMask6.setObjLabel('raised crown mask')
        self.launchProtocol(protMask6)
        self.assertIsNotNone(protMask6.outputParticles, "There was a problem with apply raised crown mask for particles")
    
    def testApplyUserMask(self):
        print "Run apply user mask for particles"
        # Created MASK
        protMask01 = self.newProtocol(XmippProtCreateMask2D,
                                     samplingRate=1.237, 
                                     size=20, 
                                     geo=0, radius=-1 )
        protMask01.setObjLabel('circular mask')
        self.launchProtocol(protMask01)
        self.assertIsNotNone(protMask01.outputMask, "There was a problem with apply user custom mask for particles")
        #Applied MASK
        protMask02 = self.newProtocol(XmippProtMaskParticles,
                                     source=1,
                                     fillType=1 )
        protMask02.inputParticles.set(self.protImport.outputParticles)
        protMask02.inputMask.set(protMask01.outputMask)
        protMask02.setObjLabel('user custom mask')
        self.launchProtocol(protMask02)
        self.assertIsNotNone(protMask02.outputParticles, "There was a problem with apply user custom mask for particles")
    

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
        oldSize = 500.
        newSize = 128
        protCropResize = self.newProtocol(XmippProtCropResizeParticles, 
                                         doResize=True, resizeOption=1, resizeDim=newSize, 
                                         doWindow=True, windowOperation=0)
        input = self.protImport.outputParticles
        protCropResize.inputParticles.set(input)
        self.launchProtocol(protCropResize)
        
        output = protCropResize.outputParticles
        self.assertIsNotNone(output, "There was a problem with resize/crop the particles")
        self.validateAcquisition(output)
        self.assertEquals(newSize, output.getDimensions()[0])#, "Output particles dimension should be equal to %d" % newSize)
        self.assertAlmostEquals(output.getSamplingRate(), input.getSamplingRate()*(oldSize/newSize))
    
    def test_cropResizePart2(self):
        print "Run crop/resize particles v2"
        protCropResize = self.newProtocol(XmippProtCropResizeParticles,  
                                         doResize=True, resizeOption=2, resizeLevel=0.5, 
                                         doWindow=True, windowOperation=1, windowSize=500)
        input = self.protImport.outputParticles
        protCropResize.inputParticles.set(input)
        self.launchProtocol(protCropResize)
        output = protCropResize.outputParticles
        
        self.assertIsNotNone(output, "There was a problem with resize/crop v2 the particles")
        self.validateAcquisition(output)
        # Since the images where downsampled twice (factor 0.5)
        # the sampling rate should be the double
        self.assertAlmostEquals(output.getSamplingRate(), input.getSamplingRate()*2)
        # Since we have done windowing operation, the dimensions should be the same
        self.assertEquals(input.getDim(), output.getDim())

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
