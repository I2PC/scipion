# **************************************************************************
# *
# * Authors:    Laura del Cano (ldelcano@cnb.csic.es)
# *             Josue Gomez Blanco (jgomez@cnb.csic.es)
# *             Jose Gutierrez (jose.gutierrez@cnb.csic.es)
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

import sys
import unittest

from pyworkflow.utils import redStr, greenStr, magentaStr
from pyworkflow.tests import *
from pyworkflow.em import *
from pyworkflow.em.packages.xmipp3 import *
from pyworkflow.em.packages.xmipp3 import XmippFilterHelper as xfh
from pyworkflow.em.packages.xmipp3 import XmippResizeHelper as xrh
from pyworkflow.em.packages.xmipp3.protocol_preprocess import (
    OP_COLUNM, OP_DOTPRODUCT, OP_MULTIPLY, OP_SQRT, OP_RADIAL, OP_ROW)


# Some utility functions to import particles that are used
# in several tests.
class TestXmippBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='xmipp_tutorial'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.particlesFn = cls.dataset.getFile('particles')
        cls.particlesDir = cls.dataset.getFile('particlesDir')
        cls.volumesFn = cls.dataset.getFile('volumes')
        cls.volumesDir = cls.dataset.getFile('volumesDir')
        cls.averagesFn = cls.dataset.getFile('averages')
        cls.averagesDir = cls.dataset.getFile('averagesDir')
    
    @classmethod
    def runImportParticles(cls, pattern, samplingRate, checkStack=False,
                           phaseFlip=False):
        """ Run an Import particles protocol. """
        cls.protImport = cls.newProtocol(ProtImportParticles, 
                                         filesPath=pattern,
                                         samplingRate=samplingRate,
                                         checkStack=checkStack,
                                         haveDataBeenPhaseFlipped=phaseFlip)
        print '_label: ', cls.protImport._label
        cls.launchProtocol(cls.protImport)
        # check that input images have been imported (a better way to do this?)
        if cls.protImport.outputParticles is None:
            raise Exception('Import of images: %s, failed. outputParticles is None.' % pattern)
        return cls.protImport
    
    @classmethod
    def runImportAverages(cls, pattern, samplingRate, checkStack=False):
        """ Run an Import particles protocol. """
        cls.protImportAvg = cls.newProtocol(ProtImportAverages,
                                            filesPath=pattern,
                                            samplingRate=samplingRate,
                                            checkStack=checkStack)
        print '_label: ', cls.protImportAvg._label
        cls.launchProtocol(cls.protImportAvg)
        # check that input images have been imported (a better way to do this?)
        if cls.protImportAvg.outputAverages is None:
            raise Exception('Import of averages: %s, failed. outputAverages is None.' % pattern)
        return cls.protImportAvg
    
    @classmethod
    def runImportVolume(cls, pattern, samplingRate, checkStack=False):
        """ Run an Import particles protocol. """
        cls.protImport = cls.newProtocol(ProtImportVolumes, 
                                         filesPath=pattern,
                                         samplingRate=samplingRate,
                                         checkStack=checkStack)
        print '_label: ', cls.protImport._label
        cls.launchProtocol(cls.protImport)
        # check that input images have been imported (a better way to do this?)
        if cls.protImport.outputVolume is None:
            raise Exception('Import of volume: %s, failed. outputVolume is None.' % pattern)
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
                                           numberOfClasses=4, maxIters=3, doMlf=False,
                                           numberOfMpi=3, numberOfThreads=2)
        cls.ProtClassify.inputParticles.set(particles)
        cls.launchProtocol(cls.ProtClassify)
        return cls.ProtClassify

    @classmethod
    def runCreateMask(cls, samplingRate, size):
        cls.protMask = cls.newProtocol(XmippProtCreateMask2D,
                                     samplingRate = samplingRate,
                                     size= size,
                                     geo=0, radius=-1 )
        cls.protMask.setObjLabel('circular mask')
        cls.launchProtocol(cls.protMask)
        return cls.protMask


class TestXmippCreateMask2D(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData()
        cls.protImport = cls.runImportParticles(cls.particlesFn, 1.237, True)
        cls.samplingRate = cls.protImport.outputParticles.getSamplingRate()
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
        print "launching protMask2"
        self.launchProtocol(protMask2)
        print "assert...."
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
                                     geo=5, innerRadius=2, outerRadius=12, borderDecay=2 )
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
        self.assertAlmostEquals(protMask1.outputParticles.getSamplingRate(), 
                                self.protImport.outputParticles.getSamplingRate(), "There was a problem with the sampling rate value for the apply user custom mask for particles")
        self.assertIsNotNone(protMask1.outputParticles, "There was a problem with apply circular mask for particles")
    
    def testApplyBoxMask(self):
        print "Run apply box mask for particles"
        protMask2 = self.newProtocol(XmippProtMaskParticles,
                                     source=0, geo=1, boxSize=-1,
                                     fillType=1 )
        protMask2.inputParticles.set(self.protImport.outputParticles)
        protMask2.setObjLabel('box mask')
        self.launchProtocol(protMask2)
        self.assertAlmostEquals(protMask2.outputParticles.getSamplingRate(), 
                                self.protImport.outputParticles.getSamplingRate(), "There was a problem with the sampling rate value for the apply user custom mask for particles")
        self.assertIsNotNone(protMask2.outputParticles, "There was a problem with apply boxed mask for particles")
    
    def testApplyCrownMask(self):
        print "Run apply crown mask for particles"
        protMask3 = self.newProtocol(XmippProtMaskParticles,
                                     source=0, geo=2, innerRadius=2, outerRadius=12,
                                     fillType=2 )
        protMask3.inputParticles.set(self.protImport.outputParticles)
        protMask3.setObjLabel('crown mask')
        self.launchProtocol(protMask3)
        self.assertAlmostEquals(protMask3.outputParticles.getSamplingRate(), 
                                self.protImport.outputParticles.getSamplingRate(), "There was a problem with the sampling rate value for the apply user custom mask for particles")
        self.assertIsNotNone(protMask3.outputParticles, "There was a problem with apply crown mask for particles")
        
    def testApplyGaussianMask(self):
        print "Run apply gaussian mask for particles"
        protMask4 = self.newProtocol(XmippProtMaskParticles,
                                     source=0, geo=3, sigma=-1,
                                     fillType=3 )
        protMask4.inputParticles.set(self.protImport.outputParticles)
        protMask4.setObjLabel('gaussian mask')
        self.launchProtocol(protMask4)
        self.assertAlmostEquals(protMask4.outputParticles.getSamplingRate(), 
                                self.protImport.outputParticles.getSamplingRate(), "There was a problem with the sampling rate value for the apply user custom mask for particles")
        self.assertIsNotNone(protMask4.outputParticles, "There was a problem with apply gaussian mask for particles")
        
    def testApplyRaisedCosineMask(self):
        print "Run apply raised cosine mask for particles"
        protMask5 = self.newProtocol(XmippProtMaskParticles,
                                     source=0, geo=4, innerRadius=2, outerRadius=12,
                                     fillType=0, fillValue=5 )
        protMask5.inputParticles.set(self.protImport.outputParticles)
        protMask5.setObjLabel('raised cosine mask')
        self.launchProtocol(protMask5)
        self.assertAlmostEquals(protMask5.outputParticles.getSamplingRate(), 
                                self.protImport.outputParticles.getSamplingRate(), "There was a problem with the sampling rate value for the apply user custom mask for particles")
        self.assertIsNotNone(protMask5.outputParticles, "There was a problem with apply raised cosine mask for particles")
        
    def testApplyRaisedCrownMask(self):
        print "Run apply raised crown mask for particles"
        protMask6 = self.newProtocol(XmippProtMaskParticles,
                                     source=0, geo=5, innerRadius=2, outerRadius=12, borderDecay=2,
                                     fillType=1 )
        protMask6.inputParticles.set(self.protImport.outputParticles)
        protMask6.setObjLabel('raised crown mask')
        self.launchProtocol(protMask6)
        self.assertAlmostEquals(protMask6.outputParticles.getSamplingRate(), 
                                self.protImport.outputParticles.getSamplingRate(), "There was a problem with the sampling rate value for the apply user custom mask for particles")
        
        self.assertIsNotNone(protMask6.outputParticles, "There was a problem with apply raised crown mask for particles")

    def testApplyUserMask(self):
        print "Run apply user mask for particles"
        # Create MASK
        protMask01 = self.newProtocol(XmippProtCreateMask2D,
                                     samplingRate=1.237, 
                                     size=500, 
                                     geo=0, radius=225)
        protMask01.setObjLabel('circular mask')
        self.launchProtocol(protMask01)
        self.assertIsNotNone(protMask01.outputMask, "There was a problem with apply user custom mask for particles")
        # Apply MASK
        protMask02 = self.newProtocol(XmippProtMaskParticles,
                                     source=1,
                                     fillType=1 )
        protMask02.inputParticles.set(self.protImport.outputParticles)
        protMask02.inputMask.set(protMask01.outputMask)
        protMask02.setObjLabel('user custom mask')
        self.launchProtocol(protMask02)
        self.assertAlmostEquals(protMask02.outputParticles.getSamplingRate(), 
                                self.protImport.outputParticles.getSamplingRate(), "There was a problem with the sampling rate value for the apply user custom mask for particles")
        
        self.assertIsNotNone(protMask02.outputParticles, "There was a problem with apply user custom mask for particles")
        


class TestXmippScreenParticles(TestXmippBase):
    """This class check if the protocol to classify particles by their similarity to discard outliers work properly"""
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 1.237, True)
        cls.samplingRate = cls.protImport.outputParticles.getSamplingRate()
        cls.size = 20
    
    def test_screenPart(self):
        from itertools import izip
        print 'Running Screen particles test'
        xpsp = XmippProtScreenParticles  # short notation
        # First test for check I/O. Input and Output SetOfParticles must be equal sized if not rejection is selected
        print '--> Running Screen without rejection'
        protScreenNone = self.newProtocol(xpsp, autoParRejection=xpsp.REJ_NONE)
        protScreenNone.inputParticles.set(self.protImport.outputParticles)
        self.launchProtocol(protScreenNone)
        self.assertIsNotNone(protScreenNone.outputParticles, 'Output has not been produced')
        print '\t --> Output is not None'
        self.assertEqual(len(protScreenNone.outputParticles), len(self.protImport.outputParticles), "Input and Output Set Of Particles don't have same size")
        print '\t --> Input/Output sets sizes are equal (%s)' % len(protScreenNone.outputParticles)
        
        for x, y in izip(self.protImport.outputParticles, protScreenNone.outputParticles):
            print "\t      compare %s with %s" % (x, y)
            self.assertEqual(x.getObjId(), y.getObjId(), "Particles differ")
            self.assertEqual(x.getSamplingRate(), y.getSamplingRate(), "Particle sampling rate differ")
        print '\t --> Input/Output sets contain the same particles'
        
        # After this, we check for errors in method with particle rejection by ZScore
        print "--> Running Screen with rejection to maxZScore upper than 2.6"
        protScreenZScore = self.newProtocol(xpsp, autoParRejection=xpsp.REJ_MAXZSCORE,
                                      maxZscore=2.6)
        protScreenZScore.inputParticles.set(self.protImport.outputParticles)
        self.launchProtocol(protScreenZScore)
        self.assertIsNotNone(protScreenZScore.outputParticles, "Output has not been produced")
        print '\t --> Output is not None'
        self.assertEqual(len(protScreenZScore.outputParticles), 71, "Output Set Of Particles must be 71, but %s found" % len(protScreenZScore.outputParticles))
        print '\t --> Output set size is correct (%s)' % len(protScreenZScore.outputParticles)
        
        for x in protScreenZScore.outputParticles:
            self.assertLess(x._xmipp_zScore.get(), 2.6, "Particle with id (%s) has a ZScore of %s, upper than supposed threshold %s" % (x.getObjId(), x._xmipp_zScore.get(), 2.6))
        print '\t --> Output particles are below the ZScore threshold'
        
        # Finally, we check for errors in method with particle rejection by percentage
        print "--> Running Screen with rejection of the 5% particles with the lowest ZScore"
        protScreenPercentage = self.newProtocol(xpsp, autoParRejection=xpsp.REJ_MAXZSCORE,
                                      maxZscore=2.6)
        protScreenPercentage.inputParticles.set(self.protImport.outputParticles)
        self.launchProtocol(protScreenPercentage)
        self.assertIsNotNone(protScreenPercentage.outputParticles, "Output has not been produced")
        print '\t --> Output is not None'
        self.assertEqual(len(protScreenPercentage.outputParticles), 71, "Output Set Of Particles must be 71, but %s found" % len(protScreenPercentage.outputParticles))
        print '\t --> Output set size is correct (%s)' % len(protScreenPercentage.outputParticles)
        
        for x, y in izip(protScreenZScore.outputParticles, protScreenPercentage.outputParticles):
            print "\t      compare %s with %s" % (x, y)
            self.assertEqual(x.getObjId(), y.getObjId(), "Particles differ")
        print '\t --> Particles rejected using maxZScore(2.6) method and percentage(5%) one are the same'


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
        
        if self.protImport.outputParticles.hasAlignment():
            from itertools import izip
            for x, y in izip(self.protImport.outputParticles.get(),
                             protPreproc.outputParticles.get()):
                print "compare ", x , " with ", y
                self.assertEquals(x.getAlignment(), y.getAlignment(),
                                  "Alignment wrong")
                
        self.assertAlmostEquals(protPreproc.outputParticles.getSamplingRate(), 
                          self.protImport.outputParticles.getSamplingRate(),
                                "There was a problem with the sampling rate "
                                " in the preprocess particles")

        self.assertIsNotNone(protPreproc.outputParticles,
                             "There was a problem with preprocess particles")


class TestXmippCropResizeParticles(TestXmippBase):
    """Check protocol crop/resize particles from Xmipp."""
    @classmethod
    def setUpClass(cls):
        print "\n", greenStr(" Crop/Resize Set Up - Collect data ".center(75, '-'))
        setupTestProject(cls)
        TestXmippBase.setData('xmipp_tutorial')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 1.237, True)
        cls.acquisition = cls.protImport.outputParticles.getAcquisition()

    def launch(self, **kwargs):
        "Launch XmippProtCropResizeParticles and return output particles."
        print magentaStr("\n==> Crop/Resize input params: %s" % kwargs)
        prot = self.newProtocol(XmippProtCropResizeParticles, **kwargs)
        prot.inputParticles.set(self.protImport.outputParticles)
        self.launchProtocol(prot)
        self.assertTrue(
            hasattr(prot, "outputParticles") and prot.outputParticles is not None,
            "There was a problem applying resize/crop to the particles")
        self.assertAlmostEqual(prot.outputParticles.getAcquisition().getVoltage(),
                               self.acquisition.getVoltage())
        return prot.outputParticles  # for more tests

    def test_newSizeAndCrop(self):
        inP = self.protImport.outputParticles  # short notation
        newSize = 128
        outP = self.launch(doResize=True, resizeOption=xrh.RESIZE_DIMENSIONS,
                           resizeDim=newSize,
                           doWindow=True, windowOperation=xrh.WINDOW_OP_CROP)

        self.assertEqual(newSize, outP.getDim()[0],
                         "Output particles dimension should be equal to %d" % newSize)
        self.assertAlmostEqual(outP.getSamplingRate(),
                               inP.getSamplingRate() * (inP.getDim()[0] / float(newSize)))

        # All other attributes remain the same. For the set:
        self.assertTrue(outP.equalAttributes(
            inP, ignore=['_mapperPath', '_samplingRate', '_firstDim'], verbose=True))
        # And for its individual particles too:
        self.assertTrue(outP.equalItemAttributes(
            inP, ignore=['_filename', '_index', '_samplingRate'], verbose=True))

    def test_factorAndWindow(self):
        inP = self.protImport.outputParticles  # short notation
        outP = self.launch(doResize=True, resizeOption=xrh.RESIZE_FACTOR,
                           resizeFactor=0.5,
                           doWindow=True, windowOperation=xrh.WINDOW_OP_WINDOW,
                           windowSize=500)

        # Since the images were resized by a factor 0.5 (downsampled), the new
        # pixel size (painfully called "sampling rate") should be 2x.
        self.assertAlmostEqual(outP.getSamplingRate(), inP.getSamplingRate() * 2)
        # After the window operation, the dimensions should be the same.
        self.assertEqual(inP.getDim(), outP.getDim())

        # All other attributes remain the same. For the set:
        self.assertTrue(outP.equalAttributes(
            inP, ignore=['_mapperPath', '_samplingRate'], verbose=True))
        # And for its individual particles too:
        self.assertTrue(outP.equalItemAttributes(
            inP, ignore=['_filename', '_index', '_samplingRate'], verbose=True))
    
    def test_pyramid(self):
        inP = self.protImport.outputParticles  # short notation
        outP = self.launch(doResize=True, resizeOption=xrh.RESIZE_PYRAMID,
                           resizeLevel=1)

        # Since the images were expanded by 2**resizeLevel (=2) the new
        # pixel size (painfully called "sampling rate") should be 0.5x.
        self.assertAlmostEqual(outP.getSamplingRate(), inP.getSamplingRate() * 0.5)
        # We did no window operation, so the dimensions will have doubled.
        self.assertAlmostEqual(outP.getDim()[0], inP.getDim()[0] * 2)

        # All other attributes remain the same. For the set:
        self.assertTrue(outP.equalAttributes(
            inP, ignore=['_mapperPath', '_samplingRate', '_firstDim'], verbose=True))
        # And for its individual particles too:
        self.assertTrue(outP.equalItemAttributes(
            inP, ignore=['_filename', '_index', '_samplingRate'], verbose=True))
        

class TestXmippCropResizeWAngles(TestXmippBase):
    """Check protocol crop/resize particles from Xmipp."""
    @classmethod
    def setUpClass(cls):
        print "\n", greenStr(" Crop/Resize Set Up - Collect data ".center(75, '-'))
        setupTestProject(cls)
        TestXmippBase.setData('relion_tutorial')
    
    def launch(self, **kwargs):
        "Launch XmippProtCropResizeParticles and return output particles."
        print magentaStr("\n==> Crop/Resize input params: %s" % kwargs)
        prot = self.newProtocol(XmippProtCropResizeParticles, **kwargs)
#         prot.inputParticles.set(self.protImport.outputParticles)
        self.launchProtocol(prot)
        self.assertTrue(
            hasattr(prot, "outputParticles") and prot.outputParticles is not None,
            "There was a problem applying resize/crop to the particles")
        return prot.outputParticles  # for more tests

    def test_CropResizeWAngles(self):
        print "Import Set of particles with angles"
        prot1 = self.newProtocol(ProtImportParticles,
                                 objLabel='from scipion (to-reconstruct)',
                                 importFrom=ProtImportParticles.IMPORT_FROM_SCIPION,
                                 sqliteFile=self.dataset.getFile('import/case2/particles.sqlite'),
                                 magnification=10000,
                                 samplingRate=7.08
                                 )
        self.launchProtocol(prot1)
        
        inP = prot1.outputParticles  # short notation
        newSize = 30
        factor = (inP.getDim()[0] / float(newSize))
        outP = self.launch(doResize=True, resizeOption=xrh.RESIZE_DIMENSIONS,
                           resizeDim=newSize,inputParticles=inP,
                           doWindow=True, windowOperation=xrh.WINDOW_OP_CROP)

        self.assertEqual(newSize, outP.getDim()[0],
                         "Output particles dimension should be equal to %d"
                         % newSize)
        self.assertAlmostEqual(outP.getSamplingRate(),
                               inP.getSamplingRate() * factor)

        # All other attributes remain the same. For the set:
        ignoreList = ['_mapperPath', '_samplingRate', '_firstDim']
        self.assertTrue(outP.equalAttributes(inP, ignore=ignoreList,
                                             verbose=True))

        # Check the scale factor is correctly applied to coordinates and
        # transform matrix
        for inPart, outPart in izip(inP, outP):
            coordIn = inPart.getCoordinate().getX()
            coordOut = outPart.getCoordinate().getX()
            self.assertAlmostEqual(coordIn, coordOut*factor, delta=2)

            tIn = inPart.getTransform()
            tOut = outPart.getTransform()
            tOut.scaleShifts(factor)
            mIn = tIn.getMatrix()
            mOut = tOut.getMatrix()
            self.assertTrue(np.allclose(mIn, mOut),
                            msg='Matrices not equal: %s, %s' % (mIn, mOut))


class TestXmippFilterParticles(TestXmippBase):
    """Check the proper behavior of Xmipp's filter particles protocol."""

    @classmethod
    def setUpClass(cls):
        print "\n", greenStr(" Set Up - Collect data ".center(75, '-'))
        setupTestProject(cls)
        TestXmippBase.setData('xmipp_tutorial')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 1.237,
                                                True, True)

    def test_filterParticles(self):
        print "\n", greenStr(" Filter Particles ".center(75, '-'))

        def test(parts=self.protImport.outputParticles, **kwargs):
            "Launch XmippProtFilterParticles on parts and check results."
            print magentaStr("\n==> Input params: %s" % kwargs)
            prot = self.newProtocol(XmippProtFilterParticles, **kwargs)
            prot.inputParticles.set(parts)
            self.launchProtocol(prot)
            self.assertIsNotNone(prot.outputParticles,
                                 "There was a problem with filter particles")
            self.assertTrue(prot.outputParticles.equalAttributes(
                parts, ignore=['_mapperPath'], verbose=True))
            # Compare the individual particles too.
            self.assertTrue(prot.outputParticles.equalItemAttributes(
                parts, ignore=['_filename', '_index'], verbose=True))

        # Check a few different cases.
        test(filterSpace=FILTER_SPACE_FOURIER, lowFreq=0.1, highFreq=0.25)
        test(filterSpace=FILTER_SPACE_REAL, filterModeReal=xfh.FM_MEDIAN)
        # For wavelets, we need the input's size to be a power of 2
        print magentaStr("\n==> Resizing particles to 256 pixels")
        protResize = self.newProtocol(XmippProtCropResizeParticles,
                                      doResize=True,
                                      resizeOption=xrh.RESIZE_DIMENSIONS,
                                      resizeDim=256)
        protResize.inputParticles.set(self.protImport.outputParticles)
        self.launchProtocol(protResize)
        test(parts=protResize.outputParticles, filterSpace=FILTER_SPACE_WAVELET,
             filterModeWavelets=xfh.FM_DAUB12, waveletMode=xfh.FM_REMOVE_SCALE)


class TestXmippOperateParticles(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        print "\n", greenStr(" Set Up - Collect data ".center(75, '-'))
        setupTestProject(cls)
        TestXmippBase.setData('xmipp_tutorial')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 1.237,
                                                True, True)

    def launchSet(self, **kwargs):
        "Launch XmippProtImageOperateParticles and return output volumes."
        print magentaStr("\n==> Operate set of volumes input params: %s" % kwargs)
        prot = XmippProtImageOperateParticles()
        prot.operation.set(kwargs.get('operation', 1))
        prot.inputParticles.set(self.protImport.outputParticles)
        prot.setObjLabel(kwargs.get('objLabel', None))
        prot.isValue.set(kwargs.get('isValue', False))
        prot.inputParticles2.set(kwargs.get('particles2', None))
        prot.value.set(kwargs.get('value', None))
        prot.intValue.set(kwargs.get('intValue', None))
        
        self.proj.launchProtocol(prot, wait=True)
        self.assertTrue(hasattr(prot, "outputParticles") and
                        prot.outputParticles is not None,
                        "There was a problem producing the output")
        return prot.outputParticles
        
    def testMultiplyVolSets(self):
        part2 = self.protImport.outputParticles  # short notation
        prot1 = self.launchSet(operation=OP_MULTIPLY,
                               objLabel='Multiply two SetOfParticles',
                               particles2=part2)

    def testMultiplyValue(self):
        prot2 = self.launchSet(operation=OP_MULTIPLY,
                               isValue=True,
                               objLabel='Multiply by a Value',
                               value=2.5)
    
    def testDotProduct(self):
        part2 = self.protImport.outputParticles  # short notation
        prot3 = self.launchSet(operation=OP_DOTPRODUCT,
                               objLabel='Dot Product',
                               particles2=part2)

    def testSqrt(self):
        prot4 = self.launchSet(operation=OP_SQRT,
                               objLabel='Sqrt')


class TestXmippML2D(TestXmippBase):
    """ This class check if the protocol to classify with ML2D in Xmipp works
    properly."""
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
    
    def test_ml2d(self):
        print "Run ML2D"
        protML2D = self.newProtocol(XmippProtML2D, 
                                   numberOfClasses=2, maxIters=3,
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
        cls.protImportAvgs = cls.runImportAverages(cls.particlesDir + '/img00007[1-4].spi', 3.5)
    
    def test_cl2d(self):
        print "Run CL2D"
        # Run CL2D with random class and core analysis
        protCL2DRandomCore = self.newProtocol(XmippProtCL2D,
                                   numberOfClasses=2, numberOfInitialClasses=1,
                                   numberOfIterations=4, numberOfMpi=2)
        protCL2DRandomCore.inputParticles.set(self.protImport.outputParticles)
        protCL2DRandomCore.setObjLabel("CL2D with random class and core analysis")
        self.launchProtocol(protCL2DRandomCore)
        self.assertIsNotNone(protCL2DRandomCore.outputClasses, "There was a problem with CL2D with random class and core analysis")

        # Run CL2D with random class and no core analysis
        protCL2DRandomNoCore = self.newProtocol(XmippProtCL2D,
                                   numberOfClasses=2, numberOfInitialClasses=1,
                                   doCore=False, numberOfIterations=4, numberOfMpi=2)
        protCL2DRandomNoCore.inputParticles.set(self.protImport.outputParticles)
        protCL2DRandomNoCore.setObjLabel("CL2D with random class and no core analysis")
        self.launchProtocol(protCL2DRandomNoCore)
        self.assertIsNotNone(protCL2DRandomNoCore.outputClasses, "There was a problem with CL2D with random class and no core analysis")

        # Run CL2D with initial classes and core analysis
        protCL2DInitialCore = self.newProtocol(XmippProtCL2D,
                                   numberOfClasses=4, randomInitialization=False,
                                   numberOfIterations=4, numberOfMpi=2)
        protCL2DInitialCore.inputParticles.set(self.protImport.outputParticles)
        protCL2DInitialCore.initialClasses.set(self.protImportAvgs.outputAverages)
        protCL2DInitialCore.setObjLabel("CL2D with initial class and core analysis")
        self.launchProtocol(protCL2DInitialCore)
        self.assertIsNotNone(protCL2DInitialCore.outputClasses, "There was a problem with CL2D with initial class and core analysis")


class TestXmippProtCL2DAlign(TestXmippBase):
    """This class check if the protocol to align particles in Xmipp works properly."""
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
    
    def test_xmippProtCL2DAlign(self):
        print "Run Only Align"
        # Run test without image reference
        CL2DAlignNoRef = self.newProtocol(XmippProtCL2DAlign,
                                    maximumShift=5, numberOfIterations=5,
                                    numberOfMpi=4, numberOfThreads=1, useReferenceImage=False)
        CL2DAlignNoRef.setObjLabel("CL2D Align without reference")
        CL2DAlignNoRef.inputParticles.set(self.protImport.outputParticles)
        self.launchProtocol(CL2DAlignNoRef)
        # Check that output is generated
        self.assertIsNotNone(CL2DAlignNoRef.outputParticles, "There was a problem generating output particles")
        # Check that it has alignment matrix
        self.assertTrue(CL2DAlignNoRef.outputParticles.hasAlignment2D(), "Output particles do not have alignment 2D")

        CL2DAlignRef = self.newProtocol(XmippProtCL2DAlign,
                                    maximumShift=5, numberOfIterations=5,
                                    numberOfMpi=4, numberOfThreads=1, useReferenceImage=True)
        CL2DAlignRef.setObjLabel("CL2D Align with reference")
        CL2DAlignRef.inputParticles.set(self.protImport.outputParticles)
        CL2DAlignRef.referenceImage.set(CL2DAlignNoRef.outputAverage)
        self.launchProtocol(CL2DAlignRef)
        # Check that output is generated
        self.assertIsNotNone(CL2DAlignRef.outputParticles, "There was a problem generating output particles")
        # Check that it has alignment matrix
        self.assertTrue(CL2DAlignRef.outputParticles.hasAlignment2D(), "Output particles do not have alignment 2D")


class TestXmippDenoiseParticles(TestXmippBase):
    """Check protocol Denoise Particles"""
    @classmethod
    def setUpClass(cls):
        from pyworkflow.em.packages.relion import ProtRelionClassify2D, ProtRelionPreprocessParticles
        # To denoise particles we need to import the particles and the
        # classes, and particles must be aligned with classes. As this
        # is the usual situation after a CL2D, we just run that protocol.

        # Set project and data
        setupTestProject(cls)
        cls.setData('mda')

        # Import particles
        psize = 3.50  # pixel size ("sampling rate"), in A/pixel
        cls.protImport = cls.runImportParticles(cls.particlesFn, psize)

        # Normalize them
        radiusInPixel = 42
        cls.protNormalize = cls.newProtocol(ProtRelionPreprocessParticles,
                                            doNormalize=True, backRadius=radiusInPixel)
        cls.protNormalize.inputParticles.set(cls.protImport.outputParticles)
        cls.launchProtocol(cls.protNormalize)
        diameterInA = 2 * radiusInPixel * psize
        cls.protCL2D = cls.newProtocol(ProtRelionClassify2D,
                                       doCTF=False, maskDiameterA=diameterInA,
                                       numberOfMpi=4, numberOfThreads=1)
        cls.protCL2D.numberOfClasses.set(4)
        cls.protCL2D.numberOfIterations.set(3)
        cls.protCL2D.inputParticles.set(cls.protNormalize.outputParticles)
        cls.launchProtocol(cls.protCL2D)

    def test_denoiseparticles(self):
        print """
*****************************
| Note: This part of the test may last for several minutes,
|       building a PCA basis for denoising is time expensive.
*****************************
"""
        protDenoise = self.newProtocol(XmippProtDenoiseParticles)
        protDenoise.inputParticles.set(self.protImport.outputParticles)
        protDenoise.inputClasses.set(self.protCL2D.outputClasses)
        self.launchProtocol(protDenoise)
        # We check that protocol generates output
        self.assertIsNotNone(protDenoise.outputParticles,
                             "There was a problem generating output particles")


class TestXmippApplyAlignment(TestXmippBase):
    """This class checks if the protocol Apply Alignment works properly"""
    @classmethod
    def setUpClass(cls):
        # For apply alignment we need to import particles that have alignment 2D information
        setupTestProject(cls)
        TestXmippBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
        cls.align2D = cls.runCL2DAlign(cls.protImport.outputParticles)

    def test_apply_alignment(self):
        protApply = self.newProtocol(XmippProtApplyAlignment)
        protApply.inputParticles.set(self.align2D.outputParticles)
        self.launchProtocol(protApply)
        # We check that protocol generates output
        self.assertIsNotNone(protApply.outputParticles, "There was a problem generating output particles")
        # Check that output particles do not have alignment information
        self.assertFalse(protApply.outputParticles.hasAlignment(), "Output particles should not have alignment information")


#TODO: Check with JM if this test should go in here since it is not a Xmipp protocol.
class TestAlignmentAssign(TestXmippBase):
    """This class checks if the protocol Alignment Assign works properly"""
    @classmethod
    def setUpClass(cls):
        # For alignment assign we need a set of particles without alignment 2D information and other set who has alignment information
        setupTestProject(cls)
        TestXmippBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
        cls.align2D = cls.runCL2DAlign(cls.protImport.outputParticles)

    def test_alignment_assign_samesize(self):
        protAssign = self.newProtocol(ProtAlignmentAssign)
        protAssign.setObjLabel("Assign alignment of same size")
        protAssign.inputParticles.set(self.protImport.outputParticles)
        protAssign.inputAlignment.set(self.align2D.outputParticles)
        self.launchProtocol(protAssign)
        # We check that protocol generates output
        self.assertIsNotNone(protAssign.outputParticles, "There was a problem generating output particles")
        # Check that output particles do not have alignment information
        self.assertTrue(protAssign.outputParticles.hasAlignment(), "Output particles should have alignment information")

    def test_alignment_assign_othersize(self):
        protResize = self.newProtocol(XmippProtCropResizeParticles,
                                      doResize=True,
                                      resizeOption=xrh.RESIZE_DIMENSIONS,
                                      resizeDim=50)
        protResize.inputParticles.set(self.protImport.outputParticles)
        self.launchProtocol(protResize)
        protAssign = self.newProtocol(ProtAlignmentAssign)
        protAssign.setObjLabel("Assign alignment of different size")
        protAssign.inputParticles.set(protResize.outputParticles)
        protAssign.inputAlignment.set(self.align2D.outputParticles)
        self.launchProtocol(protAssign)
        # We check that protocol generates output
        self.assertIsNotNone(protAssign.outputParticles, "There was a problem generating output particles")
        #TODO: Add an assert to check that sampling rate and alignment matrix is ok


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
        xmippProtRotSpectra.inputParticles.set(self.align2D.outputParticles)
        self.launchProtocol(xmippProtRotSpectra)        
        self.assertIsNotNone(xmippProtRotSpectra.outputClasses, "There was a problem with Rotational Spectra")

    def test_rotSpectraMask(self):
        print "Run Rotational Spectra with Mask"
        protMask = self.runCreateMask(3.5, 100)
        xmippProtRotSpectra = self.newProtocol(XmippProtRotSpectra, useMask=True, SomXdim=2, SomYdim=2)
        xmippProtRotSpectra.inputParticles.set(self.align2D.outputParticles)
        xmippProtRotSpectra.useMask.set(True)
        xmippProtRotSpectra.Mask.set(protMask.outputMask)
        self.launchProtocol(xmippProtRotSpectra)
        self.assertIsNotNone(xmippProtRotSpectra.outputClasses, "There was a problem with Rotational Spectra")


class TestXmippKerdensom(TestXmippBase):
    """This class check if the protocol to calculate the kerdensom from particles in Xmipp works properly."""
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
        cls.align2D = cls.runCL2DAlign(cls.protImport.outputParticles)

    def test_kerdensom(self):
        print "Run Kerdensom"
        xmippProtKerdensom = self.newProtocol(XmippProtKerdensom, SomXdim=2, SomYdim=2)
        xmippProtKerdensom.inputParticles.set(self.align2D.outputParticles)
        self.launchProtocol(xmippProtKerdensom)
        self.assertIsNotNone(xmippProtKerdensom.outputClasses, "There was a problem with Kerdensom")

    def test_kerdensomMask(self):
        print "Run Kerdensom with a mask"
        protMask = self.runCreateMask(3.5, 100)
        xmippProtKerdensom = self.newProtocol(XmippProtKerdensom, SomXdim=2, SomYdim=2)
        xmippProtKerdensom.inputParticles.set(self.align2D.outputParticles)
        xmippProtKerdensom.useMask.set(True)
        xmippProtKerdensom.Mask.set(protMask.outputMask)
        self.launchProtocol(xmippProtKerdensom)
        self.assertIsNotNone(xmippProtKerdensom.outputClasses, "There was a problem with Kerdensom")


class TestXmippCompareReprojections(TestXmippBase):
    """This class check if the protocol compare reprojections in Xmipp works properly."""
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData('mda')
        cls.protImportPart = cls.runImportParticles(cls.particlesFn, 3.5)
        cls.protImportAvgs = cls.runImportAverages(cls.particlesFn, 3.5)
        cls.protImportVol = cls.runImportVolume(cls.volumesFn, 3.5)
        cls.protClassify = cls.runClassify(cls.protImportPart.outputParticles)
        cls.protProjMatch = cls.newProtocol(XmippProtProjMatch,
                                            doCTFCorrection=False,
                                            numberOfIterations=1,
                                            outerRadius=50,
                                            angSamplingRateDeg=5,
                                            symmetry="d6",
                                            numberOfMpi=4)
        cls.protProjMatch.inputParticles.set(cls.protImportAvgs.outputAverages)
        cls.protProjMatch.input3DReferences.set(cls.protImportVol.outputVolume)
        cls.launchProtocol(cls.protProjMatch)
    
    def test_particles1(self):
        print "Run Compare Reprojections from classes"
        prot = self.newProtocol(XmippProtCompareReprojections, 
                                        symmetryGroup="d6", numberOfMpi=5)
        prot.inputSet.set(self.protClassify.outputClasses)
        prot.inputVolume.set(self.protImportVol.outputVolume)
        self.launchProtocol(prot)
        self.assertIsNotNone(prot.outputParticles, "There was a problem with Compare Reprojections from classes")

    def test_particles2(self):
        print "Run Compare Reprojections from averages"
        prot = self.newProtocol(XmippProtCompareReprojections, 
                                        symmetryGroup="d6", numberOfMpi=5)
        prot.inputSet.set(self.protImportAvgs.outputAverages)
        prot.inputVolume.set(self.protImportVol.outputVolume)
        self.launchProtocol(prot)
        self.assertIsNotNone(prot.outputParticles, "There was a problem with Compare Reprojections from averages")

    def test_particles3(self):
        print "Run Compare Reprojections from projections with angles"
        prot = self.newProtocol(XmippProtCompareReprojections, 
                                        symmetryGroup="d6", numberOfMpi=5)
        prot.inputSet.set(self.protProjMatch.outputParticles)
        prot.inputVolume.set(self.protImportVol.outputVolume)
        self.launchProtocol(prot)
        self.assertIsNotNone(prot.outputParticles, "There was a problem with Compare Reprojections from projections with angles")


class TestXmippCreateGallery(TestXmippBase):
    """This class check if the protocol create gallery in Xmipp works properly."""
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData('mda')
        cls.protImportVol = cls.runImportVolume(cls.volumesFn, 3.5)
    
    def _createGallery(self, step, projections):
        prot = self.newProtocol(XmippProtCreateGallery,
                                symmetryGroup="d6",
                                rotStep=step, tiltStep=step)
        prot.inputVolume.set(self.protImportVol.outputVolume)
        self.launchProtocol(prot)
        outSet = getattr(prot, 'outputReprojections', None)
        self.assertIsNotNone(outSet, "There was a problem with create gallery")
        self.assertEqual(projections, outSet.getSize())

        return prot

    def test_step5(self):
        prot = self._createGallery(step=5, projections=131)

    def test_step10(self):
        prot = self._createGallery(step=10, projections=32)



class TestXmippBreakSym(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
    
    def test_AngBreakSymmetry(self):
        from tempfile import NamedTemporaryFile
        import pyworkflow.em.metadata as md
        
        fileTmp = NamedTemporaryFile(delete=False, suffix='.sqlite')
        partSet = SetOfParticles(filename=fileTmp.name)
        partSet.setAlignment(ALIGN_PROJ)
        # Populate the SetOfParticles with  images
        # taken from images.mrc file
        # and setting the previous alignment parameters
        
        m = np.array([[ 0.71461016, 0.63371837, -0.29619813, 15],
                      [ -0.61309201, 0.77128059, 0.17101008, 25],
                      [ 0.33682409, 0.059391174, 0.93969262, 35],
                      [ 0,          0,           0,           1]])
        p = Particle()
        p.setLocation(1, "kk.mrc")

        p.setTransform(Transform(m))
        partSet.append(p)
        partSet.write()
        
        print "import particles"
        protImport = self.newProtocol(ProtImportParticles, 
                                         sqliteFile=fileTmp.name, samplingRate=1, importFrom=4,
                                         checkStack=False, haveDataBeenPhaseFlipped=False)
        self.launchProtocol(protImport)
        
        print "Run AngBreakSymmetry particles"
        protBreakSym = self.newProtocol(XmippProtAngBreakSymmetry, symmetryGroup="i2")
        protBreakSym.inputParticles.set(protImport.outputParticles)
        self.launchProtocol(protBreakSym)
        os.chdir(protBreakSym._getPath())
        from pyworkflow.utils import runJob
        runJob(None, 'xmipp_angular_distance',
               "--ang1 images.xmd --ang2 input_particles.xmd --sym i2 --oroot kk",
               env=getEnviron())
        mdRober = md.MetaData("kk_vec_diff_hist.txt")
        objId = mdRober.firstObject()
        count = mdRober.getValue(md.MDL_COUNT, objId)
        
        self.assertEqual(count, 1, "There was a problem with break symmetry")
        os.unlink(fileTmp.name)


class TestXmippCorrectWiener2D(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData()
    
    def test_CorrectWiener(self):
        prot1 = self.newProtocol(ProtImportParticles,
                                 importFrom=ProtImportParticles.IMPORT_FROM_XMIPP3,
                                 mdFile=self.dataset.getFile('particles/sphere_128.xmd'),
                                 magnification=10000,
                                 samplingRate=1,
                                 haveDataBeenPhaseFlipped=False
                                 )
        self.launchProtocol(prot1)
        print "Run CTFCorrectWiener2D particles"
        protCorrect = self.newProtocol(XmippProtCTFCorrectWiener2D)
        protCorrect.inputParticles.set(prot1.outputParticles)
        self.launchProtocol(protCorrect)
        self.assertIsNotNone(protCorrect.outputParticles, "There was a problem with Wiener Correction")

        
class TestXmippSubtractProjection(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsRelion = DataSet.getDataSet('relion_tutorial')
    
    def test_subtract(self):
        protParts = self.newProtocol(ProtImportParticles,
                                     objLabel='from relion auto-refine',
                                     importFrom=ProtImportParticles.IMPORT_FROM_RELION,
                                     starFile=self.dsRelion.getFile('import/refine3d/extra/relion_it001_data.star'),
                                     magnification=10000,
                                     samplingRate=7.08,
                                     haveDataBeenPhaseFlipped=True
                                     )
        self.launchProtocol(protParts)
        self.assertEqual(60, protParts.outputParticles.getXDim())
        
        protVol = self.newProtocol(ProtImportVolumes,
                                   filesPath=self.dsRelion.getFile('volumes/reference.mrc'),
                                   samplingRate=7.08)
        self.launchProtocol(protVol)
        self.assertEqual(60, protVol.outputVolume.getDim()[0])
        
        protSubtract = self.newProtocol(XmippProtSubtractProjection)
        protSubtract.inputParticles.set(protParts.outputParticles)
        protSubtract.inputVolume.set(protVol.outputVolume)
        self.launchProtocol(protSubtract)
        self.assertIsNotNone(protSubtract.outputParticles, "There was a problem with subtract projection")
        

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
