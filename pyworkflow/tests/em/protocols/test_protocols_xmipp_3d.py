# **************************************************************************
# *
# * Authors:    Josue Gomez Blanco (jgomez@cnb.csic.es)
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
from itertools import izip

from pyworkflow.utils import redStr, greenStr, magentaStr
from pyworkflow.tests import *
from pyworkflow.em import *
from pyworkflow.em.packages.xmipp3 import *
from pyworkflow.em.packages.xmipp3 import (XmippFilterHelper as xfh,
                                           XmippResizeHelper as xrh)
from pyworkflow.em.packages.xmipp3.protocol_align_volume import (
    ALIGN_ALGORITHM_EXHAUSTIVE, ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL,
    ALIGN_ALGORITHM_LOCAL)
from pyworkflow.em.packages.xmipp3.protocol_preprocess import (
    OP_COLUNM, OP_DOTPRODUCT, OP_MULTIPLY, OP_SQRT, OP_RADIAL, OP_ROW)

class TestXmippBase(BaseTest):
    """ Some utility functions to import volumes that are used in several tests."""

    @classmethod
    def setData(cls, dataProject='xmipp_tutorial'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.volumes = cls.dataset.getFile('volumes')
        cls.vol1 = cls.dataset.getFile('vol1')
        cls.vol2 = cls.dataset.getFile('vol2')
        cls.vol3 = cls.dataset.getFile('vol3')
        cls.vol4 = cls.dataset.getFile('vol4')

    @classmethod
    def runImportVolumes(cls, pattern, samplingRate):
        """ Run an Import particles protocol. """
        cls.protImport = cls.newProtocol(ProtImportVolumes,
                                         filesPath=pattern,
                                         samplingRate=samplingRate)
        cls.launchProtocol(cls.protImport)
        return cls.protImport

    @classmethod
    def runImportMask(cls, pattern, samplingRate):
        """ Run an Import particles protocol. """
        cls.protImportMask = cls.newProtocol(ProtImportMask,
                                         maskPath=pattern,
                                             samplingRate=samplingRate)
        cls.launchProtocol(cls.protImportMask)
        return cls.protImportMask

    @classmethod
    def runImportParticles(cls, pattern, samplingRate, checkStack=False):
        """ Run an Import particles protocol. """
        cls.protImport = cls.newProtocol(ProtImportParticles,
                                         filesPath=pattern,
                                         samplingRate=samplingRate,
                                         checkStack=checkStack)
        cls.launchProtocol(cls.protImport)
        # check that input images have been imported (a better way to do this?)
        if cls.protImport.outputParticles is None:
            raise Exception('Import of images: %s, failed. outputParticles '
                            'is None.' % pattern)
        return cls.protImport

    @classmethod
    def runClassify(cls, particles):
        cls.ProtClassify = cls.newProtocol(XmippProtML2D,
                                           numberOfClasses=8, maxIters=4,
                                           doMlf=False,
                                           numberOfMpi=2, numberOfThreads=2)
        cls.ProtClassify.inputParticles.set(particles)
        cls.launchProtocol(cls.ProtClassify)
        return cls.ProtClassify


class TestXmippCreateMask3D(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData()
        cls.protImport = cls.runImportVolumes(cls.vol1, 9.896)

    def testCreateMask(self):
        print "Run create threshold mask from volume"
        protMask1 = self.newProtocol(XmippProtCreateMask3D,
                                     source=0, volumeOperation=0,
                                     threshold=0.4)
        protMask1.inputVolume.set(self.protImport.outputVolume)
        protMask1.setObjLabel('threshold mask')
        self.launchProtocol(protMask1)
        self.assertIsNotNone(protMask1.outputMask,
                             "There was a problem with create mask from volume")


        print "Run create segment mask from volume"
        protMask2 = self.newProtocol(XmippProtCreateMask3D,
                                     source=0, volumeOperation=1,
                                     segmentationType=3)
        protMask2.inputVolume.set(self.protImport.outputVolume)
        protMask2.setObjLabel('segmentation automatic')
        self.launchProtocol(protMask2)
        self.assertIsNotNone(protMask2.outputMask,
                             "There was a problem with create mask from volume")

    
        print "Run create mask from another mask"
        protMask3 = self.newProtocol(XmippProtCreateMask3D,
                                     source=0, volumeOperation=2,
                                     doMorphological=True, elementSize=3)
        protMask3.inputVolume.set(protMask1.outputMask)
        protMask3.setObjLabel('dilation mask')
        self.launchProtocol(protMask3)
        self.assertIsNotNone(protMask3.outputMask,
                             "There was a problem with mask from another mask")

        print "Run create mask from geometry"
        protMask4 = self.newProtocol(XmippProtCreateMask3D,
                                     source=1, size=64, samplingRate=9.89,
                                     geo=6, innerRadius=10, outerRadius=25,
                                     borderDecay=2)
        protMask4.setObjLabel('crown mask')
        self.launchProtocol(protMask4)
        self.assertIsNotNone(protMask4.outputMask,
                             "There was a problem with create mask from geometry")


class TestXmippApplyMask3D(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData()
        cls.protImport = cls.runImportVolumes(cls.vol1, 9.896)

    def testApplyCircularMask(self):
        print "Run apply circular mask for volumes"
        protMask1 = self.newProtocol(XmippProtMaskVolumes,
                                     source=0, geo=0, radius=-1,
                                     fillType=0, fillValue=5 )
        protMask1.inputVolumes.set(self.protImport.outputVolume)
        protMask1.setObjLabel('circular mask')
        self.launchProtocol(protMask1)
        self.assertAlmostEquals(protMask1.outputVol.getSamplingRate(), 
                                self.protImport.outputVolume.getSamplingRate(),
                                "There was a problem with the sampling rate value for the apply user custom mask for Volumes")
        self.assertIsNotNone(protMask1.outputVol,
                             "There was a problem with apply circular mask for Volumes")
    
    def testApplyBoxMask(self):
        print "Run apply box mask for Volumes"
        protMask2 = self.newProtocol(XmippProtMaskVolumes,
                                     source=0, geo=1, boxSize=-1,
                                     fillType=1 )
        protMask2.inputVolumes.set(self.protImport.outputVolume)
        protMask2.setObjLabel('box mask')
        self.launchProtocol(protMask2)
        self.assertAlmostEquals(protMask2.outputVol.getSamplingRate(), 
                                self.protImport.outputVolume.getSamplingRate(),
                                "There was a problem with the sampling rate value for the apply user custom mask for Volumes")
        self.assertIsNotNone(protMask2.outputVol,
                             "There was a problem with apply boxed mask for Volumes")
    
    def testapplyCrownMask(self):
        print "Run apply crown mask for Volumes"
        protMask3 = self.newProtocol(XmippProtMaskVolumes,
                                     source=0, geo=2, innerRadius=2, outerRadius=12,
                                     fillType=2 )
        protMask3.inputVolumes.set(self.protImport.outputVolume)
        protMask3.setObjLabel('crown mask')
        self.launchProtocol(protMask3)
        self.assertAlmostEquals(protMask3.outputVol.getSamplingRate(), 
                                self.protImport.outputVolume.getSamplingRate(),
                                "There was a problem with the sampling rate value for the apply user custom mask for Volumes")
        self.assertIsNotNone(protMask3.outputVol, "There was a problem with apply crown mask for Volumes")
        
    def testApplyGaussianMask(self):
        print "Run apply gaussian mask for Volumes"
        protMask4 = self.newProtocol(XmippProtMaskVolumes,
                                     source=0, geo=3, sigma=-1,
                                     fillType=3 )
        protMask4.inputVolumes.set(self.protImport.outputVolume)
        protMask4.setObjLabel('gaussian mask')
        self.launchProtocol(protMask4)
        self.assertAlmostEquals(protMask4.outputVol.getSamplingRate(), 
                                self.protImport.outputVolume.getSamplingRate(),
                                "There was a problem with the sampling rate value for the apply user custom mask for Volumes")
        self.assertIsNotNone(protMask4.outputVol, "There was a problem with apply gaussian mask for Volumes")
        
    def testApplyRaisedCosineMask(self):
        print "Run apply raised cosine mask for Volumes"
        protMask5 = self.newProtocol(XmippProtMaskVolumes,
                                     source=0, geo=4, innerRadius=2, outerRadius=12,
                                     fillType=0, fillValue=5 )
        protMask5.inputVolumes.set(self.protImport.outputVolume)
        protMask5.setObjLabel('raised cosine mask')
        self.launchProtocol(protMask5)
        self.assertAlmostEquals(protMask5.outputVol.getSamplingRate(), 
                                self.protImport.outputVolume.getSamplingRate(),
                                "There was a problem with the sampling rate value for the apply user custom mask for Volumes")
        self.assertIsNotNone(protMask5.outputVol, "There was a problem with apply raised cosine mask for Volumes")
        
    def testApplyRaisedCrownMask(self):
        print "Run apply raised crown mask for Volumes"
        protMask6 = self.newProtocol(XmippProtMaskVolumes,
                                     source=0, geo=5, innerRadius=2, outerRadius=12, borderDecay=2,
                                     fillType=1 )
        protMask6.inputVolumes.set(self.protImport.outputVolume)
        protMask6.setObjLabel('raised crown mask')
        self.launchProtocol(protMask6)
        self.assertAlmostEquals(protMask6.outputVol.getSamplingRate(), 
                                self.protImport.outputVolume.getSamplingRate(),
                                "There was a problem with the sampling rate value for the apply user custom mask for Volumes")
        
        self.assertIsNotNone(protMask6.outputVol, "There was a problem with apply raised crown mask for Volumes")
        
    def testApplyUserMask(self):
        print "Run apply user mask for Volumes"
        # Create MASK
        protMask01 = self.newProtocol(XmippProtCreateMask3D,
                                     source=1, size=64, samplingRate=9.89,
                                     geo=6, innerRadius=10, outerRadius=25, borderDecay=2)
        protMask01.setObjLabel('crown mask')
        self.launchProtocol(protMask01)
        self.assertIsNotNone(protMask01.outputMask,
                             "There was a problem with create mask from geometry")
        # Apply MASK
        protMask02 = self.newProtocol(XmippProtMaskVolumes,
                                     source=1,
                                     fillType=1 )
        protMask02.inputVolumes.set(self.protImport.outputVolume)
        protMask02.inputMask.set(protMask01.outputMask)
        protMask02.setObjLabel('user custom mask')
        self.launchProtocol(protMask02)
        self.assertAlmostEquals(protMask02.outputVol.getSamplingRate(), 
                                self.protImport.outputVolume.getSamplingRate(),
                                "There was a problem with the sampling rate value for the apply user custom mask for Volumes")
         
        self.assertIsNotNone(protMask02.outputVol,
                             "There was a problem with apply user custom mask for Volumes")
  

class TestXmippPreprocessVolumes(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData()
        cls.protImport1 = cls.runImportVolumes(cls.volumes, 9.896)
        cls.protImport2 = cls.runImportVolumes(cls.vol1, 9.896)
        cls.protImport2 = cls.runImportVolumes(cls.vol1, 9.896)
        #test symmetryze with mask
        dataProject='SymVirus'
        dataset = DataSet.getDataSet(dataProject)
        virusVol  = dataset.getFile('whole_vol_half')
        virusMaskCapsid = dataset.getFile('large_vol_half_th')
        virusMaskPenton = dataset.getFile('small_vol_half_th')
        cls.protImportVirus = cls.runImportVolumes(virusVol, 1)
        cls.protImportvirusMaskCapsid = cls.runImportMask(virusMaskCapsid, 1)
        cls.protImportvirusMaskPenton = cls.runImportMask(virusMaskPenton, 1)


    def testPreprocessVolumes(self):
        print "Run preprocess a volume"
        protPreprocessVol1 = XmippProtPreprocessVolumes(doChangeHand=True, doRandomize=True, doSymmetrize=True, symmetryGroup='d6',
                                                        doSegment=True, doNormalize=True, backRadius=20, doInvert=True,
                                                        doThreshold=True, thresholdType=1)
        protPreprocessVol1.inputVolumes.set(self.protImport2.outputVolume)
        self.proj.launchProtocol(protPreprocessVol1, wait=True)
        self.assertIsNotNone(protPreprocessVol1.outputVol, "There was a problem with a volume")

        print "Run preprocess a SetOfVolumes"
        protPreprocessVol2 = XmippProtPreprocessVolumes(doChangeHand=True, doRandomize=True, doSymmetrize=True, symmetryGroup='d6',
                                                        doSegment=True, doNormalize=True, backRadius=20, doInvert=True,
                                                        doThreshold=True, thresholdType=1)
        protPreprocessVol2.inputVolumes.set(self.protImport1.outputVolumes)
        self.proj.launchProtocol(protPreprocessVol2, wait=True)
        self.assertIsNotNone(protPreprocessVol2.outputVol, "There was a problem with preprocess a SetOfVolumes")

        print "Run preprocess a volume using mask_1 in the symmetrization"
        protPreprocessVol3 = XmippProtPreprocessVolumes(doChangeHand=False, doRandomize=False,
                                                        doSymmetrize=True, symmetryGroup='i3',
                                                        doSegment=False, doNormalize=False,
                                                        doInvert=False, doThreshold=False,
                                                        doVolumeMask=True
                                                        )
        protPreprocessVol3.inputVolumes.set(self.protImportVirus.outputVolume)
        protPreprocessVol3.volumeMask.set(self.protImportvirusMaskCapsid.outputMask)
        self.proj.launchProtocol(protPreprocessVol3, wait=True)
        self.assertIsNotNone(protPreprocessVol3.outputVol, "There was a problem with a volume")

        print "Run preprocess a volume using mask_2 in the symmetrization"
        protPreprocessVol4 = XmippProtPreprocessVolumes(doChangeHand=False, doRandomize=False,
                                                        doSymmetrize=True, symmetryGroup='c7',
                                                        doSegment=False, doNormalize=False,
                                                        doInvert=False, doThreshold=False,
                                                        doVolumeMask=True
                                                        )
        protPreprocessVol4.inputVolumes.set(self.protImportVirus.outputVolume)
        protPreprocessVol4.volumeMask.set(self.protImportvirusMaskPenton.outputMask)
        self.proj.launchProtocol(protPreprocessVol4, wait=True)
        self.assertIsNotNone(protPreprocessVol4.outputVol, "There was a problem with a volume")


class TestXmippResolution3D(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData()
        cls.protImport1 = cls.runImportVolumes(cls.vol2, 9.896)
        cls.protImport2 = cls.runImportVolumes(cls.vol3, 9.896)

    def testCalculateResolution(self):
        print "Run resolution 3D"
        protResol3D = XmippProtResolution3D(doSSNR=False)
        protResol3D.inputVolume.set(self.protImport1.outputVolume)
        protResol3D.referenceVolume.set(self.protImport2.outputVolume)
        self.proj.launchProtocol(protResol3D, wait=True)
        self.assertIsNotNone(protResol3D._defineFscName(), "There was a problem with fsc")
        self.assertIsNotNone(protResol3D._defineStructFactorName(), "There was a problem with structure factor")


class TestXmippFilterVolumes(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        print "\n", greenStr(" Filter Volumes Set Up - Collect data ".center(75, '-'))
        setupTestProject(cls)
        TestXmippBase.setData()
        cls.protImport1 = cls.runImportVolumes(cls.volumes, 9.896)
        cls.protImport2 = cls.runImportVolumes(cls.vol1, 9.896)

    # Tests with single volume as input.
    def launchAndTestSingle(self, **kwargs):
        "Launch XmippProtFilterVolumes on single volume and check results."
        print magentaStr("\n==> Filter singe volume input params: %s" % kwargs)
        prot = XmippProtFilterVolumes(**kwargs)
        prot.inputVolumes.set(self.protImport2.outputVolume)
        self.proj.launchProtocol(prot, wait=True)
        self.assertTrue(hasattr(prot, "outputVol") and prot.outputVol is not None,
                        "There was a problem with filter single volume")
        self.assertTrue(prot.outputVol.equalAttributes(
            self.protImport2.outputVolume, ignore=['_index', '_filename'],
            verbose=True))

    def testSingleFourier(self):
        self.launchAndTestSingle(filterSpace=FILTER_SPACE_FOURIER,
                                 lowFreq=0.1, highFreq=0.25)

    def testSingleMedian(self):
        self.launchAndTestSingle(filterSpace=FILTER_SPACE_REAL,
                                 filterModeReal=xfh.FM_MEDIAN)

    def testSingleWavelets(self):
        self.launchAndTestSingle(filterSpace=FILTER_SPACE_WAVELET,
                                 filterModeWavelets=xfh.FM_DAUB12,
                                 waveletMode=xfh.FM_REMOVE_SCALE)

    # Tests with multiple volumes as input.
    def launchAndTestSet(self, **kwargs):
        "Launch XmippProtFilterVolumes on set of volumes and check results."
        print magentaStr("\n==> Filter multiple volumes input params: %s" % kwargs)
        prot = XmippProtFilterVolumes(**kwargs)
        vIn = self.protImport1.outputVolumes  # short notation
        prot.inputVolumes.set(vIn)
        self.proj.launchProtocol(prot, wait=True)
        self.assertTrue(hasattr(prot, "outputVol") and prot.outputVol is not None,
                        "There was a problem with filter multiple volumes")
        self.assertTrue(prot.outputVol.equalAttributes(
            self.protImport1.outputVolumes, ignore=['_mapperPath'],
            verbose=True))
        # Compare the individual volumes too.
        self.assertTrue(prot.outputVol.equalItemAttributes(
            self.protImport1.outputVolumes, ignore=['_index', '_filename'],
            verbose=True))

    def testSetFourier(self):
        self.launchAndTestSet(filterSpace=FILTER_SPACE_FOURIER,
                              lowFreq=0.1, highFreq=0.25)

    def testSetMedian(self):
        self.launchAndTestSet(filterSpace=FILTER_SPACE_REAL,
                              filterModeReal=xfh.FM_MEDIAN)

    def testSetWavelets(self):
        self.launchAndTestSet(filterSpace=FILTER_SPACE_WAVELET,
                              filterModeWavelets=xfh.FM_DAUB12,
                              waveletMode=xfh.FM_REMOVE_SCALE)


class TestXmippMaskVolumes(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData()
        cls.protImport1 = cls.runImportVolumes(cls.volumes, 9.896)
        cls.protImport2 = cls.runImportVolumes(cls.vol1, 9.896)

    def testMaskVolumes(self):
        print "Run mask single volume"
        protMaskVolume = self.newProtocol(XmippProtMaskVolumes,
                                          radius=23)
        protMaskVolume.inputVolumes.set(self.protImport2.outputVolume)
        self.launchProtocol(protMaskVolume)
        self.assertIsNotNone(protMaskVolume.outputVol, "There was a problem with applying mask to a volume")

        print "Run mask SetOfVolumes"
        protMaskVolumes = self.newProtocol(XmippProtMaskVolumes,
                                           geo=MASK3D_CROWN, innerRadius=18, outerRadius=23)
        protMaskVolumes.inputVolumes.set(self.protImport1.outputVolumes)
        self.proj.launchProtocol(protMaskVolumes, wait=True)
        self.assertIsNotNone(protMaskVolumes.outputVol, "There was a problem with applying mask to SetOfVolumes")


class TestXmippCropResizeVolumes(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        print "\n", greenStr(" Crop/Resize Volumes Set Up - Collect data ".center(75, '-'))
        setupTestProject(cls)
        TestXmippBase.setData()
        cls.protImport1 = cls.runImportVolumes(cls.volumes, 9.896)
        cls.protImport2 = cls.runImportVolumes(cls.vol1, 9.896)

    # Tests with single volume as input.
    def launchSingle(self, **kwargs):
        "Launch XmippProtCropResizeVolumes and return output volume."
        print magentaStr("\n==> Crop/Resize single volume input params: %s" % kwargs)
        prot = XmippProtCropResizeVolumes(**kwargs)
        prot.inputVolumes.set(self.protImport2.outputVolume)
        self.proj.launchProtocol(prot, wait=True)
        self.assertTrue(hasattr(prot, "outputVol") and prot.outputVol is not None,
                        "There was a problem with applying resize/crop to a volume")
        return prot.outputVol

    def testSingleResizeDimensions(self):
        inV = self.protImport2.outputVolume  # short notation
        newSize = 128
        outV = self.launchSingle(doResize=True,
                                 resizeOption=xrh.RESIZE_DIMENSIONS,
                                 resizeDim=newSize, doWindow=True,
                                 windowOperation=xrh.WINDOW_OP_WINDOW,
                                 windowSize=newSize*2)

        self.assertEqual(newSize * 2, outV.getDim()[0])
        self.assertAlmostEqual(outV.getSamplingRate(),
                               inV.getSamplingRate() * (inV.getDim()[0] / float(newSize)))
        self.assertTrue(outV.equalAttributes(
            inV, ignore=['_index', '_filename', '_samplingRate'], verbose=True))

    def testSingleFactorAndCrop(self):
        inV = self.protImport2.outputVolume  # short notation
        outV = self.launchSingle(doResize=True,
                                 resizeOption=xrh.RESIZE_FACTOR,
                                 resizeFactor=0.5,
                                 doWindow=True,
                                 windowOperation=xrh.WINDOW_OP_CROP)

        self.assertEqual(inV.getDim()[0] * 0.5, outV.getDim()[0])
        self.assertAlmostEqual(outV.getSamplingRate(), inV.getSamplingRate() * 2)
        self.assertTrue(outV.equalAttributes(
            inV, ignore=['_index', '_filename', '_samplingRate'], verbose=True))

    def testSinglePyramid(self):
        inV = self.protImport2.outputVolume  # short notation
        outV = self.launchSingle(doResize=True, resizeOption=xrh.RESIZE_PYRAMID,
                                 resizeLevel=1)

        # Since the images were expanded by 2**resizeLevel (=2) the new
        # pixel size (painfully called "sampling rate") should be 0.5x.
        self.assertEqual(inV.getDim()[0] * 2, outV.getDim()[0])
        self.assertAlmostEqual(outV.getSamplingRate(), inV.getSamplingRate() * 0.5)
        self.assertTrue(outV.equalAttributes(
            inV, ignore=['_index', '_filename', '_samplingRate'], verbose=True))

    # Tests with multiple volumes as input.
    def launchSet(self, **kwargs):
        "Launch XmippProtCropResizeVolumes and return output volumes."
        print magentaStr("\n==> Crop/Resize single set of volumes input params: %s" % kwargs)
        prot = XmippProtCropResizeVolumes(**kwargs)
        prot.inputVolumes.set(self.protImport1.outputVolumes)
        self.proj.launchProtocol(prot, wait=True)
        self.assertTrue(hasattr(prot, "outputVol") and prot.outputVol is not None,
                        "There was a problem with applying resize/crop to a set of volumes")
        return prot.outputVol

    def testSetResizeDimensions(self):
        inV = self.protImport1.outputVolumes  # short notation
        newSize = 128
        outV = self.launchSet(doResize=True,
                              resizeOption=xrh.RESIZE_DIMENSIONS,
                              resizeDim=newSize, doWindow=True,
                              windowOperation=xrh.WINDOW_OP_WINDOW,
                              windowSize=newSize*2)

        self.assertEqual(newSize * 2, outV.getDim()[0])
        self.assertAlmostEqual(outV.getSamplingRate(),
                               inV.getSamplingRate() * (inV.getDim()[0] / float(newSize)))
        self.assertTrue(outV.equalAttributes(
            inV, ignore=['_mapperPath', '_samplingRate', '_firstDim'], verbose=True))
        # Compare the individual volumes too.
        self.assertTrue(outV.equalItemAttributes(
            inV, ignore=['_index', '_filename', '_samplingRate'], verbose=True))

    def testSetFactorAndCrop(self):
        inV = self.protImport1.outputVolumes  # short notation
        outV = self.launchSet(doResize=True,
                              resizeOption=xrh.RESIZE_FACTOR,
                              resizeFactor=0.5,
                              doWindow=True,
                              windowOperation=xrh.WINDOW_OP_CROP)

        self.assertEqual(inV.getDim()[0] * 0.5, outV.getDim()[0])
        self.assertAlmostEqual(outV.getSamplingRate(), inV.getSamplingRate() * 2)
        self.assertTrue(outV.equalAttributes(
            inV, ignore=['_mapperPath', '_samplingRate', '_firstDim'], verbose=True))
        # Compare the individual volumes too.
        self.assertTrue(outV.equalItemAttributes(
            inV, ignore=['_index', '_filename', '_samplingRate'], verbose=True))

    def testSetPyramid(self):
        inV = self.protImport1.outputVolumes  # short notation
        outV = self.launchSet(doResize=True, resizeOption=xrh.RESIZE_PYRAMID,
                              resizeLevel=1)

        # Since the images were expanded by 2**resizeLevel (=2) the new
        # pixel size (painfully called "sampling rate") should be 0.5x.
        self.assertEqual(inV.getDim()[0] * 2, outV.getDim()[0])
        self.assertAlmostEqual(outV.getSamplingRate(), inV.getSamplingRate() * 0.5)
        self.assertTrue(outV.equalAttributes(
            inV, ignore=['_mapperPath', '_samplingRate', '_firstDim'], verbose=True))
        # Compare the individual volumes too.
        self.assertTrue(outV.equalItemAttributes(
            inV, ignore=['_index', '_filename', '_samplingRate'], verbose=True))


class TestXmippOperateVolumes(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData()
        cls.protImport1 = cls.runImportVolumes(cls.vol1, 9.896)
        cls.protImport2 = cls.runImportVolumes(cls.vol2, 9.896)
        cls.protImport3 = cls.runImportVolumes(cls.volumes, 9.896)

    # Tests with single volume as input.
    def launchSingle(self, **kwargs):
        "Launch XmippProtImageOperateVolumes and return output volume."
        print magentaStr("\n==> Operate single volume input params: %s" % kwargs)
        prot = XmippProtImageOperateVolumes()
        prot.operation.set(kwargs.get('operation', 1))
        prot.inputVolumes.set(self.protImport1.outputVolume)
        prot.setObjLabel(kwargs.get('objLabel', None))
        prot.isValue.set(kwargs.get('isValue', False))
        prot.inputVolumes2.set(kwargs.get('volumes2', None))
        prot.value.set(kwargs.get('value', None))
        prot.intValue.set(kwargs.get('intValue', None))
        
        self.proj.launchProtocol(prot, wait=True)
        self.assertTrue(hasattr(prot, "outputVol") and prot.outputVol is not None,
                        "There was a problem producing the output")
        return prot.outputVol

    def testMultiplyVolumes(self):
        vol2 = self.protImport2.outputVolume  # short notation
        prot1 = self.launchSingle(operation=OP_MULTIPLY,
                                  objLabel='Multiply two Volumes',
                                  volumes2=vol2)

    def testMultiplyValue(self):
        prot2 = self.launchSingle(operation=OP_MULTIPLY,
                                  isValue=True,
                                  objLabel='Multiply by a Value',
                                  value=2.5)
    
    def testDotProduct(self):
        vol2 = self.protImport2.outputVolume  # short notation
        prot3 = self.launchSingle(operation=OP_DOTPRODUCT,
                                  objLabel='Dot Product',
                                  volumes2=vol2)

    def testSqrt(self):
        prot4 = self.launchSingle(operation=OP_SQRT,
                                  objLabel='Sqrt')

    def testRadial(self):
        prot5 = self.launchSingle(operation=OP_RADIAL,
                                  objLabel='Radial Average')

    def testColumn(self):
        prot6 = self.launchSingle(operation=OP_COLUNM,
                                  objLabel='Column',
                                  intValue  =7)

    def testRow(self):
        prot6 = self.launchSingle(operation=OP_ROW,
                                  objLabel='Row',
                                  intValue  =8)

#     # Tests with multiple volumes as input.
    def launchSet(self, **kwargs):
        "Launch XmippProtImageOperateVolumes and return output volumes."
        print magentaStr("\n==> Operate set of volumes input params: %s" % kwargs)
        prot = XmippProtImageOperateVolumes()
        prot.operation.set(kwargs.get('operation', 1))
        prot.inputVolumes.set(self.protImport3.outputVolumes)
        prot.setObjLabel(kwargs.get('objLabel', None))
        prot.isValue.set(kwargs.get('isValue', False))
        prot.inputVolumes2.set(kwargs.get('volumes2', None))
        prot.value.set(kwargs.get('value', None))
        prot.intValue.set(kwargs.get('intValue', None))
        
        self.proj.launchProtocol(prot, wait=True)
        self.assertTrue(hasattr(prot, "outputVol") and prot.outputVol is not None,
                        "There was a problem producing the output")
        return prot.outputVol
        
    def testMultiplyVolSets(self):
        vol2 = self.protImport3.outputVolumes  # short notation
        prot6 = self.launchSet(operation=OP_MULTIPLY,
                                  objLabel='Multiply two SetOfVolumes',
                                  volumes2=vol2)

    def testMultiplyValue2(self):
        prot7 = self.launchSet(operation=OP_MULTIPLY,
                               isValue=True,
                               objLabel='Multiply by a Value 2',
                               value=2.5)
    
    def testDotProduct2(self):
        vol2 = self.protImport3.outputVolumes  # short notation
        prot8 = self.launchSet(operation=OP_DOTPRODUCT,
                               objLabel='Dot Product 2',
                               volumes2=vol2)

    def testSqrt2(self):
        prot9 = self.launchSet(operation=OP_SQRT,
                               objLabel='Sqrt 2')


class TestXmippProtAlignVolume(TestXmippBase):
    @classmethod
    def setData(cls, dataProject='xmipp_tutorial'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.volumes = cls.dataset.getFile('volumes')
        cls.vol1 = cls.dataset.getFile('vol1')
        cls.vol2 = cls.dataset.getFile('vol2')
        cls.vol3 = cls.dataset.getFile('vol3')

    @classmethod
    def runImportVolumes(cls, pattern, samplingRate):
        """ Run an Import particles protocol. """
        return cls.protImport
    
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('relion_tutorial')

        cls.protImport1 = cls.newProtocol(ProtImportVolumes,
                                         filesPath=cls.ds.getFile('volumes/reference_rotated.vol'), 
                                         samplingRate=1.0)
        cls.launchProtocol(cls.protImport1)
        
        cls.protImport2 = cls.newProtocol(ProtImportVolumes,
                                         filesPath=cls.ds.getFile('volumes/reference.mrc'), 
                                         samplingRate=1.0)
        cls.launchProtocol(cls.protImport2)
        
        # Rotate that volume rot=90 tilt=90 to create 
        # a gold rotated volume
        os.system('')

    def testExhaustive(self):
        protAlign = self.newProtocol(XmippProtAlignVolume,
                                     inputReference=self.protImport1.outputVolume,
                                     alignmentAlgorithm=ALIGN_ALGORITHM_EXHAUSTIVE,
                                     minRotationalAngle=65, 
                                     maxRotationalAngle=100,
                                     stepRotationalAngle=10,
                                     minTiltAngle=65,
                                     maxTiltAngle=100,
                                     stepTiltAngle=10,
                                     minInplaneAngle=0,
                                     maxInplaneAngle=0,
                                     stepInplaneAngle=0,
                                     numberOfMpi=1, numberOfThreads=1                               
                                     )
        protAlign.inputVolumes.append(self.protImport2.outputVolume)
        self.launchProtocol(protAlign)
        
    def testLocal(self):
        protAlign = self.newProtocol(XmippProtAlignVolume,
                                     inputReference=self.protImport1.outputVolume,
                                     alignmentAlgorithm=ALIGN_ALGORITHM_LOCAL,
                                     initialRotAngle=0,
                                     initialTiltAngle=0,
                                     initialInplaneAngle=0,
                                     initialShiftX=0,
                                     initialShiftY=0,
                                     initialShiftZ=0,
                                     initialScale=1,
                                     optimizeScale=True,
                                     numberOfMpi=1, numberOfThreads=1                               
                                     )
        protAlign.inputVolumes.append(self.protImport2.outputVolume)
        self.launchProtocol(protAlign)
        
    def testExhaustiveLocal(self):
        protAlign = self.newProtocol(XmippProtAlignVolume,
                                     inputReference=self.protImport1.outputVolume,
                                     alignmentAlgorithm=ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL,
                                     minRotationalAngle=65,
                                     maxRotationalAngle=100,
                                     stepRotationalAngle=10,
                                     minTiltAngle=65,
                                     maxTiltAngle=100,
                                     stepTiltAngle=10,
                                     minInplaneAngle=0,
                                     maxInplaneAngle=0,
                                     stepInplaneAngle=0,
                                     numberOfMpi=1, numberOfThreads=1                               
                                     )
        protAlign.inputVolumes.append(self.protImport2.outputVolume)
        self.launchProtocol(protAlign)
        

class TestXmippConvertToPseudoatoms(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData()
        cls.protImport = cls.runImportVolumes(cls.vol1, 9.896)

    def testConvertToPseudoatoms(self):
        print "Run convert to pseudoatoms"
        prot = XmippProtConvertToPseudoAtoms(pseudoAtomTarget=15)
        prot.inputStructure.set(self.protImport.outputVolume)
        self.proj.launchProtocol(prot, wait=True)
        self.assertIsNotNone(prot.outputPdb, "There was a problem with Convert to pseudoatoms output Pdb")


class TestXmippProtHelicalParameters(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('general')
        cls.vol = cls.ds.getFile('vol_helix')
        cls.protImport = cls.runImportVolumes(cls.vol, 1.0)

    def testHelicalParameters(self):
        print "Run symmetrize helical"
        protHelical = XmippProtHelicalParameters(cylinderOuterRadius=20,dihedral=True,rot0=50,rotF=70,rotStep=5,z0=5,zF=10,zStep=0.5)
        protHelical.inputVolume.set(self.protImport.outputVolume)
        self.proj.launchProtocol(protHelical, wait=True)

        self.assertIsNotNone(protHelical.outputVolume, "There was a problem with Helical output volume")
        self.assertIsNotNone(protHelical.deltaRot.get(), "Output delta rot is None")
        self.assertIsNotNone(protHelical.deltaZ.get(), "Output delta Z is None")
        print "protHelical.deltaRot.get() ", protHelical.deltaRot.get()
        self.assertAlmostEqual(protHelical.deltaRot.get(), 59.59, delta=1, msg="Output delta rot is wrong")
        self.assertAlmostEqual(protHelical.deltaZ.get(), 6.628, delta=0.2, msg="Output delta Z is wrong")


class TestXmippRansacMda(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('mda')
        cls.averages = cls.dataset.getFile('averages')
        cls.samplingRate = 3.5
        cls.symmetryGroup = 'd6'
        cls.angularSampling = 15
        cls.nRansac = 25
        cls.numSamples = 5
        cls.dimRed = False
        cls.numVolumes = 2
        cls.maxFreq = 30

    def test_ransac(self):
        #Import a set of averages
        print "Import Set of averages"
        protImportAvg = self.newProtocol(ProtImportAverages, 
                                         filesPath=self.averages, 
                                         checkStack=True,
                                         samplingRate=self.samplingRate)
        self.launchProtocol(protImportAvg)
        self.assertIsNotNone(protImportAvg.getFiles(), "There was a problem with the import")
        
        print "Run Ransac"
        protRansac = self.newProtocol(XmippProtRansac,
                                      symmetryGroup=self.symmetryGroup, 
                                      angularSampling=self.angularSampling,
                                      nRansac=self.nRansac, 
                                      numSamples=self.numSamples, 
                                      dimRed=self.dimRed,
                                      numVolumes=self.numVolumes, 
                                      maxFreq=self.maxFreq, useAll=True, 
                                      numberOfThreads=4)
        protRansac.inputSet.set(protImportAvg.outputAverages)
        self.launchProtocol(protRansac)
        self.assertIsNotNone(protRansac.outputVolumes, "There was a problem with ransac protocol")


class TestXmippRansacGroel(TestXmippRansacMda):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('groel')
        cls.averages = cls.dataset.getFile('averages')
        cls.samplingRate = 2.1
        cls.symmetryGroup = 'd7'
        cls.angularSampling = 7
        cls.nRansac = 25
        cls.numSamples = 5
        cls.dimRed = True
        cls.numVolumes = 2
        cls.maxFreq = 12


class TestXmippRotationalSymmetry(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('xmipp_tutorial')
        cls.vol = cls.dataset.getFile('vol110')

    def test_rotsym(self):
        print "Import Volume"
        protImportVol = self.newProtocol(ProtImportVolumes,
                                         objLabel='Volume',
                                         filesPath=self.vol,
                                         samplingRate=7.08)
        self.launchProtocol(protImportVol)
        self.assertIsNotNone(protImportVol.getFiles(),
                             "There was a problem with the import")
        
        print "Run find rotational symmetry axis"
        protRotSym = self.newProtocol(XmippProtRotationalSymmetry,
                                         symOrder=2,
                                         searchMode=2,
                                         tilt0=70,
                                         tiltF=110)
        protRotSym.inputVolume.set(protImportVol.outputVolume)
        self.launchProtocol(protRotSym)
        self.assertIsNotNone(protRotSym.outputVolume,
                             "There was a problem with Rotational Symmetry")


class TestXmippProjMatching(TestXmippBase):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('relion_tutorial')
        cls.vol = cls.dataset.getFile('volume')

    def testXmippProjMatching(self):
        print "Import Particles"
        protImportParts = self.newProtocol(ProtImportParticles,
                                 objLabel='Particles from scipion',
                                 importFrom=ProtImportParticles.IMPORT_FROM_SCIPION,
                                 sqliteFile=self.dataset.getFile('import/case2/particles.sqlite'),
                                 magnification=50000,
                                 samplingRate=7.08,
                                 haveDataBeenPhaseFlipped=True
                                 )
        self.launchProtocol(protImportParts)
        self.assertIsNotNone(protImportParts.getFiles(), "There was a problem with the import")
        
        print "Get a Subset of particles"
        protSubset = self.newProtocol(ProtSubSet,
                                         objLabel='100 Particles',
                                         chooseAtRandom=True,
                                         nElements=100)
        protSubset.inputFullSet.set(protImportParts.outputParticles)
        self.launchProtocol(protSubset)
        
        print "Import Volume"
        protImportVol = self.newProtocol(ProtImportVolumes,
                                         objLabel='Volume',
                                         filesPath=self.vol,
                                         samplingRate=7.08)
        self.launchProtocol(protImportVol)
        self.assertIsNotNone(protImportVol.getFiles(), "There was a problem with the import")
        
        print "Run Projection Matching"
        protProjMatch = self.newProtocol(XmippProtProjMatch,
                                         ctfGroupMaxDiff=0.00001,
                                         mpiJobSize=10,
                                         numberOfIterations=2,
                                         numberOfThreads=2,
                                         numberOfMpi=3)
        protProjMatch.inputParticles.set(protSubset.outputParticles)
        protProjMatch.input3DReferences.set(protImportVol.outputVolume)
        self.launchProtocol(protProjMatch)
        self.assertIsNotNone(protProjMatch.outputVolume, "There was a problem with Projection Matching")


class TestPdbImport(TestXmippBase):
    
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('nma')
        cls.pdb = cls.dataset.getFile('pdb')
    
    def testImportPdbFromId(self):
        print "Run convert a pdb from database"
        protConvert = self.newProtocol(ProtImportPdb, pdbId="3j3i")
        self.launchProtocol(protConvert)
        self.assertIsNotNone(protConvert.outputPdb.getFileName(), 
                             "There was a problem with the import")
        
    def testImportPdbFromFn(self):
        print "Run convert a pdb from file"
        protConvert = self.newProtocol(ProtImportPdb, 
                                       inputPdbData=ProtImportPdb.IMPORT_FROM_FILES, 
                                       pdbFile=self.pdb)
        self.launchProtocol(protConvert)
        self.assertIsNotNone(protConvert.outputPdb.getFileName(), 
                             "There was a problem with the import")
        
class TestXmippPdbConvert(TestXmippBase):
    
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('nma')
        cls.pdb = cls.dataset.getFile('pdb')
    
    def testXmippPdbConvertFromDb(self):
        print "Run convert a pdb from database"
        protConvert = self.newProtocol(XmippProtConvertPdb, pdbId="3j3i", sampling=4, setSize=True, size=100)
        self.launchProtocol(protConvert)
        self.assertIsNotNone(protConvert.outputVolume.getFileName(), "There was a problem with the conversion")
        self.assertAlmostEqual(protConvert.outputVolume.getSamplingRate(), protConvert.sampling.get(), places=1, msg="wrong sampling rate")
        self.assertAlmostEqual(protConvert.outputVolume.getDim()[0], protConvert.size.get(), places=1, msg="wrong size")
        
    def testXmippPdbConvertFromObj(self):
        print "Run convert a pdb from import"
        protImport = self.newProtocol(ProtImportPdb, 
                                      inputPdbData=ProtImportPdb.IMPORT_FROM_FILES, 
                                      pdbFile=self.pdb)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputPdb.getFileName(), "There was a problem with the import")
        
        protConvert = self.newProtocol(XmippProtConvertPdb, 
                                       inputPdbData=XmippProtConvertPdb.IMPORT_OBJ, 
                                       sampling=3, setSize=True, size=20)
        protConvert.pdbObj.set(protImport.outputPdb)
        self.launchProtocol(protConvert)
        self.assertIsNotNone(protConvert.outputVolume.getFileName(), "There was a problem with the conversion")
        self.assertAlmostEqual(protConvert.outputVolume.getSamplingRate(), protConvert.sampling.get(), places=1, msg="wrong sampling rate")
        self.assertAlmostEqual(protConvert.outputVolume.getDim()[0], protConvert.size.get(), places=1, msg="wrong size")

    def testXmippPdbConvertFromFn(self):
        print "Run convert a pdb from file"
        protConvert = self.newProtocol(XmippProtConvertPdb,inputPdbData=2, pdbFile=self.pdb, sampling=2, setSize=True)
        self.launchProtocol(protConvert)
        self.assertIsNotNone(protConvert.outputVolume.getFileName(), "There was a problem with the conversion")
        self.assertAlmostEqual(protConvert.outputVolume.getSamplingRate(), protConvert.sampling.get(), places=1, msg="wrong sampling rate")
        self.assertAlmostEqual(protConvert.outputVolume.getDim()[0], 48, places=1, msg="wrong size")


class TEstXmippValidateNonTilt(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('relion_tutorial')
        cls.vol = cls.dataset.getFile('volume')

    def testXmippValidateNonTilt(self):
        print "Import Particles"
        protImportParts = self.newProtocol(ProtImportParticles,
                                 objLabel='Particles from scipion',
                                 importFrom=ProtImportParticles.IMPORT_FROM_SCIPION,
                                 sqliteFile=self.dataset.getFile('import/case2/particles.sqlite'),
                                 magnification=50000,
                                 samplingRate=7.08,
                                 haveDataBeenPhaseFlipped=True
                                 )
        self.launchProtocol(protImportParts)
        self.assertIsNotNone(protImportParts.getFiles(), "There was a problem with the import")
        
        print "Get a Subset of particles"
        protSubset = self.newProtocol(ProtSubSet,
                                         objLabel='100 Particles',
                                         chooseAtRandom=True,
                                         nElements=100)
        protSubset.inputFullSet.set(protImportParts.outputParticles)
        self.launchProtocol(protSubset)
        
        print "Import Volume"
        protImportVol = self.newProtocol(ProtImportVolumes,
                                         objLabel='Volume',
                                         filesPath=self.vol,
                                         samplingRate=7.08)
        self.launchProtocol(protImportVol)
        self.assertIsNotNone(protImportVol.getFiles(), "There was a problem with the import")
        
        print "Run Validate Non-Tilt significant"
        protValidate = self.newProtocol(XmippProtValidateNonTilt)
        protValidate.inputParticles.set(protSubset.outputParticles)
        protValidate.inputVolumes.set(protImportVol.outputVolume)
        self.launchProtocol(protValidate)
        self.assertIsNotNone(protValidate.outputVolumes, "There was a problem with Validate Non-Tilt")
        
        print "Run Validate Non-Tilt projection matching"
        protValidate = self.newProtocol(XmippProtValidateNonTilt, alignmentMethod=1)
        protValidate.inputParticles.set(protSubset.outputParticles)
        protValidate.inputVolumes.set(protImportVol.outputVolume)
        self.launchProtocol(protValidate)
        self.assertIsNotNone(protValidate.outputVolumes, "There was a problem with Validate Non-Tilt")


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
