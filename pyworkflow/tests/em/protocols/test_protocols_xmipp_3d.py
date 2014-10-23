# **************************************************************************
# *
# * Authors:    Josue Gomez Blanco (jgomez@cnb.csic.es)
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



class TestXmippBase(BaseTest):
    """ Some utility functions to import volumes that are used in several tests."""

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
        cls.protImport = cls.newProtocol(ProtImportVolumes,
                                         pattern=pattern, samplingRate=samplingRate)
        cls.launchProtocol(cls.protImport)
        return cls.protImport

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
    def runClassify(cls, particles):
        cls.ProtClassify = cls.newProtocol(XmippProtML2D,
                                           numberOfReferences=8, maxIters=4, doMlf=False,
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

    def testCreateMask1(self):
        print "Run create mask from volume"
        protMask1 = self.newProtocol(XmippProtCreateMask3D,
                                     source=0, threshold=0.4)
        protMask1.inputVolume.set(self.protImport.outputVolume)
        protMask1.setObjLabel('thresold mask')
        self.launchProtocol(protMask1)
        self.assertIsNotNone(protMask1.outputMask, "There was a problem with create mask from volume")

        print "Run create mask from another mask"
        protMask2 = self.newProtocol(XmippProtCreateMask3D,
                                     source=2, doMorphological=True, elementSize=3)
        protMask2.inputMask.set(protMask1.outputMask)
        protMask2.setObjLabel('dilation mask')
        self.launchProtocol(protMask2)
        self.assertIsNotNone(protMask2.outputMask, "There was a problem with mask from another mask")

        print "Apply dilation mask to imported volume"
        protMaskVolume = self.newProtocol(XmippProtMaskVolumes)
        protMaskVolume.inputVolumes.set(self.protImport.outputVolume)
        protMaskVolume.source.set(SOURCE_MASK)
        protMaskVolume.inputMask.set(protMask2.outputMask)
        self.launchProtocol(protMaskVolume)
        self.assertIsNotNone(protMaskVolume.outputVol, "There was a problem with applying mask to a volume")

        print "Run create mask from geometry"
        protMask3 = self.newProtocol(XmippProtCreateMask3D,
                                     source=1, size=64, samplingRate=9.89,
                                     geo=6, innerRadius=10, outerRadius=25, borderDecay=2)
        protMask3.setObjLabel('crown mask')
        self.launchProtocol(protMask3)
        self.assertIsNotNone(protMask3.outputMask, "There was a problem with create mask from geometry")


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


class TestXmippPreprocessVolumes(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData()
        cls.protImport1 = cls.runImportVolumes(cls.volumes, 9.896)
        cls.protImport2 = cls.runImportVolumes(cls.vol1, 9.896)

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


class TestXmippFilterVolumes(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData()
        cls.protImport1 = cls.runImportVolumes(cls.volumes, 9.896)
        cls.protImport2 = cls.runImportVolumes(cls.vol1, 9.896)

    def testFilterVolumes(self):
        print "Run filter single volume"
        protFilterVolume = XmippProtFilterVolumes(lowFreq=0.1, highFreq=0.25)
        protFilterVolume.inputVolumes.set(self.protImport2.outputVolume)
        self.proj.launchProtocol(protFilterVolume, wait=True)
        self.assertIsNotNone(protFilterVolume.outputVol, "There was a problem with filter a volume")

        print "Run filter SetOfVolumes"
        protFilterVolumes = XmippProtFilterVolumes(lowFreq=0.1, highFreq=0.25)
        protFilterVolumes.inputVolumes.set(self.protImport1.outputVolumes)
        self.proj.launchProtocol(protFilterVolumes, wait=True)
        self.assertIsNotNone(protFilterVolumes.outputVol, "There was a problem with filter SetOfVolumes")


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
        setupTestProject(cls)
        TestXmippBase.setData()
        cls.protImport1 = cls.runImportVolumes(cls.volumes, 9.896)
        cls.protImport2 = cls.runImportVolumes(cls.vol1, 9.896)

    def testCropResizeVolumes(self):
        print "Run Resize-Crop single volume"
        protCropResizeVolume = XmippProtCropResizeVolumes(doResize=True, resizeOption=1, resizeDim=128, doWindow=True,
                                                    windowOperation=1, windowSize=256)
        protCropResizeVolume.inputVolumes.set(self.protImport2.outputVolume)
        self.proj.launchProtocol(protCropResizeVolume, wait=True)
        self.assertIsNotNone(protCropResizeVolume.outputVol, "There was a problem with applying resize and crop to a volume")

        print "Run Resize-Crop SetOfVolumes"
        protCropResizeVolumes = XmippProtCropResizeVolumes(doResize=True, resizeOption=1, resizeDim=128, doWindow=True,
                                                    windowOperation=1, windowSize=256)
        protCropResizeVolumes.inputVolumes.set(self.protImport1.outputVolumes)
        self.proj.launchProtocol(protCropResizeVolumes, wait=True)
        self.assertIsNotNone(protCropResizeVolumes.outputVol, "There was a problem with applying resize and crop to SetOfVolumes")


class TestXmippCLTomo(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData('tomo')
        cls.protImport = cls.runImportVolumes(cls.volumes, 9.896)

    def testCLTomo(self):
        print "Run CLTomo"
        protCLTomo = XmippProtCLTomo(numberOfReferences=1,numberOfIterations=1)
        protCLTomo.volumelist.set(self.protImport.outputVolumes)
        self.proj.launchProtocol(protCLTomo, wait=True)

        self.assertIsNotNone(protCLTomo.outputClasses, "There was a problem with CLTomo output classes")
        self.assertIsNotNone(protCLTomo.alignedVolumes, "There was a problem with CLTomo output aligned volumes")


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
        TestXmippBase.setData()
        cls.protImport = cls.runImportVolumes(cls.vol1, 9.896)

    def testHelicalParameters(self):
        print "Run symmetrize helical"
        protHelical = XmippProtHelicalParameters(cylinderRadius=20,dihedral=False,rot0=50,rotF=70,rotStep=5,z0=5,zF=10,zStep=0.5)
        protHelical.inputVolume.set(self.protImport.outputVolume)
        self.proj.launchProtocol(protHelical, wait=True)

        self.assertIsNotNone(protHelical.outputVolume, "There was a problem with Helical output volume")


class TestXmippSimAnnealing(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
        cls.Class2D = cls.runClassify(cls.protImport.outputParticles)

    def test_simAnnealing(self):
        print "Run Simulating annealing"
        protSimAnneal = self.newProtocol(XmippProtInitVolSimAnneal,
                                         symmetryGroup='d6', numberOfSimAnnealRef=2, percentRejection=0)
        protSimAnneal.inputClasses.set(self.Class2D.outputClasses)
        self.launchProtocol(protSimAnneal)
        self.assertIsNotNone(protSimAnneal.outputVolumes, "There was a problem with simulating annealing protocol")


class TestXmippRansac(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('mda')
        cls.particlesFn = cls.dataset.getFile('particles')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
        cls.Class2D = cls.runClassify(cls.protImport.outputParticles)

    def test_ransac(self):
        print "Run Ransac"
        protRansac = self.newProtocol(XmippProtRansac,
                                      symmetryGroup='d6', angularSampling=15, nRansac=25, numSamples=5,
                                      dimRed=False, numVolumes=2, maxFreq=30, useAll=True, numberOfThreads=4)
        protRansac.inputClasses.set(self.Class2D.outputClasses)
        self.launchProtocol(protRansac)
        self.assertIsNotNone(protRansac.outputVolumes, "There was a problem with simulating annealing protocol")


class TestXmippProjMatching(TestXmippBase):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('xmipp_tutorial')
        cls.allCrdsDir = cls.dataset.getFile('posAllDir')
        cls.micsFn = cls.dataset.getFile('allMics')
        cls.vol1 = cls.dataset.getFile('vol1')

    def testXmippProjMatching(self):
        #First, import a set of micrographs
        protImport = self.newProtocol(ProtImportMicrographs, pattern=self.micsFn, samplingRate=1.237, voltage=300)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputMicrographs.getFileName(), "There was a problem with the import")
#         self.validateFiles('protImport', protImport)

        #Import a set of volumes
        print "Import Volume"
        protImportVol = self.newProtocol(ProtImportVolumes, pattern=self.vol1, samplingRate=9.896)
        self.launchProtocol(protImportVol)
        self.assertIsNotNone(protImportVol.getFiles(), "There was a problem with the import")
#        self.validateFiles('protImportVol', protImportVol)

        # Perform a downsampling on the micrographs
        print "Downsampling..."
        protDownsampling = self.newProtocol(XmippProtPreprocessMicrographs, doDownsample=True, downFactor=5, doCrop=False, runMode=1)
        protDownsampling.inputMicrographs.set(protImport.outputMicrographs)
        self.launchProtocol(protDownsampling)
        self.assertIsNotNone(protDownsampling.outputMicrographs, "There was a problem with the downsampling")
#         self.validateFiles('protDownsampling', protDownsampling)

        # Now estimate CTF on the downsampled micrographs
        print "Performing CTF..."
        protCTF = self.newProtocol(XmippProtCTFMicrographs, numberOfThreads=4, minDefocus=2.2, maxDefocus=2.5)
        protCTF.inputMicrographs.set(protDownsampling.outputMicrographs)
        self.launchProtocol(protCTF)
        self.assertIsNotNone(protCTF.outputCTF, "There was a problem with the CTF estimation")
        # After CTF estimation, the output micrograph should have CTF info
#         self.validateFiles('protCTF', protCTF)

        print "Running fake particle picking..."
        protPP = self.newProtocol(XmippProtParticlePicking, importFolder=self.allCrdsDir)
        protPP.inputMicrographs.set(protDownsampling.outputMicrographs)
        self.launchProtocol(protPP)
#         self.protDict['protPicking'] = protPP
        self.assertIsNotNone(protPP.outputCoordinates, "There was a problem with the faked picking")

        print "Run extract particles with other downsampling factor"
        protExtract = self.newProtocol(XmippProtExtractParticles, boxSize=64, downsampleType=1, doFlip=True, downFactor=8, runMode=1, doInvert=True)
        protExtract.inputCoordinates.set(protPP.outputCoordinates)
        protExtract.ctfRelations.set(protCTF.outputCTF)
        protExtract.inputMicrographs.set(protImport.outputMicrographs)
        self.launchProtocol(protExtract)
        self.assertIsNotNone(protExtract.outputParticles, "There was a problem with the extract particles")
#         self.validateFiles('protExtract', protExtract)

        print "Run Projection Matching"
        protProjMatch = self.newProtocol(XmippProtProjMatch, ctfGroupMaxDiff=0.00001)
        protProjMatch.inputParticles.set(protExtract.outputParticles)
        protProjMatch.input3DReferences.set(protImportVol.outputVolume)
        self.launchProtocol(protProjMatch)
        self.assertIsNotNone(protProjMatch.outputVolumes, "There was a problem with Projection Matching")


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
