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
from test_protocols_xmipp import TestXmippBase


    
class TestXmippCeateMask3D(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupProject(cls)
        
#         cls.protImport = cls.runImportParticles(pattern=images, samplingRate=1, checkStack=False)
#         cls.iniVol = getInputPath('ml3dData', 'icoFiltered.vol')
    
    def testCreateMask1(self):
        print "Import volume"
        protImportVol = ProtImportVolumes(pattern=getInputPath('Volumes_BPV', 'BPV_scale_filtered_windowed_64.vol'), samplingRate=9.896)
        self.proj.launchProtocol(protImportVol, wait=True)
    
        print "Run create mask from volume"
        protMask1 = XmippProtCreateMask3D(source=0, threshold=0.4)
        protMask1.volume.set(protImportVol.outputVolume)
        self.proj.launchProtocol(protMask1, wait=True)        
        
        self.assertIsNotNone(protMask1.outputMask, "There was a problem with create mask from volume")          

        print "Run create mask from geometry"
        protMask2 = XmippProtCreateMask3D(source=1, size=64, samplingRate=9.89, geo=6, innerRadius=10, outerRadius=25, borderDecay=2)
        self.proj.launchProtocol(protMask2, wait=True)        
        
        self.assertIsNotNone(protMask2.outputMask, "There was a problem with create mask from geometry")          

        print "Run create mask from another mask"
        protMask3 = XmippProtCreateMask3D(source=2, doMorphological=True, elementSize=3)
        protMask3.inputMask.set(protMask1.outputMask)
        self.proj.launchProtocol(protMask3, wait=True)        
         
        self.assertIsNotNone(protMask3.outputMask, "There was a problem with mask from another mask")          

class TestXmippResolution3D(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupProject(cls)
        
#         cls.protImport = cls.runImportParticles(pattern=images, samplingRate=1, checkStack=False)
#         cls.iniVol = getInputPath('ml3dData', 'icoFiltered.vol')
    
    def testCalculateResolution(self):
        print "Import volume 1"
        protImportVol1 = ProtImportVolumes(pattern=getInputPath('Test_Resolution_Xmipp', 'volume_1_iter_002.mrc'), samplingRate=9.896)
        self.proj.launchProtocol(protImportVol1, wait=True)

        print "Import volume 2"
        protImportVol2 = ProtImportVolumes(pattern=getInputPath('Test_Resolution_Xmipp', 'volume_2_iter_002.mrc'), samplingRate=9.896)
        self.proj.launchProtocol(protImportVol2, wait=True)
    
        print "Run resolution 3D"
        protResol3D = XmippProtResolution3D(doSSNR=False)
        protResol3D.inputVolume.set(protImportVol1.outputVolume)
        protResol3D.referenceVolume.set(protImportVol2.outputVolume)
        self.proj.launchProtocol(protResol3D, wait=True)        
         
        self.assertIsNotNone(protResol3D._defineFscName(), "There was a problem with fsc")
        self.assertIsNotNone(protResol3D._defineStructFactorName(), "There was a problem with structure factor")
        
class TestXmippPreprocessVolumes(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupProject(cls)
        
#         cls.protImport = cls.runImportParticles(pattern=images, samplingRate=1, checkStack=False)
#         cls.iniVol = getInputPath('ml3dData', 'icoFiltered.vol')
    
    def testPreprocessVolumes(self):
        print "Import volume"
        protImportVol1 = ProtImportVolumes(pattern=getInputPath('Test_Resolution_Xmipp', 'volume_1_iter_002.mrc'), samplingRate=9.896)
        self.proj.launchProtocol(protImportVol1, wait=True)

        print "Import SetOfVolumes"
        protImportVol2 = ProtImportVolumes(pattern=getInputPath('Test_Resolution_Xmipp', '*.mrc'), samplingRate=9.896)
        self.proj.launchProtocol(protImportVol2, wait=True)
    
        print "Run preprocess a volume"
        protPreprocessVol1 = XmippProtPreprocessVolumes(doWindow=True, finalSize=300, doRandomize=True)
        protPreprocessVol1.inputVolumes.set(protImportVol1.outputVolume)
        self.proj.launchProtocol(protPreprocessVol1, wait=True)        
         
        self.assertIsNotNone(protPreprocessVol1.outputVol, "There was a problem with a volume")

        print "Run preprocess a SetOfVolumes"
        protPreprocessVol2 = XmippProtPreprocessVolumes(doWindow=True, finalSize=128, doFilter=True)
        protPreprocessVol2.inputVolumes.set(protImportVol2.outputVolumes)
        self.proj.launchProtocol(protPreprocessVol2, wait=True)        
         
        self.assertIsNotNone(protPreprocessVol2.outputVol, "There was a problem with preprocess a SetOfVolumes")


class TestXmippFilterVolumes(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupProject(cls)
        
#         cls.protImport = cls.runImportParticles(pattern=images, samplingRate=1, checkStack=False)
#         cls.iniVol = getInputPath('ml3dData', 'icoFiltered.vol')
    
    def testFilterVolumes(self):
        print "Import volume"
        protImportVol1 = ProtImportVolumes(pattern=getInputPath('Test_Resolution_Xmipp', 'volume_1_iter_002.mrc'), samplingRate=9.896)
        self.proj.launchProtocol(protImportVol1, wait=True)

        print "Import SetOfVolumes"
        protImportVol2 = ProtImportVolumes(pattern=getInputPath('Test_Resolution_Xmipp', 'volume_?_iter_002.mrc'), samplingRate=9.896)
        self.proj.launchProtocol(protImportVol2, wait=True)
    
        print "Run filter single volume"
        protFilterVolume = XmippProtFilterVolumes(lowFreq=0.1, highFreq=0.25)
        protFilterVolume.inputVolumes.set(protImportVol1.outputVolume)
        self.proj.launchProtocol(protFilterVolume, wait=True)        
        
        self.assertIsNotNone(protFilterVolume.outputVol, "There was a problem with filter a volume")

        print "Run filter SetOfVolumes"
        protFilterVolumes = XmippProtFilterVolumes(lowFreq=0.1, highFreq=0.25)
        protFilterVolumes.inputVolumes.set(protImportVol2.outputVolumes)
        self.proj.launchProtocol(protFilterVolumes, wait=True)        
        
        self.assertIsNotNone(protFilterVolumes.outputVol, "There was a problem with filter SetOfVolumes")


class TestXmippMaskVolumes(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupProject(cls)
        
#         cls.protImport = cls.runImportParticles(pattern=images, samplingRate=1, checkStack=False)
#         cls.iniVol = getInputPath('ml3dData', 'icoFiltered.vol')
    
    def testMaskVolumes(self):
        print "Import volume"
        protImportVol1 = ProtImportVolumes(pattern=getInputPath('Test_Resolution_Xmipp', 'volume_1_iter_002.mrc'), samplingRate=9.896)
        self.proj.launchProtocol(protImportVol1, wait=True)

        print "Import SetOfVolumes"
        protImportVol2 = ProtImportVolumes(pattern=getInputPath('Test_Resolution_Xmipp', 'volume_?_iter_002.mrc'), samplingRate=9.896)
        self.proj.launchProtocol(protImportVol2, wait=True)
    
        print "Run mask single volume"
        protMaskVolume = XmippProtMaskVolumes(radius=23)
        protMaskVolume.inputVolumes.set(protImportVol1.outputVolume)
        self.proj.launchProtocol(protMaskVolume, wait=True)        
        
        self.assertIsNotNone(protMaskVolume.outputVol, "There was a problem with applying mask to a volume")

        print "Run mask SetOfVolumes"
        protMaskVolumes = XmippProtMaskVolumes(geo=MASK3D_CROWN, innerRadius=18, outerRadius=23)
        protMaskVolumes.inputVolumes.set(protImportVol2.outputVolumes)
        self.proj.launchProtocol(protMaskVolumes, wait=True)        
        
        self.assertIsNotNone(protMaskVolumes.outputVol, "There was a problem with applying mask to SetOfVolumes")


class TestXmippSimAnnealing(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupProject(cls)
        
#         cls.protImport = cls.runImportParticles(pattern=images, samplingRate=1, checkStack=False)
#         cls.iniVol = getInputPath('ml3dData', 'icoFiltered.vol')
    
    def testSimAnnealing(self):
        
        """ Run an Import particles protocol. """
        project = self.proj
        pattern = os.environ.get('HEMOGLOBIN', getInputPath('particlesHemoglobin', '*.spi'))
        protImport = ProtImportParticles(pattern=pattern, samplingRate=3.5)
        project.launchProtocol(protImport, wait=True)
        # check that input images have been imported (a better way to do this?)
        if protImport.outputParticles is None:
            raise Exception('Import of images: %s, failed. outputParticles is None.' % pattern)

#         print "Run CL2D"
#         protCL2D = XmippProtCL2D(numberOfReferences=8, numberOfInitialReferences=2, 
#                                  numberOfIterations=3, numberOfMpi=3)
#         protCL2D.inputImages.set(protImport.outputParticles)
#         self.proj.launchProtocol(protCL2D, wait=True)        
#         
#         self.assertIsNotNone(protCL2D.outputClasses, "There was a problem with CL2D")

        print "Run ML2D"
        protML2D = XmippProtML2D(numberOfReferences=8, maxIters=4, doMlf=False,
                                 numberOfMpi=2, numberOfThreads=1)
        protML2D.inputParticles.set(protImport.outputParticles)
        self.proj.launchProtocol(protML2D, wait=True)        
        self.assertIsNotNone(protML2D.outputClasses, "There was a problem with ML2D") 

        print "Run Simulating annealing"
        protSimAnneal = XmippProtInitVolSimAnneal(symmetryGroup='d6', numberOfSimAnnealRef=2, percentRejection=0)
        protSimAnneal.inputClasses.set(protML2D.outputClasses)
        self.proj.launchProtocol(protSimAnneal, wait=True)        
        
        self.assertIsNotNone(protSimAnneal.outputVolumes, "There was a problem with simulating annealing protocol")


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
