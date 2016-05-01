# **************************************************************************
# *
# * Authors:    Laura del Cano (ldelcano@cnb.csic.es)
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

import sys, unittest

from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.relion import *


# Some utility functions to import micrographs that are used
# in several tests.
class TestRelionBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='xmipp_tutorial'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.particlesFn = cls.dataset.getFile('particles')
        cls.vol = cls.dataset.getFile('volumes')
    
    @classmethod
    def runImportParticles(cls, pattern, samplingRate, checkStack=False):
        """ Run an Import particles protocol. """
        protImport = cls.newProtocol(ProtImportParticles, 
                                      filesPath=pattern, samplingRate=samplingRate, 
                                      checkStack=checkStack)
        cls.launchProtocol(protImport)
        # check that input images have been imported (a better way to do this?)
        if protImport.outputParticles is None:
            raise Exception('Import of images: %s, failed. outputParticles is None.' % pattern)
        return protImport
    
    @classmethod
    def runNormalizeParticles(cls, particles):
        """ Run normalize particles protocol """
        protPreproc = cls.newProtocol(ProtRelionPreprocessParticles,
                                      doNormalize=True)
        protPreproc.inputParticles.set(particles)
        cls.launchProtocol(protPreproc)
        return protPreproc
    
    @classmethod
    def runImportVolumes(cls, pattern, samplingRate):
        """ Run an Import particles protocol. """
        protImport = cls.newProtocol(ProtImportVolumes, 
                                     filesPath=pattern, samplingRate=samplingRate)
        cls.launchProtocol(protImport)
        return protImport


class TestRelionClassify2D(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestRelionBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
        cls.protNormalize = cls.runNormalizeParticles(cls.protImport.outputParticles)
    
    def testRelion2D(self):                  
        print "Run relion2D"
        prot2D = self.newProtocol(ProtRelionClassify2D,
                                  doCTF=False, maskDiameterA=340,
                                  numberOfMpi=4, numberOfThreads=1)
        prot2D.numberOfClasses.set(4)
        prot2D.numberOfIterations.set(3)
        prot2D.inputParticles.set(self.protNormalize.outputParticles)
        self.launchProtocol(prot2D)        
        
        self.assertIsNotNone(getattr(prot2D, 'outputClasses', None), 
                             "There was a problem with Relion 2D:\n" + prot2D.getErrorMessage())
        self.assertAlmostEquals(self.protNormalize.outputParticles.getSamplingRate(), 
                                prot2D.outputClasses.getImages().getSamplingRate(),
                                "There was a problem with the sampling rate of the particles")
        for class2D in prot2D.outputClasses:
            self.assertTrue(class2D.hasAlignment2D())


class TestRelionClassify3D(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestRelionBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
#         cls.protNormalize = cls.runNormalizeParticles(cls.protImport.outputParticles)
        cls.protImportVol = cls.runImportVolumes(cls.vol, 3.5)
    
    def testProtRelionClassify3D(self):
        relionNormalize = self.newProtocol(ProtRelionPreprocessParticles)
        relionNormalize.inputParticles.set(self.protImport.outputParticles)
        relionNormalize.doNormalize.set(True)
        self.launchProtocol(relionNormalize)

        print "Run ProtRelionClassify3D"
        relion3DClass = self.newProtocol(ProtRelionClassify3D, 
                                         numberOfClasses=3, numberOfIterations=4, 
                                         doCTF=False, runMode=1, maskDiameterA=320,
                                         numberOfMpi=2, numberOfThreads=2)
        relion3DClass.inputParticles.set(relionNormalize.outputParticles)
        relion3DClass.referenceVolume.set(self.protImportVol.outputVolume)
        self.launchProtocol(relion3DClass)
        
        self.assertIsNotNone(getattr(relion3DClass, 'outputClasses', None), 
                             "There was a problem with Relion 3D:\n" + relion3DClass.getErrorMessage()) 

        for class3D in relion3DClass.outputClasses:
            self.assertTrue(class3D.hasAlignment3D())


class TestRelionRefine(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestRelionBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
#         cls.protNormalize = cls.runNormalizeParticles(cls.protImport.outputParticles)
        cls.protImportVol = cls.runImportVolumes(cls.vol, 3.5)
      
    def testProtRelionRefine(self):
        relionNormalize = self.newProtocol(ProtRelionPreprocessParticles)
        relionNormalize.inputParticles.set(self.protImport.outputParticles)
        relionNormalize.doNormalize.set(True)
        self.launchProtocol(relionNormalize)
  
        print "Run ProtRelionRefine"
        relionRefine = self.newProtocol(ProtRelionRefine3D, 
                                         doCTF=False, runMode=1, memoryPreThreads=1,
                                         maskDiameterA=340, symmetryGroup="d6",
                                         numberOfMpi=3, numberOfThreads=2)
        relionRefine.inputParticles.set(relionNormalize.outputParticles)
        relionRefine.referenceVolume.set(self.protImportVol.outputVolume)
        self.launchProtocol(relionRefine)
        
        self.assertIsNotNone(getattr(relionRefine, 'outputVolume', None), 
                             "There was a problem with Relion 3D:\n" + relionRefine.getErrorMessage())
        
        relionRefine._initialize() # Load filename templates
        dataSqlite =  relionRefine._getIterData(3)
        outImgSet = em.SetOfParticles(filename=dataSqlite)
        self.assertAlmostEqual(outImgSet[1].getSamplingRate(),
                               relionNormalize.outputParticles[1].getSamplingRate(),
                               "The sampling rate is wrong", delta=0.00001)
        
        self.assertAlmostEqual(outImgSet[1].getFileName(),
                               relionNormalize.outputParticles[1].getFileName(),
                               "The particles filenames are wrong")
        
class TestRelionPreprocess(TestRelionBase):
    """ This class helps to test all different preprocessing particles options on RElion. """
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestRelionBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
            
    def test_NormalizeAndDust(self):
        """ Normalize particles.
        """
        # Test now a normalization after the imported particles   
        protocol = self.newProtocol(ProtRelionPreprocessParticles,
                                    doNormalize=True, backRadius=40,
                                    doRemoveDust=True, whiteDust=4, blackDust=8)
        protocol.setObjLabel('relion: norm-dust')
        protocol.inputParticles.set(self.protImport.outputParticles)
        self.launchProtocol(protocol)

    def test_ScaleAndInvert(self):
        """ Test all options at once.
        """
        # Test now a normalization after the imported particles   
        protocol = self.newProtocol(ProtRelionPreprocessParticles,
                                    doNormalize=False,
                                    doScale=True, scaleSize=50,
                                    doInvert=True)
        protocol.setObjLabel('relion: scale-invert')
        protocol.inputParticles.set(self.protImport.outputParticles)
        
        self.launchProtocol(protocol)       


class TestRelionSubtract(TestRelionBase):
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
        
        protSubtract = self.newProtocol(ProtRelionSubtract)
        protSubtract.inputParticles.set(protParts.outputParticles)
        protSubtract.inputVolume.set(protVol.outputVolume)
        self.launchProtocol(protSubtract)
        self.assertIsNotNone(protSubtract.outputParticles, "There was a problem with subtract projection")




# class TestPolishParticles(TestRelionBase):
#     @classmethod
#     def setUpClass(cls):
#         setupTestProject(cls)
#         cls.dataset = DataSet.getDataSet('ribo_movies')
#         cls.movies = cls.dataset.getFile('movies')
#         cls.crdsDir = cls.dataset.getFile('posAllDir')
#         cls.vol = cls.dataset.getFile('volume')
#         cls.protImportVol = cls.runImportVolumes(cls.vol, 4.74)
#     
#     def test_polish(self):
#         #First, import a set of movies
#         print "Importing a set of movies..."
#         protImport = self.newProtocol(ProtImportMovies,
#                                       filesPath=self.movies,
#                                       samplingRate=2.37, magnification=59000,
#                                       voltage=300, sphericalAberration=2.0)
#         protImport.setObjLabel('import movies')
#         self.proj.launchProtocol(protImport, wait=True)
#         self.assertIsNotNone(protImport.outputMovies, "There was a problem importing movies")
#         
#         print "Aligning the movies..."
#         protAlignMov = self.newProtocol(ProtMovieAlignment, numberOfThreads=4)
#         protAlignMov.inputMovies.set(protImport.outputMovies)
#         protAlignMov.setObjLabel('align movies')
#         self.proj.launchProtocol(protAlignMov, wait=True)
#         self.assertIsNotNone(protAlignMov.outputMicrographs, "There was a problem aligning movies")
#          
#         print "Preprocessing the micrographs..."
#         protPreprocess = self.newProtocol(XmippProtPreprocessMicrographs, doCrop=True,
#                                           cropPixels=50, doDownsample=True, downFactor=2,
#                                           numberOfThreads=4)
#         protPreprocess.inputMicrographs.set(protAlignMov.outputMicrographs)
#         protPreprocess.setObjLabel('crop and Down')
#         self.proj.launchProtocol(protPreprocess, wait=True)
#         self.assertIsNotNone(protPreprocess.outputMicrographs, "There was a problem with the crop")
#          
#         # Now estimate CTF on the micrographs
#         print "Performing CTF Micrographs..."
#         protCTF = self.newProtocol(ProtCTFFind, lowRes=0.03, highRes=0.31, numberOfThreads=4)
#         protCTF.inputMicrographs.set(protPreprocess.outputMicrographs)
#         protCTF.setObjLabel('ctf')
#         self.proj.launchProtocol(protCTF, wait=True)
#         self.assertIsNotNone(protCTF.outputCTF, "There was a problem with the CTF")
#  
#         print "Running Xmipp Import Coordinates"
#         protPP = self.newProtocol(ProtImportCoordinates,
#                                  importFrom=ProtImportCoordinates.IMPORT_FROM_XMIPP,
#                                  filesPath=self.crdsDir,
#                                  filesPattern='*.pos', boxSize=120, scale=0.5)
#         protPP.inputMicrographs.set(protPreprocess.outputMicrographs)
#         protPP.setObjLabel('Picking')
#         self.launchProtocol(protPP)
#         self.assertIsNotNone(protPP.outputCoordinates, "There was a problem with the Xmipp import coordinates")
#  
#         print "Run extract particles with <Other> option"
#         protExtract = self.newProtocol(XmippProtExtractParticles,
#                                        boxSize=110, doInvert=True,
#                                                  doFlip=False, numberOfThreads=4)
#         protExtract.inputCoordinates.set(protPP.outputCoordinates)
#         protExtract.ctfRelations.set(protCTF.outputCTF)
#         protExtract.setObjLabel('extract particles')
#         self.proj.launchProtocol(protExtract, wait=True)
#         self.assertIsNotNone(protExtract.outputParticles, "There was a problem with the extract particles")
#         
#         print "Run Relion Refine"
#         proRef = self.newProtocol(ProtRelionRefine3D,
#                                   initialLowPassFilterA=40, maskDiameterA=500,
#                                   numberOfMpi=8, numberOfThreads=2)
#         proRef.inputParticles.set(protExtract.outputParticles)
#         proRef.referenceVolume.set(self.protImportVol.outputVolume)
#         proRef.setObjLabel('relion Refine')
#         self.launchProtocol(proRef)
#         self.assertIsNotNone(proRef.outputVolume, "There was a problem with Relion Refine:\n" + (proRef.getErrorMessage() or "No error set"))
#         self.assertIsNotNone(proRef.outputParticles, "There was a problem with Relion Refine:\n" + (proRef.getErrorMessage() or "No error set"))
