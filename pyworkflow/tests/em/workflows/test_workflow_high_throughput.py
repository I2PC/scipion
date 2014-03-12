import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3 import *
from pyworkflow.em.packages.brandeis import *
from pyworkflow.em.packages.eman2 import *
from test_workflow import TestWorkflow


class HighThroughputTestDay1(TestWorkflow):
        
    @classmethod
    def setUpClass(cls):    
        # Create a new project
        setupProject(cls)
        cls.pattern = join(os.environ['HOME'], 'javi_movies', 'day1', '1??_*.mrcs')
        cls.importFolder = getInputPath('FakePick_HighThroughput')
#         cls.importVol = getInputPath('Ribosomes_Sjors', 'reference.mrc')
        
    def testWorkflow(self):
        #First, import a set of movies
        print "Importing a set of movies..."
        protImport = ProtImportMovies(pattern=self.pattern, samplingRate=1.2, magnification=60000,
                                      voltage=300, sphericalAberration=2.0)
        protImport.setObjLabel('import movies - Day1')
        self.proj.launchProtocol(protImport, wait=True)
        self.assertIsNotNone(protImport.outputMovies, "There was a problem importing movies")
        
        print "Aligning the movies..."
        protAlignMov = ProtOpticalAlignment()
        protAlignMov.inputMovies.set(protImport.outputMovies)
        protAlignMov.setObjLabel('align movies - Day1')
        self.proj.launchProtocol(protAlignMov, wait=True)
        self.assertIsNotNone(protAlignMov.outputMicrographs, "There was a problem aligning movies")
        
        print "Preprocessing the micrographs..."
        protPreprocess = XmippProtPreprocessMicrographs(doCrop=True, cropPixels=50)
        protPreprocess.inputMicrographs.set(protAlignMov.outputMicrographs)
        protPreprocess.setObjLabel('crop mics 50px')
        self.proj.launchProtocol(protPreprocess, wait=True)
        self.assertIsNotNone(protPreprocess.outputMicrographs, "There was a problem with the crop")
        
        # Now estimate CTF on the micrographs
        print "Performing CTF Micrographs..."
        protCTF = XmippProtCTFMicrographs(lowRes=0.04, highRes=0.31, minDefocus=0.3, maxDefocus=2,
                              runMode=1, numberOfMpi=1, numberOfThreads=3)
        protCTF.inputMicrographs.set(protPreprocess.outputMicrographs)
        protCTF.setObjLabel('ctf - Day1')
        self.proj.launchProtocol(protCTF, wait=True)
        self.assertIsNotNone(protCTF.outputCTF, "There was a problem with the CTF")
        
        print "Running Xmipp supervised fake particle picking..."
        protPP = XmippProtParticlePicking(importFolder=self.importFolder, runMode=1)                
        protPP.inputMicrographs.set(protPreprocess.outputMicrographs)
        protPP.setObjLabel('Picking - Day1')
        self.proj.launchProtocol(protPP, wait=True)
        self.assertIsNotNone(protPP.outputCoordinates, "There was a problem with the faked picking")
        
        print "Run extract particles with <Same as picking> option"
        protExtract = XmippProtExtractParticles(boxSize=240, downsampleType=1, doInvert=True, backRadius=110, runMode=1)
        protExtract.inputCoordinates.set(protPP.outputCoordinates)
        protExtract.ctfRelations.set(protCTF.outputCTF)
        protExtract.setObjLabel('extract particles')
        self.proj.launchProtocol(protExtract, wait=True)
        self.assertIsNotNone(protExtract.outputParticles, "There was a problem with the extract particles")
        
        print "Running Spider Filter"
        protFilter = SpiderProtFilter(lowFreq=0.07, highFreq=0.43)
        protFilter.inputParticles.set(protExtract.outputParticles)
        protFilter.setObjLabel('spi filter')
        self.proj.launchProtocol(protFilter, wait=True)
        self.assertIsNotNone(protFilter.outputParticles, "There was a problem with the Spider filter")
        
        print "Run Only Align2d"
        protOnlyAlign = XmippProtCL2DAlign(maximumShift=5, numberOfIterations=10, 
                                 numberOfMpi=4, numberOfThreads=1, useReferenceImage=False)
        protOnlyAlign.inputParticles.set(protFilter.outputParticles)
        protOnlyAlign.setObjLabel('cl2d align')
        self.proj.launchProtocol(protOnlyAlign, wait=True)        
        self.assertIsNotNone(protOnlyAlign.outputParticles, "There was a problem with Only align2d")
        
        print "Running Spider Dimension Reduction"
        protCAPCA = SpiderProtCAPCA()
        protCAPCA.inputParticles.set(protOnlyAlign.outputParticles)
        protCAPCA.setObjLabel('spi PCA')
        self.proj.launchProtocol(protCAPCA, wait=True)
        self.assertIsNotNone(protCAPCA.imcFile, "There was a problem with Spider Dimension Reduction")
        
        print "Running Spider Ward Classification"
        protWard = SpiderProtClassifyWard()
        protWard.pcaFilePointer.set(protCAPCA.imcFile)
        protWard.inputParticles.set(protOnlyAlign.outputParticles)
        protWard.setObjLabel('spi ward')
        self.proj.launchProtocol(protWard, wait=True)
        self.assertIsNotNone(protWard.outputClasses, "There was a problem with Spider Ward Classification")


if __name__ == "__main__":
    unittest.main()