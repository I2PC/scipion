import unittest, sys
from pyworkflow.em import ProtImportMovies, ProtImportCoordinates
from pyworkflow.tests import DataSet, setupTestProject
from pyworkflow.em.packages.xmipp3 import (XmippProtOFAlignment, XmippProtPreprocessMicrographs,
                                           XmippProtCTFMicrographs, XmippProtExtractParticles,
                                           XmippProtCL2DAlign)
from pyworkflow.em.packages.xmipp3.constants import OTHER
from pyworkflow.em.packages.spider import SpiderProtFilter, SpiderProtCAPCA, SpiderProtClassifyWard
from test_workflow import TestWorkflow


class HighThroughputTest(TestWorkflow):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('ribo_movies')
        cls.movies = cls.dataset.getFile('movies')
        cls.crdsDir = cls.dataset.getFile('posAllDir')
    
    def test_workflow(self):
        #First, import a set of movies
        print "Importing a set of movies..."
        protImport = ProtImportMovies(filesPath=self.movies, samplingRate=2.37, magnification=59000,
                                      voltage=300, sphericalAberration=2.0)
        protImport.setObjLabel('import movies - Day1')
        self.proj.launchProtocol(protImport, wait=True)
        self.assertIsNotNone(protImport.outputMovies, "There was a problem importing movies")
        
        print "Aligning the movies..."
        protAlignMov = XmippProtOFAlignment(useAlignment=False)
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
        protCTF = XmippProtCTFMicrographs(lowRes=0.04, highRes=0.31, runMode=1, numberOfMpi=1, numberOfThreads=3)
        protCTF.inputMicrographs.set(protPreprocess.outputMicrographs)
        protCTF.setObjLabel('ctf - Day1')
        self.proj.launchProtocol(protCTF, wait=True)
        self.assertIsNotNone(protCTF.outputCTF, "There was a problem with the CTF")

        print "Running Xmipp Import Coordinates"
        protPP = self.newProtocol(ProtImportCoordinates,
                                 importFrom=ProtImportCoordinates.IMPORT_FROM_XMIPP,
                                 filesPath=self.crdsDir,
                                 filesPattern='*.pos', boxSize=110)
        protPP.inputMicrographs.set(protPreprocess.outputMicrographs)
        protPP.setObjLabel('Picking - Day1')
        self.launchProtocol(protPP)
        self.assertIsNotNone(protPP.outputCoordinates, "There was a problem with the Xmipp import coordinates")


        print "Run extract particles with <Other> option"
        protExtract = XmippProtExtractParticles(boxSize=60,
                                                downsampleType=OTHER,
                                                doInvert=True,
                                                doFlip=True,
                                                backRadius=28, runMode=1,
                                                numberOfThreads=3)
        protExtract.inputCoordinates.set(protPP.outputCoordinates)
        protExtract.ctfRelations.set(protCTF.outputCTF)
        protExtract.inputMicrographs.set(protPreprocess.outputMicrographs)
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
                                 numberOfMpi=8, numberOfThreads=1, useReferenceImage=False)
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
        protWard.pcaFile.set(protCAPCA.imcFile)
        protWard.inputParticles.set(protOnlyAlign.outputParticles)
        protWard.setObjLabel('spi ward')
        self.proj.launchProtocol(protWard, wait=True)


if __name__ == "__main__":
    unittest.main()