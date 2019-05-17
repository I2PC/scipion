# **************************************************************************
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import unittest, sys

from pyworkflow.em import ProtImportMovies, ProtImportCoordinates
from pyworkflow.tests import DataSet, setupTestProject
from test_workflow import TestWorkflow
from pyworkflow.utils import importFromPlugin

XmippProtOFAlignment = importFromPlugin('xmipp3.protocols', 'XmippProtOFAlignment', doRaise=True)
XmippProtPreprocessMicrographs = importFromPlugin('xmipp3.protocols','XmippProtPreprocessMicrographs')
XmippProtCTFMicrographs = importFromPlugin('xmipp3.protocols', 'XmippProtCTFMicrographs')
XmippProtExtractParticles = importFromPlugin('xmipp3.protocols', 'XmippProtExtractParticles')
XmippProtCL2DAlign = importFromPlugin('xmipp3.protocols', 'XmippProtCL2DAlign')
OTHER = importFromPlugin('xmipp3.constants', 'OTHER')
SpiderProtFilter = importFromPlugin('spider.protocols', 'SpiderProtFilter', doRaise=True)
SpiderProtCAPCA = importFromPlugin('spider.protocols', 'SpiderProtCAPCA')
SpiderProtClassifyWard = importFromPlugin('spider.protocols', 'SpiderProtClassifyWard')


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
        self.assertSetSize(protImport.outputMovies,msg="There was a problem importing movies")
        
        print "Aligning the movies..."
        protAlignMov = XmippProtOFAlignment(useAlignment=False,
                                            alignFrame0=1, alignFrameN=7,
                                            doApplyDoseFilter=False)
        protAlignMov.inputMovies.set(protImport.outputMovies)
        protAlignMov.setObjLabel('align movies - Day1')
        self.proj.launchProtocol(protAlignMov, wait=True)
        self.assertSetSize(protAlignMov.outputMicrographs, msg="There was a problem aligning movies")
        
        print "Preprocessing the micrographs..."
        protPreprocess = XmippProtPreprocessMicrographs(doCrop=True, cropPixels=50)
        protPreprocess.inputMicrographs.set(protAlignMov.outputMicrographs)
        protPreprocess.setObjLabel('crop mics 50px')
        self.proj.launchProtocol(protPreprocess, wait=True)
        self.assertSetSize(protPreprocess.outputMicrographs, msg="There was a problem with the crop")
        
        # Now estimate CTF on the micrographs
        print "Performing CTF Micrographs..."
        protCTF = XmippProtCTFMicrographs(lowRes=0.04, highRes=0.31, runMode=1, numberOfMpi=1, numberOfThreads=3)
        protCTF.inputMicrographs.set(protPreprocess.outputMicrographs)
        protCTF.setObjLabel('ctf - Day1')
        self.proj.launchProtocol(protCTF, wait=True)
        self.assertSetSize(protCTF.outputCTF, msg="There was a problem with the CTF")

        print "Running Xmipp Import Coordinates"
        protPP = self.newProtocol(ProtImportCoordinates,
                                 importFrom=ProtImportCoordinates.IMPORT_FROM_XMIPP,
                                 filesPath=self.crdsDir,
                                 filesPattern='*.pos', boxSize=110)
        protPP.inputMicrographs.set(protPreprocess.outputMicrographs)
        protPP.setObjLabel('Picking - Day1')
        self.launchProtocol(protPP)
        self.assertSetSize(protPP.outputCoordinates, msg="There was a problem with the Xmipp import coordinates")


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
        self.assertSetSize(protExtract.outputParticles, msg="There was a problem with the extract particles")
        
        print "Running Spider Filter"
        protFilter = SpiderProtFilter(lowFreq=0.07, highFreq=0.43)
        protFilter.inputParticles.set(protExtract.outputParticles)
        protFilter.setObjLabel('spi filter')
        self.proj.launchProtocol(protFilter, wait=True)
        self.assertSetSize(protFilter.outputParticles, msg="There was a problem with the Spider filter")
        
        print "Run Only Align2d"
        protOnlyAlign = XmippProtCL2DAlign(maximumShift=5, numberOfIterations=10, 
                                 numberOfMpi=8, numberOfThreads=1, useReferenceImage=False)
        protOnlyAlign.inputParticles.set(protFilter.outputParticles)
        protOnlyAlign.setObjLabel('cl2d align')
        self.proj.launchProtocol(protOnlyAlign, wait=True)        
        self.assertSetSize(protOnlyAlign.outputParticles, msg="There was a problem with Only align2d")
        
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