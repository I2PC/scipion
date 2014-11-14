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

import unittest, sys
from pyworkflow.tests import *
from pyworkflow.em import *
from pyworkflow.em.packages.eman2 import *


class TestEmanBase(BaseTest):
    @classmethod
    def setData(cls, projectData='xmipp_tutorial'):
        cls.dataset = DataSet.getDataSet(projectData)
        cls.micsFn = cls.dataset.getFile('allMics')
        cls.crdsDir = cls.dataset.getFile('boxingDir')
    
    @classmethod
    def runImportMicrograph(cls, pattern, samplingRate, voltage, scannedPixelSize, magnification, sphericalAberration):
        """ Run an Import micrograph protocol. """
        # We have two options: passe the SamplingRate or the ScannedPixelSize + microscope magnification
        if not samplingRate is None:
            cls.protImport = cls.newProtocol(ProtImportMicrographs, samplingRateMode=0, pattern=pattern, samplingRate=samplingRate, magnification=magnification, 
                                                   voltage=voltage, sphericalAberration=sphericalAberration)
        else:
            cls.protImport = cls.newProtocol(ProtImportMicrographs, samplingRateMode=1, pattern=pattern, scannedPixelSize=scannedPixelSize, 
                                                   voltage=voltage, magnification=magnification, sphericalAberration=sphericalAberration)
            
        cls.proj.launchProtocol(cls.protImport, wait=True)
        # check that input micrographs have been imported (a better way to do this?)
        if cls.protImport.outputMicrographs is None:
            raise Exception('Import of micrograph: %s, failed. outputMicrographs is None.' % pattern)
        return cls.protImport
    
    @classmethod
    def runImportMicrographBPV(cls, pattern):
        """ Run an Import micrograph protocol. """
        return cls.runImportMicrograph(pattern, samplingRate=1.237, voltage=300, sphericalAberration=2, scannedPixelSize=None, magnification=56000)

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
    

class TestEmanBoxing(TestEmanBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestEmanBase.setData()
        cls.protImport = cls.runImportMicrographBPV(cls.micsFn)
    
    def testCreateOutput(self):
        print "Running Eman fake particle picking..."
        protPP = EmanProtBoxing(importFolder=self.crdsDir, runMode=1)
        protPP.inputMicrographs.set(self.protImport.outputMicrographs)
        protPP.boxSize.set(550)
        self.proj.launchProtocol(protPP, wait=True)
        self.assertIsNotNone(protPP.outputCoordinates, "There was a problem with the faked picking")

class TestEmanInitialModelMda(TestEmanBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('mda')
        cls.averages = cls.dataset.getFile('averages')
        cls.samplingRate = 3.5
        cls.symmetry = 'd6'
        cls.numberOfIterations = 5
        cls.numberOfModels = 2

    def test_initialmodel(self):
        #Import a set of averages
        print "Import Set of averages"
        protImportAvg = self.newProtocol(ProtImportAverages, pattern=self.averages, checkStack=True, samplingRate=2.1)
        self.launchProtocol(protImportAvg)
        self.assertIsNotNone(protImportAvg.getFiles(), "There was a problem with the import")

        print "Run Initial model"
        protIniModel = self.newProtocol(EmanProtInitModel,
                                      symmetry=self.symmetry, numberOfIterations=self.numberOfIterations, numberOfModels=self.numberOfModels, numberOfThreads=4)
        protIniModel.inputSet.set(protImportAvg.outputAverages)
        self.launchProtocol(protIniModel)
        self.assertIsNotNone(protIniModel.outputVolumes, "There was a problem with eman initial model protocol")

class TestEmanInitialModelGroel(TestEmanInitialModelMda):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('groel')
        cls.averages = cls.dataset.getFile('averages')
        cls.samplingRate = 2.1
        cls.symmetry = 'd7'
        cls.numberOfIterations = 10
        cls.numberOfModels = 10

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestEmanBoxing)
    unittest.TextTestRunner(verbosity=2).run(suite)