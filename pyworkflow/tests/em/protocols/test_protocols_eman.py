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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import unittest, sys
from pyworkflow.tests import *
from pyworkflow.em import *
from pyworkflow.em.packages.eman2 import *
from pyworkflow.em.protocol import ProtImportParticles


class TestEmanBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsRelion = DataSet.getDataSet('relion_tutorial')
    
    @classmethod
    def setData(cls, projectData='xmipp_tutorial'):
        cls.dataset = DataSet.getDataSet(projectData)
        cls.micsFn = cls.dataset.getFile('allMics')
        cls.crdsDir = cls.dataset.getFile('boxingDir')
        cls.particlesFn = cls.dataset.getFile('particles')
        cls.vol = cls.dataset.getFile('volumes')
    
    @classmethod
    def runImportMicrograph(cls, pattern, samplingRate, voltage, scannedPixelSize, magnification, sphericalAberration):
        """ Run an Import micrograph protocol. """
        # We have two options: passe the SamplingRate or the ScannedPixelSize + microscope magnification
        if not samplingRate is None:
            cls.protImport = cls.newProtocol(ProtImportMicrographs, samplingRateMode=0, filesPath=pattern, samplingRate=samplingRate, magnification=magnification, 
                                                   voltage=voltage, sphericalAberration=sphericalAberration)
        else:
            cls.protImport = cls.newProtocol(ProtImportMicrographs, samplingRateMode=1, filesPath=pattern, scannedPixelSize=scannedPixelSize, 
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
                                         filesPath=pattern, samplingRate=samplingRate,
                                         checkStack=checkStack)
        cls.launchProtocol(cls.protImport)
        # check that input images have been imported (a better way to do this?)
        if cls.protImport.outputParticles is None:
            raise Exception('Import of images: %s, failed. outputParticles is None.' % pattern)
        return cls.protImport

    @classmethod
    def runImportVolumes(cls, pattern, samplingRate):
        """ Run an Import particles protocol. """
        protImport = cls.newProtocol(ProtImportVolumes, 
                                     filesPath=pattern, samplingRate=samplingRate)
        cls.launchProtocol(protImport)
        return protImport

#     @classmethod
#     def runClassify(cls, particles):
#         cls.ProtClassify = cls.newProtocol(XmippProtML2D,
#                                            numberOfClasses=8, maxIters=4, doMlf=False,
#                                            numberOfMpi=2, numberOfThreads=2)
#         cls.ProtClassify.inputParticles.set(particles)
#         cls.launchProtocol(cls.ProtClassify)
#         return cls.ProtClassify
    

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
        protImportAvg = self.newProtocol(ProtImportAverages, filesPath=self.averages, checkStack=True, samplingRate=2.1)
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


class TestEmanReconstruct(TestEmanBase):
    def test_ReconstructEman(self):
        print "Import Set of particles with angles"
        prot1 = self.newProtocol(ProtImportParticles,
                                 objLabel='from scipion (to-reconstruct)',
                                 importFrom=ProtImportParticles.IMPORT_FROM_SCIPION,
                                 sqliteFile=self.dsRelion.getFile('import/case2/particles.sqlite'),
                                 magnification=10000,
                                 samplingRate=7.08
                                 )
        self.launchProtocol(prot1)
        
        print "Run Eman Reconstruct"
        protReconstruct = self.newProtocol(EmanProtReconstruct)
        protReconstruct.inputParticles.set(prot1.outputParticles)
        self.launchProtocol(protReconstruct)
        self.assertIsNotNone(protReconstruct.outputVolume, "There was a problem with eman reconstruction protocol")


class TestEmanRefine(TestEmanBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestEmanBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
        cls.protImportVol = cls.runImportVolumes(cls.vol, 3.5)
        
    def test_RefineEman(self):
        print "Run Eman Refine Easy"
        protRefine = self.newProtocol(EmanProtRefine, symmetry="d6", speed=6, numberOfIterations=1)
        protRefine.inputParticles.set(self.protImport.outputParticles)
        protRefine.input3DReference.set(self.protImportVol.outputVolume)
        self.launchProtocol(protRefine)
        self.assertIsNotNone(protRefine.outputVolume, "There was a problem with eman refine protocol")


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestEmanBoxing)
    unittest.TextTestRunner(verbosity=2).run(suite)