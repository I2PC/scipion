import unittest, sys, os
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3 import *
from test_workflow import TestWorkflow
   
       
class TestSpiderWorkflow(TestWorkflow):
    @classmethod
    def setUpClass(cls):    
        # Create a new project
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('mda')
        cls.particlesFn = cls.dataset.getFile('particles')
    
    def test_mdaWorkflow(self):
        """ Run an Import particles protocol. """
        protImport = ProtImportParticles(pattern=self.particlesFn, samplingRate=3.5)
        self.launchProtocol(protImport)
        # check that input images have been imported (a better way to do this?)
        if protImport.outputParticles is None:
            raise Exception('Import of images: %s, failed. outputParticles is None.' % pattern)
        
        protFilter = SpiderProtFilter()
        protFilter.inputParticles.set(protImport)
        protFilter.inputParticles.setExtendedAttribute('outputParticles')
        self.launchProtocol(protFilter)
        
        protAPSR = SpiderProtAlignAPSR()
        protAPSR.inputParticles.set(protFilter.outputParticles)
        self.launchProtocol(protAPSR)
        
        protPairwise = SpiderProtAlignPairwise()
        protPairwise.inputParticles.set(protFilter.outputParticles)
        self.launchProtocol(protPairwise)       
         
        protMask = SpiderProtCustomMask()
        protMask.inputImage.set(protAPSR.outputAverage)
        self.launchProtocol(protMask)       
              
        protCAPCA = SpiderProtCAPCA()
        protCAPCA.maskType.set(1)
        protCAPCA.maskImage.set(protMask.outputMask)
        protCAPCA.inputParticles.set(protAPSR.outputParticles)
        self.launchProtocol(protCAPCA)
        
        protWard = SpiderProtClassifyWard()
        protWard.pcaFile.set(protCAPCA.imcFile)
        protWard.inputParticles.set(protAPSR.outputParticles)
        self.launchProtocol(protWard)
        
        protKmeans = SpiderProtClassifyKmeans()
        protKmeans.pcaFile.set(protCAPCA.imcFile)
        protKmeans.inputParticles.set(protAPSR.outputParticles)
        protKmeans.numberOfClasses.set(4)
        self.launchProtocol(protKmeans)
        
        protDiday = SpiderProtClassifyDiday()
        protDiday.pcaFile.set(protCAPCA.imcFile)
        protDiday.inputParticles.set(protAPSR.outputParticles)
        self.launchProtocol(protDiday)               
        
    def test_mdaPairwise(self):
        protImport = ProtImportParticles(pattern=self.particlesFn, samplingRate=3.5)
        self.launchProtocol(protImport)
        # check that input images have been imported (a better way to do this?)
        if protImport.outputParticles is None:
            raise Exception('Import of images: %s, failed. outputParticles is None.' % pattern)
        
        protPairwise = SpiderProtAlignPairwise()
        protPairwise.inputParticles.set(protImport.outputParticles)
        protPairwise.innerRadius.set(48)
        self.launchProtocol(protPairwise)
        
    def test_StoredPointers(self):
        """ Test the usage of pointer extensions to store a 
        workflow dependencies before the current execution 
        of the protocols.
        """
        protImport = ProtImportParticles(pattern=self.particlesFn, samplingRate=3.5)
        self.saveProtocol(protImport)
        
        protFilter = SpiderProtFilter()
        protFilter.inputParticles.set(protImport)
        protFilter.inputParticles.setExtendedAttribute('outputParticles')
        self.saveProtocol(protFilter)
        
        protAPSR = SpiderProtAlignAPSR()
        protAPSR.inputParticles.set(protFilter)
        protAPSR.inputParticles.setExtendedAttribute('outputParticles')
        self.saveProtocol(protAPSR)
        
        return 
        # Launch all stored protocol runs
        self.launchProtocol(protImport)
        # check that input images have been imported (a better way to do this?)
        if protImport.outputParticles is None:
            raise Exception('Import of images: %s, failed. outputParticles is None.' % pattern)
        self.launchProtocol(protFilter)
        self.launchProtocol(protAPSR)       

if __name__ == "__main__":
    unittest.main()
