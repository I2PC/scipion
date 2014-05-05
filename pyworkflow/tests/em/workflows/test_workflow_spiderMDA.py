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
        self.proj.launchProtocol(protImport, wait=True)
        # check that input images have been imported (a better way to do this?)
        if protImport.outputParticles is None:
            raise Exception('Import of images: %s, failed. outputParticles is None.' % pattern)
        
        protFilter = SpiderProtFilter()
        protFilter.inputParticles.set(protImport)
        protFilter.inputParticles.setExtendedAttribute('outputParticles')
        self.proj.launchProtocol(protFilter, wait=True)
        
        protAPSR = SpiderProtAlignAPSR()
        protAPSR.inputParticles.set(protFilter.outputParticles)
        self.proj.launchProtocol(protAPSR, wait=True)
         
        protMask = SpiderProtCustomMask()
        protMask.inputImage.set(protAPSR.outputAverage)
        self.proj.launchProtocol(protMask, wait=True)       
              
        protCAPCA = SpiderProtCAPCA()
        protCAPCA.maskType.set(1)
        protCAPCA.maskImage.set(protMask.outputMask)
        protCAPCA.inputParticles.set(protAPSR.outputParticles)
        self.proj.launchProtocol(protCAPCA, wait=True)
        
        protWard = SpiderProtClassifyWard()
        protWard.pcaFilePointer.set(protCAPCA.imcFile)
        protWard.inputParticles.set(protAPSR.outputParticles)
        self.proj.launchProtocol(protWard, wait=True)


if __name__ == "__main__":
    unittest.main()
