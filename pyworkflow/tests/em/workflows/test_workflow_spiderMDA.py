import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3 import *
from test_workflow import TestWorkflow
   
       
class TestSpiderWorkflow(TestWorkflow):
    
    @classmethod
    def setUpClass(cls):    
        # Create a new project
        setupProject(cls)
        cls.pattern = getInputPath('Micrographs_BPV3', '*.mrc')        
        cls.importFolder = getInputPath('Picking_XmippBPV3_Down3')
            
    def testMSAWorkflow(self):
        """ Run an Import particles protocol. """
        project = self.proj
        pattern = getInputPath('particlesHemoglobin', '*.spi')
        protImport = ProtImportParticles(pattern=pattern, samplingRate=3.5)
        project.launchProtocol(protImport, wait=True)
        # check that input images have been imported (a better way to do this?)
        if protImport.outputParticles is None:
            raise Exception('Import of images: %s, failed. outputParticles is None.' % pattern)

        
        protFilter = SpiderProtFilter()
        protFilter.inputParticles.set(protImport.outputParticles)
        project.launchProtocol(protFilter, wait=True)
        
        protAPSR = SpiderProtAlignAPSR()
        protAPSR.inputParticles.set(protFilter.outputParticles)
        project.launchProtocol(protAPSR, wait=True)
         
        protMask = SpiderProtCustomMask()
        protMask.inputImage.set(protAPSR.outputAverage)
        project.launchProtocol(protMask, wait=True)       
              
        protCAPCA = SpiderProtCAPCA()
        protCAPCA.maskType.set(1)
        protCAPCA.maskImage.set(protMask.outputMask)
        protCAPCA.inputParticles.set(protAPSR.outputParticles)
        project.launchProtocol(protCAPCA, wait=True)
        
        protWard = SpiderProtClassifyWard()
        protWard.pcaFilePointer.set(protCAPCA.imcFile)
        protWard.inputParticles.set(protAPSR.outputParticles)
        project.launchProtocol(protWard, wait=True)
        


if __name__ == "__main__":
    unittest.main()
