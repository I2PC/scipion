
from pyworkflow.em import ProtImportMicrographsTiltPairs, ProtUserSubSet
from pyworkflow.em.packages.xmipp3 import XmippProtParticlePickingPairs, XmippProtExtractParticlesPairs, XmippProtCL2D, XmippProtRCT
from pyworkflow.em.packages.xmipp3.constants import SAME_AS_PICKING, OTHER 
from pyworkflow.tests import *
from test_workflow import TestWorkflow

# update this test when RCT workflow are implemented
class TestXmippRCTWorkflow(TestWorkflow):
    @classmethod
    def setUpClass(cls):    
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('rct')
        cls.allCrdsDir = cls.dataset.getFile('positions')
        cls.micsUFn = cls.dataset.getFile('untilted')
        cls.micsTFn = cls.dataset.getFile('tilted')
        cls.classesSqlite = cls.dataset.getFile('classes')
        
    def testXmippRCTWorkflowBasic(self):
        #First, import a set of micrographs
        protImport = self.newProtocol(ProtImportMicrographsTiltPairs, 
                                      patternUntilted=self.micsUFn,
                                      patternTilted=self.micsTFn,
                                      samplingRate=2.28, voltage=100,
                                      sphericalAberration=2.9)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputMicrographsTiltPair,
                             "There was a problem with the import")
        #self.validateFiles('protImportRCT', protImport)
                
        #Then simulate a particle picking               
        print "Running fake particle picking..."   
        protPicking = self.newProtocol(XmippProtParticlePickingPairs,
                                       importFolder=self.allCrdsDir)
                
        protPicking.inputMicrographsTiltedPair.set(protImport.outputMicrographsTiltPair)
        
        self.proj.launchProtocol(protPicking, wait=True)
            
        self.assertIsNotNone(protPicking.outputCoordinatesTiltPair,
                             "There was a problem with the faked picking")
        #self.validateFiles('protPicking', protPicking)  

        #Extract particles    
        print "Run extract particles with Same as picking"
        protExtract = self.newProtocol(XmippProtExtractParticlesPairs,
                                       downFactor=2,
                                       boxSize=60,
                                       doInvert=False,
                                       downsampleType=SAME_AS_PICKING)
        protExtract.inputCoordinatesTiltedPairs.set(protPicking.outputCoordinatesTiltPair)

        self.proj.launchProtocol(protExtract, wait=True)
        
        #self.validateFiles('protExtract', protExtract)  
        self.assertIsNotNone(protExtract.outputParticlesTiltPair,
                             "There was a problem with the extract particles")

        # Classify using Xmipp CL2D
        print "Run CL2D"
        protCL2D = self.newProtocol(XmippProtCL2D,
                                    numberOfClasses=10,
                                    numberOfInitialClasses=1,
                                    numberOfIterations=3, numberOfMpi=2)
        protCL2D.inputParticles.set(protExtract.outputParticlesTiltPair.getUntilted())
        self.launchProtocol(protCL2D)        
        self.assertIsNotNone(protCL2D.outputClasses,
                             "There was a problem with CL2D")
        #self.validateFiles('protCL2D', protCL2D) 
        
        # Random Conical Tilt
        print "Run Random Conical Tilt"
        protRCT = self.newProtocol(XmippProtRCT)
        protRCT.inputParticlesTiltPair.set(protExtract.outputParticlesTiltPair)
        protRCT.inputParticles.set(protCL2D.outputClasses)

        self.proj.launchProtocol(protRCT, wait=True)
        
        #self.validateFiles('protExtract', protExtract)  
        self.assertIsNotNone(protRCT.outputVolumes, "There was a problem with the RCT")

