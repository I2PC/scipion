
from pyworkflow.em import ProtImportMicrographsTiltPairs, ProtImportMicrographs, ProtImportCoordinates
from pyworkflow.em.packages.xmipp3 import XmippProtAssignmentTiltPair
from pyworkflow.tests import *
from test_workflow import TestWorkflow

# update this test when RCT workflow are implemented
class TestXmippAssignmentTiltPairsWorkflow(TestWorkflow):
    @classmethod
    def setUpClass(cls):    
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('rct')
        cls.allCrdsDir = cls.dataset.getFile('positions')
        cls.micsUFn = cls.dataset.getFile('untilted')
        cls.micsTFn = cls.dataset.getFile('tilted')
        
    def testXmippAssignmnetTiltPairsWorkflowBasic(self):
        #First, import a set of micrographs
        protImport = self.newProtocol(ProtImportMicrographsTiltPairs, 
                                      patternUntilted=self.micsUFn, patternTilted=self.micsTFn, 
                                      samplingRate=2.28, voltage=100, sphericalAberration=2.9)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputMicrographsTiltPair, "There was a problem with the import")
        
        protImportUntilted = self.newProtocol(ProtImportMicrographs, filesPath=self.micsUFn, samplingRate=2.28, voltage=100)
        self.launchProtocol(protImportUntilted)
        protImportTilted = self.newProtocol(ProtImportMicrographs, filesPath=self.micsTFn, samplingRate=2.28, voltage=100)
        self.launchProtocol(protImportTilted)
        
        protImportCoorU = self.newProtocol(ProtImportCoordinates,
                         importFrom=ProtImportCoordinates.IMPORT_FROM_XMIPP,
                         filesPath=self.allCrdsDir,
                         filesPattern='F_rct_u_*.pos', boxSize=100
                         )
        protImportCoorU.inputMicrographs.set(protImportUntilted.outputMicrographs)
        self.launchProtocol(protImportCoorU)
        protImportCoorT = self.newProtocol(ProtImportCoordinates,
                         importFrom=ProtImportCoordinates.IMPORT_FROM_XMIPP,
                         filesPath=self.allCrdsDir,
                         filesPattern='F_rct_t_*.pos', boxSize=100
                         )
        protImportCoorT.inputMicrographs.set(protImportTilted.outputMicrographs)
        self.launchProtocol(protImportCoorT)
    
                
        #Then simulate a particle picking               
        print "Running tilt pairs assignment..."   
        protAssigning = self.newProtocol(XmippProtAssignmentTiltPair)
                
        protAssigning.tiltpair.set(protImport.outputMicrographsTiltPair)
        print self.micsUFn
        print self.micsTFn
        protAssigning.untiltCoor.set(protImportCoorU.outputCoordinates)
        protAssigning.tiltCoor.set(protImportCoorT.outputCoordinates)

        print 'AQUI'
        self.launchProtocol(protAssigning)
        print 'AQUI'    
        self.assertIsNotNone(protAssigning.outputCoordinatesTiltPair, "There was a problem with the protocol assignment tilt pairs")
        #self.validateFiles('protPicking', protPicking)  

    def launchImportCoord(self, pathName, boxSize):
        protImportCoorT = self.newProtocol(ProtImportCoordinates,
                         importFrom=ProtImportCoordinates.IMPORT_FROM_XMIPP,
                         filesPath=spathName,
                         filesPattern='*.pos', boxSize=boxSize,
                         scale=3.,
                         invertX=False,
                         invertY=False
                         )
        protImportCoorT.inputMicrographs.set(protImportTilted.outputMicrographs)
        self.launchProtocol(protImportCoorT)
