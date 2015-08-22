
import unittest, sys
from pyworkflow.em import ProtImportMicrographsTiltPairs, ProtUserSubSet
from pyworkflow.em.packages.xmipp3 import XmippProtParticlePickingPairs, XmippProtExtractParticlesPairs, XmippProtCL2D, XmippProtRCT
from pyworkflow.em.packages.xmipp3.constants import SAME_AS_PICKING, OTHER 
from pyworkflow.em import SetOfParticles
from pyworkflow.tests import *
from test_workflow import TestWorkflow

# update this test when RCT workflow are implemented
class TestXmippRCTWorkflow(TestWorkflow):
    
#     GOLD_FILES = {
#             'protImport': ['protImport/micrograph001U.mrc',
#                     'protImport/micrographs.sqlite', 
#                      'protImport/micrograph003T.mrc',
#                     'protImport/micrograph001T.mrc',
#                     'protImport/micrograph002U.mrc',
#                     'protImport/micrograph002T.mrc',
#                     'protImport/micrograph003U.mrc'],
#               'protDownsampling': ['protDownsampling/micrograph002U.mrc',
#                     'protImport/micrograph003T.mrc',
#                     'protImport/micrograph001T.mrc',
#                     'protDownsampling/micrograph001U.mrc',
#                     'protImport/micrograph001U.mrc',
#                     'protDownsampling/micrograph001T.mrc',
#                     'protImport/micrograph003U.mrc',
#                     'protDownsampling/micrograph003T.mrc',
#                     'protDownsampling/micrograph003U.mrc',
#                     'protImport/micrographs.sqlite',
#                     'protImport/micrograph002U.mrc',
#                     'protDownsampling/micrographs.xmd',
#                     'protDownsampling/micrograph002T.mrc',
#                     'protDownsampling/log/protocol.log',
#                     'protImport/micrograph002T.mrc'],
#               'protPicking':['protDownsampling/micrograph002U.mrc',
#                     'protPicking/extra/micrograph002T.pos',
#                     'protDownsampling/micrograph003U.mrc',
#                     'protPicking/extra/families.xmd',
#                     'protDownsampling/micrograph001U.mrc',
#                     'protPicking/extra/micrograph001U.pos',
#                     'protPicking/extra/micrograph003U.pos',
#                     'protDownsampling/micrograph001T.mrc',
#                     'protPicking/extra/micrograph001T.pos',
#                     'protPicking/extra/micrograph003T.pos',
#                     'protDownsampling/micrograph003T.mrc',
#                     'protPicking/log/protocol.log',
#                     'protDownsampling/micrographs.xmd',
#                     'protDownsampling/micrograph002T.mrc',
#                     'protPicking/extra/micrograph002U.pos'],
#               'protExtract':[
#                     'protPicking/extra/micrograph002T.pos',
#                     'protExtract/tmp/micrograph001T_noDust.xmp',
#                     'protExtract/tmp/micrograph003T_noDust.xmp',
#                     'protExtract/tmp/micrograph002T_noDust.xmp',
#                     'protExtract/extra/micrograph002T.stk',
#                     'protExtract/images.xmd',
#                     'protExtract/extra/micrograph003U.stk',
#                     'protExtract/extra/micrograph002U.pos',
#                     'protExtract/extra/micrograph001T.xmd',
#                     'protExtract/tmp/micrograph002U_noDust.xmp',
#                     'protExtract/extra/micrograph002T.pos',
#                     'protExtract/extra/micrograph002U.xmd',
#                     'protPicking/extra/micrograph003T.pos',
#                     'protExtract/extra/micrograph001T.pos',
#                     'protExtract/extra/micrograph001U.xmd',
#                     'protExtract/extra/micrograph001U.pos',
#                     'protExtract/extra/micrograph003T.pos',
#                     'protExtract/extra/micrograph003U.xmd',
#                     'protExtract/extra/micrograph003T.xmd',
#                     'protExtract/extra/micrograph003U.pos',
#                     'protPicking/extra/micrograph001U.pos',
#                     'protExtract/log/protocol.log',
#                     'protPicking/extra/micrograph001T.pos',
#                     'protExtract/tmp/micrograph003U_noDust.xmp',
#                     'protExtract/extra/micrograph001U.stk',
#                     'protExtract/extra/micrograph002T.xmd',
#                     'protExtract/tmp/micrograph001U_noDust.xmp',
#                     'protExtract/extra/micrograph001T.stk',
#                     'protExtract/extra/micrograph002U.stk',
#                     'protPicking/extra/micrograph002U.pos',
#                     'protPicking/extra/micrograph003U.pos',
#                     'protExtract/extra/micrograph003T.stk',
#                     ],
#               'protML2D': [
#                    'protML2D/ml2d_extra/iter003/iter_images.xmd',
#                     'protML2D/ml2d_extra/iter004/iter_classes.xmd',
#                     'protML2D/log/protocol.log',
#                     'protExtract/extra/micrograph001U.stk',
#                     'protML2D/ml2d_extra/iter004/iter_classes.stk',
#                     'protExtract/extra/micrograph002T.stk',
#                     'protML2D/ml2d_extra/iter001/iter_classes.stk',
#                     'protExtract/images.xmd',
#                     'protExtract/extra/micrograph003U.stk',
#                     'protML2D/ml2d_extra/iter002/iter_images.xmd',
#                     'protML2D/ml2d_extra/iter003/iter_classes.xmd',
#                     'protExtract/extra/micrograph002U.stk',
#                     'protML2D/ml2d_extra/iter003/iter_classes.stk',
#                     'protML2D/ml2d_extra/iter001/iter_images.xmd',
#                     'protML2D/ml2d_extra/iter004/iter_images.xmd',
#                     'protML2D/ml2d_classes.xmd',
#                     'protML2D/ml2d_extra/iter001/iter_classes.xmd',
#                     'protML2D/ml2d_extra/iter002/iter_classes.stk',
#                     'protML2D/ml2d_extra/iter002/iter_classes.xmd',
#                     'protML2D/ml2d__images_average.xmp',
#                     'protML2D/ml2d_classes.stk',
#                     'protML2D/ml2d_images.xmd',
#                     'protExtract/extra/micrograph001T.stk',
#                     'protExtract/extra/micrograph003T.stk'],
#               'protCL2D': ['protCL2D/extra/classes_core_hierarchy.txt', 
#                     'protCL2D/extra/level_01/level_classes_core.xmd', 
#                     'protCL2D/extra/level_01/level_classes.stk', 
#                     'protCL2D/extra/level_01/classes_sorted.stk', 
#                     'protCL2D/extra/level_00/level_classes.xmd', 
#                     'protCL2D/extra/level_00/classes_sorted.xmd', 
#                     'protCTF/extra/BPV_1386/xmipp_ctf.ctfparam', 
#                     'protCTF/extra/BPV_1388/xmipp_ctf.ctfparam', 
#                     'protCL2D/extra/level_00/level_classes_core.stk', 
#                     'protCL2D/log/protocol.log', 
#                     'protExtract/extra/BPV_1388.stk', 
#                     'protCTF/extra/BPV_1387/xmipp_ctf.ctfparam', 
#                     'protCL2D/extra/level_00/classes_core_sorted.xmd', 
#                     'protCL2D/extra/level_01/classes_core_sorted.stk', 
#                     'protExtract/images.xmd', 
#                     'protCL2D/extra/level_00/classes_sorted.stk', 
#                     'protCL2D/extra/level_01/classes_sorted.xmd', 
#                     'protCL2D/extra/classes_hierarchy.txt', 
#                     'protCL2D/extra/images.xmd', 
#                     'protExtract/extra/BPV_1386.stk', 
#                     'protCL2D/extra/level_01/classes_core_sorted.xmd', 
#                     'protCL2D/extra/level_00/classes_core_sorted.stk', 
#                     'protExtract/extra/BPV_1387.stk', 
#                     'protCL2D/extra/level_00/level_classes_core.xmd', 
#                     'protCL2D/extra/level_01/level_classes.xmd', 
#                     'protCL2D/extra/level_00/level_classes.stk', 
#                     'protCL2D/extra/level_01/level_classes_core.stk']
#               }
    
    @classmethod
    def setUpClass(cls):    
        setupTestProject(cls)
#         cls.dataset = DataSet.getDataSet('xmipp_tutorial')
#         cls.allCrdsDir = cls.dataset.getFile('rctCoords')
#         cls.micsUFn = cls.dataset.getFile('rctMicsU')
#         cls.micsTFn = cls.dataset.getFile('rctMicsT')
        cls.dataset = DataSet.getDataSet('rct')
        cls.allCrdsDir = cls.dataset.getFile('positions')
        cls.micsUFn = cls.dataset.getFile('untilted')
        cls.micsTFn = cls.dataset.getFile('tilted')
        cls.classesSqlite = cls.dataset.getFile('classes')
        
    def testXmippRCTWorkflowBasic(self):
        #First, import a set of micrographs
        protImport = self.newProtocol(ProtImportMicrographsTiltPairs, 
                                      patternUntilted=self.micsUFn, patternTilted=self.micsTFn, 
                                      samplingRate=2.28, voltage=100, sphericalAberration=2.9)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputMicrographsTiltPair, "There was a problem with the import")
        #self.validateFiles('protImportRCT', protImport)
                
        #Then simulate a particle picking               
        print "Running fake particle picking..."   
        protPicking = self.newProtocol(XmippProtParticlePickingPairs, importFolder=self.allCrdsDir)
                
        protPicking.inputMicrographsTiltedPair.set(protImport.outputMicrographsTiltPair)
        
        self.proj.launchProtocol(protPicking, wait=True)
            
        self.assertIsNotNone(protPicking.outputCoordinatesTiltPair, "There was a problem with the faked picking")
        #self.validateFiles('protPicking', protPicking)  

        #Extract particles    
        print "Run extract particles with Same as picking"
        protExtract = self.newProtocol(XmippProtExtractParticlesPairs, boxSize=60, downsampleType=OTHER, downFactor=2)
        protExtract.inputCoordinatesTiltedPairs.set(protPicking.outputCoordinatesTiltPair)

        self.proj.launchProtocol(protExtract, wait=True)
        
        #self.validateFiles('protExtract', protExtract)  
        self.assertIsNotNone(protExtract.outputParticlesTiltPair, "There was a problem with the extract particles")


        # Classify using Xmipp CL2D
        print "Run CL2D"
        protCL2D = self.newProtocol(XmippProtCL2D, numberOfClasses=10, numberOfInitialClasses=1,
                                 numberOfIterations=3, numberOfMpi=2)
        protCL2D.inputParticles.set(protExtract.outputParticlesTiltPair.getUntilted())
        self.launchProtocol(protCL2D)        
        self.assertIsNotNone(protCL2D.outputClasses, "There was a problem with CL2D")
        #self.validateFiles('protCL2D', protCL2D) 
        
        # Make a subset of images with only particles from class 5
#         print "Run ProtUserSubset"
#         protSubset = self.newProtocol(ProtUserSubSet)
# 
#         protSubset.inputObject.set(protCL2D.outputClasses)
#         protSubset.inputObject.setExtended(5)
#         protSubset.outputClassName.set("SetOfParticles")
#         protSubset.sqliteFile.set("%s," % self.classesSqlite)
#         self.launchProtocol(protSubset)        
#         self.assertIsNotNone(protSubset.outputParticles, "There was a problem with ProtUSerSubset")
        #self.validateFiles('ProtUserSubSet', ProtUserSubSet) 
        
        # Random Conical Tilt    
        print "Run Random Conical Tilt"
        protRCT = self.newProtocol(XmippProtRCT)
        protRCT.inputParticlesTiltPair.set(protExtract.outputParticlesTiltPair)
        protRCT.inputParticles.set(protCL2D.outputClasses)

        self.proj.launchProtocol(protRCT, wait=True)
        
        #self.validateFiles('protExtract', protExtract)  
        self.assertIsNotNone(protRCT.outputVolumes, "There was a problem with the RCT")

        
if __name__ == "__main__":
#    suite = unittest.TestLoader().loadTestsFromTestCase(TestXmippTiltedWorkflow)    
#    unittest.TextTestRunner(verbosity=2).run(suite)
    unittest.main()