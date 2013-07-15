
import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3 import *
from test_workflow import TestWorkflow


class TestXmippTiltedWorkflow(TestWorkflow):
    
    GOLD_FILES = {
            'protImport': ['protImport/micrograph001U.mrc',
                    'protImport/micrographs.sqlite', 
                     'protImport/micrograph003T.mrc',
                    'protImport/micrograph001T.mrc',
                    'protImport/micrograph002U.mrc',
                    'protImport/micrograph002T.mrc',
                    'protImport/micrograph003U.mrc'],
              'protDownsampling': ['protDownsampling/micrograph002U.mrc',
                    'protImport/micrograph003T.mrc',
                    'protImport/micrograph001T.mrc',
                    'protDownsampling/micrograph001U.mrc',
                    'protImport/micrograph001U.mrc',
                    'protDownsampling/micrograph001T.mrc',
                    'protImport/micrograph003U.mrc',
                    'protDownsampling/micrograph003T.mrc',
                    'protDownsampling/micrograph003U.mrc',
                    'protImport/micrographs.sqlite',
                    'protImport/micrograph002U.mrc',
                    'protDownsampling/micrographs.xmd',
                    'protDownsampling/micrograph002T.mrc',
                    'protDownsampling/log/protocol.log',
                    'protImport/micrograph002T.mrc'],
              'protPicking':['protDownsampling/micrograph002U.mrc',
                    'protPicking/extra/micrograph002T.pos',
                    'protDownsampling/micrograph003U.mrc',
                    'protPicking/extra/families.xmd',
                    'protDownsampling/micrograph001U.mrc',
                    'protPicking/extra/micrograph001U.pos',
                    'protPicking/extra/micrograph003U.pos',
                    'protDownsampling/micrograph001T.mrc',
                    'protPicking/extra/micrograph001T.pos',
                    'protPicking/extra/micrograph003T.pos',
                    'protDownsampling/micrograph003T.mrc',
                    'protPicking/log/protocol.log',
                    'protDownsampling/micrographs.xmd',
                    'protDownsampling/micrograph002T.mrc',
                    'protPicking/extra/micrograph002U.pos'],
              'protExtract':[
                    'protPicking/extra/micrograph002T.pos',
                    'protExtract/tmp/micrograph001T_noDust.xmp',
                    'protExtract/tmp/micrograph003T_noDust.xmp',
                    'protExtract/tmp/micrograph002T_noDust.xmp',
                    'protExtract/extra/micrograph002T.stk',
                    'protExtract/images.xmd',
                    'protExtract/extra/micrograph003U.stk',
                    'protExtract/extra/micrograph002U.pos',
                    'protExtract/extra/micrograph001T.xmd',
                    'protExtract/tmp/micrograph002U_noDust.xmp',
                    'protExtract/extra/micrograph002T.pos',
                    'protExtract/extra/micrograph002U.xmd',
                    'protPicking/extra/micrograph003T.pos',
                    'protExtract/extra/micrograph001T.pos',
                    'protExtract/extra/micrograph001U.xmd',
                    'protExtract/extra/micrograph001U.pos',
                    'protExtract/extra/micrograph003T.pos',
                    'protExtract/extra/micrograph003U.xmd',
                    'protExtract/extra/micrograph003T.xmd',
                    'protExtract/extra/micrograph003U.pos',
                    'protPicking/extra/micrograph001U.pos',
                    'protExtract/log/protocol.log',
                    'protPicking/extra/micrograph001T.pos',
                    'protExtract/tmp/micrograph003U_noDust.xmp',
                    'protExtract/extra/micrograph001U.stk',
                    'protExtract/extra/micrograph002T.xmd',
                    'protExtract/tmp/micrograph001U_noDust.xmp',
                    'protExtract/extra/micrograph001T.stk',
                    'protExtract/extra/micrograph002U.stk',
                    'protPicking/extra/micrograph002U.pos',
                    'protPicking/extra/micrograph003U.pos',
                    'protExtract/extra/micrograph003T.stk',
                    ],
              'protML2D': [
                   'protML2D/ml2d_extra/iter003/iter_images.xmd',
                    'protML2D/ml2d_extra/iter004/iter_classes.xmd',
                    'protML2D/log/protocol.log',
                    'protExtract/extra/micrograph001U.stk',
                    'protML2D/ml2d_extra/iter004/iter_classes.stk',
                    'protExtract/extra/micrograph002T.stk',
                    'protML2D/ml2d_extra/iter001/iter_classes.stk',
                    'protExtract/images.xmd',
                    'protExtract/extra/micrograph003U.stk',
                    'protML2D/ml2d_extra/iter002/iter_images.xmd',
                    'protML2D/ml2d_extra/iter003/iter_classes.xmd',
                    'protExtract/extra/micrograph002U.stk',
                    'protML2D/ml2d_extra/iter003/iter_classes.stk',
                    'protML2D/ml2d_extra/iter001/iter_images.xmd',
                    'protML2D/ml2d_extra/iter004/iter_images.xmd',
                    'protML2D/ml2d_classes.xmd',
                    'protML2D/ml2d_extra/iter001/iter_classes.xmd',
                    'protML2D/ml2d_extra/iter002/iter_classes.stk',
                    'protML2D/ml2d_extra/iter002/iter_classes.xmd',
                    'protML2D/ml2d__images_average.xmp',
                    'protML2D/ml2d_classes.stk',
                    'protML2D/ml2d_images.xmd',
                    'protExtract/extra/micrograph001T.stk',
                    'protExtract/extra/micrograph003T.stk'],
              'protCL2D': ['protCL2D/extra/classes_core_hierarchy.txt', 
                    'protCL2D/extra/level_01/level_classes_core.xmd', 
                    'protCL2D/extra/level_01/level_classes.stk', 
                    'protCL2D/extra/level_01/classes_sorted.stk', 
                    'protCL2D/extra/level_00/level_classes.xmd', 
                    'protCL2D/extra/level_00/classes_sorted.xmd', 
                    'protCTF/extra/BPV_1386/xmipp_ctf.ctfparam', 
                    'protCTF/extra/BPV_1388/xmipp_ctf.ctfparam', 
                    'protCL2D/extra/level_00/level_classes_core.stk', 
                    'protCL2D/log/protocol.log', 
                    'protExtract/extra/BPV_1388.stk', 
                    'protCTF/extra/BPV_1387/xmipp_ctf.ctfparam', 
                    'protCL2D/extra/level_00/classes_core_sorted.xmd', 
                    'protCL2D/extra/level_01/classes_core_sorted.stk', 
                    'protExtract/images.xmd', 
                    'protCL2D/extra/level_00/classes_sorted.stk', 
                    'protCL2D/extra/level_01/classes_sorted.xmd', 
                    'protCL2D/extra/classes_hierarchy.txt', 
                    'protCL2D/extra/images.xmd', 
                    'protExtract/extra/BPV_1386.stk', 
                    'protCL2D/extra/level_01/classes_core_sorted.xmd', 
                    'protCL2D/extra/level_00/classes_core_sorted.stk', 
                    'protExtract/extra/BPV_1387.stk', 
                    'protCL2D/extra/level_00/level_classes_core.xmd', 
                    'protCL2D/extra/level_01/level_classes.xmd', 
                    'protCL2D/extra/level_00/level_classes.stk', 
                    'protCL2D/extra/level_01/level_classes_core.stk']
              }
    
    @classmethod
    def setUpClass(cls):    
        # Create a new project
        setupProject(cls)
        cls.pattern = getInputPath('Micrographs_TiltedPhantom', '*.mrc')        
        cls.importFolder = getInputPath('Picking_TiltedPhantom')

    def atestXmippTiltedWorkflow(self):
        #First, import a set of micrographs
        protImport = ProtImportMicrographs(pattern=self.pattern, samplingRate=1.237, 
                                           voltage=300, tiltPairs=True)
        self.proj.launchProtocol(protImport, wait=True)
        
        self.printSet('current', protImport.getFiles())
        self.printSet('reference', self.getProtocolFiles('protImport'))
                
        self.assertIsNotNone(protImport.outputMicrographs, "There was a problem with the import")
        self.validateFiles('protImport', protImport)  
                
        # Perform a downsampling on the micrographs
        print "Downsampling..."
        protDownsampling = XmippProtPreprocessMicrographs(doDownsample=True, downFactor=2)
        protDownsampling.inputMicrographs.set(protImport.outputMicrographs)
        self.proj.launchProtocol(protDownsampling, wait=True)
          
        self.assertIsNotNone(protDownsampling.outputMicrographs, "There was a problem with the downsampling")
        self.validateFiles('protDownsampling', protDownsampling)  
                
        print "Running fake particle picking..."   
        protPicking = XmippProtParticlePicking(importFolder=self.importFolder)
                
        protPicking.inputMicrographs.set(protDownsampling.outputMicrographs)
        
        self.proj.launchProtocol(protPicking, wait=True)
            
        self.assertIsNotNone(protPicking.outputCoordinates, "There was a problem with the faked picking")
        self.validateFiles('protPicking', protPicking)  

            
        print "Run extract particles with Same as picking"
        protExtract = XmippProtExtractParticles(boxSize=20, downsampleType=1, doFlip=False)
        protExtract.inputCoordinates.set(protPicking.outputCoordinates)

        self.proj.launchProtocol(protExtract, wait=True)
        
        self.validateFiles('protExtract', protExtract)  
        self.assertIsNotNone(protExtract.outputParticles, "There was a problem with the extract particles")
        
        print "Run ML2D"
        protML2D = XmippProtML2D(numberOfReferences=1, maxIters=4, numberOfMpi=1)
        protML2D.inputImages.set(protExtract.outputParticles)
        self.proj.launchProtocol(protML2D, wait=True)        
        
        self.validateFiles('protML2D', protML2D)  
        self.assertIsNotNone(protML2D.outputClassification, "There was a problem with ML2D")  
        
        print "Run CL2D"
        protCL2D = XmippProtCL2D(numberOfReferences=2, numberOfInitialReferences=1, 
                                 numberOfIterations=4, numberOfMpi=1)
        protCL2D.inputImages.set(protExtract.outputParticles)
        self.proj.launchProtocol(protCL2D, wait=True)        
        
        self.printSet('current', protCL2D.getFiles())
        
        self.validateFiles('protCL2D', protCL2D)  
        self.assertIsNotNone(protCL2D.outputClassification, "There was a problem with CL2D")         

        
if __name__ == "__main__":
#    suite = unittest.TestLoader().loadTestsFromTestCase(TestXmippTiltedWorkflow)    
#    unittest.TextTestRunner(verbosity=2).run(suite)
    unittest.main()