import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3 import *
from pyworkflow.em.packages.brandeis import *
from pyworkflow.em.packages.eman2 import *
from test_workflow import TestWorkflow
    
    
class TestMixedWorkflow_1(TestWorkflow):

    GOLD_FILES = {'protImport': [
                    'protImport/BPV_1386.mrc',
                    'protImport/micrographs.sqlite'],
              'protDownsampling': [
                    'protDownsampling/BPV_1386.mrc', 
                    'protImport/BPV_1386.mrc',
                    'protDownsampling/micrographs.xmd', 
                    'protImport/micrographs.sqlite', 
                    'protDownsampling/logs/run.log',
                    'protDownsampling/logs/run.db'
                    ],
              'protCTF': [
                    'protCTF/extra/BPV_1386/ctffind_psd.mrc', 
                    'protCTF/extra/BPV_1386/ctffind.out', 
                    'protCTF/micrographs.sqlite',
                    'protDownsampling/micrographs.xmd', 
                    'protDownsampling/BPV_1386.mrc',
                    'protCTF/logs/run.log', 
                    'protCTF/logs/run.db'],
              'protExtract':[
                    'protPicking/extra/BPV_1386.pos', 
                    'protExtract/tmp/BPV_1386_noDust.xmp', 
                    'protExtract/extra/BPV_1386.xmd', 
                    'protExtract/images.xmd', 
                    'protExtract/extra/BPV_1386.stk', 
                    'protExtract/extra/BPV_1386.pos', 
                    'protExtract/tmp/BPV_1386_flipped.xmp',
                    'protExtract/logs/run.log',
                    'protExtract/logs/run.db'],
              }
    
    @classmethod
    def setUpClass(cls):    
        # Create a new project
        setupProject(cls)
        cls.pattern = getInputPath('Micrographs_BPV1', '*.mrc')        
        cls.importFolder = getInputPath('Picking_XmippBPV1')
        
    def testWorkflow(self):
        #First, import a set of micrographs
        protImport = ProtImportMicrographs(pattern=self.pattern, samplingRate=1.237, voltage=300)
        self.proj.launchProtocol(protImport, wait=True)
        
        self.assertIsNotNone(protImport.outputMicrographs, "There was a problem with the import")
        self.validateFiles('protImport', protImport) 
        
        # Perform a downsampling on the micrographs

        print "Downsampling..."
        protDownsampling = XmippProtPreprocessMicrographs(doDownsample=True, downFactor=3, doCrop=False)
        protDownsampling.inputMicrographs.set(protImport.outputMicrographs)
        self.proj.launchProtocol(protDownsampling, wait=True)
          
        self.assertIsNotNone(protDownsampling.outputMicrographs, "There was a problem with the downsampling")
        self.validateFiles('protDownsampling', protDownsampling) 


        # Now estimate CTF on the downsampled micrographs 
        print "Performing CTFfind..."   
        protCTF = ProtCTFFind(runMode=1)         
        protCTF.inputMicrographs.set(protDownsampling.outputMicrographs)        
        self.proj.launchProtocol(protCTF, wait=True)
        
        self.validateFiles('protCTF', protCTF) 
        
        print "Running fake particle picking..."   
        protPicking = XmippProtParticlePicking(importFolder=self.importFolder, runMode=1)                
        protPicking.inputMicrographs.set(protCTF.outputMicrographs)        
        self.proj.launchProtocol(protPicking, wait=True)
        self.protDict['protPicking'] = protPicking
            
        self.assertIsNotNone(protPicking.outputCoordinates, "There was a problem with the faked picking")
            
        print "Run extract particles with Same as picking"
        protExtract = XmippProtExtractParticles(boxSize=171, downsampleType=1, runMode=1)
        protExtract.inputCoordinates.set(protPicking.outputCoordinates)
        #protExtract.inputMicrographs.set(protDownsampling.outputMicrographs)
        self.proj.launchProtocol(protExtract, wait=True)
        
        self.assertIsNotNone(protExtract.outputImages, "There was a problem with the extract particles")
        self.validateFiles('protExtract', protExtract)
        
class TestMixedWorkflow_2(TestWorkflow):

    GOLD_FILES = {'protImport': [
                    'protImport/BPV_1388.mrc',
                    'protImport/micrographs.sqlite', 
                    'protImport/BPV_1387.mrc',
                    'protImport/BPV_1386.mrc'],
              'protDownsampling': ['protDownsampling/BPV_1388.mrc', 
                    'protDownsampling/BPV_1387.mrc', 
                    'protImport/BPV_1386.mrc', 
                    'protDownsampling/micrographs.xmd', 
                    'protImport/BPV_1388.mrc', 
                    'protImport/micrographs.sqlite', 
                    'protDownsampling/BPV_1386.mrc', 
                    'protImport/BPV_1387.mrc',
                    'protDownsampling/logs/run.log',
                    'protDownsampling/logs/run.db'],
              'protCTF': [
                    'protCTF/extra/BPV_1387/xmipp_ctf_ctfmodel_quadrant.xmp', 
                    'protCTF/extra/BPV_1387/xmipp_ctf.psd', 
                    'protCTF/extra/BPV_1386/xmipp_ctf_enhanced_psd.xmp', 
                    'protCTF/extra/BPV_1386/xmipp_ctf.ctfparam', 
                    'protCTF/extra/BPV_1388/xmipp_ctf.ctfparam', 
                    'protCTF/extra/BPV_1388/xmipp_ctf.psd', 
                    'protCTF/extra/BPV_1388/xmipp_ctf_ctfmodel_halfplane.xmp', 
                    'protCTF/extra/BPV_1387/xmipp_ctf.ctfparam', 
                    'protCTF/extra/BPV_1388/xmipp_ctf_enhanced_psd.xmp', 
                    'protCTF/tmp/micrographs.xmd', 
                    'protCTF/extra/BPV_1387/xmipp_ctf_ctfmodel_halfplane.xmp', 
                    'protCTF/extra/BPV_1387/xmipp_ctf_enhanced_psd.xmp', 
                    'protCTF/extra/BPV_1386/xmipp_ctf_ctfmodel_quadrant.xmp', 
                    'protDownsampling/BPV_1388.mrc', 
                    'protDownsampling/BPV_1387.mrc', 
                    'protDownsampling/micrographs.xmd', 
                    'protCTF/extra/BPV_1386/xmipp_ctf.psd', 
                    'protCTF/extra/BPV_1386/xmipp_ctf_ctfmodel_halfplane.xmp', 
                    'protCTF/micrographs.xmd', 
                    'protDownsampling/BPV_1386.mrc', 
                    'protCTF/extra/BPV_1388/xmipp_ctf_ctfmodel_quadrant.xmp',
                    'protCTF/logs/run.log', 
                    'protCTF/logs/run.db'],
              'protExtract':[
                    'protPicking/BPV_1387.box', 
                    'protPicking/BPV_1386.box', 
                    'protPicking/BPV_1388.box', 
                    'protExtract/extra/BPV_1387.pos', 
                    'protExtract/tmp/BPV_1388_flipped.xmp', 
                    'protExtract/tmp/BPV_1387_flipped.xmp', 
                    'protExtract/tmp/BPV_1386_noDust.xmp', 
                    'protExtract/extra/BPV_1388.pos', 
                    'protExtract/extra/BPV_1386.xmd', 
                    'protExtract/extra/BPV_1388.stk', 
                    'protExtract/images.xmd', 
                    'protExtract/extra/BPV_1386.stk', 
                    'protExtract/extra/BPV_1388.xmd', 
                    'protExtract/extra/BPV_1386.pos', 
                    'protExtract/extra/BPV_1387.stk', 
                    'protExtract/tmp/BPV_1387_noDust.xmp', 
                    'protExtract/tmp/BPV_1388_noDust.xmp', 
                    'protExtract/extra/BPV_1387.xmd', 
                    'protExtract/tmp/BPV_1386_flipped.xmp',
                    'protExtract/logs/run.log',
                    'protExtract/logs/run.db',
                    ],
                'protML2D': [
                    'protExtract/extra/BPV_1386.stk', 
                    'protExtract/extra/BPV_1387.stk', 
                    'protExtract/extra/BPV_1388.stk', 
                    'protExtract/images.xmd', 
                    'protCTF/extra/BPV_1386/xmipp_ctf.ctfparam',
                    'protCTF/extra/BPV_1387/xmipp_ctf.ctfparam',
                    'protCTF/extra/BPV_1388/xmipp_ctf.ctfparam',
                    'protML2D/ml2d_extra/iter002/iter_classes.xmd', 
                    'protML2D/ml2d_extra/iter001/iter_classes.xmd', 
                    'protML2D/ml2d_extra/iter002/iter_images.xmd', 
                    'protML2D/ml2d_classes.stk', 
                    'protML2D/ml2d_extra/iter003/iter_images.xmd', 
                    'protML2D/ml2d_images.xmd', 
                    'protML2D/ml2d_extra/iter004/iter_classes.stk', 
                    'protML2D/ml2d__images_average.xmp', 
                    'protML2D/ml2d_extra/iter001/iter_images.xmd', 
                    'protML2D/ml2d_extra/iter002/iter_classes.stk', 
                    'protML2D/ml2d_extra/iter004/iter_classes.xmd', 
                    'protML2D/ml2d_extra/iter004/iter_images.xmd', 
                    'protML2D/ml2d_extra/iter003/iter_classes.xmd',  
                    'protML2D/logs/run.log',
                    'protML2D/logs/run.db',
                    'protML2D/ml2d_classes.xmd', 
                    'protML2D/ml2d_extra/iter001/iter_classes.stk', 
                    'protML2D/ml2d_extra/iter003/iter_classes.stk'],
              }
        
    @classmethod
    def setUpClass(cls):    
        # Create a new project
        setupProject(cls)
        cls.pattern = getInputPath('Micrographs_BPV3', '*.mrc')        
        cls.importFolder = getInputPath('EmanTestProject')
        
    def testWorkflow(self):
        #First, import a set of micrographs
        protImport = ProtImportMicrographs(pattern=self.pattern, samplingRate=1.237, voltage=300)
        self.proj.launchProtocol(protImport, wait=True)
        
        self.assertIsNotNone(protImport.outputMicrographs, "There was a problem with the import")
        self.validateFiles('protImport', protImport)
        
        # Perform a downsampling on the micrographs
        print "Downsampling..."
        protDownsampling = XmippProtPreprocessMicrographs(doDownsample=True, downFactor=5, doCrop=False)
        protDownsampling.inputMicrographs.set(protImport.outputMicrographs)
        self.proj.launchProtocol(protDownsampling, wait=True)
          
        self.assertIsNotNone(protDownsampling.outputMicrographs, "There was a problem with the downsampling")
        self.validateFiles('protDownsampling', protDownsampling) 

        # Now estimate CTF on the downsampled micrographs 
        print "Performing CTF estimation..."   
        protCTF = XmippProtCTFMicrographs(numberOfThreads=3, runMode=1)         
        protCTF.inputMicrographs.set(protDownsampling.outputMicrographs)        
        self.proj.launchProtocol(protCTF, wait=True)
        
        self.validateFiles('protCTF', protCTF) 
        
        print "Running Eman fake particle picking..."   
        protPP = EmanProtBoxing(importFolder=self.importFolder, runMode=1)                
#        protPP.inputMicrographs.set(protCTF.outputMicrographs)        
        protPP.inputMicrographs.set(protCTF.outputMicrographs)
        self.proj.launchProtocol(protPP, wait=True)
        self.protDict['protPicking'] = protPP
            
        self.assertIsNotNone(protPP.outputCoordinates, "There was a problem with the faked picking")
            
        print "Run extract particles with Same as picking"
        protExtract = XmippProtExtractParticles(boxSize=110, downsampleType=1, doFlip=True, runMode=1)
        protExtract.inputCoordinates.set(protPP.outputCoordinates)
        #protExtract.inputMicrographs.set(protDownsampling.outputMicrographs)
        self.proj.launchProtocol(protExtract, wait=True)
        
        self.assertIsNotNone(protExtract.outputImages, "There was a problem with the extract particles")
        self.validateFiles('protExtract', protExtract)
        
        print "Run ML2D"
        protML2D = XmippProtML2D(numberOfReferences=1, maxIters=4, 
                                 numberOfMpi=2, numberOfThreads=1)
#        protML2D.inputImages.set(protExtract.outputImages)
        protML2D.inputImages.set(protExtract.outputImages)
        self.proj.launchProtocol(protML2D, wait=True)        
        
        self.assertIsNotNone(protML2D.outputClassification, "There was a problem with ML2D")  
        self.validateFiles('protML2D', protML2D)
               
               
        
if __name__ == "__main__":
    unittest.main()