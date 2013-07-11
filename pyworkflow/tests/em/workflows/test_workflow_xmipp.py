import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3 import *
from test_workflow import TestWorkflow
   
       
class TestXmippWorkflow(TestWorkflow):
    
    GOLD_FILES = {'protImport': ['protImport/BPV_1388.mrc',
                    'protImport/micrographs.sqlite', 
                    'protImport/BPV_1387.mrc',
                    'protImport/BPV_1386.mrc'],
              'protDownsampling': [
                    'protImport/BPV_1386.mrc', 
                    'protImport/BPV_1388.mrc', 
                    'protImport/BPV_1387.mrc',
                    'protImport/micrographs.sqlite', 
                    'protDownsampling/BPV_1386.mrc', 
                    'protDownsampling/BPV_1387.mrc', 
                    'protDownsampling/BPV_1388.mrc', 
                    'protDownsampling/micrographs.xmd', 
                    'protDownsampling/logs/run.log',
                    'protDownsampling/logs/run.db', 
                    ],
              'protCTF': ['protCTF/extra/BPV_1387/xmipp_ctf_ctfmodel_quadrant.xmp', 
                    'protCTF/logs/run.log', 
                    'protCTF/logs/run.db', 
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
                    'protCTF/extra/BPV_1388/xmipp_ctf_ctfmodel_quadrant.xmp'],
              'protExtract':[
                    'protImport/BPV_1386.mrc',
                    'protImport/BPV_1388.mrc',
                    'protImport/BPV_1387.mrc',
                    'protImport/micrographs.sqlite',
                    'protPicking/extra/BPV_1387.pos', 
                    'protPicking/extra/BPV_1386.pos', 
                    'protPicking/extra/BPV_1388.pos', 
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
                    'protExtract/tmp/BPV_1387_downsampled.xmp',
                    'protExtract/tmp/BPV_1386_downsampled.xmp',
                    'protExtract/tmp/BPV_1388_downsampled.xmp',
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
                    'protML2D/mlf2d_extra/iter002/iter_classes.xmd', 
                    'protML2D/mlf2d_extra/iter001/iter_classes.xmd', 
                    'protML2D/mlf2d_extra/iter002/iter_images.xmd', 
                    'protML2D/mlf2d_classes.stk', 
                    'protML2D/mlf2d_extra/iter003/iter_images.xmd', 
                    'protML2D/mlf2d_images.xmd', 
                    'protML2D/mlf2d_extra/iter004/iter_classes.stk', 
                    'protML2D/mlf2d_extra/iter001/iter_images.xmd', 
                    'protML2D/mlf2d_extra/iter002/iter_classes.stk', 
                    'protML2D/mlf2d_extra/iter004/iter_classes.xmd', 
                    'protML2D/mlf2d_extra/iter004/iter_images.xmd', 
                    'protML2D/mlf2d_extra/iter003/iter_classes.xmd', 
                    'protML2D/mlf2d_extra/iter001/iter_noise.xmd',
                    'protML2D/mlf2d_extra/iter002/iter_noise.xmd',
                    'protML2D/mlf2d_extra/iter002/iter_ssnr.xmd',
                    'protML2D/mlf2d_extra/iter000/iter_classes.xmd',
                    'protML2D/mlf2d_extra/iter001/iter_ssnr.xmd',
                    'protML2D/mlf2d_extra/iter003/iter_noise.xmd',
                    'protML2D/mlf2d_extra/iter000/iter_noise.xmd',
                    'protML2D/mlf2d_extra/iter004/iter_ssnr.xmd',
                    'protML2D/mlf2d_noise.xmd',
                    'protML2D/mlf2d_extra/iter000/iter_classes.stk',
                    'protML2D/mlf2d_extra/iter003/iter_ssnr.xmd',
                    'protML2D/mlf2d_extra/cref_classes.stk',
                    'protML2D/mlf2d_extra/iter000/iter_ssnr.xmd',
                    'protML2D/mlf2d_extra/cref_classes.xmd',
                    'protML2D/mlf2d_extra/iter004/iter_noise.xmd',
                    'protML2D/logs/run.log', 
                    'protML2D/logs/run.db',
                    'protML2D/mlf2d_classes.xmd', 
                    'protML2D/mlf2d_extra/iter001/iter_classes.stk', 
                    'protML2D/mlf2d_extra/iter003/iter_classes.stk'],
              'protCL2D': ['protCL2D/extra/classes_core_hierarchy.txt', 
                    'protCL2D/extra/level_01/level_classes_core.xmd', 
                    'protCL2D/extra/level_01/level_classes.stk', 
                    'protCL2D/extra/level_01/classes_sorted.stk', 
                    'protCL2D/extra/level_00/level_classes.xmd', 
                    'protCL2D/extra/level_00/classes_sorted.xmd', 
                    'protCTF/extra/BPV_1386/xmipp_ctf.ctfparam', 
                    'protCTF/extra/BPV_1388/xmipp_ctf.ctfparam', 
                    'protCL2D/extra/level_00/level_classes_core.stk', 
                    'protCL2D/logs/run.log', 
                    'protCL2D/logs/run.db', 
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
        cls.pattern = getInputPath('Micrographs_BPV3', '*.mrc')        
        cls.importFolder = getInputPath('Picking_XmippBPV3')
            
    def testXmippWorkflow(self):
        #First, import a set of micrographs
        protImport = ProtImportMicrographs(pattern=self.pattern, samplingRate=1.237, voltage=300)
        self.proj.launchProtocol(protImport, wait=True)
        
        self.assertIsNotNone(protImport.outputMicrographs, "There was a problem with the import")
        self.validateFiles('protImport', protImport)        
        
        # Perform a downsampling on the micrographs

        print "Downsampling..."
        protDownsampling = XmippProtPreprocessMicrographs(doDownsample=True, downFactor=3, doCrop=False, runMode=1)
        protDownsampling.inputMicrographs.set(protImport.outputMicrographs)
        self.proj.launchProtocol(protDownsampling, wait=True)
          
        self.assertIsNotNone(protDownsampling.outputMicrographs, "There was a problem with the downsampling")
        self.validateFiles('protDownsampling', protDownsampling)
          
        # Now estimate CTF on the downsampled micrographs 
        print "Performing CTF..."   
        protCTF = XmippProtCTFMicrographs(numberOfThreads=3)                
        protCTF.inputMicrographs.set(protDownsampling.outputMicrographs)        
        self.proj.launchProtocol(protCTF, wait=True)
        
        # After CTF estimation, the output micrograph should have CTF info
        self.assertTrue(protCTF.outputMicrographs.hasCTF())
        self.validateFiles('protCTF', protCTF)
        
        print "Running fake particle picking..."   
        protPP = XmippProtParticlePicking(importFolder=self.importFolder)                
        protPP.inputMicrographs.set(protCTF.outputMicrographs)        
        self.proj.launchProtocol(protPP, wait=True)
        self.protDict['protPicking'] = protPP
            
        self.assertIsNotNone(protPP.outputCoordinates, "There was a problem with the faked picking")
            
        print "Run extract particles with other downsampling factor"
        protExtract = XmippProtExtractParticles(boxSize=64, downsampleType=2, downFactor=8, runMode=1)
        protExtract.inputCoordinates.set(protPP.outputCoordinates)
        protExtract.inputMicrographs.set(protImport.outputMicrographs)
        self.proj.launchProtocol(protExtract, wait=True)
        
        self.assertIsNotNone(protExtract.outputImages, "There was a problem with the extract particles")
        self.validateFiles('protExtract', protExtract)
        
        print "Run ML2D"
        protML2D = XmippProtML2D(numberOfReferences=1, maxIters=4, doMlf=True,
                                 numberOfMpi=2, numberOfThreads=1)
        protML2D.inputImages.set(protExtract.outputImages)
        self.proj.launchProtocol(protML2D, wait=True)        
        
        self.assertIsNotNone(protML2D.outputClassification, "There was a problem with ML2D") 
        # Check that images related to each class have ctf model
        for class2D in protML2D.outputClassification:
            for imgCA in class2D:
                xmippImg = imgCA.getImage()
                self.assertTrue(imgCA.getImage().hasCTF(), "Image class has not CTF information.")
             
        self.validateFiles('protML2D', protML2D)
        
        print "Run CL2D"
        protCL2D = XmippProtCL2D(numberOfReferences=2, numberOfInitialReferences=1, 
                                 numberOfIterations=4, numberOfMpi=2)
        protCL2D.inputImages.set(protExtract.outputImages)
        self.proj.launchProtocol(protCL2D, wait=True)        
        
        self.assertIsNotNone(protCL2D.outputClassification, "There was a problem with CL2D")
        # Check that images related to each class have ctf model
        for class2D in protCL2D.outputClassification:
            for imgCA in class2D:
                xmippImg = imgCA.getImage()
                self.assertTrue(imgCA.getImage().hasCTF(), "Image class has not CTF information.")
        self.validateFiles('protCL2D', protCL2D) 

        print "Run Only Align2d"
        protOnlyalign = XmippProtCL2DAlign(maximumShift=5, numberOfIterations=5, 
                                 numberOfMpi=2, numberOfThreads=1, useReferenceImage=False)

        protOnlyalign.inputImages.set(protExtract.outputImages)
        self.proj.launchProtocol(protOnlyalign, wait=True)        
        
        self.assertIsNotNone(protOnlyalign.outputClassification, "There was a problem with Only align2d")  
        self.validateFiles('protOnlyalign', protOnlyalign)

        print "Run kerdensom"
        XmippProtKerdensom = XmippProtKerdensom()

        protOnlyalign.inputImages.set(protExtract.outputImages)
        self.proj.launchProtocol(XmippProtKerdensom, wait=True)        
        
        self.assertIsNotNone(XmippProtKerdensom.outputClassification, "There was a problem with kerdensom")  
        self.validateFiles('XmippProtKerdensom', XmippProtKerdensom)

if __name__ == "__main__":
    unittest.main()
