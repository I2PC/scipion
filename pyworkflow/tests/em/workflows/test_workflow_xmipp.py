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
                    'protDownsampling/micrographs.sqlite', 
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
                    'protCTF/extra/BPV_1387/xmipp_ctf_ctfmodel_halfplane.xmp', 
                    'protCTF/extra/BPV_1387/xmipp_ctf_enhanced_psd.xmp', 
                    'protCTF/extra/BPV_1386/xmipp_ctf_ctfmodel_quadrant.xmp', 
                    'protDownsampling/BPV_1388.mrc', 
                    'protDownsampling/BPV_1387.mrc', 
                    'protDownsampling/micrographs.sqlite', 
                    'protCTF/extra/BPV_1386/xmipp_ctf.psd', 
                    'protCTF/extra/BPV_1386/xmipp_ctf_ctfmodel_halfplane.xmp', 
                    'protCTF/micrographs.sqlite',
                    'protCTF/micrographs.xmd',  
                    'protDownsampling/BPV_1386.mrc', 
                    'protCTF/extra/BPV_1388/xmipp_ctf_ctfmodel_quadrant.xmp'],
              'protExtract':[
                    'protImport/BPV_1386.mrc',
                    'protImport/BPV_1388.mrc',
                    'protImport/BPV_1387.mrc',
                    'protImport/micrographs.sqlite',
                    'protPicking/coordinates.sqlite',
                    'protExtract/tmp/BPV_1388_flipped.xmp', 
                    'protExtract/tmp/BPV_1387_flipped.xmp', 
                    'protExtract/tmp/BPV_1386_noDust.xmp', 
                    'protExtract/extra/BPV_1386.xmd', 
                    'protExtract/extra/BPV_1388.stk', 
                    'protExtract/images.xmd', 
                    'protExtract/particles.sqlite',
                    'protExtract/extra/BPV_1386.stk', 
                    'protExtract/extra/BPV_1388.xmd', 
                    'protExtract/extra/BPV_1387.stk', 
                    'protExtract/tmp/BPV_1387_noDust.xmp', 
                    'protExtract/tmp/BPV_1388_noDust.xmp', 
                    'protExtract/extra/BPV_1387.xmd',
                    'protExtract/extra/BPV_1386.pos',
                    'protExtract/extra/BPV_1387.pos',
                    'protExtract/extra/BPV_1388.pos',  
                    'protExtract/tmp/BPV_1386_flipped.xmp',
                    'protExtract/tmp/BPV_1387_downsampled.xmp',
                    'protExtract/tmp/BPV_1386_downsampled.xmp',
                    'protExtract/tmp/BPV_1388_downsampled.xmp',
                    'protExtract/tmp/BPV_1386.ctfParam',
                    'protExtract/tmp/BPV_1387.ctfParam',
                    'protExtract/tmp/BPV_1388.ctfParam',
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
                    'protCL2D/extra/level_01/level_classes_core.stk'],
              'protOnlyAlign': ['protOnlyAlign/extra/images.xmd', 
                    'protOnlyAlign/extra/level_00/class_classes.stk', 
                    'protOnlyAlign/extra/level_00/class_classes.xmd', 
                    'protCTF/extra/BPV_1386/xmipp_ctf.ctfparam', 
                    'protCTF/extra/BPV_1388/xmipp_ctf.ctfparam', 
                    'protOnlyAlign/logs/run.log', 
                    'protOnlyAlign/logs/run.db', 
                    'protExtract/extra/BPV_1388.stk', 
                    'protCTF/extra/BPV_1387/xmipp_ctf.ctfparam', 
                    'protExtract/images.xmd', 
                    'protExtract/extra/BPV_1386.stk', 
                    'protExtract/extra/BPV_1387.stk'],
                'protML3D': ['protML3D/GeneratedReferences/vol001mlf2dextra/iter001/iter_classes.stk',
                    'protML3D/extra/generated_volumes.stk',
                    'protML3D/extra/iter000/vol000002.vol',
                    'protML3D/CorrectGreyscale/vol001/proj_match_weight.doc',
                    'protML3D/CorrectGreyscale/vol001/projections.doc',
                    'protML3D/mlf2dextra/iter002/iter_classes.xmd',
                    'protML3D/GeneratedReferences/vol002mlf2dclasses.xmd',
                    'protML3D/GeneratedReferences/vol002extra/projections.xmd',
                    'protML3D/GeneratedReferences/vol002extra/iter000/vol000001.vol',
                    'protML3D/CorrectGreyscale/vol001/projections.stk',
                    'protML3D/mlf2dextra/iter001/iter_classes.stk',
                    'protML3D/mlf2dextra/iter001/iter_classes.xmd',
                    'protCTF/extra/BPV_1386/xmipp_ctf.ctfparam',
                    'protML3D/mlf2dextra/iter001/iter_images.xmd',
                    'protML3D/classes.xmd',
                    'protML3D/extra/initial_volumes.stk',
                    'protML3D/GeneratedReferences/vol002mlf2dimages.xmd',
                    'protML3D/extra/projections.stk',
                    'protExtract/extra/BPV_1386.stk',
                    'protML3D/extra/iter001/iter_volumes.xmd',
                    'protML3D/GeneratedReferences/vol001extra/projections.xmd',
                    'protML3D/extra/iter002/vol000001.vol',
                    'protML3D/GeneratedReferences/vol001extra/iter001/vol000001.vol',
                    'protML3D/GeneratedReferences/vol002extra/iter001/vol000001.vol',
                    'protML3D/extra/iter002/vol000001.projections.xmd',
                    'protCTF/extra/BPV_1388/xmipp_ctf.ctfparam',
                    'protML3D/CorrectGreyscale/vol001/projections_angles.doc',
                    'protML3D/GeneratedReferences/vol001extra/iter001/vol000001.projections.xmd',
                    'protML3D/extra/iter002/vol000002.projections.xmd',
                    'protML3D/logs/run.log',
                    'protML3D/GeneratedReferences/vol002extra/projections.stk',
                    'protML3D/GeneratedReferences/vol001extra/iter001/vol',
                    'protML3D/GeneratedReferences/vol002mlf2dclasses.stk',
                    'protML3D/GeneratedReferences/vol002extra/iter001/vol000001.projections.xmd',
                    'protML3D/extra/filtered_volumes.stk',
                    'protML3D/extra/iter001/vol000002.vol',
                    'protML3D/GeneratedReferences/vol002mlf2dextra/iter001/iter_images.xmd',
                    'protML3D/extra/iter002/vol',
                    'protML3D/GeneratedReferences/images000002.xmd',
                    'protML3D/GeneratedReferences/vol002extra/iter001/iter_volumes.xmd',
                    'protExtract/extra/BPV_1388.stk',
                    'protML3D/extra/iter001/vol000002.projections.xmd',
                    'protML3D/extra/iter002/iter_volumes.xmd',
                    'protML3D/GeneratedReferences/vol002mlf2dextra/iter001/iter_classes.stk',
                    'protML3D/CorrectGreyscale/vol001/corrected_refs_Ref3D_001.xmd',
                    'protML3D/GeneratedReferences/vol001mlf2dimages.xmd',
                    'protCTF/extra/BPV_1387/xmipp_ctf.ctfparam',
                    'protML3D/CorrectGreyscale/vol001/corrected_refs_Ref3D_001.stk',
                    'protML3D/GeneratedReferences/vol001mlf2dextra/iter001/iter_images.xmd',
                    'protML3D/logs/run.db',
                    'protML3D/classes.stk',
                    'protML3D/GeneratedReferences/vol001extra/projections.stk',
                    'protML3D/GeneratedReferences/vol001extra/iter000/vol000001.vol',
                    'protExtract/images.xmd',
                    'protML3D/GeneratedReferences/vol001mlf2dextra/iter001/iter_classes.xmd',
                    'protML3D/CorrectGreyscale/vol001/projections_sampling.xmd',
                    'protML3D/CorrectGreyscale/vol001/corrected_refs_discarded.xmd',
                    'protML3D/extra/corrected_volumes.stk',
                    'protML3D/mlf2dextra/iter002/iter_classes.stk',
                    'protML3D/extra/iter002/vol000002.vol',
                    'protML3D/GeneratedReferences/images000001.xmd',
                    'protML3D/GeneratedReferences/vol001mlf2dclasses.xmd',
                    'protML3D/GeneratedReferences/vol001extra/iter001/iter_volumes.xmd',
                    'protML3D/extra/projections.xmd',
                    'protExtract/extra/BPV_1387.stk',
                    'protML3D/GeneratedReferences/vol002extra/iter001/vol',
                    'protML3D/GeneratedReferences/vol002mlf2dextra/iter001/iter_classes.xmd',
                    'protML3D/mlf2dextra/iter002/iter_images.xmd',
                    'protML3D/GeneratedReferences/vol002extra/noise_images.xmd',
                    'protML3D/GeneratedReferences/vol001extra/noise_vol000001.vol',
                    'protML3D/GeneratedReferences/vol002extra/noise_vol000001.projections.xmd',
                    'protML3D/GeneratedReferences/vol001extra/cref_vol',
                    'protML3D/GeneratedReferences/vol002extra/noise_vol',
                    'protML3D/GeneratedReferences/vol002extra/cref_vol',
                    'protML3D/extra/iter002/iter_3dssnr.log',
                    'protML3D/extra/iter001/vol',
                    'protML3D/extra/iter001/vol000001.vol',
                    'protML3D/extra/iter000/vol000001.vol',
                    'protML3D/extra/iter001/vol000001.projections.xmd',
                    'protML3D/GeneratedReferences/vol001mlf2dclasses.stk',
                    'protML3D/CorrectGreyscale/vol001/proj_match.doc',
                    'protML3D/GeneratedReferences/vol002mlf2dextra/iter000/iter_noise.xmd',
                    'protML3D/mlf2dextra/cref_classes.xmd',
                    'protML3D/mlf2dextra/iter002/iter_noise.xmd',
                    'protML3D/GeneratedReferences/vol001extra/noise_vol000001.projections.xmd',
                    'protML3D/mlf2dextra/iter000/iter_ssnr.xmd',
                    'protML3D/GeneratedReferences/vol001mlf2dextra/cref_classes.xmd',
                    'protML3D/extra/cref_vol000001.projections.xmd',
                    'protML3D/mlf2dextra/cref_classes.stk',
                    'protML3D/GeneratedReferences/vol002extra/noise_images.stk',
                    'protML3D/extra/cref_vol000002.projections.xmd',
                    'protML3D/GeneratedReferences/vol001mlf2dextra/iter001/iter_noise.xmd',
                    'protML3D/mlf2dextra/iter001/iter_noise.xmd',
                    'protML3D/images.xmd',
                    'protML3D/GeneratedReferences/vol001extra/noise_images.xmd',
                    'protML3D/extra/noise_vol000002.vol',
                    'protML3D/GeneratedReferences/vol001mlf2dextra/iter001/iter_ssnr.xmd',
                    'protML3D/extra/noise_images.xmd',
                    'protML3D/mlf2dextra/iter001/iter_ssnr.xmd',
                    'protML3D/GeneratedReferences/vol001mlf2dnoise.xmd',
                    'protML3D/GeneratedReferences/vol002extra/cref_vol000001.vol',
                    'protML3D/GeneratedReferences/vol002mlf2dextra/cref_classes.stk',
                    'protML3D/extra/noise_vol000001.vol',
                    'protML3D/GeneratedReferences/vol002mlf2dnoise.xmd',
                    'protML3D/extra/cref_vol000002.vol',
                    'protML3D/GeneratedReferences/vol002extra/cref_vol000001.projections.xmd',
                    'protML3D/GeneratedReferences/vol002extra/iter001/iter_3dssnr.log',
                    'protML3D/GeneratedReferences/vol001extra/cref_vol000001.projections.xmd',
                    'protML3D/GeneratedReferences/vol002mlf2dextra/iter000/iter_ssnr.xmd',
                    'protML3D/extra/noise_vol000001.projections.xmd',
                    'protML3D/extra/noise_vol',
                    'protML3D/GeneratedReferences/vol002extra/noise_vol000001.vol',
                    'protML3D/GeneratedReferences/vol001extra/noise_vol',
                    'protML3D/GeneratedReferences/vol002mlf2dextra/iter001/iter_ssnr.xmd',
                    'protML3D/mlf2dnoise.xmd',
                    'protML3D/GeneratedReferences/vol001mlf2dextra/cref_classes.stk',
                    'protML3D/GeneratedReferences/vol002mlf2dextra/cref_classes.xmd',
                    'protML3D/GeneratedReferences/vol001extra/iter001/iter_3dssnr.log',
                    'protML3D/extra/iter001/iter_3dssnr.log',
                    'protML3D/mlf2dextra/iter002/iter_ssnr.xmd',
                    'protML3D/mlf2dextra/iter000/iter_noise.xmd',
                    'protML3D/GeneratedReferences/vol001mlf2dextra/iter000/iter_noise.xmd',
                    'protML3D/extra/cref_vol',
                    'protML3D/GeneratedReferences/vol001extra/noise_images.stk',
                    'protML3D/GeneratedReferences/vol001mlf2dextra/iter000/iter_ssnr.xmd',
                    'protML3D/GeneratedReferences/vol002mlf2dextra/iter001/iter_noise.xmd',
                    'protML3D/extra/cref_vol000001.vol',
                    'protML3D/extra/noise_images.stk',
                    'protML3D/GeneratedReferences/vol001extra/cref_vol000001.vol',
                    'protML3D/extra/noise_vol000002.projections.xmd']            
              }

    @classmethod
    def setUpClass(cls):    
        # Create a new project
        setupProject(cls)
        cls.pattern = getInputPath('Micrographs_BPV3', '*.mrc')        
        cls.importFolder = getInputPath('Picking_XmippBPV3_Down3')
            
    def testXmippWorkflow(self):
        #First, import a set of micrographs
        protImport = ProtImportMicrographs(pattern=self.pattern, samplingRate=1.237, voltage=300)
        self.proj.launchProtocol(protImport, wait=True)
        
        self.assertIsNotNone(protImport.outputMicrographs.getFileName(), "There was a problem with the import")
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
        protExtract = XmippProtExtractParticles(boxSize=64, downsampleType=2, downFactor=8, runMode=1, doInvert=True)
        protExtract.inputCoordinates.set(protPP.outputCoordinates)
        protExtract.inputMicrographs.set(protImport.outputMicrographs)
        self.proj.launchProtocol(protExtract, wait=True)
        
        self.assertIsNotNone(protExtract.outputParticles, "There was a problem with the extract particles")
        self.validateFiles('protExtract', protExtract)
        
#        print "Run ML2D"
#        protML2D = XmippProtML2D(numberOfReferences=1, maxIters=4, doMlf=False,#True,
#                                 numberOfMpi=2, numberOfThreads=1)
#        protML2D.inputImages.set(protExtract.outputParticles)
#        self.proj.launchProtocol(protML2D, wait=True)        
#        
#        self.assertIsNotNone(protML2D.outputClasses, "There was a problem with ML2D") 
        # Check that images related to each class have ctf model
#        for class2D in protML2D.outputClasses:
#            for imgCA in class2D:
#                xmippImg = imgCA.getImage()
#                self.assertTrue(imgCA.getImage().hasCTF(), "Image class has not CTF information.")
             
#        self.validateFiles('protML2D', protML2D)
        
        print "Run CL2D"
        protCL2D = XmippProtCL2D(numberOfReferences=2, numberOfInitialReferences=1, 
                                 numberOfIterations=4, numberOfMpi=2)
        protCL2D.inputImages.set(protExtract.outputParticles)
        self.proj.launchProtocol(protCL2D, wait=True)        
        
        self.assertIsNotNone(protCL2D.outputClasses, "There was a problem with CL2D")
        # Check that images related to each class have ctf model
#        for class2D in protCL2D.outputClasses:
#            for imgCA in class2D:
#                xmippImg = imgCA.getImage()
#                self.assertTrue(imgCA.getImage().hasCTF(), "Image class has not CTF information.")
        self.validateFiles('protCL2D', protCL2D) 

        print "Run Only Align2d"
        protOnlyAlign = XmippProtCL2DAlign(maximumShift=5, numberOfIterations=2, 
                                 numberOfMpi=2, numberOfThreads=1, useReferenceImage=False)

        protOnlyAlign.inputImages.set(protExtract.outputParticles)
        self.proj.launchProtocol(protOnlyAlign, wait=True)        
        
        self.assertIsNotNone(protOnlyAlign.outputParticles, "There was a problem with Only align2d")  
        self.validateFiles('protOnlyAlign', protOnlyAlign)

        print "Run kerdensom"
        ProtKerdensom = XmippProtKerdensom(useMask=False, SomXdim=2, SomYdim=2,
                                 SomReg0=800, SomReg1=400, SomSteps=2)

        ProtKerdensom.inputImages.set(protOnlyAlign.outputParticles)
        self.proj.launchProtocol(ProtKerdensom, wait=True)        
        
        self.assertIsNotNone(ProtKerdensom.outputClasses, "There was a problem with kerdensom")  
        #self.validateFiles('ProtKerdensom', ProtKerdensom)
        
        print "Run Rotational Spectra"
        xmippProtRotSpectra = XmippProtRotSpectra(SomXdim=2, SomYdim=2)
        xmippProtRotSpectra.inputImages.set(protOnlyAlign.outputParticles)
        self.proj.launchProtocol(xmippProtRotSpectra, wait=True)        
        
        self.assertIsNotNone(xmippProtRotSpectra.outputClasses, "There was a problem with Rotational Spectra")

        # The ML3D test is taking too long now
        # skipping until revision
        #return 
        
#        print "ML3D"
#        protML3D = XmippProtML3D(angularSampling=15, numberOfIterations=2, runMode=1, numberOfMpi=2, numberOfThreads=2)
#        protML3D.inputImages.set(protExtract.outputParticles)
#        protML3D.ini3DrefVolumes.set(getInputPath('Volumes_BPV', 'BPV_scale_filtered_windowed_64.vol'))
#        protML3D.doCorrectGreyScale.set(True)
#        protML3D.doMlf.set(True)
#        protML3D.numberOfSeedsPerRef.set(2)
#
#        self.proj.launchProtocol(protML3D, wait=True)        
#        
#        self.assertIsNotNone(protML3D.outputVolumes, "There was a problem with ML3D")
#        self.validateFiles('protML3D', protML3D)

if __name__ == "__main__":
    unittest.main()
