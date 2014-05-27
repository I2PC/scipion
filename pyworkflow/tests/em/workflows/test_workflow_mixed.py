import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3 import *
from pyworkflow.em.packages.brandeis import *
from pyworkflow.em.packages.eman2 import *
from pyworkflow.em.packages.relion import *
from test_workflow import TestWorkflow


class TestMixedBPV(TestWorkflow):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('xmipp_tutorial')
        cls.crdsDir = cls.dataset.getFile('boxingDir')
        cls.micsFn = cls.dataset.getFile('allMics')
        cls.vol1 = cls.dataset.getFile('vol1')
        
    def test_workflow(self):
        #First, import a set of micrographs
        protImport = ProtImportMicrographs(pattern=self.micsFn, samplingRate=1.237, voltage=300)
        self.proj.launchProtocol(protImport, wait=True)
        self.assertIsNotNone(protImport.outputMicrographs, "There was a problem with the import")
#         self.validateFiles('protImport', protImport) 
        
        #Import a set of volumes        
        print "Import Volume"
        protImportVol = ProtImportVolumes(pattern=self.vol1, samplingRate=9.896)
        self.proj.launchProtocol(protImportVol, wait=True)
        self.assertIsNotNone(protImportVol.getFiles(), "There was a problem with the import")
#        self.validateFiles('protImportVol', protImportVol)
        
        # Perform a downsampling on the micrographs
        print "Downsampling..."
        protDownsampling = XmippProtPreprocessMicrographs(doDownsample=True, downFactor=5, doCrop=False, runMode=1)
        protDownsampling.inputMicrographs.set(protImport.outputMicrographs)
        self.proj.launchProtocol(protDownsampling, wait=True)
        self.assertIsNotNone(protDownsampling.outputMicrographs, "There was a problem with the downsampling")
#         self.validateFiles('protDownsampling', protDownsampling)
        
        # Estimate CTF on the downsampled micrographs
        print "Performing CTFfind..."   
        protCTF = ProtCTFFind(numberOfThreads=4)
        protCTF.inputMicrographs.set(protDownsampling.outputMicrographs)        
        self.proj.launchProtocol(protCTF, wait=True)
        self.assertIsNotNone(protCTF.outputCTF, "There was a problem with the CTF estimation")
#         self.validateFiles('protCTF', protCTF)
        
        print "Running Eman fake particle picking..."
        protPP = EmanProtBoxing(importFolder=self.crdsDir, runMode=1) 
        protPP.inputMicrographs.set(protDownsampling.outputMicrographs)
        protPP.boxSize.set(110)
        self.proj.launchProtocol(protPP, wait=True)
        self.assertIsNotNone(protPP.outputCoordinates, "There was a problem with the faked picking")
#         self.protDict['protPP'] = protPP
        
        # Extract the SetOfParticles.
        print "Run extract particles with other downsampling factor"
        protExtract = XmippProtExtractParticles(boxSize=64, downsampleType=2, doFlip=False, downFactor=5, runMode=1, doInvert=False)
        protExtract.inputCoordinates.set(protPP.outputCoordinates)
        protExtract.ctfRelations.set(protCTF.outputCTF)
        protExtract.inputMicrographs.set(protImport.outputMicrographs)
        self.proj.launchProtocol(protExtract, wait=True)
        self.assertIsNotNone(protExtract.outputParticles, "There was a problem with the extract particles")
#         self.validateFiles('protExtract', protExtract)
        
        # Refine the SetOfParticles and reconstruct a refined volume.
        print "Running Frealign..."
        protFrealign = ProtFrealign(angStepSize=15, numberOfIterations=2, mode=1, doExtraRealSpaceSym=True,
                                    innerRadius=100, outerRadius=320, symmetry='I2', PhaseResidual=30,molMass=19400,
                                    score=5, resolution=20, runMode=1, numberOfMpi=1, numberOfThreads=4)
        protFrealign.inputParticles.set(protExtract.outputParticles)
        protFrealign.input3DReference.set(protImportVol.outputVolume)
        self.proj.launchProtocol(protFrealign, wait=True)
        self.assertIsNotNone(protFrealign.outputVolume, "There was a problem with Frealign")
#         self.validateFiles('protFrealign', protFrealign)


class TestMixedBPV2(TestWorkflow):

#     GOLD_FILES = {'protImport': [
#                     'protImport/BPV_1388.mrc',
#                     'protImport/micrographs.sqlite', 
#                     'protImport/BPV_1387.mrc',
#                     'protImport/BPV_1386.mrc'],
#               'protDownsampling': ['protDownsampling/BPV_1388.mrc', 
#                     'protDownsampling/BPV_1387.mrc', 
#                     'protImport/BPV_1386.mrc', 
#                     'protDownsampling/micrographs.xmd', 
#                     'protImport/BPV_1388.mrc', 
#                     'protImport/micrographs.sqlite', 
#                     'protDownsampling/BPV_1386.mrc', 
#                     'protImport/BPV_1387.mrc',
#                     'protDownsampling/logs/run.log',
#                     'protDownsampling/logs/run.db'],
#               'protCTF': [
#                     'protCTF/extra/BPV_1387/xmipp_ctf_ctfmodel_quadrant.xmp', 
#                     'protCTF/extra/BPV_1387/xmipp_ctf.psd', 
#                     'protCTF/extra/BPV_1386/xmipp_ctf_enhanced_psd.xmp', 
#                     'protCTF/extra/BPV_1386/xmipp_ctf.ctfparam', 
#                     'protCTF/extra/BPV_1388/xmipp_ctf.ctfparam', 
#                     'protCTF/extra/BPV_1388/xmipp_ctf.psd', 
#                     'protCTF/extra/BPV_1388/xmipp_ctf_ctfmodel_halfplane.xmp', 
#                     'protCTF/extra/BPV_1387/xmipp_ctf.ctfparam', 
#                     'protCTF/extra/BPV_1388/xmipp_ctf_enhanced_psd.xmp', 
#                     'protCTF/tmp/micrographs.xmd', 
#                     'protCTF/extra/BPV_1387/xmipp_ctf_ctfmodel_halfplane.xmp', 
#                     'protCTF/extra/BPV_1387/xmipp_ctf_enhanced_psd.xmp', 
#                     'protCTF/extra/BPV_1386/xmipp_ctf_ctfmodel_quadrant.xmp', 
#                     'protDownsampling/BPV_1388.mrc', 
#                     'protDownsampling/BPV_1387.mrc', 
#                     'protDownsampling/micrographs.xmd', 
#                     'protCTF/extra/BPV_1386/xmipp_ctf.psd', 
#                     'protCTF/extra/BPV_1386/xmipp_ctf_ctfmodel_halfplane.xmp', 
#                     'protCTF/micrographs.xmd', 
#                     'protDownsampling/BPV_1386.mrc', 
#                     'protCTF/extra/BPV_1388/xmipp_ctf_ctfmodel_quadrant.xmp',
#                     'protCTF/logs/run.log', 
#                     'protCTF/logs/run.db'],
#               'protExtract':[
#                     'protPP/info/BPV_1386_info.json',
#                     'protPP/info/BPV_1387_info.json',
#                     'protPP/info/BPV_1388_info.json',
#                     'protExtract/extra/scipion_micrographs_coordinates.xmd',
#                     'protExtract/images.xmd', 
#                     'protExtract/extra/BPV_1386.pos', 
#                     'protExtract/extra/BPV_1387.pos', 
#                     'protExtract/extra/BPV_1388.pos', 
#                     'protExtract/tmp/BPV_1388_flipped.xmp', 
#                     'protExtract/tmp/BPV_1387_flipped.xmp', 
#                     'protExtract/tmp/BPV_1386_flipped.xmp',
#                     'protExtract/tmp/BPV_1386_noDust.xmp', 
#                     'protExtract/tmp/BPV_1387_noDust.xmp', 
#                     'protExtract/tmp/BPV_1388_noDust.xmp', 
#                     'protExtract/extra/BPV_1386.xmd', 
#                     'protExtract/extra/BPV_1387.xmd', 
#                     'protExtract/extra/BPV_1388.xmd', 
#                     'protExtract/extra/BPV_1388.stk', 
#                     'protExtract/extra/BPV_1386.stk', 
#                     'protExtract/extra/BPV_1387.stk', 
#                     'protExtract/logs/run.log',
#                     'protExtract/logs/run.db',
#                     ],
#                 'protML2D': [
#                     'protExtract/extra/BPV_1386.stk', 
#                     'protExtract/extra/BPV_1387.stk', 
#                     'protExtract/extra/BPV_1388.stk', 
#                     'protExtract/images.xmd', 
#                     'protCTF/extra/BPV_1386/xmipp_ctf.ctfparam',
#                     'protCTF/extra/BPV_1387/xmipp_ctf.ctfparam',
#                     'protCTF/extra/BPV_1388/xmipp_ctf.ctfparam',
#                     'protML2D/ml2d_extra/iter000/iter_images.xmd', 
#                     'protML2D/ml2d_extra/iter000/iter_classes.xmd', 
#                     'protML2D/ml2d_extra/iter000/iter_classes.stk', 
#                     'protML2D/ml2d_extra/iter001/iter_images.xmd', 
#                     'protML2D/ml2d_extra/iter001/iter_classes.xmd', 
#                     'protML2D/ml2d_extra/iter001/iter_classes.stk', 
#                     'protML2D/ml2d_extra/iter002/iter_images.xmd', 
#                     'protML2D/ml2d_extra/iter002/iter_classes.xmd', 
#                     'protML2D/ml2d_extra/iter002/iter_classes.stk', 
# #                     'protML2D/ml2d_extra/iter003/iter_images.xmd', 
# #                     'protML2D/ml2d_extra/iter003/iter_classes.xmd',  
# #                     'protML2D/ml2d_extra/iter003/iter_classes.stk',
# #                     'protML2D/ml2d_extra/iter004/iter_images.xmd', 
# #                     'protML2D/ml2d_extra/iter004/iter_classes.xmd',  
# #                     'protML2D/ml2d_extra/iter004/iter_classes.stk',
#                     'protML2D/ml2d_classes.stk', 
#                     'protML2D/ml2d_images.xmd', 
#                     'protML2D/ml2d__images_average.xmp', 
#                     'protML2D/logs/run.log',
#                     'protML2D/logs/run.db',
#                     'protML2D/ml2d_classes.xmd',
#                     ],
#                 'protIniModel': [
#                     'protML2D/ml2d_classes.stk',
#                     'protML2D/ml2d_classes.xmd',
#                     'protIniModel/initial_models/model_00_01.hdf',
#                     'protIniModel/initial_models/model_00_02.hdf',
#                     'protIniModel/initial_models/model_00_03.hdf',
#                     'protIniModel/initial_models/model_00_04.hdf',
#                     'protIniModel/initial_models/model_00_01_init.hdf',
#                     'protIniModel/initial_models/model_00_02_init.hdf',
#                     'protIniModel/initial_models/model_00_03_init.hdf',
#                     'protIniModel/initial_models/model_00_04_init.hdf',
#                     'protIniModel/initial_models/model_00_01_proj.hdf',
#                     'protIniModel/initial_models/model_00_02_proj.hdf',
#                     'protIniModel/initial_models/model_00_03_proj.hdf',
#                     'protIniModel/initial_models/model_00_04_proj.hdf',
#                     'protIniModel/initial_models/model_00_01_aptcl.hdf',
#                     'protIniModel/initial_models/model_00_02_aptcl.hdf',
#                     'protIniModel/initial_models/model_00_03_aptcl.hdf',
#                     'protIniModel/initial_models/model_00_04_aptcl.hdf',
#                     'protIniModel/initial_models/particles_00.hdf',
#                     'protIniModel/tasks_did2name.json',
#                     'protIniModel/.eman2log.txt',
#                     'protIniModel/logs/run.log',
#                     'protIniModel/logs/run.db',
#                     'protIniModel/tasks_active.json',
#                     'protIniModel/tasks_complete.txt',
#                     'protIniModel/tasks_name2did.json',
#                     'protIniModel/precache_files.json',
#                     'protIniModel/scipion_volumes.json'],
# 
#                 }

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('xmipp_tutorial')
        cls.crdsDir = cls.dataset.getFile('boxingDir')
        cls.micsFn = cls.dataset.getFile('allMics')
        cls.vol1 = cls.dataset.getFile('vol1')
        
    def test_workflow(self):
        #First, import a set of micrographs
        protImport = ProtImportMicrographs(pattern=self.micsFn, samplingRate=1.237, voltage=300)
        self.proj.launchProtocol(protImport, wait=True)
        self.assertIsNotNone(protImport.outputMicrographs, "There was a problem with the import")
#         self.validateFiles('protImport', protImport) 
        
        #Import a set of volumes        
        print "Import Volume"
        protImportVol = ProtImportVolumes(pattern=self.vol1, samplingRate=9.896)
        self.proj.launchProtocol(protImportVol, wait=True)
        self.assertIsNotNone(protImportVol.getFiles(), "There was a problem with the import")
#        self.validateFiles('protImportVol', protImportVol)
        
        # Perform a downsampling on the micrographs
        print "Downsampling..."
        protDownsampling = XmippProtPreprocessMicrographs(doDownsample=True, downFactor=2, doCrop=False, runMode=1)
        protDownsampling.inputMicrographs.set(protImport.outputMicrographs)
        self.proj.launchProtocol(protDownsampling, wait=True)
        self.assertIsNotNone(protDownsampling.outputMicrographs, "There was a problem with the downsampling")
#         self.validateFiles('protDownsampling', protDownsampling)
        
        # Estimate CTF on the downsampled micrographs
        print "Performing CTFfind..."   
        protCTF = ProtCTFFind(numberOfThreads=3)
        protCTF.inputMicrographs.set(protDownsampling.outputMicrographs)        
        self.proj.launchProtocol(protCTF, wait=True)
        self.assertIsNotNone(protCTF.outputCTF, "There was a problem with the CTF estimation")
#         self.validateFiles('protCTF', protCTF)
        
        print "Running Eman fake particle picking..."
        protPP = EmanProtBoxing(importFolder=self.crdsDir, runMode=1) 
        protPP.inputMicrographs.set(protImport.outputMicrographs)
        protPP.boxSize.set(550)
        self.proj.launchProtocol(protPP, wait=True)
        self.assertIsNotNone(protPP.outputCoordinates, "There was a problem with the faked picking")
#         self.protDict['protPP'] = protPP
        
        print "<Run extract particles with Same as picking>"
        protExtract = XmippProtExtractParticles(boxSize=550, downsampleType=1, doFlip=True, doInvert=True, runMode=1)
        protExtract.inputCoordinates.set(protPP.outputCoordinates)
        protExtract.ctfRelations.set(protCTF.outputCTF)
        self.proj.launchProtocol(protExtract, wait=True)
        self.assertIsNotNone(protExtract.outputParticles, "There was a problem with the extract particles")
        #self.validateFiles('protExtract', protExtract)
        
        print "Run Preprocess particles"
        protCropResize = XmippProtCropResizeParticles(doResize=True, resizeOption=1, resizeDim=110)
        protCropResize.inputParticles.set(protExtract.outputParticles)
        self.proj.launchProtocol(protCropResize, wait=True)
        
        self.assertIsNotNone(protCropResize.outputParticles, "There was a problem with resize/crop the particles")
        
        print "Run ML2D"
        protML2D = XmippProtML2D(numberOfReferences=8, maxIters=2, 
                                 numberOfMpi=2, numberOfThreads=2)
        protML2D.inputParticles.set(protCropResize.outputParticles)
        self.proj.launchProtocol(protML2D, wait=True)        
        self.assertIsNotNone(protML2D.outputClasses, "There was a problem with ML2D")  
        #self.validateFiles('protML2D', protML2D)
        
        print "Run Initial Model"
        protIniModel = EmanProtInitModel(numberOfIterations=1, numberOfModels=2,
                                 shrink=1, symmetry='icos', numberOfThreads=3)
        protIniModel.inputClasses.set(protML2D.outputClasses)
        self.proj.launchProtocol(protIniModel, wait=True)        
        self.assertIsNotNone(protIniModel.outputVolumes, "There was a problem with Initial Model")  
        #self.validateFiles('protIniModel', protIniModel)


class TestMixedRelionTutorial(TestWorkflow):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('relion_tutorial')
        cls.crdsEmanDir = cls.dataset.getFile('boxingDir')
        cls.crdsXmippDir = cls.dataset.getFile('posAllDir')
        cls.micsFn = cls.dataset.getFile('allMics')
        cls.vol = cls.dataset.getFile('volume')
    
    def test_workflow(self):
        #First, import a set of micrographs
        print "Importing a set of micrographs..."
        protImport = ProtImportMicrographs(pattern=self.micsFn, samplingRateMode=1, magnification=79096,
                                           scannedPixelSize=56, voltage=300, sphericalAberration=2.0)
        protImport.setObjLabel('import 20 mics')
        self.proj.launchProtocol(protImport, wait=True)
        self.assertIsNotNone(protImport.outputMicrographs, "There was a problem with the import")
        
        print "Importing a volume..."
        protImportVol = ProtImportVolumes(pattern=self.vol, samplingRate=7.08)
        protImportVol.setObjLabel('import single vol')
        self.proj.launchProtocol(protImportVol, wait=True)
        self.assertIsNotNone(protImportVol.outputVolume, "There was a problem with the import")
        
        print "Preprocessing the micrographs..."
        protPreprocess = XmippProtPreprocessMicrographs(doCrop=True, cropPixels=50)
        protPreprocess.inputMicrographs.set(protImport.outputMicrographs)
        protPreprocess.setObjLabel('crop 50px')
        self.proj.launchProtocol(protPreprocess, wait=True)
        self.assertIsNotNone(protPreprocess.outputMicrographs, "There was a problem with the downsampling")

        # Now estimate CTF on the micrographs with ctffind 
        print "Performing CTFfind..."   
        protCTF = ProtCTFFind(lowRes=0.04, highRes=0.45, minDefocus=1.2, maxDefocus=3,
                              runMode=1, numberOfMpi=1, numberOfThreads=16)         
        protCTF.inputMicrographs.set(protPreprocess.outputMicrographs)
        protCTF.setObjLabel('CTF ctffind')
        self.proj.launchProtocol(protCTF, wait=True)
        
        print "Running Eman fake particle picking..."
        protPP = EmanProtBoxing(importFolder=self.crdsEmanDir, runMode=1)                
        protPP.inputMicrographs.set(protPreprocess.outputMicrographs)  
        protPP.boxSize.set(60)
        protPP.setObjLabel('Eman boxing') 
        self.proj.launchProtocol(protPP, wait=True)
        self.assertIsNotNone(protPP.outputCoordinates, "There was a problem with the Eman faked picking")
        
        print "Run extract particles with <Same as picking> option"
        protExtract = XmippProtExtractParticles(boxSize=60, downsampleType=1, doRemoveDust=False,
                                                doFlip=False, backRadius=28, runMode=1)
        protExtract.inputCoordinates.set(protPP.outputCoordinates)
        protExtract.ctfRelations.set(protCTF.outputCTF)
        protExtract.setObjLabel('Extract particles')
        self.proj.launchProtocol(protExtract, wait=True)
        self.assertIsNotNone(protExtract.outputParticles, "There was a problem with the extract particles")
        
        print "Run CL2D"
        protCL2D = XmippProtCL2D(numberOfReferences=32, numberOfInitialReferences=4, 
                                 numberOfIterations=2, numberOfMpi=16)
        protCL2D.inputImages.set(protExtract.outputParticles)
        protCL2D.setObjLabel('CL2D')
        self.proj.launchProtocol(protCL2D, wait=True)   
        self.assertIsNotNone(protCL2D.outputClasses, "There was a problem with CL2D")
        
        # Refine the SetOfParticles and reconstruct a refined volume.
        print "Running Frealign..."
        protFrealign = ProtFrealign(angStepSize=20, numberOfIterations=2, mode=1, doExtraRealSpaceSym=True,
                                    outerRadius=180, PhaseResidual=65, lowResolRefine=300, highResolRefine=15,
                                    resolution=15, runMode=1, numberOfMpi=1, numberOfThreads=16)
        protFrealign.inputParticles.set(protExtract.outputParticles)
        protFrealign.input3DReference.set(protImportVol.outputVolume)
        protFrealign.setObjLabel('Frealign')
        self.proj.launchProtocol(protFrealign, wait=True)        
        self.assertIsNotNone(protFrealign.outputVolume, "There was a problem with Frealign")
        
        # Now estimate CTF on the micrographs with xmipp
        print "Performing Xmipp CTF..."   
        protCTF2 = XmippProtCTFMicrographs(lowRes=0.04, highRes=0.45, minDefocus=1.2, maxDefocus=3,
                              runMode=1, numberOfMpi=1, numberOfThreads=16)         
        protCTF2.inputMicrographs.set(protPreprocess.outputMicrographs)
        protCTF2.setObjLabel('CTF xmipp')
        self.proj.launchProtocol(protCTF2, wait=True)
        
        print "Running Xmipp fake particle picking..."
        protPP2 = XmippProtParticlePicking(importFolder=self.crdsXmippDir, runMode=1)                
        protPP2.inputMicrographs.set(protPreprocess.outputMicrographs)  
        protPP2.setObjLabel('Xmipp Picking') 
        self.proj.launchProtocol(protPP2, wait=True)
        self.assertIsNotNone(protPP2.outputCoordinates, "There was a problem with the Xmipp faked picking")
        
        print "Run extract particles with <Same as picking> option"
        protExtract2 = XmippProtExtractParticles(boxSize=60, downsampleType=1, doRemoveDust=False, doInvert=True,
                                                doFlip=False, backRadius=28, runMode=1)
        protExtract2.inputCoordinates.set(protPP2.outputCoordinates)
        protExtract2.ctfRelations.set(protCTF2.outputCTF)
        protExtract2.setObjLabel('Extract particles')
        self.proj.launchProtocol(protExtract2, wait=True)
        self.assertIsNotNone(protExtract2.outputParticles, "There was a problem with the extract particles")
        
        print "Run Relion Classification2d"
        prot2D = ProtRelionClassify2D(regularisationParamT=2, numberOfMpi=4, numberOfThreads=4)
        prot2D.numberOfClasses.set(50)
        prot2D.numberOfIterations.set(25)
        prot2D.inputParticles.set(protExtract2.outputParticles)
        prot2D.setObjLabel('relion 2D')
        self.proj.launchProtocol(prot2D, wait=True)        
        self.assertIsNotNone(prot2D.outputClasses, "There was a problem with Relion 2D:\n" + (prot2D.getErrorMessage() or "No error set"))


class TestMixedFrealignClassify(TestWorkflow):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('relion_tutorial')
        cls.crdsEmanDir = cls.dataset.getFile('boxingDir')
        cls.crdsXmippDir = cls.dataset.getFile('posAllDir')
        cls.micsFn = cls.dataset.getFile('allMics')
        cls.vol = cls.dataset.getFile('volume')
    
    def test_workflow(self):
        #First, import a set of micrographs
        print "Importing a set of micrographs..."
        protImport = ProtImportMicrographs(pattern=self.micsFn, samplingRateMode=1, magnification=79096,
                                           scannedPixelSize=56, voltage=300, sphericalAberration=2.0)
        protImport.setObjLabel('import 20 mics')
        self.proj.launchProtocol(protImport, wait=True)
        self.assertIsNotNone(protImport.outputMicrographs, "There was a problem with the import")
        
        print "Importing a volume..."
        protImportVol = ProtImportVolumes(pattern=self.vol, samplingRate=7.08)
        protImportVol.setObjLabel('import single vol')
        self.proj.launchProtocol(protImportVol, wait=True)
        self.assertIsNotNone(protImportVol.outputVolume, "There was a problem with the import")
        
        print "Preprocessing the micrographs..."
        protPreprocess = XmippProtPreprocessMicrographs(doCrop=True, cropPixels=50)
        protPreprocess.inputMicrographs.set(protImport.outputMicrographs)
        protPreprocess.setObjLabel('crop 50px')
        self.proj.launchProtocol(protPreprocess, wait=True)
        self.assertIsNotNone(protPreprocess.outputMicrographs, "There was a problem with the downsampling")

        # Now estimate CTF on the micrographs with ctffind 
        print "Performing CTFfind..."   
        protCTF = ProtCTFFind(lowRes=0.04, highRes=0.45, minDefocus=1.2, maxDefocus=3,
                              runMode=1, numberOfMpi=1, numberOfThreads=16)         
        protCTF.inputMicrographs.set(protPreprocess.outputMicrographs)
        protCTF.setObjLabel('CTF ctffind')
        self.proj.launchProtocol(protCTF, wait=True)
        
        print "Running Eman fake particle picking..."
        protPP = EmanProtBoxing(importFolder=self.crdsEmanDir, runMode=1)                
        protPP.inputMicrographs.set(protPreprocess.outputMicrographs)  
        protPP.boxSize.set(60)
        protPP.setObjLabel('Eman boxing') 
        self.proj.launchProtocol(protPP, wait=True)
        self.assertIsNotNone(protPP.outputCoordinates, "There was a problem with the Eman faked picking")
        
        print "Run extract particles with <Same as picking> option"
        protExtract = XmippProtExtractParticles(boxSize=60, downsampleType=1, doRemoveDust=False,
                                                doFlip=False, backRadius=28, runMode=1)
        protExtract.inputCoordinates.set(protPP.outputCoordinates)
        protExtract.ctfRelations.set(protCTF.outputCTF)
        protExtract.setObjLabel('Extract particles')
        self.proj.launchProtocol(protExtract, wait=True)
        self.assertIsNotNone(protExtract.outputParticles, "There was a problem with the extract particles")
        
        # Classify the SetOfParticles.
        print "Running Frealign Classification..."
        protFrealign = ProtFrealignClassify(numberOfClasses=3, itRefineAngles=2, itRefineShifts=3, angStepSize=20, numberOfIterations=6, mode=1, doExtraRealSpaceSym=True,
                                    outerRadius=180, PhaseResidual=65, lowResolRefine=300, highResolRefine=15,
                                    resolution=15, runMode=1, numberOfThreads=16)
        protFrealign.inputParticles.set(protExtract.outputParticles)
        protFrealign.input3DReference.set(protImportVol.outputVolume)
        protFrealign.setObjLabel('Frealign')
        self.proj.launchProtocol(protFrealign, wait=True)        
        self.assertIsNotNone(protFrealign.outputClasses, "There was a problem with Frealign")


if __name__ == "__main__":
    unittest.main()
