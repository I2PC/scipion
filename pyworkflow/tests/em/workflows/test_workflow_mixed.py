from pyworkflow.tests import DataSet, setupTestProject
from pyworkflow.em.packages.xmipp3 import (XmippProtPreprocessMicrographs,
                                           XmippProtExtractParticles,
                                           XmippProtCropResizeParticles, XmippProtML2D)
from pyworkflow.em.packages.xmipp3.constants import SAME_AS_PICKING, OTHER 
from pyworkflow.em.packages.brandeis import ProtCTFFind, ProtFrealign
from pyworkflow.em.packages.eman2 import EmanProtBoxing, EmanProtInitModel
from pyworkflow.em.packages.relion import ProtImportMicrographs, ProtImportVolumes
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
        protImport = self.newProtocol(ProtImportMicrographs, pattern=self.micsFn, 
                                      samplingRate=1.237, voltage=300)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputMicrographs, "There was a problem with the import")
#         self.validateFiles('protImport', protImport) 
        
        #Import a set of volumes        
        print "Import Volume"
        protImportVol = self.newProtocol(ProtImportVolumes, pattern=self.vol1, samplingRate=9.896)
        self.launchProtocol(protImportVol)
        self.assertIsNotNone(protImportVol.getFiles(), "There was a problem with the import")
#        self.validateFiles('protImportVol', protImportVol)
        
        # Perform a downsampling on the micrographs
        print "Downsampling..."
        protDownsampling = self.newProtocol(XmippProtPreprocessMicrographs, doDownsample=True, downFactor=5, doCrop=False, runMode=1)
        protDownsampling.inputMicrographs.set(protImport.outputMicrographs)
        self.launchProtocol(protDownsampling)
        self.assertIsNotNone(protDownsampling.outputMicrographs, "There was a problem with the downsampling")
#         self.validateFiles('protDownsampling', protDownsampling)
        
        # Estimate CTF on the downsampled micrographs
        print "Performing CTFfind..."   
        protCTF = self.newProtocol(ProtCTFFind, numberOfThreads=4, minDefocus=2.2, maxDefocus=2.5)
        protCTF.inputMicrographs.set(protDownsampling.outputMicrographs)        
        self.launchProtocol(protCTF)
        self.assertIsNotNone(protCTF.outputCTF, "There was a problem with the CTF estimation")
#         self.validateFiles('protCTF', protCTF)
        
        print "Running Eman fake particle picking..."
        protPP = self.newProtocol(EmanProtBoxing, importFolder=self.crdsDir, runMode=1) 
        protPP.inputMicrographs.set(protDownsampling.outputMicrographs)
        protPP.boxSize.set(110)
        self.launchProtocol(protPP)
        self.assertIsNotNone(protPP.outputCoordinates, "There was a problem with the faked picking")
#         self.protDict['protPP'] = protPP
        
        # Extract the SetOfParticles.
        print "Run extract particles with other downsampling factor"
        protExtract = self.newProtocol(XmippProtExtractParticles, boxSize=64, downsampleType=OTHER, doFlip=False, downFactor=5, runMode=1, doInvert=False)
        protExtract.inputCoordinates.set(protPP.outputCoordinates)
        protExtract.ctfRelations.set(protCTF.outputCTF)
        protExtract.inputMicrographs.set(protImport.outputMicrographs)
        self.launchProtocol(protExtract)
        self.assertIsNotNone(protExtract.outputParticles, "There was a problem with the extract particles")
#         self.validateFiles('protExtract', protExtract)
        
        # Refine the SetOfParticles and reconstruct a refined volume.
        print "Running Frealign..."
        protFrealign = self.newProtocol(ProtFrealign, angStepSize=15, numberOfIterations=2, mode=1, doExtraRealSpaceSym=True,
                                    innerRadius=130, outerRadius=300, symmetry='I1', PhaseResidual=30,molMass=19400,
                                    score=5, resolution=20, runMode=1, numberOfMpi=5)
        protFrealign.inputParticles.set(protExtract.outputParticles)
        protFrealign.input3DReference.set(protImportVol.outputVolume)
        self.launchProtocol(protFrealign)
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
        protImport = self.newProtocol(ProtImportMicrographs, pattern=self.micsFn, samplingRate=1.237, voltage=300)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputMicrographs, "There was a problem with the import")
#         self.validateFiles('protImport', protImport) 
        
        #Import a set of volumes        
        print "Import Volume"
        protImportVol = self.newProtocol(ProtImportVolumes, pattern=self.vol1, samplingRate=9.896)
        self.launchProtocol(protImportVol)
        self.assertIsNotNone(protImportVol.getFiles(), "There was a problem with the import")
#        self.validateFiles('protImportVol', protImportVol)
        
        # Perform a downsampling on the micrographs
        print "Downsampling..."
        protDownsampling = self.newProtocol(XmippProtPreprocessMicrographs, doDownsample=True, downFactor=5, doCrop=False, runMode=1)
        protDownsampling.inputMicrographs.set(protImport.outputMicrographs)
        self.launchProtocol(protDownsampling)
        self.assertIsNotNone(protDownsampling.outputMicrographs, "There was a problem with the downsampling")
#         self.validateFiles('protDownsampling', protDownsampling)
        
        # Estimate CTF on the downsampled micrographs
        print "Performing CTFfind..."   
        protCTF = self.newProtocol(ProtCTFFind, numberOfThreads=4, minDefocus=2.2, maxDefocus=2.5)
        protCTF.inputMicrographs.set(protImport.outputMicrographs)        
        self.launchProtocol(protCTF)
        self.assertIsNotNone(protCTF.outputCTF, "There was a problem with the CTF estimation")
#         self.validateFiles('protCTF', protCTF)
        
        print "Running Eman fake particle picking..."
        protPP = self.newProtocol(EmanProtBoxing, importFolder=self.crdsDir, runMode=1) 
        protPP.inputMicrographs.set(protDownsampling.outputMicrographs)
        protPP.boxSize.set(550)
        self.launchProtocol(protPP)
        self.assertIsNotNone(protPP.outputCoordinates, "There was a problem with the faked picking")
#         self.protDict['protPP'] = protPP
        
        print "<Run extract particles with Same as picking>"
        protExtract = self.newProtocol(XmippProtExtractParticles, boxSize=110, downsampleType=SAME_AS_PICKING, doFlip=True, doInvert=True, runMode=1)
        protExtract.inputCoordinates.set(protPP.outputCoordinates)
        protExtract.ctfRelations.set(protCTF.outputCTF)
        self.launchProtocol(protExtract)
        self.assertIsNotNone(protExtract.outputParticles, "There was a problem with the extract particles")
        #self.validateFiles('protExtract', protExtract)
        
        print "Run Preprocess particles"
        protCropResize = self.newProtocol(XmippProtCropResizeParticles, doResize=True, resizeOption=1, resizeDim=110)
        protCropResize.inputParticles.set(protExtract.outputParticles)
        self.launchProtocol(protCropResize)
        
        self.assertIsNotNone(protCropResize.outputParticles, "There was a problem with resize/crop the particles")
        
        print "Run ML2D"
        protML2D = self.newProtocol(XmippProtML2D, numberOfReferences=8, maxIters=2, 
                                 numberOfMpi=2, numberOfThreads=2)
        protML2D.inputParticles.set(protCropResize.outputParticles)
        self.launchProtocol(protML2D)        
        self.assertIsNotNone(protML2D.outputClasses, "There was a problem with ML2D")  
        #self.validateFiles('protML2D', protML2D)
        
#         #FIXME: Check the initial model with EMAn and restore the next step
#         return
        
        print "Run Initial Model"
        protIniModel = self.newProtocol(EmanProtInitModel, numberOfIterations=1, numberOfModels=2,
                                 shrink=1, symmetry='icos', numberOfThreads=3)
        protIniModel.inputClasses.set(protML2D.outputClasses)
        self.launchProtocol(protIniModel)        
        self.assertIsNotNone(protIniModel.outputVolumes, "There was a problem with Initial Model")  


