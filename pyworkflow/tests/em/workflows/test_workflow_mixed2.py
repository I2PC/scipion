import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3 import *
from pyworkflow.em.packages.brandeis import *
from pyworkflow.em.packages.eman2 import *
from test_workflow import TestWorkflow
    
    
#class TestMixedWorkflow_1(TestWorkflow):
#
#    GOLD_FILES = {'protImport': [
#                    'protImport/BPV_1386.mrc',
#                    'protImport/micrographs.sqlite'],
#              'protDownsampling': [
#                    'protDownsampling/BPV_1386.mrc', 
#                    'protImport/BPV_1386.mrc',
#                    'protDownsampling/micrographs.xmd', 
#                    'protImport/micrographs.sqlite', 
#                    'protDownsampling/logs/run.log',
#                    'protDownsampling/logs/run.db'
#                    ],
#              'protCTF': [
#                    'protCTF/extra/BPV_1386/ctffind_psd.mrc', 
#                    'protCTF/extra/BPV_1386/ctffind.out', 
#                    'protCTF/micrographs.sqlite',
#                    'protDownsampling/micrographs.xmd', 
#                    'protDownsampling/BPV_1386.mrc',
#                    'protCTF/logs/run.log', 
#                    'protCTF/logs/run.db'],
#              'protExtract':[
#                    'protPicking/extra/BPV_1386.pos', 
#                    'protExtract/tmp/BPV_1386_noDust.xmp', 
#                    'protExtract/extra/BPV_1386.xmd', 
#                    'protExtract/images.xmd', 
#                    'protExtract/extra/BPV_1386.stk', 
#                    'protExtract/extra/BPV_1386.pos', 
#                    'protExtract/tmp/BPV_1386_flipped.xmp',
#                    'protExtract/logs/run.log',
#                    'protExtract/logs/run.db'],
#              }
#    
#    @classmethod
#    def setUpClass(cls):    
#        # Create a new project
#        setupProject(cls)
#        cls.pattern = getInputPath('Micrographs_BPV1', '*.mrc')        
#        cls.importFolder = getInputPath('Picking_XmippBPV1')
#        
#    def testWorkflow(self):
#        #First, import a set of micrographs
#        protImport = ProtImportMicrographs(pattern=self.pattern, samplingRate=1.237, voltage=300)
#        self.proj.launchProtocol(protImport, wait=True)
#        
#        self.assertIsNotNone(protImport.outputMicrographs, "There was a problem with the import")
#        self.validateFiles('protImport', protImport) 
#        
#        # Perform a downsampling on the micrographs
#
#        print "Downsampling..."
#        protDownsampling = XmippProtPreprocessMicrographs(doDownsample=True, downFactor=3, doCrop=False,
#                                                          numberOfMpi=1, numberOfThreads=3)
#        protDownsampling.inputMicrographs.set(protImport.outputMicrographs)
#        self.proj.launchProtocol(protDownsampling, wait=True)
#          
#        self.assertIsNotNone(protDownsampling.outputMicrographs, "There was a problem with the downsampling")
#        self.validateFiles('protDownsampling', protDownsampling) 
#
#
#        # Now estimate CTF on the downsampled micrographs 
#        print "Performing CTFfind..."   
#        protCTF = ProtCTFFind(runMode=1, numberOfMpi=1, numberOfThreads=3)         
#        protCTF.inputMicrographs.set(protDownsampling.outputMicrographs)        
#        self.proj.launchProtocol(protCTF, wait=True)
#        
#        self.validateFiles('protCTF', protCTF) 
#        
#        print "Running fake particle picking..."   
#        protPicking = XmippProtParticlePicking(importFolder=self.importFolder, runMode=1)                
#        protPicking.inputMicrographs.set(protCTF.outputMicrographs)        
#        self.proj.launchProtocol(protPicking, wait=True)
#        self.protDict['protPicking'] = protPicking
#            
#        self.assertIsNotNone(protPicking.outputCoordinates, "There was a problem with the faked picking")
#            
#        print "Run extract particles with Same as picking"
#        protExtract = XmippProtExtractParticles(boxSize=171, downsampleType=1, runMode=1)
#        protExtract.inputCoordinates.set(protPicking.outputCoordinates)
#        #protExtract.inputMicrographs.set(protDownsampling.outputMicrographs)
#        self.proj.launchProtocol(protExtract, wait=True)
#        
#        self.assertIsNotNone(protExtract.outputParticles, "There was a problem with the extract particles")
#        self.validateFiles('protExtract', protExtract)
        
        
class ScipionMixedWorkflow(TestWorkflow):
        
    @classmethod
    def setUpClass(cls):    
        # Create a new project
        setupProject(cls)
        cls.pattern = getInputPath('Ribosomes_Sjors', 'Mics', '*.mrc')
        cls.importFolder = getInputPath('Ribosomes_Sjors', 'EmanBoxing')
        cls.importVol = getInputPath('Ribosomes_Sjors', 'reference.mrc')
        
    def testWorkflow(self):
        #First, import a set of micrographs
        print "Importing a set of micrographs..."
        protImport = ProtImportMicrographs(pattern=self.pattern, samplingRateMode=1, magnification=79096,
                                           scannedPixelSize=56, voltage=300, sphericalAberration=2.0)
        self.proj.launchProtocol(protImport, wait=True)
        self.assertIsNotNone(protImport.outputMicrographs, "There was a problem with the import")
        
        print "Importing a volume..."
        protImportVol = ProtImportVolumes(pattern=self.importVol, samplingRate=7.08)
        self.proj.launchProtocol(protImportVol, wait=True)
        self.assertIsNotNone(protImportVol.outputVolume, "There was a problem with the import")
        
        print "Preprocessing the micrographs..."
        protPreprocess = XmippProtPreprocessMicrographs(doCrop=True, cropPixels=50)
        protPreprocess.inputMicrographs.set(protImport.outputMicrographs)
        self.proj.launchProtocol(protPreprocess, wait=True)
        self.assertIsNotNone(protPreprocess.outputMicrographs, "There was a problem with the downsampling")

        # Now estimate CTF on the micrographs 
        print "Performing CTFfind..."   
        protCTF = ProtCTFFind(lowRes=0.04, highRes=0.45, minDefocus=1.2, maxDefocus=3,
                              runMode=1, numberOfMpi=1, numberOfThreads=4)         
        protCTF.inputMicrographs.set(protPreprocess.outputMicrographs)        
        self.proj.launchProtocol(protCTF, wait=True)
        
        print "Running Eman fake particle picking..."
        protPP = EmanProtBoxing(importFolder=self.importFolder, runMode=1)                
        protPP.inputMicrographs.set(protPreprocess.outputMicrographs)  
        protPP.boxSize.set(60)      
        self.proj.launchProtocol(protPP, wait=True)
        self.assertIsNotNone(protPP.outputCoordinates, "There was a problem with the faked picking")
        #self.protDict['protPP'] = protPP
        
        print "Run extract particles with <Same as picking> option"
        protExtract = XmippProtExtractParticles(boxSize=60, downsampleType=1, doRemoveDust=False,
                                                doFlip=False, backRadius=28, runMode=1)
        protExtract.inputCoordinates.set(protPP.outputCoordinates)
        protExtract.ctfRelations.set(protCTF.outputCTF)
        self.proj.launchProtocol(protExtract, wait=True)
        self.assertIsNotNone(protExtract.outputParticles, "There was a problem with the extract particles")
        
        print "Run CL2D"
        protCL2D = XmippProtCL2D(numberOfReferences=8, numberOfInitialReferences=2, 
                                 numberOfIterations=2, numberOfMpi=4)
        protCL2D.inputImages.set(protExtract.outputParticles)
        self.proj.launchProtocol(protCL2D, wait=True)   
        self.assertIsNotNone(protCL2D.outputClasses, "There was a problem with CL2D")
        
        # Refine the SetOfParticles and reconstruct a refined volume.    
        print "Running Frealign..."
        protFrealign = ProtFrealign(angStepSize=10, numberOfIterations=2, mode=1, doExtraRealSpaceSym=True,
                                    outerRadius=180, PhaseResidual=65, lowResolRefine=300, highResolRefine=15,
                                    resolution=15, runMode=1, numberOfMpi=1, numberOfThreads=4)
        protFrealign.inputParticles.set(protExtract.outputParticles)
        protFrealign.input3DReferences.set(protImportVol.outputVolume)
        self.proj.launchProtocol(protFrealign, wait=True)        
        self.assertIsNotNone(protFrealign.outputVolume, "There was a problem with Frealign")

if __name__ == "__main__":
    unittest.main()