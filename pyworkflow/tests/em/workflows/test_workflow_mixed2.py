import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3 import *
from pyworkflow.em.packages.brandeis import *
from pyworkflow.em.packages.eman2 import *
from pyworkflow.em.packages.relion import *
from test_workflow import TestWorkflow
    
    
class ScipionMixedWorkflow(TestWorkflow):
        
    @classmethod
    def setUpClass(cls):    
        # Create a new project
        setupProject(cls)
        cls.pattern = getInputPath('Ribosomes_Sjors', 'Mics', '*.mrc')
        cls.importFolder1 = getInputPath('Ribosomes_Sjors', 'EmanBoxing')
        cls.importFolder2 = getInputPath('Ribosomes_Sjors', 'XmippPicking')
        cls.importVol = getInputPath('Ribosomes_Sjors', 'reference.mrc')
        
    def testWorkflow(self):
        #First, import a set of micrographs
        print "Importing a set of micrographs..."
        protImport = ProtImportMicrographs(pattern=self.pattern, samplingRateMode=1, magnification=79096,
                                           scannedPixelSize=56, voltage=300, sphericalAberration=2.0)
        protImport.setObjLabel('import 20 mics')
        self.proj.launchProtocol(protImport, wait=True)
        self.assertIsNotNone(protImport.outputMicrographs, "There was a problem with the import")
        
        print "Importing a volume..."
        protImportVol = ProtImportVolumes(pattern=self.importVol, samplingRate=7.08)
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
        protPP = EmanProtBoxing(importFolder=self.importFolder1, runMode=1)                
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
        protFrealign.input3DReferences.set(protImportVol.outputVolume)
        protFrealign.setObjLabel('Frealign')
        self.proj.launchProtocol(protFrealign, wait=True)        
        self.assertIsNotNone(protFrealign.outputVolume, "There was a problem with Frealign")
        
        print "Run Preprocess particles"
        protPreproc = XmippProtPreprocessParticles(doRemoveDust=False, doNormalize=True, doInvert=True, doThreshold=False)
        protPreproc.inputParticles.set(protExtract.outputParticles)
        protPreproc.setObjLabel('Invert')
        self.proj.launchProtocol(protPreproc, wait=True)
        self.assertIsNotNone(protPreproc.outputParticles, "There was a problem with preprocess particles")
        
        print "Run Relion Classification2d"
        prot2D = ProtRelionClassify2D(regularisationParamT=2, numberOfMpi=4, numberOfThreads=4)
        prot2D.numberOfClasses.set(50)
        prot2D.numberOfIterations.set(5)
        prot2D.inputParticles.set(protPreproc.outputParticles)
        prot2D.setObjLabel('relion 2D')
        self.proj.launchProtocol(prot2D, wait=True)        
        self.assertIsNotNone(prot2D.outputClasses, "There was a problem with Relion 2D:\n" + (prot2D.getErrorMessage() or "No error set"))
        
        # Now estimate CTF on the micrographs with xmipp
        print "Performing Xmipp CTF..."   
        protCTF2 = XmippProtCTFMicrographs(lowRes=0.04, highRes=0.45, minDefocus=1.2, maxDefocus=3,
                              runMode=1, numberOfMpi=1, numberOfThreads=16)         
        protCTF2.inputMicrographs.set(protPreprocess.outputMicrographs)
        protCTF2.setObjLabel('CTF xmipp')
        self.proj.launchProtocol(protCTF2, wait=True)
        
        print "Running Xmipp fake particle picking..."
        protPP2 = XmippProtParticlePicking(importFolder=self.importFolder2, runMode=1)                
        protPP2.inputMicrographs.set(protPreprocess.outputMicrographs)  
        protPP2.setObjLabel('Xmipp Picking') 
        self.proj.launchProtocol(protPP2, wait=True)
        self.assertIsNotNone(protPP2.outputCoordinates, "There was a problem with the Xmipp faked picking")
        
        print "Run extract particles with <Same as picking> option"
        protExtract2 = XmippProtExtractParticles(boxSize=60, downsampleType=1, doRemoveDust=False,
                                                doFlip=False, backRadius=28, runMode=1)
        protExtract2.inputCoordinates.set(protPP2.outputCoordinates)
        protExtract2.ctfRelations.set(protCTF2.outputCTF)
        protExtract2.setObjLabel('Extract particles')
        self.proj.launchProtocol(protExtract2, wait=True)
        self.assertIsNotNone(protExtract2.outputParticles, "There was a problem with the extract particles")

if __name__ == "__main__":
    unittest.main()