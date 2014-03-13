import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3 import *
from pyworkflow.em.packages.brandeis import *
from pyworkflow.em.packages.eman2 import *
from test_workflow import TestWorkflow
    
    
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

        # Now estimate CTF on the micrographs 
        print "Performing CTFfind..."   
        protCTF = ProtCTFFind(lowRes=0.04, highRes=0.45, minDefocus=1.2, maxDefocus=3,
                              runMode=1, numberOfMpi=1, numberOfThreads=16)         
        protCTF.inputMicrographs.set(protPreprocess.outputMicrographs)
        protCTF.setObjLabel('ctf estimation')
        self.proj.launchProtocol(protCTF, wait=True)
        
        print "Running Eman fake particle picking..."
        protPP = EmanProtBoxing(importFolder=self.importFolder, runMode=1)                
        protPP.inputMicrographs.set(protPreprocess.outputMicrographs)  
        protPP.boxSize.set(60)
        protPP.setObjLabel('Eman boxing') 
        self.proj.launchProtocol(protPP, wait=True)
        self.assertIsNotNone(protPP.outputCoordinates, "There was a problem with the faked picking")
        #self.protDict['protPP'] = protPP
        
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

if __name__ == "__main__":
    unittest.main()