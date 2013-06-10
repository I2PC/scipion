import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3 import *
    
    
class TestXmippWorkflow(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):    
        # Create a new project
        setupProject(cls)
        cls.pattern = getInputPath('Micrographs_BPV3', '*.mrc')        
        cls.importFolder = getInputPath('Picking_XmippBPV3')

    def validateFiles(self, prot, filesSet):
        """ Validate if the produced files are the expected ones.
        Params:
            prot: the protocol to validate. 
            filesSet: the known files that should be produced (set)
        """
        self.assertEqual(prot.getFiles(), filesSet)
        
    def testXmippWorkflow(self):
        #First, import a set of micrographs
        protImport = ProtImportMicrographs(pattern=self.pattern, samplingRate=1.237, voltage=300)
        self.proj.launchProtocol(protImport, wait=True)
        
        print("******** IMPORT MICROGRAPHS ***********")
        if (protImport.getFiles() is not None):
            print(str(protImport.getFiles()))
        else:
            print ("NONEEEEEEEEEEEEEEEEE")
        
        self.assertIsNotNone(protImport.outputMicrographs, "There was a problem with the import")
        
        # Perform a downsampling on the micrographs

        print "Downsampling..."
        protDownsampling = XmippProtPreprocessMicrographs(doDownsample=True, downFactor=3, doCrop=False)
        protDownsampling.inputMicrographs.set(protImport.outputMicrographs)
        self.proj.launchProtocol(protDownsampling, wait=True)
          
        print("******** DOWNSAMPLING ***********")
        if (protDownsampling.getFiles() is not None):
            print(str(protDownsampling.getFiles())) 
        else:
            print ("NONEEEEEEEEEEEEEEEEE") 
          
        self.assertIsNotNone(protDownsampling.outputMicrographs, "There was a problem with the downsampling")
          
        # Now estimate CTF on the downsampled micrographs 
        print "Performing CTF..."   
        protCTF = XmippProtCTFMicrographs()                
        protCTF.inputMicrographs.set(protDownsampling.outputMicrographs)        
        self.proj.launchProtocol(protCTF, wait=True)
        
        print("******** CTF ***********")
        if (protCTF.getFiles() is not None):
            print(str(protCTF.getFiles()))  
        else:
            print ("NONEEEEEEEEEEEEEEEEE") 
        
        print "Running fake particle picking..."   
        protPP = XmippProtParticlePicking(importFolder=self.importFolder)                
        protPP.inputMicrographs.set(protCTF.outputMicrographs)        
        self.proj.launchProtocol(protPP, wait=True)
            
        self.assertIsNotNone(protPP.outputCoordinates, "There was a problem with the faked picking")
            
        print "Run extract particles with Same as picking"
        protExtract = XmippProtExtractParticles(boxSize=258, downsampleType=1)
        protExtract.inputCoordinates.set(protPP.outputCoordinates)
        #protExtract.inputMicrographs.set(protDownsampling.outputMicrographs)
        self.proj.launchProtocol(protExtract, wait=True)
        
        print("******** EXTRACT ***********")
        if (protExtract.getFiles() is not None):
            print(str(protExtract.getFiles()))  
        else:
            print ("NONEEEEEEEEEEEEEEEEE") 
        
        self.assertIsNotNone(protExtract.outputImages, "There was a problem with the extract particles")
        
        print "Run ML2D"
        protML2D = XmippProtML2D(numberOfReferences=1, maxIters=4, numberOfMpi=1)
        protML2D.inputImages.set(protExtract.outputImages)
        self.proj.launchProtocol(protML2D, wait=True)        
        
        print("******** ML2D ***********")
        if (protML2D.getFiles() is not None):
            print(str(protML2D.getFiles()))  
        else:
            print ("NONEEEEEEEEEEEEEEEEE") 
        
        self.assertIsNotNone(protML2D.outputClassification, "There was a problem with ML2D")  
        
        print "Run CL2D"
        protCL2D = XmippProtCL2D(numberOfReferences=2, numberOfInitialReferences=1, 
                                 numberOfIterations=4, numberOfMpi=1)
        protCL2D.inputImages.set(protExtract.outputImages)
        self.proj.launchProtocol(protCL2D, wait=True)        
        
        print("******** CL2D ***********")
        if (protCL2D.getFiles() is not None):
            print(str(protCL2D.getFiles()))  
        else:
            print ("NONEEEEEEEEEEEEEEEEE") 
        
        self.assertIsNotNone(protCL2D.outputClassification, "There was a problem with CL2D")         
        
if __name__ == "__main__":
    unittest.main()