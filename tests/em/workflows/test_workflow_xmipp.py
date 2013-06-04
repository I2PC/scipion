import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3 import *
    
    
class TestXmippWorkflow(unittest.TestCase):
    
    def setUp(self):    
        # Create a new project
        setupProject(self, 'TestProject_Xmipp')
        self.pattern = getInputPath('Micrographs_BPV3', '*.mrc')        
        self.importFolder = getInputPath('Picking_XmippBPV3')

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
        
        self.assertTrue(protImport.outputMicrographs is not None, "There was a problem with the import")
        
        # Perform a downsampling on the micrographs

        print "Downsampling..."
        protDownsampling = XmippProtPreprocessMicrographs(doDownsample=True, downFactor=3, doCrop=False)
        protDownsampling.inputMicrographs.set(protImport.outputMicrographs)
        self.proj.launchProtocol(protDownsampling, wait=True)
          
        self.assertTrue(protDownsampling.outputMicrographs is not None, "There was a problem with the downsampling")
          
        # Now estimate CTF on the downsampled micrographs    

        print "Performing CTF..."   
        protCTF = XmippProtCTFMicrographs()                
        protCTF.inputMicrographs.set(protDownsampling.outputMicrographs)        
        self.proj.launchProtocol(protCTF, wait=True)
        
        print "Running fake particle picking..."   
        protPP = XmippProtParticlePicking(importFolder=self.importFolder)                
        protPP.inputMicrographs.set(protCTF.outputMicrographs)        
        self.proj.launchProtocol(protPP, wait=True)
            
        self.assertTrue(protPP.outputCoordinates is not None, "There was a problem with the faked picking")
            
        print "Run extract particles with Same as picking"
        protExtract = XmippProtExtractParticles(boxSize=258, downsampleType=1)
        protExtract.inputCoordinates.set(protPP.outputCoordinates)
        #protExtract.inputMicrographs.set(protDownsampling.outputMicrographs)
        self.proj.launchProtocol(protExtract, wait=True)
        
       
        
if __name__ == "__main__":
    unittest.main()