import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3 import *
    
"""
The main goal of this project is to test 
execution of remote protocols. 
"""
 
class TestXmippWorkflow(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('xmipp_tutorial')
        cls.allCrdsDir = cls.dataset.getFile('posAllDir')
        cls.micsFn = cls.dataset.getFile('allMics')
        cls.vol1 = cls.dataset.getFile('vol1')
    
    def testXmippWorkflow(self):
        #First, import a set of micrographs
        protImport = ProtImportMicrographs(pattern=self.micsFn, samplingRate=1.237, voltage=300)
        self.proj.launchProtocol(protImport, wait=True)
        self.assertIsNotNone(protImport.outputMicrographs.getFileName(), "There was a problem with the import")
        self.validateFiles('protImport', protImport)      

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
        self.validateFiles('protDownsampling', protDownsampling)
     
        # Now estimate CTF on the downsampled micrographs 
        print "Performing CTF..."   
        protCTF = XmippProtCTFMicrographs(numberOfThreads=3)                
        protCTF.inputMicrographs.set(protDownsampling.outputMicrographs)        
        self.proj.launchProtocol(protCTF, wait=True)
        self.assertIsNotNone(protCTF.outputCTF, "There was a problem with the CTF estimation")
        # After CTF estimation, the output micrograph should have CTF info
        self.validateFiles('protCTF', protCTF)
        
        print "Running fake particle picking..."   
        protPP = XmippProtParticlePicking(importFolder=self.allCrdsDir)                
        protPP.inputMicrographs.set(protDownsampling.outputMicrographs)        
        self.proj.launchProtocol(protPP, wait=True)
        self.protDict['protPicking'] = protPP
        self.assertIsNotNone(protPP.outputCoordinates, "There was a problem with the faked picking")
            
        print "Run extract particles with other downsampling factor"
        protExtract = XmippProtExtractParticles(boxSize=64, downsampleType=2, doFlip=True, downFactor=8, runMode=1, doInvert=True)
        protExtract.inputCoordinates.set(protPP.outputCoordinates)
        protExtract.ctfRelations.set(protCTF.outputCTF)
        protExtract.inputMicrographs.set(protImport.outputMicrographs)
        self.proj.launchProtocol(protExtract, wait=True)
        self.assertIsNotNone(protExtract.outputParticles, "There was a problem with the extract particles")
        self.validateFiles('protExtract', protExtract)


if __name__ == "__main__":
    unittest.main()