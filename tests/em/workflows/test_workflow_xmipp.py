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
        
        self.assertIsNotNone(protImport.outputMicrographs, "There was a problem with the import")
        self.assertTrue(protImport.getFiles() == getImportMicrographFiles())
        
        # Perform a downsampling on the micrographs

        print "Downsampling..."
        protDownsampling = XmippProtPreprocessMicrographs(doDownsample=True, downFactor=3, doCrop=False)
        protDownsampling.inputMicrographs.set(protImport.outputMicrographs)
        self.proj.launchProtocol(protDownsampling, wait=True)
          
        self.assertIsNotNone(protDownsampling.outputMicrographs, "There was a problem with the downsampling")
        self.assertTrue(protDownsampling.getFiles() == getDownsamplingFiles())
          
        # Now estimate CTF on the downsampled micrographs 
        print "Performing CTF..."   
        protCTF = XmippProtCTFMicrographs(numberOfThreads=3)                
        protCTF.inputMicrographs.set(protDownsampling.outputMicrographs)        
        self.proj.launchProtocol(protCTF, wait=True)
        
        # After CTF estimation, the output micrograph should have CTF info
        self.assertTrue(protCTF.outputMicrographs.hasCTF())
        self.assertTrue(protCTF.getFiles() == getCTFFiles())
        
        print "Running fake particle picking..."   
        protPP = XmippProtParticlePicking(importFolder=self.importFolder)                
        protPP.inputMicrographs.set(protCTF.outputMicrographs)        
        self.proj.launchProtocol(protPP, wait=True)
            
        self.assertIsNotNone(protPP.outputCoordinates, "There was a problem with the faked picking")
            
        print "Run extract particles with Same as picking"
        protExtract = XmippProtExtractParticles(boxSize=64, downsampleType=1)
        protExtract.inputCoordinates.set(protPP.outputCoordinates)
        #protExtract.inputMicrographs.set(protDownsampling.outputMicrographs)
        self.proj.launchProtocol(protExtract, wait=True)
        
        self.assertIsNotNone(protExtract.outputImages, "There was a problem with the extract particles")
        self.assertTrue(protExtract.getFiles() == getExtractPrticlesFiles())
        
        print "Run ML2D"
        protML2D = XmippProtML2D(numberOfReferences=1, maxIters=4, 
                                 numberOfMpi=2, numberOfThreads=1)
        protML2D.inputImages.set(protExtract.outputImages)
        self.proj.launchProtocol(protML2D, wait=True)        
        
        self.assertIsNotNone(protML2D.outputClassification, "There was a problem with ML2D")  
        self.assertTrue(protML2D.getFiles() == getML2DFiles())
        
        print "Run CL2D"
        protCL2D = XmippProtCL2D(numberOfReferences=2, numberOfInitialReferences=1, 
                                 numberOfIterations=4, numberOfMpi=2)
        protCL2D.inputImages.set(protExtract.outputImages)
        self.proj.launchProtocol(protCL2D, wait=True)        
        
        self.assertIsNotNone(protCL2D.outputClassification, "There was a problem with CL2D")
        self.assertTrue(protCL2D.getFiles() == getCL2DFiles())         

def getImportMicrographFiles():
    return set(['Runs/ProtImportMicrographs1/BPV_1388.mrc', 
                'Runs/ProtImportMicrographs1/micrographs.sqlite', 
                'Runs/ProtImportMicrographs1/BPV_1387.mrc', 
                'Runs/ProtImportMicrographs1/BPV_1386.mrc'])

def getDownsamplingFiles():
    return set(['Runs/XmippProtPreprocessMicrographs45/BPV_1388.mrc', 
                'Runs/XmippProtPreprocessMicrographs45/BPV_1387.mrc', 
                'Runs/ProtImportMicrographs1/BPV_1386.mrc', 
                'Runs/XmippProtPreprocessMicrographs45/micrographs.xmd', 
                'Runs/ProtImportMicrographs1/BPV_1388.mrc', 
                'Runs/ProtImportMicrographs1/micrographs.sqlite', 
                'Runs/XmippProtPreprocessMicrographs45/log/protocol.log', 
                'Runs/XmippProtPreprocessMicrographs45/BPV_1386.mrc', 
                'Runs/ProtImportMicrographs1/BPV_1387.mrc'])

def getCTFFiles():
    return set(['Runs/XmippProtCTFMicrographs119/extra/BPV_1387/xmipp_ctf_ctfmodel_quadrant.xmp', 
                'Runs/XmippProtCTFMicrographs119/log/protocol.log', 
                'Runs/XmippProtCTFMicrographs119/extra/BPV_1387/xmipp_ctf.psd', 
                'Runs/XmippProtCTFMicrographs119/extra/BPV_1386/xmipp_ctf_enhanced_psd.xmp', 
                'Runs/XmippProtCTFMicrographs119/extra/BPV_1386/xmipp_ctf.ctfparam', 
                'Runs/XmippProtCTFMicrographs119/extra/BPV_1388/xmipp_ctf.ctfparam', 
                'Runs/XmippProtCTFMicrographs119/extra/BPV_1388/xmipp_ctf.psd', 
                'Runs/XmippProtCTFMicrographs119/extra/BPV_1388/xmipp_ctf_ctfmodel_halfplane.xmp', 
                'Runs/XmippProtCTFMicrographs119/extra/BPV_1387/xmipp_ctf.ctfparam', 
                'Runs/XmippProtCTFMicrographs119/extra/BPV_1388/xmipp_ctf_enhanced_psd.xmp', 
                'Runs/XmippProtCTFMicrographs119/tmp/micrographs.xmd', 
                'Runs/XmippProtCTFMicrographs119/extra/BPV_1387/xmipp_ctf_ctfmodel_halfplane.xmp', 
                'Runs/XmippProtCTFMicrographs119/extra/BPV_1387/xmipp_ctf_enhanced_psd.xmp', 
                'Runs/XmippProtCTFMicrographs119/extra/BPV_1386/xmipp_ctf_ctfmodel_quadrant.xmp', 
                'Runs/XmippProtPreprocessMicrographs45/BPV_1388.mrc', 
                'Runs/XmippProtPreprocessMicrographs45/BPV_1387.mrc', 
                'Runs/XmippProtPreprocessMicrographs45/micrographs.xmd', 
                'Runs/XmippProtCTFMicrographs119/extra/BPV_1386/xmipp_ctf.psd', 
                'Runs/XmippProtCTFMicrographs119/extra/BPV_1386/xmipp_ctf_ctfmodel_halfplane.xmp', 
                'Runs/XmippProtCTFMicrographs119/micrographs.xmd', 
                'Runs/XmippProtPreprocessMicrographs45/BPV_1386.mrc', 
                'Runs/XmippProtCTFMicrographs119/extra/BPV_1388/xmipp_ctf_ctfmodel_quadrant.xmp'])
    
def getExtractPrticlesFiles():
    return set(['Runs/XmippProtExtractParticles230/extra/BPV_1387.pos', 
                'Runs/XmippProtExtractParticles230/tmp/BPV_1388_flipped.xmp', 
                'Runs/XmippProtExtractParticles230/tmp/BPV_1387_flipped.xmp', 
                'Runs/XmippProtExtractParticles230/tmp/BPV_1386_noDust.xmp', 
                'Runs/XmippProtExtractParticles230/extra/BPV_1388.pos', 
                'Runs/XmippProtParticlePicking189/extra/BPV_1387.pos', 
                'Runs/XmippProtParticlePicking189/extra/BPV_1386.pos', 
                'Runs/XmippProtExtractParticles230/extra/BPV_1386.xmd', 
                'Runs/XmippProtExtractParticles230/extra/BPV_1388.stk', 
                'Runs/XmippProtParticlePicking189/extra/BPV_1388.pos', 
                'Runs/XmippProtExtractParticles230/images.xmd', 
                'Runs/XmippProtExtractParticles230/extra/BPV_1386.stk', 
                'Runs/XmippProtExtractParticles230/extra/BPV_1388.xmd', 
                'Runs/XmippProtExtractParticles230/extra/BPV_1386.pos', 
                'Runs/XmippProtExtractParticles230/extra/BPV_1387.stk', 
                'Runs/XmippProtExtractParticles230/tmp/BPV_1387_noDust.xmp', 
                'Runs/XmippProtExtractParticles230/tmp/BPV_1388_noDust.xmp', 
                'Runs/XmippProtExtractParticles230/log/protocol.log', 
                'Runs/XmippProtExtractParticles230/extra/BPV_1387.xmd', 
                'Runs/XmippProtExtractParticles230/tmp/BPV_1386_flipped.xmp'])

def getML2DFiles():
    return set(['Runs/XmippProtML2D366/ml2d_extra/iter002/iter_classes.xmd', 
                'Runs/XmippProtML2D366/ml2d_extra/iter001/iter_classes.xmd', 
                'Runs/XmippProtML2D366/ml2d_extra/iter002/iter_images.xmd', 
                'Runs/XmippProtML2D366/ml2d_classes.stk', 
                'Runs/XmippProtML2D366/ml2d_extra/iter003/iter_images.xmd', 
                'Runs/XmippProtML2D366/ml2d_images.xmd', 
                'Runs/XmippProtCTFMicrographs119/extra/BPV_1386/xmipp_ctf.ctfparam', 
                'Runs/XmippProtML2D366/ml2d_extra/iter004/iter_classes.stk', 
                'Runs/XmippProtCTFMicrographs119/extra/BPV_1388/xmipp_ctf.ctfparam', 
                'Runs/XmippProtML2D366/ml2d__images_average.xmp', 
                'Runs/XmippProtExtractParticles230/extra/BPV_1388.stk', 
                'Runs/XmippProtCTFMicrographs119/extra/BPV_1387/xmipp_ctf.ctfparam', 
                'Runs/XmippProtML2D366/ml2d_extra/iter001/iter_images.xmd', 
                'Runs/XmippProtML2D366/ml2d_extra/iter002/iter_classes.stk', 
                'Runs/XmippProtExtractParticles230/images.xmd', 
                'Runs/XmippProtML2D366/ml2d_extra/iter004/iter_classes.xmd', 
                'Runs/XmippProtExtractParticles230/extra/BPV_1386.stk', 
                'Runs/XmippProtML2D366/ml2d_extra/iter004/iter_images.xmd', 
                'Runs/XmippProtML2D366/ml2d_extra/iter003/iter_classes.xmd', 
                'Runs/XmippProtML2D366/log/protocol.log', 
                'Runs/XmippProtExtractParticles230/extra/BPV_1387.stk', 
                'Runs/XmippProtML2D366/ml2d_classes.xmd', 
                'Runs/XmippProtML2D366/ml2d_extra/iter001/iter_classes.stk', 
                'Runs/XmippProtML2D366/ml2d_extra/iter003/iter_classes.stk'])

def getCL2DFiles():
    return set(['Runs/XmippProtCL2D420/extra/classes_core_hierarchy.txt', 
                'Runs/XmippProtCL2D420/extra/level_01/level_classes_core.xmd', 
                'Runs/XmippProtCL2D420/extra/level_01/level_classes.stk', 
                'Runs/XmippProtCL2D420/extra/level_01/classes_sorted.stk', 
                'Runs/XmippProtCL2D420/extra/level_00/level_classes.xmd', 
                'Runs/XmippProtCL2D420/extra/level_00/classes_sorted.xmd', 
                'Runs/XmippProtCTFMicrographs119/extra/BPV_1386/xmipp_ctf.ctfparam', 
                'Runs/XmippProtCTFMicrographs119/extra/BPV_1388/xmipp_ctf.ctfparam', 
                'Runs/XmippProtCL2D420/extra/level_00/level_classes_core.stk', 
                'Runs/XmippProtCL2D420/log/protocol.log', 
                'Runs/XmippProtExtractParticles230/extra/BPV_1388.stk', 
                'Runs/XmippProtCTFMicrographs119/extra/BPV_1387/xmipp_ctf.ctfparam', 
                'Runs/XmippProtCL2D420/extra/level_00/classes_core_sorted.xmd', 
                'Runs/XmippProtCL2D420/extra/level_01/classes_core_sorted.stk', 
                'Runs/XmippProtExtractParticles230/images.xmd', 
                'Runs/XmippProtCL2D420/extra/level_00/classes_sorted.stk', 
                'Runs/XmippProtCL2D420/extra/level_01/classes_sorted.xmd', 
                'Runs/XmippProtCL2D420/extra/classes_hierarchy.txt', 
                'Runs/XmippProtCL2D420/extra/images.xmd', 
                'Runs/XmippProtExtractParticles230/extra/BPV_1386.stk', 
                'Runs/XmippProtCL2D420/extra/level_01/classes_core_sorted.xmd', 
                'Runs/XmippProtCL2D420/extra/level_00/classes_core_sorted.stk', 
                'Runs/XmippProtExtractParticles230/extra/BPV_1387.stk', 
                'Runs/XmippProtCL2D420/extra/level_00/level_classes_core.xmd', 
                'Runs/XmippProtCL2D420/extra/level_01/level_classes.xmd', 
                'Runs/XmippProtCL2D420/extra/level_00/level_classes.stk', 
                'Runs/XmippProtCL2D420/extra/level_01/level_classes_core.stk'])


if __name__ == "__main__":
    unittest.main()
