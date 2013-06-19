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

    def validateFiles(self, key, prot):
        """ Validate if the produced files are the expected ones.
        Params:
            prot: the protocol to validate. 
            filesSet: the known files that should be produced (set)
        """
        self.protDict[key] = prot
        filesSet = self.getProtocolFiles(key)
        self.assertEqual(prot.getFiles(), filesSet)
        
    def printSet(self, msg, s):
        print "============= %s ==========" % msg
        for i in s:
            print i
            
    def testXmippWorkflow(self):
        self.protDict = {}
        #First, import a set of micrographs
        protImport = ProtImportMicrographs(pattern=self.pattern, samplingRate=1.237, voltage=300)
        self.proj.launchProtocol(protImport, wait=True)
        
        self.assertIsNotNone(protImport.outputMicrographs, "There was a problem with the import")
        self.validateFiles('protImport', protImport)        
        
        # Perform a downsampling on the micrographs

        print "Downsampling..."
        protDownsampling = XmippProtPreprocessMicrographs(doDownsample=True, downFactor=3, doCrop=False)
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
        protExtract = XmippProtExtractParticles(boxSize=64, downsampleType=2, downFactor=8, runMode=1)
        protExtract.inputCoordinates.set(protPP.outputCoordinates)
        protExtract.inputMicrographs.set(protImport.outputMicrographs)
        self.proj.launchProtocol(protExtract, wait=True)
        
        self.assertIsNotNone(protExtract.outputImages, "There was a problem with the extract particles")
        self.printSet('reference', self.getProtocolFiles('protExtract'))
        self.printSet('current', protExtract.getFiles())
        self.validateFiles('protExtract', protExtract)
        
        print "Run ML2D"
        protML2D = XmippProtML2D(numberOfReferences=1, maxIters=4, 
                                 numberOfMpi=2, numberOfThreads=1)
        protML2D.inputImages.set(protExtract.outputImages)
        self.proj.launchProtocol(protML2D, wait=True)        
        
        self.assertIsNotNone(protML2D.outputClassification, "There was a problem with ML2D")  
        self.validateFiles('protML2D', protML2D)
        
        print "Run CL2D"
        protCL2D = XmippProtCL2D(numberOfReferences=2, numberOfInitialReferences=1, 
                                 numberOfIterations=4, numberOfMpi=2)
        protCL2D.inputImages.set(protExtract.outputImages)
        self.proj.launchProtocol(protCL2D, wait=True)        
        
        self.assertIsNotNone(protCL2D.outputClassification, "There was a problem with CL2D")
        self.validateFiles('protCL2D', protCL2D) 

    def getProtocolFiles(self, key):
        fileList = GOLD_FILES[key]
        fileSet = set([self.__replaceFilename(f) for f in fileList])
        
        return fileSet
    
    def __replaceFilename(self, filename):
        """ Convert list to set and replace the key
        in the filename by the protocol working dir. 
        """
        for k, v in self.protDict.iteritems():
            if filename.startswith(k):
                return filename.replace(k, v.getWorkingDir())
        return filename
                
    
GOLD_FILES = {'protImport': ['protImport/BPV_1388.mrc',
                    'protImport/micrographs.sqlite', 
                    'protImport/BPV_1387.mrc',
                    'protImport/BPV_1386.mrc'],
              'protDownsampling': ['protDownsampling/BPV_1388.mrc', 
                    'protDownsampling/BPV_1387.mrc', 
                    'protImport/BPV_1386.mrc', 
                    'protDownsampling/micrographs.xmd', 
                    'protImport/BPV_1388.mrc', 
                    'protImport/micrographs.sqlite', 
                    'protDownsampling/log/protocol.log', 
                    'protDownsampling/BPV_1386.mrc', 
                    'protImport/BPV_1387.mrc'],
              'protCTF': ['protCTF/extra/BPV_1387/xmipp_ctf_ctfmodel_quadrant.xmp', 
                    'protCTF/log/protocol.log', 
                    'protCTF/extra/BPV_1387/xmipp_ctf.psd', 
                    'protCTF/extra/BPV_1386/xmipp_ctf_enhanced_psd.xmp', 
                    'protCTF/extra/BPV_1386/xmipp_ctf.ctfparam', 
                    'protCTF/extra/BPV_1388/xmipp_ctf.ctfparam', 
                    'protCTF/extra/BPV_1388/xmipp_ctf.psd', 
                    'protCTF/extra/BPV_1388/xmipp_ctf_ctfmodel_halfplane.xmp', 
                    'protCTF/extra/BPV_1387/xmipp_ctf.ctfparam', 
                    'protCTF/extra/BPV_1388/xmipp_ctf_enhanced_psd.xmp', 
                    'protCTF/tmp/micrographs.xmd', 
                    'protCTF/extra/BPV_1387/xmipp_ctf_ctfmodel_halfplane.xmp', 
                    'protCTF/extra/BPV_1387/xmipp_ctf_enhanced_psd.xmp', 
                    'protCTF/extra/BPV_1386/xmipp_ctf_ctfmodel_quadrant.xmp', 
                    'protDownsampling/BPV_1388.mrc', 
                    'protDownsampling/BPV_1387.mrc', 
                    'protDownsampling/micrographs.xmd', 
                    'protCTF/extra/BPV_1386/xmipp_ctf.psd', 
                    'protCTF/extra/BPV_1386/xmipp_ctf_ctfmodel_halfplane.xmp', 
                    'protCTF/micrographs.xmd', 
                    'protDownsampling/BPV_1386.mrc', 
                    'protCTF/extra/BPV_1388/xmipp_ctf_ctfmodel_quadrant.xmp'],
              'protExtract':[
                    'protImport/BPV_1386.mrc',
                    'protImport/BPV_1388.mrc',
                    'protImport/BPV_1387.mrc',
                    'protImport/micrographs.sqlite',
                    'protPicking/extra/BPV_1387.pos', 
                    'protPicking/extra/BPV_1386.pos', 
                    'protPicking/extra/BPV_1388.pos', 
                    'protExtract/extra/BPV_1387.pos', 
                    'protExtract/tmp/BPV_1388_flipped.xmp', 
                    'protExtract/tmp/BPV_1387_flipped.xmp', 
                    'protExtract/tmp/BPV_1386_noDust.xmp', 
                    'protExtract/extra/BPV_1388.pos', 
                    'protExtract/extra/BPV_1386.xmd', 
                    'protExtract/extra/BPV_1388.stk', 
                    'protExtract/images.xmd', 
                    'protExtract/extra/BPV_1386.stk', 
                    'protExtract/extra/BPV_1388.xmd', 
                    'protExtract/extra/BPV_1386.pos', 
                    'protExtract/extra/BPV_1387.stk', 
                    'protExtract/tmp/BPV_1387_noDust.xmp', 
                    'protExtract/tmp/BPV_1388_noDust.xmp', 
                    'protExtract/extra/BPV_1387.xmd', 
                    'protExtract/tmp/BPV_1386_flipped.xmp',
                    'protExtract/tmp/BPV_1387_downsampled.xmp',
                    'protExtract/tmp/BPV_1386_downsampled.xmp',
                    'protExtract/tmp/BPV_1388_downsampled.xmp',
                    ],
              'protML2D': [
                    'protExtract/extra/BPV_1386.stk', 
                    'protExtract/extra/BPV_1387.stk', 
                    'protExtract/extra/BPV_1388.stk', 
                    'protExtract/images.xmd', 
                    'protCTF/extra/BPV_1386/xmipp_ctf.ctfparam',
                    'protCTF/extra/BPV_1387/xmipp_ctf.ctfparam',
                    'protCTF/extra/BPV_1388/xmipp_ctf.ctfparam',
                    'protML2D/ml2d_extra/iter002/iter_classes.xmd', 
                    'protML2D/ml2d_extra/iter001/iter_classes.xmd', 
                    'protML2D/ml2d_extra/iter002/iter_images.xmd', 
                    'protML2D/ml2d_classes.stk', 
                    'protML2D/ml2d_extra/iter003/iter_images.xmd', 
                    'protML2D/ml2d_images.xmd', 
                    'protML2D/ml2d_extra/iter004/iter_classes.stk', 
                    'protML2D/ml2d__images_average.xmp', 
                    'protML2D/ml2d_extra/iter001/iter_images.xmd', 
                    'protML2D/ml2d_extra/iter002/iter_classes.stk', 
                    'protML2D/ml2d_extra/iter004/iter_classes.xmd', 
                    'protML2D/ml2d_extra/iter004/iter_images.xmd', 
                    'protML2D/ml2d_extra/iter003/iter_classes.xmd', 
                    'protML2D/log/protocol.log', 
                    'protML2D/ml2d_classes.xmd', 
                    'protML2D/ml2d_extra/iter001/iter_classes.stk', 
                    'protML2D/ml2d_extra/iter003/iter_classes.stk'],
              'protCL2D': ['protCL2D/extra/classes_core_hierarchy.txt', 
                    'protCL2D/extra/level_01/level_classes_core.xmd', 
                    'protCL2D/extra/level_01/level_classes.stk', 
                    'protCL2D/extra/level_01/classes_sorted.stk', 
                    'protCL2D/extra/level_00/level_classes.xmd', 
                    'protCL2D/extra/level_00/classes_sorted.xmd', 
                    'protCTF/extra/BPV_1386/xmipp_ctf.ctfparam', 
                    'protCTF/extra/BPV_1388/xmipp_ctf.ctfparam', 
                    'protCL2D/extra/level_00/level_classes_core.stk', 
                    'protCL2D/log/protocol.log', 
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
                    'protCL2D/extra/level_01/level_classes_core.stk']
              }

if __name__ == "__main__":
    unittest.main()
