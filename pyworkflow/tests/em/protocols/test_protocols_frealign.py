import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.brandeis import *
from pyworkflow.em.packages.xmipp3 import *
from pyworkflow.em.packages.brandeis import *

class TestXmippBase(unittest.TestCase):
    
    @classmethod
    def runImportMicrograph(cls, pattern, samplingRate, voltage, scannedPixelSize, magnification, sphericalAberration):
        """ Run an Import micrograph protocol. """
        # We have two options: passe the SamplingRate or the ScannedPixelSize + microscope magnification
        if not samplingRate is None:
            cls.protImport = ProtImportMicrographs(samplingRateMode=0, pattern=pattern, samplingRate=samplingRate, magnification=magnification, 
                                                   voltage=voltage, sphericalAberration=sphericalAberration)
        else:
            cls.protImport = ProtImportMicrographs(samplingRateMode=1, pattern=pattern, scannedPixelSize=scannedPixelSize, 
                                                   voltage=voltage, magnification=magnification, sphericalAberration=sphericalAberration)
            
        cls.proj.launchProtocol(cls.protImport, wait=True)
        # check that input micrographs have been imported (a better way to do this?)
        if cls.protImport.outputMicrographs is None:
            raise Exception('Import of micrograph: %s, failed. outputMicrographs is None.' % pattern)
        return cls.protImport
        
    @classmethod
    def runImportMicrographBPV3(cls):
        """ Run an Import micrograph protocol. """
        pattern = getInputPath('Micrographs_BPV3', '*.mrc')
        return cls.runImportMicrograph(pattern, samplingRate=1.237, voltage=300, sphericalAberration=2, scannedPixelSize=None, magnification=56000)


class TestFrealignProtocol(TestXmippBase):

    GOLD_FILES = {'protDownsampling': ['protDownsampling/BPV_1388.mrc', 
                    'protDownsampling/BPV_1387.mrc', 
                    'protImport/BPV_1386.mrc', 
                    'protDownsampling/micrographs.xmd', 
                    'protImport/BPV_1388.mrc', 
                    'protImport/micrographs.sqlite', 
                    'protDownsampling/BPV_1386.mrc', 
                    'protImport/BPV_1387.mrc',
                    'protDownsampling/logs/run.log',
                    'protDownsampling/logs/run.db'],
            'protCTF': [
                  'protCTF/extra/BPV_1386/ctffind_psd.mrc', 
                  'protCTF/extra/BPV_1386/ctffind.out', 
                  'protDownsampling/BPV_1386.mrc',
                  'protCTF/extra/BPV_1387/ctffind_psd.mrc', 
                  'protCTF/extra/BPV_1387/ctffind.out', 
                  'protDownsampling/BPV_1387.mrc',
                  'protCTF/extra/BPV_1388/ctffind_psd.mrc', 
                  'protCTF/extra/BPV_1388/ctffind.out', 
                  'protDownsampling/BPV_1388.mrc',
                  'protCTF/micrographs.sqlite',
                  'protDownsampling/micrographs.xmd', 
                  'protCTF/logs/run.log', 
                  'protCTF/logs/run.db'],
              'protExtract':[
                    'protPP/info/BPV_1386_info.json',
                    'protPP/info/BPV_1387_info.json',
                    'protPP/info/BPV_1388_info.json',
                    'protExtract/extra/scipion_micrographs_coordinates.xmd',
                    'protExtract/images.xmd', 
                    'protExtract/extra/BPV_1386.pos', 
                    'protExtract/extra/BPV_1387.pos', 
                    'protExtract/extra/BPV_1388.pos', 
                    'protExtract/tmp/BPV_1388_flipped.xmp', 
                    'protExtract/tmp/BPV_1387_flipped.xmp', 
                    'protExtract/tmp/BPV_1386_flipped.xmp',
                    'protExtract/tmp/BPV_1386_noDust.xmp', 
                    'protExtract/tmp/BPV_1387_noDust.xmp', 
                    'protExtract/tmp/BPV_1388_noDust.xmp', 
                    'protExtract/extra/BPV_1386.xmd', 
                    'protExtract/extra/BPV_1387.xmd', 
                    'protExtract/extra/BPV_1388.xmd', 
                    'protExtract/extra/BPV_1388.stk', 
                    'protExtract/extra/BPV_1386.stk', 
                    'protExtract/extra/BPV_1387.stk', 
                    'protExtract/logs/run.log',
                    'protExtract/logs/run.db',
                    ],
                }
        
    @classmethod
    def setUpClass(cls):    
        # Create a new project
        setupProject(cls)
        cls.protImport = cls.runImportMicrographBPV3()
        
    def testXmippWorkflow(self):
        #Import a set of volumes        
        print "Import Volume"
        protImportVol = ProtImportVolumes(pattern=getInputPath('Volumes_BPV', 'BPV_scale_filtered_windowed_64.vol'), samplingRate=9.896)
        self.proj.launchProtocol(protImportVol, wait=True)
        
        self.assertIsNotNone(protImportVol.getFiles(), "There was a problem with the import")
#        self.validateFiles('protImportVol', protImportVol)        

        # Perform a downsampling on the micrographs
        print "Downsampling..."
        protDownsampling = XmippProtPreprocessMicrographs(doDownsample=True, downFactor=3, doCrop=False, runMode=1)
        protDownsampling.inputMicrographs.set(self.protImport.outputMicrographs)
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
     
        
        print "Running fake particle picking..."
        protPP = XmippProtParticlePicking(importFolder=getInputPath('Picking_XmippBPV3_Down3'))                
        protPP.inputMicrographs.set(protDownsampling.outputMicrographs)        
        self.proj.launchProtocol(protPP, wait=True)
#        self.protDict['protPicking'] = protPP
            
        self.assertIsNotNone(protPP.outputCoordinates, "There was a problem with the faked picking")
        
        # Extract the SetOfParticles.    
        print "Run extract particles with other downsampling factor"
        protExtract = XmippProtExtractParticles(boxSize=64, downsampleType=2, doFlip=False, downFactor=8, runMode=1, doInvert=False)
        protExtract.inputCoordinates.set(protPP.outputCoordinates)
        protExtract.ctfRelations.set(protCTF.outputCTF)
        protExtract.inputMicrographs.set(self.protImport.outputMicrographs)
        
        self.proj.launchProtocol(protExtract, wait=True)
        
        self.assertIsNotNone(protExtract.outputParticles, "There was a problem with the extract particles")
#         self.validateFiles('protExtract', protExtract)

        # Refine the SetOfParticles and reconstruct a refined volume.    
        print "Running Frealign..."
        protFrealign = ProtFrealign(angularSampling=15, numberOfIterations=2, mode=1, doExtraRealSpaceSym=True,
                                    innerRadius=150, outerRadius=315, symmetry='I2', PhaseResidual=70,
                                    resolution=20, runMode=1, numberOfMpi=1, numberOfThreads=4)
        protFrealign.inputParticles.set(protExtract.outputParticles)
        protFrealign.input3DReferences.set(protImportVol.outputVolumes)

        self.proj.launchProtocol(protFrealign, wait=True)        
        
        self.assertIsNotNone(protFrealign.outputVolume, "There was a problem with Frealign")
#         self.validateFiles('protFrealign', protFrealign)
        
if __name__ == "__main__":
    unittest.main()