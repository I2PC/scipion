import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3 import *
from pyworkflow.em.packages.xmipp3.constants import SAME_AS_PICKING 
from pyworkflow.em.packages.grigoriefflab import *
from pyworkflow.em.packages.eman2 import *
from pyworkflow.em.packages.relion import *
from test_workflow import TestWorkflow



class TestMixedRelionTutorial(TestWorkflow):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('relion_tutorial')
        cls.crdsEmanDir = cls.dataset.getFile('boxingDir')
        cls.crdsXmippDir = cls.dataset.getFile('posAllDir')
        cls.micsFn = cls.dataset.getFile('allMics')
        cls.vol = cls.dataset.getFile('volume')
    
    def test_workflow(self):
        #First, import a set of micrographs
        print "Importing a set of micrographs..."
        protImport = self.newProtocol(ProtImportMicrographs, filesPath=self.micsFn, samplingRateMode=1, magnification=79096,
                                           scannedPixelSize=56, voltage=300, sphericalAberration=2.0)
        protImport.setObjLabel('import 20 mics')
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputMicrographs, "There was a problem with the import")
        
        print "Importing a volume..."
        protImportVol = self.newProtocol(ProtImportVolumes, filesPath=self.vol, samplingRate=7.08)
        protImportVol.setObjLabel('import single vol')
        self.launchProtocol(protImportVol)
        self.assertIsNotNone(protImportVol.outputVolume, "There was a problem with the import")
        
        print "Preprocessing the micrographs..."
        protPreprocess = self.newProtocol(XmippProtPreprocessMicrographs, doCrop=True, cropPixels=25)
        protPreprocess.inputMicrographs.set(protImport.outputMicrographs)
        protPreprocess.setObjLabel('crop 50px')
        self.launchProtocol(protPreprocess)
        self.assertIsNotNone(protPreprocess.outputMicrographs, "There was a problem with the downsampling")
        
        print "Running Eman import coordinates..."
        protPP = self.newProtocol(ProtImportCoordinates,
                                 importFrom=ProtImportCoordinates.IMPORT_FROM_EMAN,
                                 filesPath=self.crdsEmanDir,
                                 filesPattern='info/*_info.json', boxSize=60)
        protPP.inputMicrographs.set(protPreprocess.outputMicrographs)
        protPP.setObjLabel('import from Eman boxing')
        self.launchProtocol(protPP)
        self.assertIsNotNone(protPP.outputCoordinates, "There was a problem with the Eman import coordinates")

        # Now estimate CTF on the micrographs with ctffind 
        print "Performing CTFfind..."   
        protCTF = self.newProtocol(ProtCTFFind, lowRes=0.04, highRes=0.45,
                                   minDefocus=1.2, maxDefocus=3,
                                   runMode=1, numberOfMpi=1, numberOfThreads=16)
        protCTF.inputMicrographs.set(protPreprocess.outputMicrographs)
        protCTF.setObjLabel('CTF ctffind')
        self.launchProtocol(protCTF)

        print "Run extract particles with <Same as picking> option"
        protExtract = self.newProtocol(XmippProtExtractParticles, boxSize=60,
                                       downsampleType=SAME_AS_PICKING,
                                       doRemoveDust=False,
                                       doInvert=False,
                                       doFlip=False, backRadius=28, runMode=1)
        protExtract.inputCoordinates.set(protPP.outputCoordinates)
        protExtract.ctfRelations.set(protCTF.outputCTF)
        protExtract.setObjLabel('Extract particles')
        self.launchProtocol(protExtract)
        self.assertIsNotNone(protExtract.outputParticles,
                             "There was a problem with the extract particles")
        
        print "Run CL2D"
        protCL2D = XmippProtCL2D(numberOfClasses=32, numberOfInitialClasses=4,
                                 numberOfIterations=2, numberOfMpi=16)
        protCL2D.inputParticles.set(protExtract.outputParticles)
        protCL2D.setObjLabel('CL2D')
        self.launchProtocol(protCL2D)   
        self.assertIsNotNone(protCL2D.outputClasses,
                             "There was a problem with CL2D")
        
        # Now estimate CTF on the micrographs with xmipp
        print "Performing Xmipp CTF..."   
        protCTF2 = self.newProtocol(XmippProtCTFMicrographs, lowRes=0.04,
                                    highRes=0.45, minDefocus=1.2, maxDefocus=3,
                                    runMode=1, numberOfMpi=1, numberOfThreads=16)
        protCTF2.inputMicrographs.set(protPreprocess.outputMicrographs)
        protCTF2.setObjLabel('CTF xmipp')
        self.launchProtocol(protCTF2)
        
        print "Running Xmipp Import Coordinates"
        protPP2 = self.newProtocol(ProtImportCoordinates,
                                 importFrom=ProtImportCoordinates.IMPORT_FROM_XMIPP,
                                 filesPath=self.crdsXmippDir,
                                 filesPattern='*.pos', boxSize=60)
        protPP2.inputMicrographs.set(protPreprocess.outputMicrographs)  
        protPP2.setObjLabel('Xmipp Import Coords')
        self.launchProtocol(protPP2)
        self.assertIsNotNone(protPP2.outputCoordinates,
                             "There was a problem with the Xmipp import coordinates")
        
        print "Run extract particles with <Same as picking> option"
        protExtract2 = self.newProtocol(XmippProtExtractParticles, boxSize=60,
                                        downsampleType=SAME_AS_PICKING,
                                        doRemoveDust=False, doInvert=True,
                                        doFlip=False, backRadius=28, runMode=1)
        protExtract2.inputCoordinates.set(protPP2.outputCoordinates)
        protExtract2.ctfRelations.set(protCTF2.outputCTF)
        protExtract2.setObjLabel('Extract particles')
        self.launchProtocol(protExtract2)
        self.assertIsNotNone(protExtract2.outputParticles,
                             "There was a problem with the extract particles")
        
        print "Run Relion Classification2d"
        prot2D = ProtRelionClassify2D(regularisationParamT=2,
                                      numberOfMpi=4, numberOfThreads=4)
        prot2D.numberOfClasses.set(50)
        prot2D.numberOfIterations.set(25)
        prot2D.inputParticles.set(protExtract2.outputParticles)
        prot2D.setObjLabel('relion 2D')
        self.launchProtocol(prot2D)        
        self.assertIsNotNone(prot2D.outputClasses, "There was a problem with Relion 2D:\n" + (prot2D.getErrorMessage() or "No error set"))
        
        print "Run Relion Refine"
        proRef = ProtRelionRefine3D(numberOfMpi=4, numberOfThreads=4)
        proRef.inputParticles.set(protExtract2.outputParticles)
        proRef.referenceVolume.set(protImportVol.outputVolume)
        proRef.setObjLabel('relion Refine')
        self.launchProtocol(proRef)        
        self.assertIsNotNone(proRef.outputVolume, "There was a problem with Relion Refine:\n" + (proRef.getErrorMessage() or "No error set"))
        self.assertIsNotNone(proRef.outputParticles, "There was a problem with Relion Refine:\n" + (proRef.getErrorMessage() or "No error set"))


class TestMixedFrealignClassify(TestWorkflow):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('relion_tutorial')
        cls.crdsEmanDir = cls.dataset.getFile('boxingDir')
        cls.crdsXmippDir = cls.dataset.getFile('posAllDir')
        cls.micsFn = cls.dataset.getFile('allMics')
        cls.vol = cls.dataset.getFile('volume')
    
    def test_workflow(self):
        #First, import a set of micrographs
        print "Importing a set of micrographs..."
        #TODO missing contrst amplitude ROB
        protImport = self.newProtocol(ProtImportMicrographs,
                                      filesPath=self.micsFn,
                                      samplingRateMode=1, magnification=79096,
                                      scannedPixelSize=56, voltage=300,
                                      sphericalAberration=2.0)
        protImport.setObjLabel('import 20 mics')
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputMicrographs,
                             "There was a problem with the import")
        
        print "Importing a volume..."
        protImportVol = self.newProtocol(ProtImportVolumes,
                                         filesPath=self.vol, samplingRate=7.08)
        protImportVol.setObjLabel('import single vol')
        self.launchProtocol(protImportVol)
        self.assertIsNotNone(protImportVol.outputVolume,
                             "There was a problem with the import")
        
        print "Preprocessing the micrographs..."
        protPreprocess = self.newProtocol(XmippProtPreprocessMicrographs,
                                          doCrop=True, cropPixels=25)
        protPreprocess.inputMicrographs.set(protImport.outputMicrographs)
        protPreprocess.setObjLabel('crop 50px')
        self.launchProtocol(protPreprocess)
        self.assertIsNotNone(protPreprocess.outputMicrographs,
                             "There was a problem with the downsampling")

        print "Running Eman import coordinates..."
        protPP = self.newProtocol(ProtImportCoordinates,
                                 importFrom=ProtImportCoordinates.IMPORT_FROM_EMAN,
                                 filesPath=self.crdsEmanDir,
                                 filesPattern='info/*_info.json', boxSize=60)
        protPP.inputMicrographs.set(protPreprocess.outputMicrographs)
        protPP.setObjLabel('import from Eman boxing')
        self.launchProtocol(protPP)
        self.assertIsNotNone(protPP.outputCoordinates,
                             "There was a problem with the Eman import coordinates")

        # Now estimate CTF on the micrographs with ctffind 
        print "Performing CTFfind..."   
        protCTF = self.newProtocol(ProtCTFFind, lowRes=0.04, highRes=0.45,
                                   minDefocus=1.2, maxDefocus=3,
                                   runMode=1, numberOfMpi=1, numberOfThreads=16)
        protCTF.inputMicrographs.set(protPreprocess.outputMicrographs)
        protCTF.setObjLabel('CTF ctffind')
        self.launchProtocol(protCTF)
        
        print "Run extract particles with <Same as picking> option"
        protExtract = self.newProtocol(XmippProtExtractParticles, boxSize=60,
                                       downsampleType=SAME_AS_PICKING,
                                       doRemoveDust=False,
                                       doInvert=False,
                                       doFlip=False, backRadius=28, runMode=1)
        protExtract.inputCoordinates.set(protPP.outputCoordinates)
        protExtract.ctfRelations.set(protCTF.outputCTF)
        protExtract.setObjLabel('Extract particles')
        self.launchProtocol(protExtract)
        self.assertIsNotNone(protExtract.outputParticles,
                             "There was a problem with the extract particles")
        
        # Classify the SetOfParticles.
        print "Running Frealign Classification..."
        protFrealign = self.newProtocol(ProtFrealignClassify, doInvert=False,
                                        numberOfClasses=3, itRefineAngles=2,
                                        itRefineShifts=3, angStepSize=20,
                                        numberOfIterations=6, mode=1,
                                        doExtraRealSpaceSym=True,
                                        outerRadius=180, PhaseResidual=65,
                                        lowResolRefine=300, highResolRefine=15,
                                        resolution=15, runMode=1,
                                        numberOfMpi=16)
        protFrealign.inputParticles.set(protExtract.outputParticles)
        protFrealign.input3DReference.set(protImportVol.outputVolume)
        protFrealign.setObjLabel('Frealign')
        self.launchProtocol(protFrealign)        
        self.assertIsNotNone(protFrealign.outputClasses,
                             "There was a problem with Frealign")

