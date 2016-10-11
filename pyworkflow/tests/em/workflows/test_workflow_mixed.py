from pyworkflow.tests import DataSet, setupTestProject
from pyworkflow.em.packages.xmipp3 import (XmippProtPreprocessMicrographs,
                                           XmippProtExtractParticles,
                                           XmippProtCropResizeParticles, XmippProtML2D)
from pyworkflow.em.packages.xmipp3.constants import SAME_AS_PICKING, OTHER 
from pyworkflow.em.packages.grigoriefflab import ProtCTFFind, ProtFrealign
from pyworkflow.em.packages.eman2 import EmanProtInitModel
from pyworkflow.em.packages.relion import ProtImportMicrographs, ProtImportVolumes
from test_workflow import TestWorkflow
from pyworkflow.em.protocol.protocol_import import ProtImportCoordinates


class TestMixedBPV(TestWorkflow):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('xmipp_tutorial')
        cls.crdsDir = cls.dataset.getFile('boxingDir')
        cls.micsFn = cls.dataset.getFile('allMics')
        cls.vol1 = cls.dataset.getFile('vol1')
        
    def test_workflow(self):
        from itertools import izip
        #First, import a set of micrographs
        protImport = self.newProtocol(ProtImportMicrographs,
                                      filesPath=self.micsFn,
                                      samplingRate=1.237, voltage=300)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputMicrographs,
                             "There was a problem with the import")
#         self.validateFiles('protImport', protImport) 
        
        #Import a set of volumes        
        print "Import Volume"
        protImportVol = self.newProtocol(ProtImportVolumes, filesPath=self.vol1,
                                         samplingRate=9.896)
        self.launchProtocol(protImportVol)
        self.assertIsNotNone(protImportVol.getFiles(),
                             "There was a problem with the import")
#        self.validateFiles('protImportVol', protImportVol)
        
        # Perform a downsampling on the micrographs
        print "Downsampling..."
        protDownsampling = self.newProtocol(XmippProtPreprocessMicrographs,
                                            doDownsample=True, downFactor=5,
                                            doCrop=False, runMode=1)
        protDownsampling.inputMicrographs.set(protImport.outputMicrographs)
        self.launchProtocol(protDownsampling)
        self.assertIsNotNone(protDownsampling.outputMicrographs,
                             "There was a problem with the downsampling")
#         self.validateFiles('protDownsampling', protDownsampling)
        
        # Estimate CTF on the downsampled micrographs
        print "Performing CTFfind..."   
        protCTF = self.newProtocol(ProtCTFFind, numberOfThreads=4,
                                   minDefocus=2.2, maxDefocus=2.5)
        protCTF.inputMicrographs.set(protDownsampling.outputMicrographs)        
        self.launchProtocol(protCTF)
        self.assertIsNotNone(protCTF.outputCTF,
                             "There was a problem with the CTF estimation")
        
        valuesList = [[24000, 24000], [22548, 22518], [23058, 23029]]
        for ctfModel, values in izip(protCTF.outputCTF, valuesList):
            self.assertAlmostEquals(ctfModel.getDefocusU(),values[0], delta=1000)
            self.assertAlmostEquals(ctfModel.getDefocusV(),values[1], delta=1000)
            self.assertAlmostEquals(ctfModel.getMicrograph().getSamplingRate(),
                                    6.185, delta=0.001)

#         self.validateFiles('protCTF', protCTF)

        print "Running Eman import coordinates..."
        protPP = self.newProtocol(ProtImportCoordinates,
                                 importFrom=ProtImportCoordinates.IMPORT_FROM_EMAN,
                                 filesPath=self.crdsDir,
                                 filesPattern='*_info.json', boxSize=110)
        protPP.inputMicrographs.set(protDownsampling.outputMicrographs)
        self.launchProtocol(protPP)
        self.assertIsNotNone(protPP.outputCoordinates,
                             "There was a problem with the Eman import coordinates")


        # Extract the SetOfParticles.
        print "Run extract particles with other downsampling factor"
        protExtract = self.newProtocol(XmippProtExtractParticles,
                                       boxSize=64,
                                       downsampleType=OTHER,
                                       downFactor=8.0,
                                       doFlip=False, runMode=1, doInvert=False)
        protExtract.inputCoordinates.set(protPP.outputCoordinates)
        protExtract.ctfRelations.set(protCTF.outputCTF)
        protExtract.inputMicrographs.set(protImport.outputMicrographs)
        self.launchProtocol(protExtract)
        self.assertIsNotNone(protExtract.outputParticles,
                             "There was a problem with the extract particles")
#         self.validateFiles('protExtract', protExtract)
        
        # Refine the SetOfParticles and reconstruct a refined volume.
        print "Running Frealign..."
        protFrealign = self.newProtocol(ProtFrealign, doInvert=False,
                                        angStepSize=15,
                                        numberOfIterations=2, mode=1,
                                        doExtraRealSpaceSym=True,
                                    innerRadius=130, outerRadius=300,
                                    symmetry='I1', PhaseResidual=30,
                                    molMass=19400,
                                    score=5, resolution=20,
                                    runMode=1, numberOfThreads=4)
        protFrealign.inputParticles.set(protExtract.outputParticles)
        protFrealign.input3DReference.set(protImportVol.outputVolume)
        self.launchProtocol(protFrealign)
        self.assertIsNotNone(protFrealign.outputVolume,
                             "There was a problem with Frealign")
#         self.validateFiles('protFrealign', protFrealign)


class TestMixedBPV2(TestWorkflow):



    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('xmipp_tutorial')
        cls.crdsDir = cls.dataset.getFile('boxingDir')
        cls.micsFn = cls.dataset.getFile('allMics')
        cls.vol1 = cls.dataset.getFile('vol1')
        
    def test_workflow(self):
        #First, import a set of micrographs
        protImport = self.newProtocol(ProtImportMicrographs,
                                      filesPath=self.micsFn,
                                      samplingRate=1.237, voltage=300)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputMicrographs,
                             "There was a problem with the import")
#         self.validateFiles('protImport', protImport) 
        
        #Import a set of volumes        
        print "Import Volume"
        protImportVol = self.newProtocol(ProtImportVolumes, filesPath=self.vol1,
                                         samplingRate=9.896)
        self.launchProtocol(protImportVol)
        self.assertIsNotNone(protImportVol.getFiles(),
                             "There was a problem with the import")
#        self.validateFiles('protImportVol', protImportVol)
        
        # Perform a downsampling on the micrographs
        print "Downsampling..."
        protDownsampling = self.newProtocol(XmippProtPreprocessMicrographs,
                                            doDownsample=True, downFactor=5,
                                            doCrop=False, runMode=1)
        protDownsampling.inputMicrographs.set(protImport.outputMicrographs)
        self.launchProtocol(protDownsampling)
        self.assertIsNotNone(protDownsampling.outputMicrographs,
                             "There was a problem with the downsampling")
#         self.validateFiles('protDownsampling', protDownsampling)
        
        # Estimate CTF on the downsampled micrographs
        print "Performing CTFfind..."   
        protCTF = self.newProtocol(ProtCTFFind, numberOfThreads=4,
                                   minDefocus=2.2, maxDefocus=2.5)
        protCTF.inputMicrographs.set(protImport.outputMicrographs)        
        self.launchProtocol(protCTF)
        self.assertIsNotNone(protCTF.outputCTF,
                             "There was a problem with the CTF estimation")
#         self.validateFiles('protCTF', protCTF)

        print "Running Eman import coordinates..."
        protPP = self.newProtocol(ProtImportCoordinates,
                                 importFrom=ProtImportCoordinates.IMPORT_FROM_EMAN,
                                 filesPath=self.crdsDir,
                                 filesPattern='*_info.json', boxSize=110)
        protPP.inputMicrographs.set(protDownsampling.outputMicrographs)
        self.launchProtocol(protPP)
        self.assertIsNotNone(protPP.outputCoordinates,
                             "There was a problem with the Eman import coordinates")


        print "<Run extract particles with Same as picking>"
        protExtract = self.newProtocol(XmippProtExtractParticles, boxSize=110,
                                       downsampleType=SAME_AS_PICKING,
                                       doFlip=True, doInvert=True, runMode=1)
        protExtract.inputCoordinates.set(protPP.outputCoordinates)
        protExtract.ctfRelations.set(protCTF.outputCTF)
        self.launchProtocol(protExtract)
        self.assertIsNotNone(protExtract.outputParticles,
                             "There was a problem with the extract particles")
        #self.validateFiles('protExtract', protExtract)
        
        print "Run Preprocess particles"
        protCropResize = self.newProtocol(XmippProtCropResizeParticles,
                                          doResize=True, resizeOption=1,
                                          resizeDim=110)
        protCropResize.inputParticles.set(protExtract.outputParticles)
        self.launchProtocol(protCropResize)
        
        self.assertIsNotNone(protCropResize.outputParticles,
                             "There was a problem with resize/crop the particles")
        
        print "Run ML2D"
        protML2D = self.newProtocol(XmippProtML2D, numberOfClasses=8, maxIters=2,
                                 numberOfMpi=2, numberOfThreads=2)
        protML2D.inputParticles.set(protCropResize.outputParticles)
        self.launchProtocol(protML2D)        
        self.assertIsNotNone(protML2D.outputClasses,
                             "There was a problem with ML2D")
        #self.validateFiles('protML2D', protML2D)
        
#         #FIXME: Check the initial model with EMAn and restore the next step
#         return
        
        print "Run Initial Model"
        protIniModel = self.newProtocol(EmanProtInitModel, numberOfIterations=1,
                                        numberOfModels=2,
                                        shrink=5, symmetry='icos',
                                        numberOfThreads=3)
        protIniModel.inputSet.set(protML2D.outputClasses)
        self.launchProtocol(protIniModel)        
        self.assertIsNotNone(protIniModel.outputVolumes,
                             "There was a problem with Initial Model")


