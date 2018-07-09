# **************************************************************************
# *
# * Authors:    Laura del Cano (ldelcano@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
from pyworkflow.tests import *
from pyworkflow.em.packages.relion import *
from pyworkflow.em.protocol import ImageHandler


class TestRelionBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='mda'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.particlesFn = cls.dataset.getFile('particles')
        cls.vol = cls.dataset.getFile('volumes')

    def checkOutput(self, prot, outputName, conditions=[]):
        """ Check that an output was generated and
        the condition is valid.
        """
        o = getattr(prot, outputName, None)
        locals()[outputName] = o
        self.assertIsNotNone(o, "Output: %s is None" % outputName)
        for cond in conditions:
            self.assertTrue(eval(cond), 'Condition failed: ' + cond)
    
    @classmethod
    def runImportParticles(cls, pattern, samplingRate, checkStack=False):
        """ Run an Import particles protocol. """
        protImport = cls.newProtocol(ProtImportParticles, 
                                     filesPath=pattern,
                                     samplingRate=samplingRate,
                                     checkStack=checkStack)
        cls.launchProtocol(protImport)
        # check that input images have been imported (a better way to do this?)
        if protImport.outputParticles is None:
            raise Exception('Import of images: %s, failed. outputParticles '
                            'is None.' % pattern)
        return protImport

    @classmethod
    def runImportParticlesStar(cls, partStar, mag, samplingRate):
        """ Import particles from Relion star file. """
        protImport = cls.newProtocol(ProtImportParticles,
                                     importFrom=ProtImportParticles.IMPORT_FROM_RELION,
                                     starFile=partStar,
                                     magnification=mag,
                                     samplingRate=samplingRate,
                                     haveDataBeenPhaseFlipped=True
                                     )
        cls.launchProtocol(protImport)
        # check that input images have been imported (a better way to do this?)
        if protImport.outputParticles is None:
            raise Exception('Import of images: %s, failed. outputParticles '
                            'is None.' % partStar)
        return protImport

    @classmethod
    def runNormalizeParticles(cls, particles):
        """ Run normalize particles protocol """
        protPreproc = cls.newProtocol(ProtRelionPreprocessParticles,
                                      doNormalize=True)
        protPreproc.inputParticles.set(particles)
        cls.launchProtocol(protPreproc)
        return protPreproc
    
    @classmethod
    def runImportVolumes(cls, pattern, samplingRate):
        """ Run an Import particles protocol. """
        protImport = cls.newProtocol(ProtImportVolumes, 
                                     filesPath=pattern,
                                     samplingRate=samplingRate)
        cls.launchProtocol(protImport)
        return protImport

    @classmethod
    def runImportMovies(cls, pattern, mag, samplingRate, dose):
        """ Run an Import movies protocol. """
        protImport = cls.newProtocol(ProtImportMovies,
                                     filesPath=pattern,
                                     magnification=mag,
                                     samplingRate=samplingRate,
                                     dosePerFrame=dose)
        cls.launchProtocol(protImport)
        return protImport


class TestRelionClassify2D(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestRelionBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
        cls.protNormalize = cls.runNormalizeParticles(cls.protImport.outputParticles)
    
    def testRelion2D(self):
        def _runRelionClassify2D(doGpu=False, label=''):
            prot2D = self.newProtocol(ProtRelionClassify2D,
                                      doCTF=False, maskDiameterA=340,
                                      numberOfMpi=4, numberOfThreads=1)
            prot2D.numberOfClasses.set(4)
            prot2D.numberOfIterations.set(3)
            prot2D.inputParticles.set(self.protNormalize.outputParticles)
            prot2D.setObjLabel(label)

            if isVersion2():
                prot2D.doGpu.set(doGpu)

            self.launchProtocol(prot2D)
            return prot2D

        def _checkAsserts(relionProt):

            self.assertIsNotNone(relionProt.outputClasses,
                                 "There was a problem with Relion 2D classify")
        
            partsPixSize = self.protNormalize.outputParticles.getSamplingRate()
            classsesPixSize = relionProt.outputClasses.getImages().getSamplingRate()
            self.assertAlmostEquals(partsPixSize,classsesPixSize,
                                    "There was a problem with the sampling rate "
                                    "of the particles")
            for class2D in relionProt.outputClasses:
                self.assertTrue(class2D.hasAlignment2D())

        if isVersion2():
            relionNoGpu = _runRelionClassify2D(False, "Relion classify2D No GPU")
            _checkAsserts(relionNoGpu)

            environ = Environ(os.environ)
            cudaPath = environ.getFirst(('RELION_CUDA_LIB', 'CUDA_LIB'))

            if cudaPath is not None and os.path.exists(cudaPath):
                relionGpu = _runRelionClassify2D(True, "Relion classify2D GPU")
                _checkAsserts(relionGpu)
        else:
            relionProt = _runRelionClassify2D(label="Run Relion classify2D")
            _checkAsserts(relionProt)


class TestRelionClassify3D(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestRelionBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
        cls.protImportVol = cls.runImportVolumes(cls.vol, 3.5)
    
    def testProtRelionClassify3D(self):
        relionNormalize = self.newProtocol(ProtRelionPreprocessParticles)
        relionNormalize.inputParticles.set(self.protImport.outputParticles)
        relionNormalize.doNormalize.set(True)
        self.launchProtocol(relionNormalize)

        def _runRelionClassify3D(doGpu=False, label=''):

            print label
            relion3DClass = self.newProtocol(ProtRelionClassify3D,
                                             numberOfClasses=3,
                                             numberOfIterations=4,
                                             doCTF=False, runMode=1,
                                             maskDiameterA=320,
                                             numberOfMpi=2, numberOfThreads=2)

            relion3DClass.setObjLabel(label)
            relion3DClass.inputParticles.set(relionNormalize.outputParticles)
            relion3DClass.referenceVolume.set(self.protImportVol.outputVolume)


            if isVersion2():
                relion3DClass.doGpu.set(doGpu)

            self.launchProtocol(relion3DClass)
            return relion3DClass

        def _checkAsserts(relionProt):
            self.assertIsNotNone(relionProt.outputClasses, "There was a "
                                                           "problem with "
                                                           "Relion 3D classify")

            for class3D in relionProt.outputClasses:
                self.assertTrue(class3D.hasAlignmentProj())

        if isVersion2():
            relionNoGpu = _runRelionClassify3D(False, "Relion classify3D No GPU")
            _checkAsserts(relionNoGpu)

            environ = Environ(os.environ)
            cudaPath = environ.getFirst(('RELION_CUDA_LIB', 'CUDA_LIB'))

            if cudaPath is not None and os.path.exists(cudaPath):
                relionGpu = _runRelionClassify3D(True, "Relion classify3D GPU")
                _checkAsserts(relionGpu)
        else:
            relionProt = _runRelionClassify3D(label="Run Relion classify3D")
            _checkAsserts(relionProt)


class TestRelionRefine(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestRelionBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
        cls.protImportVol = cls.runImportVolumes(cls.vol, 3.5)
    
    def testProtRelionRefine(self):
        relNorm = self.newProtocol(ProtRelionPreprocessParticles)
        relNorm.inputParticles.set(self.protImport.outputParticles)
        relNorm.doNormalize.set(True)
        self.launchProtocol(relNorm)
        
        def _runRelionRefine(doGpu=False, label=''):
            relionRefine = self.newProtocol(ProtRelionRefine3D,
                                            doCTF=False, runMode=1,
                                            memoryPreThreads=1,
                                            maskDiameterA=340,
                                            symmetryGroup="d6",
                                            numberOfMpi=3, numberOfThreads=2)
            relionRefine.setObjLabel(label)
            relionRefine.inputParticles.set(relNorm.outputParticles)
            relionRefine.referenceVolume.set(self.protImportVol.outputVolume)
            
            if isVersion2():
                relionRefine.doGpu.set(doGpu)
            
            self.launchProtocol(relionRefine)
            return relionRefine
        
        def _checkAsserts(relionRefine):
            relionRefine._initialize()  # Load filename templates
            dataSqlite = relionRefine._getIterData(3)
            outImgSet = em.SetOfParticles(filename=dataSqlite)
            
            self.assertIsNotNone(relionRefine.outputVolume,
                                 "There was a problem with Relion autorefine")
            self.assertAlmostEqual(outImgSet[1].getSamplingRate(),
                                   relNorm.outputParticles[1].getSamplingRate(),
                                   "The sampling rate is wrong", delta=0.00001)
            
            self.assertAlmostEqual(outImgSet[1].getFileName(),
                                   relNorm.outputParticles[1].getFileName(),
                                   "The particles filenames are wrong")
        
        if isVersion2():
            environ = Environ(os.environ)
            cudaPath = environ.getFirst(('RELION_CUDA_LIB', 'CUDA_LIB'))

            hasCuda = (cudaPath is not None and
                       all(os.path.exists(p) for p in cudaPath.split(os.pathsep)))

            relionRefine = _runRelionRefine(hasCuda,
                                            "Relion auto-refine %sGPU"
                                            % ('' if hasCuda else 'NO-'))
            _checkAsserts(relionRefine)
        else:
            relionProt = _runRelionRefine(label="Run Relion auto-refine")
            _checkAsserts(relionProt)


class TestRelionInitialModel(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('relion_tutorial')
        cls.partFn = cls.dataset.getFile('import/classify2d/extra/relion_it015_data.star')
        cls.protImport = cls.runImportParticlesStar(cls.partFn, 50000, 7.08)

    def testProtRelionIniModel(self):
        if getVersion() in [V1_3, V1_4, V2_0]:
            raise Exception('Initial model protocol exists only for Relion v2.1 or higher!')

        def _runRelionIniModel(doGpu=True, label=''):
            print label
            relionIniModel = self.newProtocol(ProtRelionInitialModel,
                                              doCTF=False, doGpu=doGpu,
                                              maskDiameterA=340,
                                              numberOfIterations=2,
                                              symmetryGroup="C1",
                                              numberOfMpi=3, numberOfThreads=2)
            relionIniModel.setObjLabel(label)
            relionIniModel.inputParticles.set(self.protImport.outputParticles)
            self.launchProtocol(relionIniModel)

            return relionIniModel

        def _checkAsserts(relionProt):
            relionProt._initialize()  # Load filename templates
            dataSqlite = relionProt._getIterData(2)
            outImgSet = em.SetOfParticles(filename=dataSqlite)

            self.assertIsNotNone(relionProt.outputVolume,
                                 "There was a problem with Relion initial model")
            self.assertAlmostEqual(outImgSet[1].getSamplingRate(),
                                   self.protImport.outputParticles[1].getSamplingRate(),
                                   "The sampling rate is wrong", delta=0.00001)

        environ = Environ(os.environ)
        cudaPath = environ.getFirst(('RELION_CUDA_LIB', 'CUDA_LIB'))

        if cudaPath is not None and os.path.exists(cudaPath):
            relionGpu = _runRelionIniModel(doGpu=True, label="Relion initial model GPU")
            _checkAsserts(relionGpu)
        else:
            print "Warning: running this test on CPU might take a lot of time!"
            relion = _runRelionIniModel(doGpu=False, label="Relion initial model CPU")
            _checkAsserts(relion)

        
class TestRelionPreprocess(TestRelionBase):
    """ This class helps to test all different preprocessing particles options
    on Relion. """
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestRelionBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)

    def _validations(self, imgSet, dims, pxSize):
        self.assertIsNotNone(imgSet, "There was a problem with preprocess "
                                     "particles")
        xDim = imgSet.getXDim()
        sr = imgSet.getSamplingRate()
        self.assertEqual(xDim, dims, "The dimension of your particles are %d x "
                                     "%d and must be  %d x %d" % (xDim, xDim,
                                                                  dims, dims))
        self.assertAlmostEqual(sr, pxSize, 0.0001,
                               "Pixel size of your particles are  %0.2f and"
                               " must be %0.2f" % (sr, pxSize))

    def test_NormalizeAndDust(self):
        """ Normalize particles.
        """
        # Test now a normalization after the imported particles
        protocol = self.newProtocol(ProtRelionPreprocessParticles,
                                    doNormalize=True, backRadius=40,
                                    doRemoveDust=True, whiteDust=4, blackDust=8)
        protocol.setObjLabel('relion: norm-dust')
        protocol.inputParticles.set(self.protImport.outputParticles)
        self.launchProtocol(protocol)
        
        self._validations(protocol.outputParticles, 100, 3.5)
        
    def test_ScaleAndInvert(self):
        """ Test all options at once.
        """
        # Test now a normalization after the imported particles   
        protocol = self.newProtocol(ProtRelionPreprocessParticles,
                                    doNormalize=False,
                                    doScale=True, scaleSize=50,
                                    doInvert=True)
        protocol.setObjLabel('relion: scale-invert')
        protocol.inputParticles.set(self.protImport.outputParticles)
        
        self.launchProtocol(protocol)
        self._validations(protocol.outputParticles, 50, 7.0)


class TestRelionSubtract(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsRelion = DataSet.getDataSet('relion_tutorial')
    
    def test_subtract(self):
        protParts = self.newProtocol(ProtImportParticles,
                                     objLabel='from relion auto-refine',
                                     importFrom=ProtImportParticles.IMPORT_FROM_RELION,
                                     starFile=self.dsRelion.getFile('import/refine3d/extra/relion_it001_data.star'),
                                     magnification=10000,
                                     samplingRate=7.08,
                                     haveDataBeenPhaseFlipped=True
                                     )
        self.launchProtocol(protParts)
        self.assertEqual(60, protParts.outputParticles.getXDim())
        
        protVol = self.newProtocol(ProtImportVolumes,
                                   filesPath=self.dsRelion.getFile('volumes/reference.mrc'),
                                   samplingRate=7.08)
        self.launchProtocol(protVol)
        self.assertEqual(60, protVol.outputVolume.getDim()[0])
        
        protSubtract = self.newProtocol(ProtRelionSubtract)
        protSubtract.inputParticles.set(protParts.outputParticles)
        protSubtract.inputVolume.set(protVol.outputVolume)
        self.launchProtocol(protSubtract)
        self.assertIsNotNone(protSubtract.outputParticles,
                             "There was a problem with subtract projection")


class TestRelionSortParticles(TestRelionBase):
    """ This class helps to test sort particles protocol from Relion. """

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('relion_tutorial')

        cls.dataset = DataSet.getDataSet('relion_tutorial')
        cls.partFn = cls.dataset.getFile('import/particles.sqlite')
        cls.partAvg = cls.dataset.getFile('import/averages.mrcs')
        cls.partCl2dFn = cls.dataset.getFile('import/classify2d/extra/relion_it015_data.star')
        cls.partCl3dFn = cls.dataset.getFile('import/classify3d/extra/relion_it015_data.star')
        cls.partRef3dFn = cls.dataset.getFile('import/refine3d/extra/relion_data.star')
        cls.volFn = cls.dataset.getFile('import/refine3d/extra/relion_class001.mrc')

    def importParticles(self, partStar):
        """ Import particles from Relion star file. """
        protPart = self.newProtocol(ProtImportParticles,
                                    importFrom=ProtImportParticles.IMPORT_FROM_RELION,
                                    starFile=partStar,
                                    magnification=10000,
                                    samplingRate=7.08,
                                    haveDataBeenPhaseFlipped=True
                                    )
        self.launchProtocol(protPart)
        return protPart

    def importParticlesFromScipion(self):
        partFn = self.ds.getFile('import/particles.sqlite')
        protPart = self.newProtocol(ProtImportParticles,
                                    objLabel='from xmipp extract (after relion auto-picking)',
                                    importFrom=ProtImportParticles.IMPORT_FROM_SCIPION,
                                    sqliteFile=partFn,
                                    magnification=10000,
                                    samplingRate=7.08,
                                    haveDataBeenPhaseFlipped=True
                                    )

        self.launchProtocol(protPart)
        return protPart

    def importAverages(self):
        """ Import averages used for relion autopicking. """
        partAvg = self.ds.getFile('import/averages.mrcs')
        protAvg = self.newProtocol(ProtImportAverages,
                                   importFrom=ProtImportParticles.IMPORT_FROM_FILES,
                                   filesPath=partAvg,
                                   samplingRate=7.08
                                   )
        self.launchProtocol(protAvg)
        return protAvg

    def importVolume(self):
        volFn = self.ds.getFile('import/refine3d/extra/relion_class001.mrc')
        protVol = self.newProtocol(ProtImportVolumes,
                                   objLabel='import volume',
                                   filesPath=volFn,
                                   samplingRate=7.08)
        self.launchProtocol(protVol)
        return protVol

    def test_afterCl2D(self):
        partCl2dFn = self.ds.getFile(
            'import/classify2d/extra/relion_it015_data.star')
        importRun = self.importParticles(partCl2dFn)

        prot = self.newProtocol(ProtRelionSortParticles)
        prot.setObjLabel('relion - sort after cl2d')
        prot.inputSet.set(importRun.outputClasses)
        self.launchProtocol(prot)

    def test_afterCl3D(self):
        prot = self.newProtocol(ProtRelionSortParticles)
        prot.setObjLabel('relion - sort after cl3d')
        partCl3dFn = self.ds.getFile(
            'import/classify3d/extra/relion_it015_data.star')
        importRun = self.importParticles(partCl3dFn)
        prot.inputSet.set(importRun.outputClasses)
        self.launchProtocol(prot)

    def test_after3DRefinement(self):
        prot = self.newProtocol(ProtRelionSortParticles)
        prot.setObjLabel('relion - sort after ref3d')
        partRef3dFn = self.ds.getFile('import/refine3d/extra/relion_data.star')
        importRun = self.importParticles(partRef3dFn)
        prot.inputSet.set(importRun.outputParticles)
        prot.referenceVolume.set(self.importVolume().outputVolume)
        self.launchProtocol(prot)

    def test_afterPicking(self):
        prot = self.newProtocol(ProtRelionSortParticles)
        prot.setObjLabel('relion - sort after picking')
        prot.inputSet.set(self.importParticlesFromScipion().outputParticles)
        prot.referenceAverages.set(self.importAverages().outputAverages)
        self.launchProtocol(prot)


class TestRelionPostprocess(TestRelionBase):
    """ This class helps to test postprocess protocol from Relion. """

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('relion_tutorial')
    
    def importVolume(self):
        volFn = self.ds.getFile('import/refine3d/extra/relion_class001.mrc')
        protVol = self.newProtocol(ProtImportVolumes,
                                   objLabel='import volume',
                                   filesPath=volFn,
                                   samplingRate=3)
        self.launchProtocol(protVol)
        return protVol

    def importPartsFromScipion(self):
        partFn = self.ds.getFile('import/particles.sqlite')
        protPart = self.newProtocol(ProtImportParticles,
                                    objLabel='Import Particles',
                                    importFrom=ProtImportParticles.IMPORT_FROM_SCIPION,
                                    sqliteFile=partFn,
                                    magnification=10000,
                                    samplingRate=3,
                                    haveDataBeenPhaseFlipped=True
                                    )
        self.launchProtocol(protPart)
        return protPart
    
    def _createRef3DProtBox(self, label, protocol, volPatt="kk.mrc",
                            storeIter=False, iterN=0):
        from pyworkflow.protocol.constants import STATUS_FINISHED
        
        prot = self.newProtocol(protocol)
        self.saveProtocol(prot)
        
        prot.setObjLabel(label)
        project = prot.getProject()
        makePath(prot._getPath())
        makePath(prot._getExtraPath())
        makePath(prot._getTmpPath())

        prot.inputParticles.set(self.importPartsFromScipion().outputParticles)
        
        protClassName = prot.getClassName()
        if protClassName.startswith('ProtRelionRefine3D'):
            prot.referenceVolume.set(self.importVolume().outputVolume)
        elif protClassName.startswith('ProtFrealign'):
            prot.input3DReference.set(self.importVolume().outputVolume)
        elif protClassName.startswith('XmippProtProjMatch'):
            prot.input3DReferences.set(self.importVolume().outputVolume)
        elif protClassName.startswith('EmanProtRefine'):
            pass
        
        volume = em.Volume()
        volume.setFileName(prot._getExtraPath(volPatt))
        pxSize = prot.inputParticles.get().getSamplingRate()
        volume.setSamplingRate(pxSize)
        if storeIter:
            prot._setLastIter(iterN)
        prot._defineOutputs(outputVolume=volume)

        prot.setStatus(STATUS_FINISHED)
        project._storeProtocol(prot)
        return prot
    
    def _validations(self, vol, dims, pxSize, prot=""):
        self.assertIsNotNone(vol, "There was a problem with postprocess "
                                  "protocol, using %s protocol as input" % prot)
        xDim = vol.getXDim()
        sr = vol.getSamplingRate()
        self.assertEqual(xDim, dims, "The dimension of your volume is (%d)^3 "
                                     "and must be (%d)^3" % (xDim, dims))

        self.assertAlmostEqual(sr, pxSize, 0.0001,
                               "Pixel size of your volume is %0.2f and"
                               " must be %0.2f" % (sr, pxSize))
    
    def test_postProcess_from_autorefine(self):
        pathFns = 'import/refine3d/extra'
        vol = self.ds.getFile(join(pathFns, 'relion_class001.mrc'))
        half1 = self.ds.getFile(join(pathFns,
                                     'relion_it025_half1_class001.mrc'))
        half2 = self.ds.getFile(join(pathFns,
                                     'relion_it025_half2_class001.mrc'))
        volPatt = 'relion_class001.mrc'
        
        protRef = self._createRef3DProtBox("auto-refine",
                                           ProtRelionRefine3D,
                                           volPatt)
        
        copyFile(vol, protRef._getExtraPath(volPatt))
        copyFile(half1,
                 protRef._getExtraPath('relion_half1_class001_unfil.mrc'))
        copyFile(half2,
                 protRef._getExtraPath('relion_half2_class001_unfil.mrc'))

        postProt = self.newProtocol(ProtRelionPostprocess,
                                    protRefine=protRef)
        postProt.setObjLabel('post process Auto-refine')
        
        self.launchProtocol(postProt)
        self._validations(postProt.outputVolume, 60, 3, "Relion auto-refine")
        
    def test_postProcess_from_frealign(self):
        from pyworkflow.em.packages.grigoriefflab import ProtFrealign
        
        pathFns = 'import/refine3d/extra'
        vol = self.ds.getFile(join(pathFns, 'relion_class001.mrc'))
        half1 = self.ds.getFile(join(pathFns,
                                     'relion_it025_half1_class001.mrc'))
        half2 = self.ds.getFile(join(pathFns,
                                     'relion_it025_half2_class001.mrc'))
        volPatt = 'iter_002/volume_iter_002.mrc'

        protRef = self._createRef3DProtBox("frealign", ProtFrealign, volPatt,
                                           storeIter=True, iterN=2)
        makePath(join(protRef._getExtraPath(), 'iter_002'))
        
        copyFile(vol, protRef._getExtraPath(volPatt))
        copyFile(half1,
                 protRef._getExtraPath('iter_002/volume_1_iter_002.mrc'))
        copyFile(half2,
                 protRef._getExtraPath('iter_002/volume_2_iter_002.mrc'))

        postProt = self.newProtocol(ProtRelionPostprocess,
                                    protRefine=protRef)
        postProt.setObjLabel('post process frealign')
        self.launchProtocol(postProt)
        self._validations(postProt.outputVolume, 60, 3, "Frealign")

    def test_postProcess_from_projMatch(self):
        from pyworkflow.em.packages.xmipp3 import XmippProtProjMatch
        
        pathFns = 'import/refine3d/extra'
        vol = self.ds.getFile(join(pathFns, 'relion_class001.mrc'))
        half1 = self.ds.getFile(join(pathFns,
                                     'relion_it025_half1_class001.mrc'))
        half2 = self.ds.getFile(join(pathFns,
                                     'relion_it025_half2_class001.mrc'))
        
        volPatt = 'iter_002/reconstruction_Ref3D_001.vol'
        protRef = self._createRef3DProtBox("Proj Match", XmippProtProjMatch,
                                           volPatt, storeIter=True, iterN=2)
        
        makePath(join(protRef._getExtraPath(), 'iter_002'))
        
        protRef._initialize()
        
        volXmipp = protRef._getFileName('reconstructedFileNamesIters',
                                        iter=2, ref=1)
        half1Xmipp = protRef._getFileName('reconstructedFileNamesItersSplit1',
                                          iter=2, ref=1)
        half2Xmipp = protRef._getFileName('reconstructedFileNamesItersSplit2',
                                          iter=2, ref=1)
        
        ih = ImageHandler()
        ih.convert(ih.getVolFileName(vol), volXmipp)
        ih.convert(ih.getVolFileName(half1), half1Xmipp)
        ih.convert(ih.getVolFileName(half2), half2Xmipp)

        postProt = self.newProtocol(ProtRelionPostprocess,
                                    protRefine=protRef)
        postProt.setObjLabel('post process Xmipp Projection Matching')
        self.launchProtocol(postProt)
        self._validations(postProt.outputVolume, 60, 3, "Projection Matching")
    
    def test_postProcess_from_eman_refineEasy(self):
        from pyworkflow.em.packages.eman2 import EmanProtRefine
        from pyworkflow.em.packages.eman2.convert import convertImage
        
        pathFns = 'import/refine3d/extra'
        vol = self.ds.getFile(join(pathFns, 'relion_class001.mrc'))
        half1 = self.ds.getFile(join(pathFns,
                                     'relion_it025_half1_class001.mrc'))
        half2 = self.ds.getFile(join(pathFns,
                                     'relion_it025_half2_class001.mrc'))

        volPatt = 'refine_01/threed_02.hdf'
        protRef = self._createRef3DProtBox("Eman refine Easy",
                                           EmanProtRefine, volPatt)
        makePath(join(protRef._getExtraPath(), 'refine_01'))
        protRef._createFilenameTemplates()
        protRef._createFilenameTemplates()

        volEman = protRef._getFileName("mapFull", run=1, iter=2)
        half1Eman = protRef._getFileName("mapEvenUnmasked", run=1)
        half2Eman = protRef._getFileName("mapOddUnmasked", run=1)
        
        convertImage(vol, volEman)
        convertImage(half1, half1Eman)
        convertImage(half2, half2Eman)

        postProt = self.newProtocol(ProtRelionPostprocess,
                                    protRefine=protRef)
        postProt.setObjLabel('post process Eman2 refine-easy')
        self.launchProtocol(postProt)
        self._validations(postProt.outputVolume, 60, 3, "Eman refine easy")


class TestRelionLocalRes(TestRelionBase):
    """ This class helps to test local resolution protocol from Relion. """

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('relion_tutorial')
        cls.partFn = cls.ds.getFile('import/refine3d/extra/relion_it025_data.star')
        cls.protImport = cls.runImportParticlesStar(cls.partFn, 50000, 7.08)
        cls.volFn = cls.ds.getFile('import/refine3d/extra/relion_class001.mrc')

    def importVolume(self):
        protVol = self.newProtocol(ProtImportVolumes,
                                   objLabel='import volume',
                                   filesPath=self.volFn,
                                   samplingRate=7.08)
        self.launchProtocol(protVol)
        return protVol

    def _createRef3DProtBox(self, label):
        from pyworkflow.protocol.constants import STATUS_FINISHED

        prot = self.newProtocol(ProtRelionRefine3D)
        self.saveProtocol(prot)
        prot.setObjLabel(label)
        project = prot.getProject()
        makePath(prot._getPath())
        makePath(prot._getExtraPath())
        makePath(prot._getTmpPath())

        prot.inputParticles.set(self.protImport.outputParticles)
        prot.referenceVolume.set(self.importVolume().outputVolume)

        volume = em.Volume()
        volume.setFileName(prot._getExtraPath("relion_class001.mrc"))
        pxSize = prot.inputParticles.get().getSamplingRate()
        volume.setSamplingRate(pxSize)

        prot._defineOutputs(outputVolume=volume)
        prot.setStatus(STATUS_FINISHED)
        project._storeProtocol(prot)
        return prot

    def _validations(self, vol, dims, pxSize):
        self.assertIsNotNone(vol, "There was a problem with localres protocol ")
        xDim = vol.getXDim()
        sr = vol.getSamplingRate()
        self.assertEqual(xDim, dims, "The dimension of your volume is (%d)^3 "
                                     "and must be (%d)^3" % (xDim, dims))
        self.assertAlmostEqual(sr, pxSize, 0.0001,
                               "Pixel size of your volume is %0.2f and"
                               " must be %0.2f" % (sr, pxSize))

    def test_runRelionLocalRes(self):
        if not isVersion2():
            raise Exception('Local resolution protocol exists only for Relion v2.0 or higher!')

        pathFns = 'import/refine3d/extra'
        vol = self.ds.getFile(join(pathFns, 'relion_class001.mrc'))
        half1 = self.ds.getFile(join(pathFns, 'relion_it025_half1_class001.mrc'))
        half2 = self.ds.getFile(join(pathFns, 'relion_it025_half2_class001.mrc'))
        volPatt = 'relion_class001.mrc'
        modelFn = self.ds.getFile(join(pathFns, 'relion_model.star'))
        protRef = self._createRef3DProtBox("auto-refine")

        copyFile(vol, protRef._getExtraPath(volPatt))
        copyFile(half1,
                 protRef._getExtraPath('relion_half1_class001_unfil.mrc'))
        copyFile(half2,
                 protRef._getExtraPath('relion_half2_class001_unfil.mrc'))
        copyFile(modelFn, protRef._getExtraPath('relion_model.star'))

        postProt = self.newProtocol(ProtRelionLocalRes, protRefine=protRef)
        postProt.setObjLabel('Relion local resolution')

        self.launchProtocol(postProt)
        self._validations(postProt.outputVolume, 60, 7.08)


class TestRelionExpandSymmetry(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('relion_tutorial')
        cls.partRef3dFn = cls.dataset.getFile('import/refine3d/extra/relion_data.star')

    def importParticles(self, partStar):
        """ Import particles from Relion star file. """
        protPart = self.newProtocol(ProtImportParticles,
                                    importFrom=ProtImportParticles.IMPORT_FROM_RELION,
                                    starFile=partStar,
                                    magnification=10000,
                                    samplingRate=7.08,
                                    haveDataBeenPhaseFlipped=True
                                    )
        self.launchProtocol(protPart)
        return protPart

    def test_ExpandSymmetry(self):
        if not isVersion2():
            raise Exception('Expand symmetry protocol exists only for Relion v2.0 or higher!')

        prot = self.newProtocol(ProtRelionExpandSymmetry)
        print "Import particles"
        importRun = self.importParticles(self.partRef3dFn)
        prot.inputParticles.set(importRun.outputParticles)
        prot.symmetryGroup.set("D2")
        print "Run expand symmetry"
        self.launchProtocol(prot)

        self.assertIsNotNone(prot.outputParticles,
                             "There was a problem with expand symmetry protocol")
        sizeIn = importRun.outputParticles.getSize()
        sizeOut = prot.outputParticles.getSize()
        self.assertAlmostEqual(sizeIn * 4, sizeOut, 0.0001,
                               "Number of output particles is %d and"
                               " must be %d" % (sizeOut, sizeIn * 4))


class TestRelionExtractMovieParticles(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('relion_tutorial')
        cls.partRef3dFn = cls.ds.getFile('import/refine3d_case2/relion_data.star')
        cls.movies = DataSet.getDataSet('movies').getFile('ribo/*.mrcs')
        cls.protImportMovies = cls.runImportMovies(cls.movies, 50000, 3.54, 1.0)
        cls.protImportParts = cls.runImportParticlesStar(cls.partRef3dFn,
                                                         50000, 3.54)

    def _runExtract(self, boxsize, frame0, frameN, avgFrames, doInvert):
        if not isVersion2():
            raise Exception('Extract movie particles protocol exists only for Relion v2.0 or higher!')

        prot = self.newProtocol(ProtRelionExtractMovieParticles,
                                boxSize=boxsize, frame0=frame0, frameN=frameN,
                                avgFrames=avgFrames, doInvert=doInvert)
        prot.inputMovies.set(self.protImportMovies.outputMovies)
        prot.inputParticles.set(self.protImportParts.outputParticles)
        self.launchProtocol(prot)

        self.assertIsNotNone(prot.outputParticles,
                             "There was a problem with extract movie particles protocol")
        sizeIn = prot.inputParticles.get().getSize()
        sizeOut = prot.outputParticles.getSize()
        factor = (frameN - frame0 + 1) / avgFrames
        if factor % avgFrames > 0:
            factor += 1
        self.assertEqual(sizeIn * factor, sizeOut,
                         "Number of output movie particles is %d and must be "
                         "%d" % (sizeOut, sizeIn * factor))

    def test_extractMovieParticlesAvg1(self):
        self._runExtract(100, 1, 16, 1, True)

    def test_extractMovieParticlesAvg3(self):
        self._runExtract(100, 1, 16, 3, True)

        
class TestRelionCreate3dMask(TestRelionBase):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('relion_tutorial')

    def importVolume(self):
        volFn = self.ds.getFile('import/refine3d/extra/relion_class001.mrc')
        protVol = self.newProtocol(ProtImportVolumes,
                                   objLabel='import volume',
                                   filesPath=volFn,
                                   samplingRate=3)
        self.launchProtocol(protVol)
        return protVol

    def _validations(self, mask, dims, pxSize, prot):
        self.assertIsNotNone(mask, "There was a problem with mask 3d protocol, "
                                  "using %s protocol as input" % prot)
        xDim = mask.getXDim()
        sr = mask.getSamplingRate()
        self.assertEqual(xDim, dims, "The dimension of your volume is (%d)^3 "
                                     "and must be (%d)^3" % (xDim, dims))

        self.assertAlmostEqual(sr, pxSize, 0.0001,
                               "Pixel size of your volume is %0.2f and"
                               " must be %0.2f" % (sr, pxSize))

    def test_createMask(self):
        importProt = self.importVolume()

        maskProt = self.newProtocol(ProtRelionCreateMask3D,
                                    initialLowPassFilterA=10)  # filter at 10 A
        vol = importProt.outputVolume
        maskProt.inputVolume.set(vol)
        self.launchProtocol(maskProt)

        self._validations(maskProt.outputMask, vol.getXDim(),
                          vol.getSamplingRate(), maskProt)

        ih = ImageHandler()
        img = ih.read(maskProt.outputMask)
        mean, std, _min, _max = img.computeStats()
        # Check the mask is non empty and between 0 and 1
        self.assertAlmostEqual(_min, 0)
        self.assertAlmostEqual(_max, 1)
        self.assertTrue(mean > 0)


class TestRelionExtractParticles(TestRelionBase):
    """This class check if the protocol to extract particles
    in Relion works properly.
    """

    @classmethod
    def runImportMicrograph(cls, pattern, samplingRate, voltage,
                            scannedPixelSize, magnification,
                            sphericalAberration):
        """ Run an Import micrograph protocol. """

        # We have two options:
        # 1) pass the SamplingRate or
        # 2) the ScannedPixelSize + microscope magnification
        if not samplingRate is None:
            cls.protImport = cls.newProtocol(ProtImportMicrographs,
                                             samplingRateMode=0,
                                             filesPath=pattern,
                                             samplingRate=samplingRate,
                                             magnification=magnification,
                                             voltage=voltage,
                                             sphericalAberration=sphericalAberration)
        else:
            cls.protImport = cls.newProtocol(ProtImportMicrographs,
                                             samplingRateMode=1,
                                             filesPath=pattern,
                                             scannedPixelSize=scannedPixelSize,
                                             voltage=voltage,
                                             magnification=magnification,
                                             sphericalAberration=sphericalAberration)

        cls.protImport.setObjLabel('import mics')
        cls.launchProtocol(cls.protImport)
        if cls.protImport.isFailed():
            raise Exception("Protocol has failed. Error: ", cls.protImport.getErrorMessage())
        # check that input micrographs have been imported (a better way to do this?)
        if cls.protImport.outputMicrographs is None:
            raise Exception('Import of micrograph: %s, failed. outputMicrographs is None.' % pattern)
        return cls.protImport

    @classmethod
    def runImportMicrographBPV(cls, pattern):
        """ Run an Import micrograph protocol. """
        return cls.runImportMicrograph(pattern, samplingRate=1.237,
                                       voltage=300, sphericalAberration=2,
                                       scannedPixelSize=None, magnification=56000)

    @classmethod
    def runDownsamplingMicrographs(cls, mics, downFactorValue, threads=1):
        # test downsampling a set of micrographs
        from pyworkflow.em.packages.xmipp3 import XmippProtPreprocessMicrographs
        cls.protDown = XmippProtPreprocessMicrographs(doDownsample=True,
                                                      downFactor=downFactorValue,
                                                      numberOfThreads=threads)
        cls.protDown.inputMicrographs.set(mics)
        cls.proj.launchProtocol(cls.protDown, wait=True)
        return cls.protDown

    @classmethod
    def runFakedPicking(cls, mics, pattern):
        """ Run a faked particle picking. Coordinates already existing. """
        from pyworkflow.em.packages.xmipp3 import XmippProtParticlePicking
        cls.protPP = XmippProtParticlePicking(importFolder=pattern, runMode=1)
        cls.protPP.inputMicrographs.set(mics)
        cls.proj.launchProtocol(cls.protPP, wait=True)
        # check that faked picking has run ok
        if cls.protPP.outputCoordinates is None:
            raise Exception('Faked particle picking failed. outputCoordinates is None.')
        return cls.protPP

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('xmipp_tutorial')
        cls.micFn = cls.dataset.getFile('mic1')
        cls.micsFn = cls.dataset.getFile('allMics')
        cls.coordsDir = cls.dataset.getFile('posSupervisedDir')
        cls.allCrdsDir = cls.dataset.getFile('posAllDir')

        cls.DOWNSAMPLING = 5.0
        cls.protImport = cls.runImportMicrographBPV(cls.micsFn)
        cls.protDown = cls.runDownsamplingMicrographs(cls.protImport.outputMicrographs,
                                                      cls.DOWNSAMPLING)

        cls.protCTF = cls.newProtocol(ProtImportCTF,
                                      importFrom=ProtImportCTF.IMPORT_FROM_XMIPP3,
                                      filesPath=cls.dataset.getFile('ctfsDir'),
                                      filesPattern='*.ctfparam')
        cls.protCTF.inputMicrographs.set(cls.protImport.outputMicrographs)
        cls.proj.launchProtocol(cls.protCTF, wait=True)

        cls.protPP = cls.runFakedPicking(cls.protDown.outputMicrographs, cls.allCrdsDir)

    def _checkSamplingConsistency(self, outputSet):
        """ Check that the set sampling is the same as item sampling. """
        first = outputSet.getFirstItem()

        self.assertAlmostEqual(outputSet.getSamplingRate(),
                               first.getSamplingRate())

    def testExtractSameAsPicking(self):
        print "Run extract particles from same micrographs as picking"
        protExtract = self.newProtocol(ProtRelionExtractParticles,
                                       boxSize=110,
                                       doInvert=False)
        protExtract.setObjLabel("extract-same as picking")
        protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
        self.launchProtocol(protExtract)

        inputCoords = protExtract.inputCoordinates.get()
        outputParts = protExtract.outputParticles
        micSampling = protExtract.inputCoordinates.get().getMicrographs().getSamplingRate()

        self.assertIsNotNone(outputParts,
                             "There was a problem generating the output.")
        self.assertAlmostEqual(outputParts.getSamplingRate() / micSampling,
                               1, 1,
                               "There was a problem generating the output.")
        self._checkSamplingConsistency(outputParts)

        def compare(objId, delta=0.001):
            cx, cy = inputCoords[objId].getPosition()
            px, py = outputParts[objId].getCoordinate().getPosition()
            micNameCoord = inputCoords[objId].getMicName()
            micNamePart = outputParts[objId].getCoordinate().getMicName()
            self.assertAlmostEquals(cx, px, delta=delta)
            self.assertAlmostEquals(cy, py, delta=delta)
            self.assertEqual(micNameCoord, micNamePart,
                             "The micName should be %s and its %s"
                             % (micNameCoord, micNamePart))

        compare(83)
        compare(228)

    def testExtractOriginal(self):
        print "Run extract particles from the original micrographs"
        protExtract = self.newProtocol(ProtRelionExtractParticles,
                                       boxSize=550,
                                       downsampleType=OTHER,
                                       doInvert=False)
        protExtract.setObjLabel("extract-original")
        protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
        protExtract.inputMicrographs.set(self.protImport.outputMicrographs)
        self.launchProtocol(protExtract)

        inputCoords = protExtract.inputCoordinates.get()
        outputParts = protExtract.outputParticles
        samplingCoords = self.protPP.outputCoordinates.getMicrographs().getSamplingRate()
        samplingFinal = self.protImport.outputMicrographs.getSamplingRate()
        samplingMics = protExtract.inputMicrographs.get().getSamplingRate()
        factor = samplingFinal / samplingCoords

        def compare(objId, delta=1.0):
            cx, cy = inputCoords[objId].getPosition()
            px, py = outputParts[objId].getCoordinate().getPosition()
            micNameCoord = inputCoords[objId].getMicName()
            micNamePart = outputParts[objId].getCoordinate().getMicName()
            self.assertAlmostEquals(cx / factor, px, delta=delta)
            self.assertAlmostEquals(cy / factor, py, delta=delta)
            self.assertEqual(micNameCoord, micNamePart,
                             "The micName should be %s and its %s"
                             % (micNameCoord, micNamePart))

        compare(111)
        compare(7)

        self.assertIsNotNone(outputParts,
                             "There was a problem generating the output.")
        self.assertEqual(outputParts.getSamplingRate(), samplingMics,
                         "Output sampling rate should be equal to input "
                         "sampling rate.")
        self._checkSamplingConsistency(outputParts)

    def testExtractOther(self):
        print "Run extract particles from original micrographs, with downsampling"
        downFactor = 2.989
        protExtract = self.newProtocol(ProtRelionExtractParticles,
                                       boxSize=550, downsampleType=OTHER,
                                       doRescale=True,
                                       rescaledSize=184,
                                       doInvert=False,
                                       doFlip=False)
        # Get all the micrographs ids to validate that all particles
        # has the micId properly set
        micsId = [mic.getObjId() for mic in
                  self.protPP.outputCoordinates.getMicrographs()]

        protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
        protExtract.inputMicrographs.set(self.protImport.outputMicrographs)
        protExtract.setObjLabel("extract-other")
        self.launchProtocol(protExtract)

        inputCoords = protExtract.inputCoordinates.get()
        outputParts = protExtract.outputParticles
        samplingCoords = self.protPP.outputCoordinates.getMicrographs().getSamplingRate()
        samplingFinal = self.protImport.outputMicrographs.getSamplingRate() * downFactor
        samplingMics = protExtract.inputMicrographs.get().getSamplingRate()
        factor = samplingFinal / samplingCoords
        self.assertIsNotNone(outputParts,
                             "There was a problem generating the output.")

        def compare(objId, delta=2.0):
            cx, cy = inputCoords[objId].getPosition()
            px, py = outputParts[objId].getCoordinate().getPosition()
            micNameCoord = inputCoords[objId].getMicName()
            micNamePart = outputParts[objId].getCoordinate().getMicName()
            self.assertAlmostEquals(cx / factor, px, delta=delta)
            self.assertAlmostEquals(cy / factor, py, delta=delta)
            self.assertEqual(micNameCoord, micNamePart,
                             "The micName should be %s and its %s"
                             % (micNameCoord, micNamePart))

        compare(45)
        compare(229)

        outputSampling = outputParts.getSamplingRate()
        self.assertAlmostEqual(outputSampling / samplingMics,
                               downFactor, 1,
                               "There was a problem generating the output.")
        for particle in outputParts:
            self.assertTrue(particle.getCoordinate().getMicId() in micsId)
            self.assertAlmostEqual(outputSampling, particle.getSamplingRate())

    def testExtractCTF(self):
        print "Run extract particles with CTF"
        protExtract = self.newProtocol(ProtRelionExtractParticles,
                                       boxSize=110,
                                       downsampleType=SAME_AS_PICKING,
                                       doInvert=False,
                                       doFlip=True)
        protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
        protExtract.ctfRelations.set(self.protCTF.outputCTF)
        protExtract.setObjLabel("extract-ctf")
        self.launchProtocol(protExtract)

        inputCoords = protExtract.inputCoordinates.get()
        outputParts = protExtract.outputParticles

        def compare(objId, delta=0.001):
            cx, cy = inputCoords[objId].getPosition()
            px, py = outputParts[objId].getCoordinate().getPosition()
            micNameCoord = inputCoords[objId].getMicName()
            micNamePart = outputParts[objId].getCoordinate().getMicName()
            self.assertAlmostEquals(cx, px, delta=delta)
            self.assertAlmostEquals(cy, py, delta=delta)
            self.assertEqual(micNameCoord, micNamePart,
                             "The micName should be %s and its %s"
                             % (micNameCoord, micNamePart))

        compare(228)
        compare(83)

        def compareCTF(partId, ctfId):
            partDefU = outputParts[partId].getCTF().getDefocusU()
            defU = protExtract.ctfRelations.get()[ctfId].getDefocusU()
            self.assertAlmostEquals(partDefU, defU, delta=1)

        compareCTF(1, 1)
        compareCTF(150, 2)
        compareCTF(300, 3)

        self.assertIsNotNone(outputParts,
                             "There was a problem generating the output.")
        self.assertTrue(outputParts.hasCTF(), "Output does not have CTF.")
        self._checkSamplingConsistency(outputParts)


class TestRelionCenterAverages(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('mda')

    def test_basic(self):
        """ Run an Import particles protocol. """
        protImport = self.newProtocol(ProtImportAverages,
                                     filesPath=self.ds.getFile('averages/averages.stk'),
                                     samplingRate=5.04)
        self.launchProtocol(protImport)
        inputAvgs = protImport.outputAverages
        protCenter = self.newProtocol(ProtRelionCenterAverages)
        protCenter.inputAverages.set(inputAvgs)
        self.launchProtocol(protCenter)

        conditions = ['outputAverages.getSize()==%d' % inputAvgs.getSize(),
                      'outputAverages.getSamplingRate() - %f < 0.00001'
                      % inputAvgs.getSamplingRate()]
        self.checkOutput(protCenter, 'outputAverages', conditions)


class TestRelionExportParticles(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('xmipp_tutorial')
        cls.particlesFn = cls.ds.getFile('particles')
        cls.runImportParticles(cls.particlesFn, 1.237, True)

    @classmethod
    def runImportParticles(cls, pattern, samplingRate, checkStack=False,
                           phaseFlip=False):
        """ Run an Import particles protocol. """
        cls.protImport = cls.newProtocol(ProtImportParticles,
                                         filesPath=pattern,
                                         samplingRate=samplingRate,
                                         checkStack=checkStack,
                                         haveDataBeenPhaseFlipped=phaseFlip)
        print('_label: ', cls.protImport._label)
        cls.launchProtocol(cls.protImport)
        # check that input images have been imported (a better way to do this?)
        if cls.protImport.outputParticles is None:
            raise Exception('Import of images: %s, failed. outputParticles is None.' % pattern)
        return cls.protImport

    def test_basic(self):
        """ Run an Import particles protocol. """

        inputParts = self.protImport.outputParticles

        stackNames = set(pwutils.removeBaseExt(p.getFileName()) for p in inputParts)

        paramsList = [
            {'stackType': 0, 'useAlignment': True},
            {'stackType': 0, 'useAlignment': False},
            {'stackType': 1, 'useAlignment': True},
            {'stackType': 1, 'useAlignment': False},
            {'stackType': 2, 'useAlignment': True},
            {'stackType': 2, 'useAlignment': False}
        ]

        def _checkProt(prot, params):
            stackFiles = glob(prot._getPath('Particles', '*mrcs'))
            n = len(stackFiles)
            if params['stackType'] == 0:
                self.assertEqual(n, 0)
            elif params['stackType'] == 1:
                self.assertGreaterEqual(n, 1)
            else:
                self.assertEqual(n, 1)

        for i, params in enumerate(paramsList):
            exportProt = self.newProtocol(ProtRelionExportParticles,
                                          objectLabel='export %d' % (i+1),
                                          **params)
            exportProt.inputParticles.set(inputParts)
            self.launchProtocol(exportProt)
            _checkProt(exportProt, params)

