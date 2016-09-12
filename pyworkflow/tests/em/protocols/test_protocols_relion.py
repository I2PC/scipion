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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

import sys, unittest

from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.relion import *


class TestRelionBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='hemoglobin_mda'):
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
                                      filesPath=pattern, samplingRate=samplingRate, 
                                      checkStack=checkStack)
        cls.launchProtocol(protImport)
        # check that input images have been imported (a better way to do this?)
        if protImport.outputParticles is None:
            raise Exception('Import of images: %s, failed. outputParticles is None.' % pattern)
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
                                     filesPath=pattern, samplingRate=samplingRate)
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
        print "Run relion2D"
        prot2D = self.newProtocol(ProtRelionClassify2D,
                                  doCTF=False, maskDiameterA=340,
                                  numberOfMpi=4, numberOfThreads=1)
        prot2D.numberOfClasses.set(4)
        prot2D.numberOfIterations.set(3)
        prot2D.inputParticles.set(self.protNormalize.outputParticles)
        self.launchProtocol(prot2D)        
        
        self.assertIsNotNone(getattr(prot2D, 'outputClasses', None), 
                             "There was a problem with Relion 2D:\n" + prot2D.getErrorMessage())
        self.assertAlmostEquals(self.protNormalize.outputParticles.getSamplingRate(), 
                                prot2D.outputClasses.getImages().getSamplingRate(),
                                "There was a problem with the sampling rate of the particles")
        for class2D in prot2D.outputClasses:
            self.assertTrue(class2D.hasAlignment2D())


class TestRelionClassify3D(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestRelionBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
#         cls.protNormalize = cls.runNormalizeParticles(cls.protImport.outputParticles)
        cls.protImportVol = cls.runImportVolumes(cls.vol, 3.5)
    
    def testProtRelionClassify3D(self):
        relionNormalize = self.newProtocol(ProtRelionPreprocessParticles)
        relionNormalize.inputParticles.set(self.protImport.outputParticles)
        relionNormalize.doNormalize.set(True)
        self.launchProtocol(relionNormalize)

        print "Run ProtRelionClassify3D"
        relion3DClass = self.newProtocol(ProtRelionClassify3D, 
                                         numberOfClasses=3, numberOfIterations=4, 
                                         doCTF=False, runMode=1, maskDiameterA=320,
                                         numberOfMpi=2, numberOfThreads=2)
        relion3DClass.inputParticles.set(relionNormalize.outputParticles)
        relion3DClass.referenceVolume.set(self.protImportVol.outputVolume)
        self.launchProtocol(relion3DClass)
        
        self.assertIsNotNone(getattr(relion3DClass, 'outputClasses', None), 
                             "There was a problem with Relion 3D:\n" + relion3DClass.getErrorMessage()) 

        for class3D in relion3DClass.outputClasses:
            self.assertTrue(class3D.hasAlignment3D())


class TestRelionRefine(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestRelionBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
#         cls.protNormalize = cls.runNormalizeParticles(cls.protImport.outputParticles)
        cls.protImportVol = cls.runImportVolumes(cls.vol, 3.5)
      
    def testProtRelionRefine(self):
        relionNormalize = self.newProtocol(ProtRelionPreprocessParticles)
        relionNormalize.inputParticles.set(self.protImport.outputParticles)
        relionNormalize.doNormalize.set(True)
        self.launchProtocol(relionNormalize)
  
        print "Run ProtRelionRefine"
        relionRefine = self.newProtocol(ProtRelionRefine3D, 
                                         doCTF=False, runMode=1, memoryPreThreads=1,
                                         maskDiameterA=340, symmetryGroup="d6",
                                         numberOfMpi=3, numberOfThreads=2)
        relionRefine.inputParticles.set(relionNormalize.outputParticles)
        relionRefine.referenceVolume.set(self.protImportVol.outputVolume)
        self.launchProtocol(relionRefine)
        
        self.assertIsNotNone(getattr(relionRefine, 'outputVolume', None), 
                             "There was a problem with Relion 3D:\n" + relionRefine.getErrorMessage())
        
        relionRefine._initialize() # Load filename templates
        dataSqlite =  relionRefine._getIterData(3)
        outImgSet = em.SetOfParticles(filename=dataSqlite)
        self.assertAlmostEqual(outImgSet[1].getSamplingRate(),
                               relionNormalize.outputParticles[1].getSamplingRate(),
                               "The sampling rate is wrong", delta=0.00001)
        
        self.assertAlmostEqual(outImgSet[1].getFileName(),
                               relionNormalize.outputParticles[1].getFileName(),
                               "The particles filenames are wrong")
        
class TestRelionPreprocess(TestRelionBase):
    """ This class helps to test all different preprocessing particles options on RElion. """
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestRelionBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
            
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
        self.assertIsNotNone(protSubtract.outputParticles, "There was a problem with subtract projection")


class TestRelionSortParticles(TestRelionBase):
    """ This class helps to test sort particles protocol from Relion. """

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('relion_tutorial')

        # FIXME
        cls.dataset = DataSet.getDataSet('relion_tutorial')
        cls.partFn = cls.dataset.getFile('import/particles.sqlite')
        cls.partAvg = cls.dataset.getFile('import/averages.mrcs')
        cls.partCl2dFn = cls.dataset.getFile('import/classify2d/extra/relion_it015_data.star')
        cls.partCl3dFn = cls.dataset.getFile('import/classify3d/extra/relion_it015_data.star')
        #FIXME: import from completed relion 3d refine run is not working
        #cls.partRef3dFn = cls.dataset.getFile('import/refine3d/extra/relion_data.star')
        cls.partRef3dFn = cls.dataset.getFile('import/refine3d/extra/relion_it025_data.star')
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
        # FIXME: import from completed relion 3d refine run is not working
        # partRef3dFn = self.ds.getFile('import/refine3d/extra/relion_data.star')
        partRef3dFn = self.ds.getFile(
            'import/refine3d/extra/relion_it025_data.star')
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

