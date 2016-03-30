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

import unittest, sys
import numpy as np

from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.grigoriefflab import *
from pyworkflow.em.protocol import ProtImportParticles, ProtImportVolumes


class TestBrandeisBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='xmipp_tutorial'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.micFn = cls.dataset.getFile('allMics')
        cls.volFn = cls.dataset.getFile('vol2')
        #cls.parFn = cls.dataset.getFile('aligned_particles')

    @classmethod
    def runImportMicrograph(cls, pattern, samplingRate, voltage, scannedPixelSize, magnification, sphericalAberration):
        """ Run an Import micrograph protocol. """
        # We have two options: passe the SamplingRate or the ScannedPixelSize + microscope magnification
        if not samplingRate is None:
            cls.protImport = ProtImportMicrographs(samplingRateMode=0, filesPath=pattern, samplingRate=samplingRate, magnification=magnification, 
                                                   voltage=voltage, sphericalAberration=sphericalAberration)
        else:
            cls.protImport = ProtImportMicrographs(samplingRateMode=1, filesPath=pattern, scannedPixelSize=scannedPixelSize, 
                                                   voltage=voltage, magnification=magnification, sphericalAberration=sphericalAberration)
        
        cls.proj.launchProtocol(cls.protImport, wait=True)
        # check that input micrographs have been imported (a better way to do this?)
        if cls.protImport.outputMicrographs is None:
            raise Exception('Import of micrograph: %s, failed. outputMicrographs is None.' % pattern)
        return cls.protImport

    @classmethod
    def runImportVolumes(cls, pattern, samplingRate,
                         importFrom = ProtImportParticles.IMPORT_FROM_FILES):
        """ Run an Import particles protocol. """
        cls.protImport = cls.newProtocol(ProtImportVolumes,
                                         filesPath=pattern,
                                         samplingRate=samplingRate
                                        )
                                         #not yet implemented
                                         #importFrom=importFrom
                                         #)
        cls.launchProtocol(cls.protImport)
        return cls.protImport

    @classmethod
    def runImportParticles(cls, pattern, samplingRate, checkStack=False,
                           importFrom = ProtImportParticles.IMPORT_FROM_FILES):
        """ Run an Import particles protocol. """
        if importFrom == ProtImportParticles.IMPORT_FROM_SCIPION:
            objLabel = 'from scipion (particles)'
        elif importFrom == ProtImportParticles.IMPORT_FROM_FILES:
            objLabel = 'from file (particles)'


        cls.protImport = cls.newProtocol(ProtImportParticles,
                                         objLabel=objLabel,
                                         filesPath=pattern,
                                         sqliteFile=pattern,
                                         samplingRate=samplingRate,
                                         checkStack=checkStack,
                                         importFrom=importFrom)

        cls.launchProtocol(cls.protImport)
        # check that input images have been imported (a better way to do this?)
        if cls.protImport.outputParticles is None:
            raise Exception('Import of images: %s, failed. outputParticles is None.' % pattern)
        return cls.protImport


    @classmethod
    def runImportMicrographBPV(cls, pattern):
        """ Run an Import micrograph protocol. """
        return cls.runImportMicrograph(pattern,
                                       samplingRate=1.237,
                                       voltage=300,
                                       sphericalAberration=2,
                                       scannedPixelSize=None,
                                       magnification=56000)

    @classmethod
    def runImportParticleGrigorieff(cls, pattern):
        """ Run an Import micrograph protocol. """
        return cls.runImportParticles(pattern,
                                      samplingRate=4.,
                                      checkStack=True,
                            importFrom=ProtImportParticles.IMPORT_FROM_SCIPION)
    @classmethod
    def runImportVolumesGrigorieff(cls, pattern):
        """ Run an Import micrograph protocol. """
        return cls.runImportVolumes(pattern,
                                    samplingRate=4.,
                                    importFrom=ProtImportParticles.IMPORT_FROM_FILES)


class TestImportParticles(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('grigorieff')

    def test_import(self):
        parFile = self.dataset.getFile('particles/particles_iter_002.par')
        stackFile = self.dataset.getFile('particles/particles.mrc')
        
        protImport = self.newProtocol(ProtImportParticles,
                                         objLabel='import parfile & stk',
                                         parFile=parFile,
                                         stackFile=stackFile,
                                         samplingRate=9.90,
                                         importFrom=ProtImportParticles.IMPORT_FROM_FREALIGN)

        self.launchProtocol(protImport)      
        # check that input images have been imported (a better way to do this?)
        if protImport.outputParticles is None:
            raise Exception('Import failed. Par file: %s' % parFile)
        
        self.assertTrue(protImport.outputParticles.getSize() == 180)
        
        goldFile = self.dataset.getFile('particles/particles.sqlite')
        goldSet = SetOfParticles(filename=goldFile)
        
        for p1, p2 in izip(goldSet, protImport.outputParticles.iterItems(orderBy=['_micId', 'id'],
                                                                     direction='ASC')):
            m1 = p1.getTransform().getMatrix()
            m2 = p2.getTransform().getMatrix()
            self.assertTrue(np.allclose(m1, m2, atol=0.01))
            
        self.assertTrue(protImport.outputParticles.hasCTF())
        self.assertTrue(protImport.outputParticles.hasAlignmentProj())

    
class TestBrandeisCtffind(TestBrandeisBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestBrandeisBase.setData()
        cls.protImport = cls.runImportMicrographBPV(cls.micFn)
    
    def testCtffind(self):
        protCTF = ProtCTFFind(useCftfind4=False)
        protCTF.inputMicrographs.set(self.protImport.outputMicrographs)
        protCTF.ctfDownFactor.set(2)
        protCTF.numberOfThreads.set(4)
        self.proj.launchProtocol(protCTF, wait=True)
        self.assertIsNotNone(protCTF.outputCTF, "SetOfCTF has not been produced.")
        
        valuesList = [[23861, 23664, 56], [22383, 22153, 52.6], [22716, 22526, 59.1]]
        for ctfModel, values in izip(protCTF.outputCTF, valuesList):
            self.assertAlmostEquals(ctfModel.getDefocusU(),values[0], delta=1000)
            self.assertAlmostEquals(ctfModel.getDefocusV(),values[1], delta=1000)
            self.assertAlmostEquals(ctfModel.getDefocusAngle(),values[2], delta=5)
            self.assertAlmostEquals(ctfModel.getMicrograph().getSamplingRate(), 2.474, delta=0.001)

    def testCtffind2(self):
        protCTF = ProtCTFFind(useCftfind4=False)
        protCTF.inputMicrographs.set(self.protImport.outputMicrographs)
        protCTF.numberOfThreads.set(4)
        self.proj.launchProtocol(protCTF, wait=True)
        self.assertIsNotNone(protCTF.outputCTF, "SetOfCTF has not been produced.")
        
        valuesList = [[23863, 23640, 64], [22159, 21983, 50.6], [22394, 22269, 45]]
        for ctfModel, values in izip(protCTF.outputCTF, valuesList):
            self.assertAlmostEquals(ctfModel.getDefocusU(),values[0], delta=1000)
            self.assertAlmostEquals(ctfModel.getDefocusV(),values[1], delta=1000)
            self.assertAlmostEquals(ctfModel.getDefocusAngle(),values[2], delta=5)
            self.assertAlmostEquals(ctfModel.getMicrograph().getSamplingRate(), 1.237, delta=0.001)


class TestBrandeisCtffind4(TestBrandeisBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestBrandeisBase.setData()
        cls.protImport = cls.runImportMicrographBPV(cls.micFn)
    
    def testCtffind4V1(self):
        protCTF = ProtCTFFind()
        protCTF.inputMicrographs.set(self.protImport.outputMicrographs)
        protCTF.ctfDownFactor.set(2)
        protCTF.numberOfThreads.set(4)
        self.proj.launchProtocol(protCTF, wait=True)
        self.assertIsNotNone(protCTF.outputCTF, "SetOfCTF has not been produced.")
        
        valuesList = [[23861, 23664, 53], [22383, 22153, 48.5], [22716, 22526, 54.3]]
        for ctfModel, values in izip(protCTF.outputCTF, valuesList):
            self.assertAlmostEquals(ctfModel.getDefocusU(),values[0], delta=1000)
            self.assertAlmostEquals(ctfModel.getDefocusV(),values[1], delta=1000)
            self.assertAlmostEquals(ctfModel.getDefocusAngle(),values[2], delta=5)
            self.assertAlmostEquals(ctfModel.getMicrograph().getSamplingRate(), 2.474, delta=0.001)

    def testCtffind4V22(self):
        protCTF = ProtCTFFind()
        protCTF.inputMicrographs.set(self.protImport.outputMicrographs)
        protCTF.numberOfThreads.set(4)
        self.proj.launchProtocol(protCTF, wait=True)
        self.assertIsNotNone(protCTF.outputCTF, "SetOfCTF has not been produced.")
        
        valuesList = [[23863, 23640, 54], [22159, 21983, 45.8], [22394, 22269, 171]]
        for ctfModel, values in izip(protCTF.outputCTF, valuesList):
            self.assertAlmostEquals(ctfModel.getDefocusU(),values[0], delta=1000)
            self.assertAlmostEquals(ctfModel.getDefocusV(),values[1], delta=1000)
            self.assertAlmostEquals(ctfModel.getDefocusAngle(),values[2], delta=5)
            self.assertAlmostEquals(ctfModel.getMicrograph().getSamplingRate(), 1.237, delta=0.001)


class TestFrealignRefine(TestBrandeisBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        dataProject='grigorieff'
        dataset = DataSet.getDataSet(dataProject)
        TestBrandeisBase.setData()
        particlesPattern   = dataset.getFile('particles.sqlite')
        volFn              = dataset.getFile('ref_volume.vol')
        cls.protImportPart = cls.runImportParticleGrigorieff(particlesPattern)
        cls.protImportVol  = cls.runImportVolumesGrigorieff(volFn)

    def testFrealign(self):
        frealign = self.newProtocol(ProtFrealign,
                                    inputParticles = self.protImportPart.outputParticles,
                                    doInvert=False,
                                    input3DReference = self.protImportVol.outputVolume,
                                    useInitialAngles=True,
                                    mode=MOD_REFINEMENT,
                                    innerRadius=0.,
                                    outerRadius=241.,
                                    symmetry='C1',
                                    numberOfThreads=4,
                                    numberOfIterations=2,
                                    doWienerFilter=False,
                                    resolution=2.,
                                    highResolRefine=2.,
                                    resolClass=2.,
                                    writeMatchProjections=False,
                                    score=0,
                                    )
        frealign.inputParticles.set(self.protImportPart.outputParticles)
        frealign.input3DReference.set(self.protImportVol.outputVolume)
        self.launchProtocol(frealign)
        print "Assert is missing: testFrealign"


class TestFrealignClassify(TestBrandeisBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        dataProject='grigorieff'
        dataset = DataSet.getDataSet(dataProject)
        TestBrandeisBase.setData()
        particlesPattern   = dataset.getFile('particles.sqlite')
        volFn              = dataset.getFile('ref_volume.vol')
        cls.protImportPart = cls.runImportParticleGrigorieff(particlesPattern)
        cls.protImportVol  = cls.runImportVolumesGrigorieff(volFn)

    def testFrealignClassify(self):
        frealign = self.newProtocol(ProtFrealignClassify,
                                    inputParticles = self.protImportPart.outputParticles,
                                    doInvert=False,
                                    input3DReference = self.protImportVol.outputVolume,
                                    numberOfIterations=3,
                                    itRefineAngles = 2,
                                    itRefineShifts = 3,
                                    numberOfClasses = 2,
                                    useInitialAngles=True,
                                    mode=MOD_REFINEMENT,
                                    innerRadius=0.,
                                    outerRadius=241.,
                                    symmetry='C1',
                                    numberOfThreads=4,
                                    doWienerFilter=False,
                                    resolution=2.,
                                    highResolRefine=2.,
                                    resolClass=2.,
                                    writeMatchProjections=False,
                                    score=0,
                                    )
        frealign.inputParticles.set(self.protImportPart.outputParticles)
        frealign.input3DReference.set(self.protImportVol.outputVolume)
        self.launchProtocol(frealign)
        print "Assert is missing: testFrealign"


class TestSummovie(BaseTest):
    @classmethod
    def setData(cls):
        cls.ds = DataSet.getDataSet('movies')

    @classmethod
    def runImportMovies(cls, pattern, **kwargs):
        """ Run an Import micrograph protocol. """
        params = {'samplingRate': 1.14,
                  'voltage': 300,
                  'sphericalAberration': 2.7,
                  'magnification': 50000,
                  'scannedPixelSize': None,
                  'filesPattern': pattern
                  }
        if 'samplingRate' not in kwargs:
            del params['samplingRate']
            params['samplingRateMode'] = 0
        else:
            params['samplingRateMode'] = 1

        params.update(kwargs)

        protImport = cls.newProtocol(ProtImportMovies, **params)
        cls.launchProtocol(protImport)
        return protImport
    
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.setData()
        cls.protImport1 = cls.runImportMovies(cls.ds.getFile('qbeta/qbeta.mrc'),
                                              magnification=50000)
        cls.protImport2 = cls.runImportMovies(cls.ds.getFile('cct/cct_1.em'),
                                              magnification=61000)

    def _checkMicrographs(self, protocol):
        self.assertIsNotNone(getattr(protocol, 'outputMicrographs', None),
                             "Output SetOfMicrographs were not created.")

#     def _checkAlignment(self, movie, goldRange, goldRoi):
#         alignment = movie.getAlignment()
#         range = alignment.getRange()
#         msgRange = "Alignment range must be %s (%s) and it is %s (%s)"
#         self.assertEqual(goldRange, range, msgRange % (goldRange, range, type(goldRange), type(range)))
#         roi = alignment.getRoi()
#         msgRoi = "Alignment ROI must be %s (%s) and it is %s (%s)"
#         self.assertEqual(goldRoi, roi, msgRoi % (goldRoi, roi, type(goldRoi), type(roi)))

    def test_qbeta(self):
        prot = self.newProtocol(ProtSummovie,
                                exposurePerFrame=1.4)
        prot.inputMovies.set(self.protImport1.outputMovies)
        self.launchProtocol(prot)

        self._checkMicrographs(prot)
#         self._checkAlignment(prot.outputMovies[1],
#                              (1,7), [0, 0, 0, 0])

    def test_cct(self):
        prot = self.newProtocol(ProtSummovie,
                                exposurePerFrame=1.4)
        prot.inputMovies.set(self.protImport2.outputMovies)
        self.launchProtocol(prot)

        self._checkMicrographs(prot)
#         self._checkAlignment(prot.outputMovies[1],
#                              (1,7), [0, 0, 0, 0])

    def test_qbeta_SkipCrop(self):
        from pyworkflow.em.packages.motioncorr import ProtMotionCorr
        
        prot = self.newProtocol(ProtMotionCorr,
                                alignFrame0=2, alignFrameN=2,
                                cropOffsetX=10, cropOffsetY=10)
        prot.inputMovies.set(self.protImport1.outputMovies)
        self.launchProtocol(prot)
        
        protSum = self.newProtocol(ProtSummovie,
                                sumFrame0=2, sumFrameN=2,
                                exposurePerFrame=1.4)
        prot.inputMovies.set(prot.outputMovies)
        self.launchProtocol(protSum)

        self._checkMicrographs(protSum)
#         self._checkAlignment(prot.outputMovies[1],
#                              (3,5), [10, 10, 0, 0])

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestBrandeisCtffind)
    suite = unittest.TestLoader().loadTestsFromName('test_protocols_brandeis.TestBrandeisCtffind.testCtffind')
    unittest.TextTestRunner(verbosity=2).run(suite)