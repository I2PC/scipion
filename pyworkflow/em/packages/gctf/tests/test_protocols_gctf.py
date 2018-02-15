# **************************************************************************
# *
# * Authors:    Grigory Sharov (sharov@igbmc.fr)
# *
# * L'Institut de genetique et de biologie moleculaire et cellulaire (IGBMC)
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

import unittest
import sys

from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.gctf import *


class TestGctfBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='xmipp_tutorial'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.micFn = cls.dataset.getFile('allMics')
        cls.partFn1 = cls.dataset.getFile('particles2')
        cls.partFn2 = cls.dataset.getFile('particles3')
        cls.ctfFn = cls.dataset.getFile('ctf')

    @classmethod
    def runImportMicrograph(cls, pattern, samplingRate, voltage,
                            scannedPixelSize, magnification,
                            sphericalAberration):
        """ Run an Import micrograph protocol. """
        # We have two options: passe the SamplingRate or the
        # ScannedPixelSize + microscope magnification
        if samplingRate is not None:
            cls.protImport = ProtImportMicrographs(
                objLabel='import mics',
                samplingRateMode=0, filesPath=pattern,
                samplingRate=samplingRate, magnification=magnification,
                voltage=voltage, sphericalAberration=sphericalAberration)
        else:
            cls.protImport = ProtImportMicrographs(
                objLabel='import mics',
                samplingRateMode=1, filesPath=pattern,
                scannedPixelSize=scannedPixelSize,
                voltage=voltage, magnification=magnification,
                sphericalAberration=sphericalAberration)

        cls.proj.launchProtocol(cls.protImport, wait=True)
        # check that input micrographs have been imported
        # (a better way to do this?)
        if cls.protImport.outputMicrographs is None:
            raise Exception('Import of micrograph: %s, failed. '
                            'outputMicrographs is None.' % pattern)
        return cls.protImport

    @classmethod
    def runImportParticles(cls, pattern, label, samplingRate, voltage, magnification):
        """ Run an Import particles protocol. """
        cls.protImport = cls.newProtocol(ProtImportParticles,
                                         objLabel=label,
                                         sqliteFile=pattern,
                                         samplingRate=samplingRate,
                                         voltage=voltage,
                                         magnification=magnification,
                                         importFrom=ProtImportParticles.IMPORT_FROM_SCIPION)

        cls.proj.launchProtocol(cls.protImport, wait=True)
        # check that input particles have been imported
        # (a better way to do this?)
        if cls.protImport.outputParticles is None:
            raise Exception('Import of particles: %s, failed. '
                            'outputParticles is None.' % pattern)
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
    def runImportParticlesBPV(cls, pattern, label):
        """ Run an Import particles protocol. """
        return cls.runImportParticles(pattern,
                                      label=label,
                                      samplingRate=4.95,
                                      voltage=300,
                                      magnification=56000)


class TestGctf(TestGctfBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestGctfBase.setData()
        cls.protImport = cls.runImportMicrographBPV(cls.micFn)

    def testRunGctf(self):
        protCTF = ProtGctf()
        protCTF.inputMicrographs.set(self.protImport.outputMicrographs)
        protCTF.ctfDownFactor.set(2)
        self.proj.launchProtocol(protCTF, wait=True)
        self.assertIsNotNone(protCTF.outputCTF,
                             "SetOfCTF has not been produced.")

        valuesList = [[23956, 23537, 52, 5, 0.3],
                      [22319, 21973, 52, 5, 0.3],
                      [22657, 22362, 56, 5, 0.3]]
        for ctfModel, values in izip(protCTF.outputCTF, valuesList):
            self.assertAlmostEquals(
                ctfModel.getDefocusU(), values[0], delta=1000)
            self.assertAlmostEquals(
                ctfModel.getDefocusV(), values[1], delta=1000)
            self.assertAlmostEquals(
                ctfModel.getDefocusAngle(), values[2], delta=5)
            self.assertAlmostEquals(
                ctfModel.getMicrograph().getSamplingRate(),
                2.474, delta=0.001)
            self.assertAlmostEquals(
                ctfModel.getResolution(), values[3], delta=1)
            self.assertAlmostEquals(
                ctfModel.getFitQuality(), values[4], delta=.1)

    def testRunGctf2(self):
        protCTF = ProtGctf()
        protCTF.inputMicrographs.set(self.protImport.outputMicrographs)
        self.proj.launchProtocol(protCTF, wait=True)
        self.assertIsNotNone(protCTF.outputCTF,
                             "SetOfCTF has not been produced.")

        valuesList = [[23885, 23552, 56, 3, 0.38],
                      [22250, 21892, 52, 4, 0.38],
                      [22482, 22268, 39, 4, 0.38]]
        for ctfModel, values in izip(protCTF.outputCTF, valuesList):
            self.assertAlmostEquals(
                ctfModel.getDefocusU(), values[0], delta=1000)
            self.assertAlmostEquals(
                ctfModel.getDefocusV(), values[1], delta=1000)
            self.assertAlmostEquals(
                ctfModel.getDefocusAngle(), values[2], delta=5)
            self.assertAlmostEquals(
                ctfModel.getMicrograph().getSamplingRate(),
                1.237, delta=0.001)
            self.assertAlmostEquals(
                ctfModel.getResolution(), values[3], delta=1)
            self.assertAlmostEquals(
                ctfModel.getFitQuality(), values[4], delta=.1)


class TestGctfRefine(TestGctfBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestGctfBase.setData()
        cls.protImport1 = cls.runImportParticlesBPV(cls.partFn1,
                                                    label='import particles (no alignment)')
        cls.protImport2 = cls.runImportParticlesBPV(cls.partFn2,
                                                    label='import particles (with alignment)')
        cls.protImportMics = cls.runImportMicrographBPV(cls.micFn)

    def testRunGctf1(self):
        protCTF = ProtGctfRefine(objLabel='gCTF local refinement')
        protCTF.inputParticles.set(self.protImport1.outputParticles)
        protCTF.inputMicrographs.set(self.protImportMics.outputMicrographs)
        protCTF.ctfDownFactor.set(2)

        self.proj.launchProtocol(protCTF, wait=True)
        self.assertIsNotNone(protCTF.outputParticles,
                             "SetOfParticles has not been produced.")
        self.assertEqual(protCTF.inputParticles.get().getSize(),
                         protCTF.outputParticles.getSize())

    def testRunGctf2(self):
        protCTF = ProtGctfRefine(objLabel='gCTF local refinement (with input CTFs)')
        protCTF.inputParticles.set(self.protImport2.outputParticles)
        protCTF.ctfDownFactor.set(2)
        protCTF.useInputCtf.set(True)
        protCTF.applyShifts.set(True)
        protCTF.inputMicrographs.set(self.protImportMics.outputMicrographs)

        # run import CTF protocol
        protImportCTF = self.newProtocol(ProtImportCTF,
                                         objLabel='import CTF from scipion',
                                         importFrom=ProtImportCTF.IMPORT_FROM_SCIPION,
                                         filesPath=self.ctfFn)
        protImportCTF.inputMicrographs.set(self.protImportMics.outputMicrographs)
        self.launchProtocol(protImportCTF)
        self.assertIsNotNone(protImportCTF.outputCTF,
                             "There was a problem when importing ctfs.")

        protCTF.ctfRelations.set(protImportCTF.outputCTF)

        self.proj.launchProtocol(protCTF, wait=True)
        self.assertIsNotNone(protCTF.outputParticles,
                             "SetOfParticles has not been produced.")
        self.assertEqual(protCTF.inputParticles.get().getSize(),
                         protCTF.outputParticles.getSize())


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestGctf)
    suite = unittest.TestLoader().loadTestsFromName(
        'test_protocols_gctf.TestGctf.testRunGctf')
    unittest.TextTestRunner(verbosity=2).run(suite)
