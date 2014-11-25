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
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.brandeis import *


class TestBrandeisBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='xmipp_tutorial'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.micFn = cls.dataset.getFile('micrographs/BPV_1386.mrc')
        cls.volFn = cls.dataset.getFile('volumes/volume_1_iter_002.mrc')
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
    def runImportVolumes(cls, pattern, samplingRate):
        """ Run an Import particles protocol. """
        cls.protImport = cls.newProtocol(ProtImportVolumes,
                                         filesPath=pattern, samplingRate=samplingRate)
        cls.launchProtocol(cls.protImport)
        return cls.protImport

    @classmethod
    def runImportParticles(cls, pattern, samplingRate, checkStack=False):
        """ Run an Import particles protocol. """
        cls.protImport = cls.newProtocol(ProtImportParticles,
                                         filesPath=pattern, samplingRate=samplingRate,
                                         checkStack=checkStack)
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
    def runImportParticleBPV(cls, pattern):
        """ Run an Import micrograph protocol. """
        return cls.runImportParticles(pattern,
                                      samplingRate=1.237,
                                      checkStack=True)
    @classmethod
    def runImportVolumesBPV(cls, pattern):
        """ Run an Import micrograph protocol. """
        return cls.runImportVolumes(pattern,
                                    samplingRate=1.237)


class TestBrandeisCtffind(TestBrandeisBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestBrandeisBase.setData()
        cls.protImport = cls.runImportMicrographBPV(cls.micFn)
    
    def testCtffind(self):
        protCTF = ProtCTFFind()
        protCTF.inputMicrographs.set(self.protImport.outputMicrographs)
        self.proj.launchProtocol(protCTF, wait=True)
        self.assertIsNotNone(protCTF.outputCTF, "SetOfCTF has not been produced.")
        ctfModel = protCTF.outputCTF.getFirstItem()
        self.assertAlmostEquals(ctfModel.getDefocusU(),23873.5, places=1)
        self.assertAlmostEquals(ctfModel.getDefocusV(),23640.28, places=1)
        self.assertAlmostEquals(ctfModel.getDefocusAngle(),64.08, places=2)

class TestBrandeisFrealign(TestBrandeisBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestBrandeisBase.setData()
        particlesPattern   = cls.dataset.getFile('particles/BPV_####_ptcls_64.spi')
        cls.protImportPart = cls.runImportParticleBPV(particlesPattern)
        cls.protImportVol  = cls.runImportVolumesBPV(cls.volFn)

    def testFrealign(self):
        frealign = self.newProtocol(ProtFrealign,
                                    useInitialAngles=True,
                                    mode=MOD_RECONSTRUCTION,
                                    innerRadius=0,
                                    outerRadius=32.,
                                    symmetry='i1'
                                    )
        frealign.inputParticles.set(self.protImportPart.outputParticles)
        frealign.input3DReference.set(self.protImportVol.outputVolume)
        self.launchProtocol(frealign)

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestBrandeisCtffind)
    suite = unittest.TestLoader().loadTestsFromName('test_protocols_brandeis.TestBrandeisCtffind.testCtffind')
    unittest.TextTestRunner(verbosity=2).run(suite)