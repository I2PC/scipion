# **************************************************************************
# *
# * Authors:    Jose Luis Vilas (jlvilas@cnb.csic.es)
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

import unittest, sys

from pyworkflow.em import exists
from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pyworkflow.em.packages.xmipp3 import XmippProtAddNoiseVolumes, XmippProtAddNoiseParticles
from pyworkflow.em.protocol import ProtImportVolumes, ProtImportParticles


class TestAddNoiseBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='resmap'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.map3D = cls.dataset.getFile('betagal')
        cls.setVols = cls.dataset.getFile('*.mrc')
        cls.dsParticles = DataSet.getDataSet('xmipp_tutorial')

    @classmethod
    def runImportVolumes(cls, pattern, samplingRate):
        """ Run an Import volumes protocol. """
        cls.protImport = cls.newProtocol(ProtImportVolumes,
                                         filesPath=pattern,
                                         samplingRate=samplingRate
                                        )
        cls.launchProtocol(cls.protImport)
        return cls.protImport
    
    @classmethod
    def runImportParticles(cls):
        """ Import Particles.
        """
        args = {'importFrom': ProtImportParticles.IMPORT_FROM_FILES,
                'filesPath': cls.dsParticles.getFile('particles/'),
                'filesPattern': 'BPV_????_ptcls.hdf',
                'amplitudConstrast': 0.1,
                'sphericalAberration': 2.,
                'voltage': 100,
                'samplingRate': 2.1,
                'haveDataBeenPhaseFlipped': True
                }

        # Id's should be set increasing from 1 if ### is not in the 
        # pattern
        protMicImport = cls.newProtocol(ProtImportParticles, **args)
        protMicImport.setObjLabel('import particles')
        cls.launchProtocol(protMicImport)
        return protMicImport


class TestAddNoise(TestAddNoiseBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestAddNoiseBase.setData()
        cls.protImportVol  = cls.runImportVolumes(cls.map3D, 3.54)
        cls.protImportVols = cls.runImportVolumes(cls.setVols, 3.54)
        cls.protImportParts = cls.runImportParticles()

    def testAddNoise1(self):
        addnoise = self.newProtocol(XmippProtAddNoiseVolumes,
                                    input = self.protImportVol.outputVolume,
                                    noiseType = XmippProtAddNoiseVolumes.GAUSSIAN_NOISE,
                                    gaussianStd = 0.08,
                                    gaussianMean = 0
                                    )
        self.launchProtocol(addnoise)
        self.assertTrue(exists(
                        addnoise._getNoisyOutputPath(
                        self.protImportVol.outputVolume.getFileName())),
                         "AddNoise with gaussian noise has failed")

    def testAddNoise2(self):
        addnoise = self.newProtocol(XmippProtAddNoiseVolumes,
                                    input = self.protImportVols.outputVolumes,
                                    noiseType = XmippProtAddNoiseVolumes.STUDENT_NOISE,
                                    studentDf = 1,
                                    studentStd = 0.08,
                                    studentMean = 0
                                  )
        self.launchProtocol(addnoise)
        fnVols = self.protImportVols.outputVolumes.getFiles() 
        for fnvol in  fnVols:
            fn = addnoise._getNoisyOutputPath(fnvol)
            self.assertTrue(exists(fn), "AddNoise with student noise has failed %s" % (fn))


    def testAddNoise3(self):
        addnoise = self.newProtocol(XmippProtAddNoiseParticles,
                                    input = self.protImportParts.outputParticles,
                                    noiseType = XmippProtAddNoiseParticles.UNIFORM_NOISE,
                                    uniformMax = 1,
                                    uniformMin = 0
                                  )
        self.launchProtocol(addnoise)
        self.assertTrue(exists(addnoise._getExtraPath('Noisy.stk')), 
                        "AddNoise with uniform noise has failed")
