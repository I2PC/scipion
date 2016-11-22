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
from pyworkflow.em.packages.xmipp3 import XmippProtAddNoise
from pyworkflow.em.protocol import ProtImportVolumes


class TestAddNoiseBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='resmap'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.map3D = cls.dataset.getFile('betagal')

    @classmethod
    def runImportVolumes(cls, pattern, samplingRate):
        """ Run an Import volumes protocol. """
        cls.protImport = cls.newProtocol(ProtImportVolumes,
                                         filesPath=pattern,
                                         samplingRate=samplingRate
                                        )
        cls.launchProtocol(cls.protImport)
        return cls.protImport


class TestAddNoise(TestAddNoiseBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestAddNoiseBase.setData()
        cls.protImportVol  = cls.runImportVolumes(cls.map3D, 3.54)

    def testAddNoise1(self):
        addnoise = self.newProtocol(XmippProtAddNoise,
                                    input = self.protImportVol.outputVolume,
                                    noiseType = 0,
                                    gaussianStd = 0.08,
                                    gaussianMean = 0
                                    )
        self.launchProtocol(addnoise)
        self.assertTrue(exists(addnoise._getExtraPath('Noisy.vol')),
                         "AddNoise with gaussian noise has failed")

    def testAddNoise2(self):
        addnoise = self.newProtocol(XmippProtAddNoise,
                                    input = self.protImportVol.outputVolume,
                                    noiseType = 1,
                                    studentDf = 1,
                                    studentStd = 0.08,
                                    studentMean = 0
                                  )
        self.launchProtocol(addnoise)
        self.assertTrue(exists(addnoise._getExtraPath('Noisy.vol')), 
                        "AddNoise with student noise has failed")

    def testAddNoise3(self):
        addnoise = self.newProtocol(XmippProtAddNoise,
                                    input = self.protImportVol.outputVolume,
                                    noiseType = 2,
                                    uniformMax = 1,
                                    uniformMin = 0
                                  )
        self.launchProtocol(addnoise)
        self.assertTrue(exists(addnoise._getExtraPath('Noisy.vol')), 
                        "AddNoise with uniform noise has failed")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        className = sys.argv[1]
        cls = globals().get(className, None)
        if cls:
            suite = unittest.TestLoader().loadTestsFromTestCase(cls)
            unittest.TextTestRunner(verbosity=2).run(suite)
        else:
            print "Test: '%s' not found." % className
    else:
        unittest.main()
