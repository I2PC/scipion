# **************************************************************************
# *
# * Authors:    Josue Gomez Blanco (jgomez@cnb.csic.es)
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
# import numpy as np

from pyworkflow.em import *
from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pyworkflow.em.packages.resmap import ProtResMap
from pyworkflow.em.protocol import ProtImportVolumes, ProtImportMask


class TestResMapBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='resmap'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.map3D = cls.dataset.getFile('betagal')
        cls.half1 = cls.dataset.getFile('betagal_half1')
        cls.half2 = cls.dataset.getFile('betagal_half2')
        cls.mask = cls.dataset.getFile('betagal_mask')

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
    def runImportMask(cls, pattern, samplingRate):
        """ Run an Import volumes protocol. """
        cls.protImport = cls.newProtocol(ProtImportMask,
                                         maskPath=pattern,
                                         samplingRate=samplingRate
                                        )
        cls.launchProtocol(cls.protImport)
        return cls.protImport


class TestResMap(TestResMapBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestResMapBase.setData()
        cls.protImportVol  = cls.runImportVolumes(cls.map3D, 3.54)
        cls.protImportHalf1  = cls.runImportVolumes(cls.half1, 3.54)
        cls.protImportHalf2  = cls.runImportVolumes(cls.half2, 3.54)
        cls.protImportMask  = cls.runImportMask(cls.mask, 3.54)

    def testResmap1(self):
        resMap = self.newProtocol(ProtResMap,
                                    inputVolume = self.protImportVol.outputVolume,
                                    prewhitenAng = 23.77,
                                    prewhitenRamp = 1,
                                    stepRes = 0.5,
                                    minRes = 7.5,
                                    maxRes = 20
                                    )
        self.launchProtocol(resMap)
        self.assertTrue(exists(resMap._getExtraPath("histogram.png")), "resmap (no split and no mask) has failed")

    def testResmap2(self):
        resMap = self.newProtocol(ProtResMap,
                                  useSplitVolume=True,
                                  volumeHalf1 = self.protImportHalf1.outputVolume,
                                  volumeHalf2 = self.protImportHalf2.outputVolume,
                                  prewhitenAng = 23.77,
                                  prewhitenRamp = 1,
                                  stepRes = 0.5,
                                  minRes = 7.5,
                                  maxRes = 20
                                  )
        self.launchProtocol(resMap)
        self.assertTrue(exists(resMap._getExtraPath("histogram.png")), "resmap (split and no mask) has failed")

    def testResmap3(self):
        resMap = self.newProtocol(ProtResMap,
                                  inputVolume = self.protImportVol.outputVolume,
                                  applyMask = True,
                                  maskVolume = self.protImportMask.outputMask,
                                  prewhitenAng = 23.77,
                                  prewhitenRamp = 1,
                                  stepRes = 0.5,
                                  minRes = 7.5,
                                  maxRes = 20
                                  )
        self.launchProtocol(resMap)
        self.assertTrue(exists(resMap._getExtraPath("histogram.png")), "resmap (no split and no mask) has failed")

    def testResmap4(self):
        resMap = self.newProtocol(ProtResMap,
                                  useSplitVolume=True,
                                  volumeHalf1 = self.protImportHalf1.outputVolume,
                                  volumeHalf2 = self.protImportHalf2.outputVolume,
                                  applyMask = True,
                                  maskVolume = self.protImportMask.outputMask,
                                  prewhitenAng = 23.77,
                                  prewhitenRamp = 1,
                                  stepRes = 0.5,
                                  minRes = 7.5,
                                  maxRes = 20
                                  )
        self.launchProtocol(resMap)
        self.assertTrue(exists(resMap._getExtraPath("histogram.png")), "resmap (split and no mask) has failed")


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
