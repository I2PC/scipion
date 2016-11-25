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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

import unittest, sys
# import numpy as np
from pyworkflow.em import exists 
from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pyworkflow.em.packages.xmipp3 import XmippProtMonoRes
from pyworkflow.em.protocol import ProtImportVolumes, ProtImportMask


class TestMonoResBase(BaseTest):
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


class TestMonoRes(TestMonoResBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestMonoResBase.setData()
        cls.protImportVol  = cls.runImportVolumes(cls.map3D, 3.54)
        cls.protImportHalf1  = cls.runImportVolumes(cls.half1, 3.54)
        cls.protImportHalf2  = cls.runImportVolumes(cls.half2, 3.54)
        cls.protImportMask  = cls.runImportMask(cls.mask, 3.54)

    def testMonoRes1(self):
        MonoRes = self.newProtocol(XmippProtMonoRes,
                                   halfVolums = False,
                                   inputVolume = self.protImportVol.outputVolume,
                                   Mask = self.protImportMask.outputVolume,
                                   symmetry = 'd2',
                                   minRes = 1,
                                   maxRes = 100,
                                   significance = 0.95,
                                   exact = False,
                                   filterInput = False,
                                   trimming = False,
                                    )
        self.launchProtocol(MonoRes)
        self.assertTrue(exists(MonoRes._getExtraPath("MGresolution.vol")), 
			"MonoRes (no split and no mask) has failed")

    def testMonoRes2(self):
        MonoRes = self.newProtocol(XmippProtMonoRes,
                                  halfVolums=True,
                                  inputVolume = self.protImportHalf1.outputVolume,
                                  inputVolume2 = self.protImportHalf2.outputVolume,
                                  provideMaskInHalves = True,
                                  Mask = self.protImportMask.outputVolume,
                                  symmetry = 'd2',
                                  minRes = 1,
                                  maxRes = 100,
                                  significance = 0.95,
                                  exact = False,
                                  filterInput = False,
                                  trimming = False,
                                  )
        self.launchProtocol(MonoRes)
        self.assertTrue(exists(MonoRes._getExtraPath("MGresolution.vol")), 
			"MonoRes (split and no mask) has failed")

    def testMonoRes3(self):
        MonoRes = self.newProtocol(XmippProtMonoRes,
                                  halfVolums = True,
                                  inputVolume = self.protImportVol.outputVolume,
                                  provideMaskInHalves = False,
                                  Mask = '',
                                  symmetry = 'd2',
                                  minRes = 1,
                                  maxRes = 100,
                                  significance = 0.95,
                                  exact = False,
                                  filterInput = False,
                                  trimming = False,
                                  )
        self.launchProtocol(MonoRes)
        self.assertTrue(exists(MonoRes._getExtraPath("MGresolution.vol")), 
			"MonoRes (no split and no mask) has failed")

    def testMonoRes4(self):
        MonoRes = self.newProtocol(XmippProtMonoRes,
                                   halfVolums = False,
                                   inputVolume = self.protImportVol.outputVolume,
                                   Mask = self.protImportMask.outputVolume,
                                   symmetry = 'd2',
                                   minRes = 1,
                                   maxRes = 100,
                                   significance = 0.95,
                                   exact = True,
                                   filterInput = True,
                                   trimming = True,
                                   kValue = 5
                                  )
        self.launchProtocol(MonoRes)
        self.assertTrue(exists(MonoRes._getExtraPath("MGresolution.vol")), 
			"MonoRes (split and no mask) has failed")


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
