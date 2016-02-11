# **************************************************************************
# *
# * Authors:    J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

from pyworkflow.em import *
from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pyworkflow.em.packages.resmap import ProtResMap
from pyworkflow.em.protocol import ProtImportVolumes, ProtImportMask
from pyworkflow.em.packages.xmipp3 import XmippProtMultipleFSCs



class TestMultipleFSCsBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='resmap'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.map3D = cls.dataset.getFile('betagal')
        cls.half1 = cls.dataset.getFile('betagal_half1')
        cls.half2 = cls.dataset.getFile('betagal_half2')
        cls.mask = cls.dataset.getFile('betagal_mask')

    @classmethod
    def runImportVolumes(cls, pattern, samplingRate, label):
        """ Run an Import volumes protocol. """
        cls.protImport = cls.newProtocol(ProtImportVolumes,
                                         objLabel=label,
                                         filesPath=pattern,
                                         samplingRate=samplingRate
                                        )
        cls.launchProtocol(cls.protImport)
        return cls.protImport
    
    @classmethod
    def runImportMask(cls, pattern, samplingRate, label):
        """ Run an Import volumes protocol. """
        cls.protImport = cls.newProtocol(ProtImportMask,
                                         objLabel=label,
                                         maskPath=pattern,
                                         samplingRate=samplingRate
                                        )
        cls.launchProtocol(cls.protImport)
        return cls.protImport


class TestMultipleFSCs(TestMultipleFSCsBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestMultipleFSCsBase.setData()
        cls.protImportVol  = cls.runImportVolumes(cls.map3D, 3.54,
                                                  'import vol')
        cls.protImportHalf1  = cls.runImportVolumes(cls.half1, 3.54,
                                                    'import half1')
        cls.protImportHalf2  = cls.runImportVolumes(cls.half2, 3.54,
                                                    'import half2')
        cls.protImportMask  = cls.runImportMask(cls.mask, 3.54,
                                                'import mask')

    def _runFSC(self, useMask):
        prot = self.newProtocol(XmippProtMultipleFSCs,
                                referenceVolume=self.protImportVol.outputVolume)
        prot.inputVolumes.append(self.protImportHalf1.outputVolume)
        prot.inputVolumes.append(self.protImportHalf2.outputVolume)

        if useMask:
            prot.mask.set(self.protImportMask.outputMask)
        self.launchProtocol(prot)

    def test_case1(self):
        """  Compute the FSC without mask       """
        self._runFSC(useMask=True)

    def test_case2(self):
        """  Compute the FSC without mask       """
        self._runFSC(useMask=False)
