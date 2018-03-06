# **************************************************************************
# *
# * Authors:    Jose Luis Vilas Prieto (jlvilas@cnb.csic.es)
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

from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pyworkflow.em.packages.bsoft.protocol_blocres import BsoftProtBlocres
from pyworkflow.em.protocol import ProtImportVolumes, ProtImportMask



class TestBsoftBlocresBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='resmap'):
        cls.dataset = DataSet.getDataSet(dataProject)
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


class TestBsoftBlocres(TestBsoftBlocresBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestBsoftBlocresBase.setData()
        cls.protImportHalf1 = cls.runImportVolumes(cls.half1, 3.54)
        cls.protImportHalf2 = cls.runImportVolumes(cls.half2, 3.54)
        cls.protCreateMask = cls.runImportMask(cls.mask, 3.54)


    def testBlocres1(self):
        blocres = self.newProtocol(BsoftProtBlocres,
                                   objLabel='blocres',
                                   inputVolume=self.protImportHalf1.outputVolume,
                                   inputVolume2=self.protImportHalf2.outputVolume,
                                   mask=self.protCreateMask.outputMask,
                                   method=True,
                                   box=20,
                                   resolutionCriterion=0,
                                   cutoff=0.5,
                                   step=1,
                                   maxresolution=2,
                                   fill=0,
                                   pad=1,
                                   symmetry='',
                                   smooth=True)
        self.launchProtocol(blocres)
        self.assertIsNotNone(blocres.resolution_Volume,
                        "blocres has failed")
