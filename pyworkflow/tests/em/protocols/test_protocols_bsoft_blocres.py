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

import unittest, sys
# import numpy as np
from pyworkflow.em import exists
from pyworkflow.em.packages.bsoft.protocol_blocres import (OUTPUT_RESOLUTION_FILE, CRITERION_FSC,
                                                           CRITERION_DPR, CRITERION_SSNR, CRITERION_RAB,
                                                           PAD_NONE, PAD_BOX, PAD_SHELL)
from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pyworkflow.em.packages.xmipp3 import XmippProtCreateMask3D
from pyworkflow.em.packages.bsoft import BsoftProtBlocres
from pyworkflow.em.protocol import ProtImportVolumes



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
    def runCreateMask(cls, pattern, thr):
        """ Create a volume mask. """
        cls.msk = cls.newProtocol(XmippProtCreateMask3D,
                                  inputVolume=pattern,
                                  volumeOperation=0,  # OPERATION_THRESHOLD,
                                  threshold=thr,
                                  doSmall=True,
                                  doBig=True,
                                  doSymmetrize=False,
                                  doMorphological=False,
                                  doInvert=False,
                                  doSmooth=False,
                                  sigmaConvolution=2
                                  )
        cls.launchProtocol(cls.msk)
        return cls.msk


class TestBsoftBlocres(TestBsoftBlocresBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestBsoftBlocresBase.setData()
        cls.protImportHalf1 = cls.runImportVolumes(cls.half1, 3.54)
        cls.protImportHalf2 = cls.runImportVolumes(cls.half2, 3.54)
        cls.protCreateMask = cls.runCreateMask(cls.protImportHalf1.outputVolume, 0.02)


    def testBlocres1(self):
        blocres = self.newProtocol(BsoftProtBlocres,
                                   objLabel='blocres',
                                   inputVolume=self.protImportHalf1.outputVolume,
                                   inputVolume2=self.protImportHalf2.outputVolume,
                                   mask=self.protCreateMask.outputMask,
                                   box=20,
                                   shell=20,
                                   resolutionCriterion=CRITERION_FSC,
                                   cutoff=0.5,
                                   step=1,
                                   maxresolution=2,
                                   fill=0,
                                   pad=PAD_BOX,
                                   symmetry='',
                                   smooth=True)
        self.launchProtocol(blocres)
        self.assertTrue(exists(blocres._getExtraPath(OUTPUT_RESOLUTION_FILE)),
                        "blocres has failed")
