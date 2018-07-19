# **************************************************************************
# *
# * Authors:    Erney ramirez-Aportela (jlvilas@cnb.csic.es)
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
#from pyworkflow.em.packages.xmipp3.protocol_resolution_monogenic_signal import OUTPUT_RESOLUTION_FILE
from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pyworkflow.em.packages.xmipp3 import XmippProtMonoRes, XmippProtCreateMask3D, XmippProtLocSharp
from pyworkflow.em.protocol import ProtImportVolumes


class TestLocalDeblurBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='resmap'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.map3D = cls.dataset.getFile('betagal')
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
                                  doSmall=False,
                                  smallSize=False,
                                  doBig=False,
                                  doSymmetrize=False,
                                  doMorphological=False,
                                  doInvert=False,
                                  doSmooth=False,
                                  sigmaConvolution=2
                                  )
        cls.launchProtocol(cls.msk)
        return cls.msk


class TestLocalDeblur(TestLocalDeblurBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestLocalDeblurBase.setData()
        cls.protImportVol = cls.runImportVolumes(cls.map3D, 3.54)
        cls.protCreateMask = cls.runCreateMask(cls.protImportVol.outputVolume, 0.02)

    def testLocalDeblur(self):
        MonoRes = self.newProtocol(XmippProtMonoRes,
                                   objLabel='single volume monores',
                                   halfVolumes=False,
                                   inputVolumes=self.protImportVol.outputVolume,
                                   Mask=self.protCreateMask.outputMask,
                                   symmetry='c1',
                                   minRes=1,
                                   maxRes=25,
                                   isPremasked = False,
                                   filterInput=False,
                                   )
        self.launchProtocol(MonoRes)
        self.assertTrue(exists(MonoRes._getExtraPath('mgresolution.vol')),
                        "MonoRes (no split, no premasked) has failed")

        LocalDeblur = self.newProtocol(XmippProtLocSharp,
                                   objLabel='sharpening localdeblur',    
                                   inputVolume=self.protImportVol.outputVolume,
                                   resolutionVolume=MonoRes.resolution_Volume,
                                   const=1,
                                   )
        self.launchProtocol(LocalDeblur)
        self.assertTrue(exists(LocalDeblur._getExtraPath('sharpenedMap.mrc')),
                        "LocalDeblur  has failed")
 
 
