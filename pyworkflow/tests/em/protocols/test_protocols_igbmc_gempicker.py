# **************************************************************************
# *
# * Authors:    Grigory Sharov (sharov@igbmc.fr)
# *             J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
from os.path import join, basename

from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.igbmc.protocol_gempicker import *

import pyworkflow.utils as pwutils



class TestIgbmcBase(BaseTest):
    @classmethod
    def setData(cls):
        cls.ds = DataSet.getDataSet('igbmc_gempicker')
    
    @classmethod
    def runImportAverages(cls):
        """ Run an Import averages protocol. """
        cls.protImportAvg = cls.newProtocol(ProtImportAverages, 
                                            objLabel='import averages (klh)',
                                            filesPath=cls.ds.getFile('templates/*.mrc'), 
                                            samplingRate=4.4)

        cls.launchProtocol(cls.protImportAvg)
        return cls.protImportAvg
 
    @classmethod
    def runImportMaskTop(cls):
        """ Run an Import mask protocol. """
        cls.protImportMskTop = cls.newProtocol(ProtImportMask, 
                                               objLabel='import mask top view (klh)',
                                               maskPath=cls.ds.getFile('masks/mask_topview.tif'), 
                                               samplingRate=4.4)

        cls.launchProtocol(cls.protImportMskTop)
        return cls.protImportMskTop

    @classmethod
    def runImportMaskSide(cls):
        """ Run an Import mask protocol. """
        cls.protImportMskSide = cls.newProtocol(ProtImportMask, 
                                               objLabel='import mask side view (klh)',
                                               maskPath=cls.ds.getFile('masks/mask_sideview.tif'), 
                                               samplingRate=4.4)

        cls.launchProtocol(cls.protImportMskSide)
        return cls.protImportMskSide

    @classmethod
    def runImportMicrograph(cls, pattern, samplingRate, voltage, magnification, sphericalAberration):
        """ Run an Import micrograph protocol. """
        cls.protImport = cls.newProtocol(ProtImportMicrographs,
                                         objLabel='import mics (klh)', 
                                         samplingRateMode=0, 
                                         filesPath=pattern, 
                                         samplingRate=samplingRate, 
                                         magnification=magnification, 
                                         voltage=voltage, 
                                         sphericalAberration=sphericalAberration)
            
        cls.launchProtocol(cls.protImport)
        return cls.protImport
    
    @classmethod
    def runImportMicrographKLH(cls):
        """ Run an Import micrograph protocol. """
        return cls.runImportMicrograph(cls.ds.getFile('micrographs/*.mrc'), 
                                       samplingRate=4.4, 
                                       voltage=120, sphericalAberration=2,
                                       magnification=66000)
    
    @classmethod
    def runPickingMask1(cls):
        """ Run a particle picking. """
        protGP = ProtGemPicker(objLabel='Gempicker with a single circular mask (klh)',
                               refsHaveInvertedContrast=True,
                               thresholdLow=0.2,
                               thresholdHigh=0.9,
                               boxSize=100,
                               boxDist=10,
                               maskType=MASK_CIRCULAR,
                               useGPU=False)                
        protGP.inputMicrographs.set(cls.protImportMics.outputMicrographs)
        protGP.inputReferences.set(cls.protImportAvgs.outputAverages)
        cls.launchProtocol(protGP)
        return protGP

    @classmethod
    def runPickingMask2(cls):
        """ Run a particle picking. """
        masks = PointerList()
        masks.append(cls.protImportMskSide.outputMask)
        masks.append(cls.protImportMskTop.outputMask)
        protGP = ProtGemPicker(objLabel='Gempicker with a individual masks (klh)',
                               refsHaveInvertedContrast=True,
                               thresholdLow=0.2,
                               thresholdHigh=0.9,
                               boxSize=100,
                               boxDist=10,
                               maskType=MASK_OBJECT,
                               inputMasks=masks,
                               useGPU=False)                
        protGP.inputMicrographs.set(cls.protImportMics.outputMicrographs)
        protGP.inputReferences.set(cls.protImportAvgs.outputAverages)
        cls.launchProtocol(protGP)
        return protGP

class TestGempickerAutomaticPicking(TestIgbmcBase):
    """This class check if the protocol to pick the micrographs automatically by gempicker works properly."""
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestIgbmcBase.setData()
        cls.protImportMics = cls.runImportMicrographKLH()
        cls.protImportMaskSide = cls.runImportMaskSide()
        cls.protImportMaskTop = cls.runImportMaskTop()
        cls.protImportAvgs = cls.runImportAverages()
    
    def testAutomaticPicking(self):
        print "Run automatic particle picking"
        self.runPickingMask1()
        self.runPickingMask2()

