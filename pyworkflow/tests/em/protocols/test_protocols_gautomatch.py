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
from pyworkflow.em.packages.gautomatch.protocol_gautomatch import *

import pyworkflow.utils as pwutils


class TestGautomatchBase(BaseTest):
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
    def runImportCoords(cls):
        """ Run an Import coords protocol. """
        cls.protImportCoords = cls.newProtocol(ProtImportCoordinates,
                                               importFrom=ProtImportCoordinates.IMPORT_FROM_XMIPP,
                                               objLabel='import bad coords',
                                               filesPath=cls.ds.getFile('coords/'),
                                               filesPattern='*.pos',
                                               boxSize=100)
        cls.protImportCoords.inputMicrographs.set(cls.protImportMics.outputMicrographs)
        cls.launchProtocol(cls.protImportCoords)
        return cls.protImportCoords

    @classmethod
    def runPicking1(cls):
        """ Run a particle picking. """
        protGM = ProtGautomatch(objLabel='Gautomatch auto-picking (klh)',
                               invertTemplatesContrast=True,
                               threshold=0.18,
                               particleSize=250,
                               advanced='False',
                               boxSize=150,
                               localSigmaCutoff=2.0)
        protGM.inputMicrographs.set(cls.protImportMics.outputMicrographs)
        protGM.inputReferences.set(cls.protImportAvgs.outputAverages)
        cls.launchProtocol(protGM)
        return protGM

    @classmethod
    def runPicking2(cls):
        """ Run a particle picking with excludsive options. """
        protGM2 = ProtGautomatch(objLabel='Gautomatch auto-picking 2 (klh)',
                                invertTemplatesContrast=True,
                                threshold=0.18,
                                particleSize=250,
                                advanced='False',
                                boxSize=150,
                                localSigmaCutoff=2.0,
                                exclusive=True)
        protGM2.inputMicrographs.set(cls.protImportMics.outputMicrographs)
        protGM2.inputReferences.set(cls.protImportAvgs.outputAverages)
        protGM2.inputBadCoords.set(cls.protImportCoords.outputCoordinates)
        cls.launchProtocol(protGM2)
        return protGM2

class TestGautomatchAutomaticPicking(TestGautomatchBase):
    """This class check if the protocol to pick the micrographs automatically by gautomatch works properly."""
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestGautomatchBase.setData()
        cls.protImportMics = cls.runImportMicrographKLH()
        cls.protImportAvgs = cls.runImportAverages()
        cls.protImportCoords = cls.runImportCoords()
    
    def testAutomaticPicking(self):
        self.runPicking1()
        self.runPicking2()
