# **************************************************************************
# *
# * Authors:    Laura del Cano (ldelcano@cnb.csic.es)
# *             Josue Gomez Blanco (jgomez@cnb.csic.es)
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
from os.path import join, basename

from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3 import *
from pyworkflow.em.packages.igbmc import *
import pyworkflow.utils as pwutils

# Some utility functions to import micrographs that are used
# in several tests.


class TestIgbmcBase(BaseTest):
    @classmethod
    def setData(cls):
        cls.dataset = DataSet.getDataSet('igbmc_gempicker')
        cls.micFn = cls.dataset.getFile('mic1')
        cls.micsFn = cls.dataset.getFile('allMics')
        cls.maskDir = cls.dataset.getFile('mask')
        cls.templateDir = cls.dataset.getFile('templates')
    
    @classmethod
    def runImportAverages(cls, templateDir):
        """ Run an Import averages protocol. """
        cls.protImportAvg = cls.newProtocol(ProtImportAverages, 
                                      objLabel='import averages (klh)',
                                      filesPath=templateDir, 
                                      samplingRate=4.4)

        cls.launchProtocol(cls.protImportAvg)
        if cls.protImportAvg.isFailed():
            raise Exception("Protocol has failed. Error: ", cls.protImportAvg.getErrorMessage())
        # check that input averages have been imported (a better way to do this?)
        if cls.protImportAvg.outputAverages is None:
            raise Exception('Import of averages: %s, failed. outputAverages is None.' % templateDir)
        return cls.protImportAvg
 
    @classmethod
    def runImportMask(cls, maskDir):
        """ Run an Import mask protocol. """
        cls.protImportMsk = cls.newProtocol(ProtImportMask, 
                                      objLabel='import mask (klh)',
                                      maskPath=maskDir, 
                                      samplingRate=4.4)

        cls.launchProtocol(cls.protImportMsk)
        if cls.protImportMsk.isFailed():
            raise Exception("Protocol has failed. Error: ", cls.protImportMsk.getErrorMessage())
        # check that input masks have been imported (a better way to do this?)
        if cls.protImportMsk.outputMask is None:
            raise Exception('Import of mask: %s, failed. outputMask is None.' % maskDir)
        return cls.protImportMsk

    @classmethod
    def runImportMicrograph(cls, pattern, samplingRate, voltage, scannedPixelSize, magnification, sphericalAberration):
        """ Run an Import micrograph protocol. """
        
        # We have two options: pass the SamplingRate or the ScannedPixelSize + microscope magnification
        if not samplingRate is None:
            cls.protImport = cls.newProtocol(ProtImportMicrographs, 
                                             samplingRateMode=0, 
                                             filesPath=pattern, 
                                             samplingRate=samplingRate, 
                                             magnification=magnification, 
                                             voltage=voltage, 
                                             sphericalAberration=sphericalAberration)
        else:
            cls.protImport = cls.newProtocol(ProtImportMicrographs, 
                                             samplingRateMode=1, 
                                             filesPath=pattern, 
                                             scannedPixelSize=scannedPixelSize, 
                                             voltage=voltage, 
                                             magnification=magnification, 
                                             sphericalAberration=sphericalAberration)
            
        cls.protImport.setObjLabel('import mics (klh)')
        cls.launchProtocol(cls.protImport)
        if cls.protImport.isFailed():
            raise Exception("Protocol has failed. Error: ", cls.protImport.getErrorMessage())
        # check that input micrographs have been imported (a better way to do this?)
        if cls.protImport.outputMicrographs is None:
            raise Exception('Import of micrograph: %s, failed. outputMicrographs is None.' % pattern)
        return cls.protImport
    
    @classmethod
    def runImportMicrographKLH(cls, pattern):
        """ Run an Import micrograph protocol. """
        pattern = cls.micsFn
        return cls.runImportMicrograph(pattern, samplingRate=2.2, 
                                       voltage=120, sphericalAberration=2, 
                                       scannedPixelSize=None, magnification=66000)
    
    @classmethod
    def runDownsamplingMicrographs(cls, mics, downFactorValue, threads=1):
        # test downsampling a set of micrographs
        mics = cls.protImport.outputMicrographs
        cls.protDown = XmippProtPreprocessMicrographs(doDownsample=True, 
                                                      downFactor=downFactorValue, 
                                                      numberOfThreads=threads)
        cls.protDown.inputMicrographs.set(mics)
        cls.proj.launchProtocol(cls.protDown, wait=True)
        return cls.protDown
    
    @classmethod
    def runPicking(cls, micsdown, mask, refs):
        """ Run a particle picking. """
        micsdown = cls.protDown.outputMicrographs
        mask = cls.protImportMsk.outputMask
        refs = cls.protImportAvg.outputAverages
        cls.protPP = ProtGemPicker(inputReferences=refs,
                                   refsHaveInvertedContrast=True,
                                   maskType=1,
                                   inputMasks=mask)                
        cls.protPP.inputMicrographs.set(micsdown)               
        cls.proj.launchProtocol(cls.protPP, wait=True)
        # check that picking has run ok
        if cls.protPP.outputCoordinates is None:
            raise Exception('Particle picking failed. outputCoordinates is None.')
        return cls.protPP


class TestGempickerAutomaticPicking(TestIgbmcBase):
    """This class check if the protocol to pick the micrographs automatically by gempicker works properly."""
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestIgbmcBase.setData()
        cls.protImport1 = cls.runImportMicrographKLH(cls.micsFn)
        cls.protImport2 = cls.runImportMask(cls.maskDir)
        cls.protImport3 = cls.runImportAverages(cls.templateDir)
#        cls.protImport2 = cls.runImportMicrographKLH(cls.micFn)
        cls.protDown1 = cls.runDownsamplingMicrographs(cls.protImport1.outputMicrographs, 2)
#        cls.protDown2 = cls.runDownsamplingMicrographs(cls.protImport2.outputMicrographs, 2)
        cls.protPP = cls.runPicking(cls.protDown1.outputMicrographs, cls.protImport2.outputMask, cls.protImport3.outputAverages)
    
    def testAutomaticPicking(cls):
        print "Run automatic particle picking"
        protAutomaticPP = ProtGemPicker()
        protAutomaticPP.set(cls.protPP)
        cls.proj.launchProtocol(protAutomaticPP, wait=True)
        if cls.protAutomaticPP.outputCoordinates is None: 
            raise Exception('There was a problem with the automatic particle picking.')
    
#    def testAutomaticPickingOther(self):
#        print "Run automatic particle picking"
#        protAutomaticPP = ProtGemPicker()
#        protAutomaticPP.set(self.protPP)
#        protAutomaticPP.inputMicrographs.set(self.protDown2.outputMicrographs)
#        protAutomaticPP.micsToPick.set(1)
#        self.proj.launchProtocol(protAutomaticPP, wait=True)
#        self.assertIsNotNone(protAutomaticPP.outputCoordinates, 
#                             "There was a problem with the automatic particle picking")

