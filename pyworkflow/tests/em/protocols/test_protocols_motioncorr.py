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
import pyworkflow.em as em
from pyworkflow.tests import *
from pyworkflow.em.packages.motioncorr import ProtMotionCorr


# Some utility functions to import movies that are used
# in several tests.
class TestMotioncorrBase(BaseTest):
    @classmethod
    def setData(cls):
        cls.dataset = DataSet.getDataSet('movies')
        cls.movie1 = cls.dataset.getFile('qbeta/qbeta.mrc')
        cls.movie2 = cls.dataset.getFile('cct/cct_1.em')
    
    @classmethod
    def runImportMovie(cls, pattern, samplingRate, voltage, scannedPixelSize, magnification, sphericalAberration):
        """ Run an Import micrograph protocol. """
        # We have two options: passe the SamplingRate or the ScannedPixelSize + microscope magnification
        if not samplingRate is None:
            cls.protImport = em.ProtImportMovies(samplingRateMode=0, filesPath=pattern, samplingRate=samplingRate, magnification=magnification, 
                                                   voltage=voltage, sphericalAberration=sphericalAberration)
        else:
            cls.protImport = em.ProtImportMovies(samplingRateMode=1, filesPath=pattern, scannedPixelSize=scannedPixelSize, 
                                                   voltage=voltage, magnification=magnification, sphericalAberration=sphericalAberration)
        
        cls.proj.launchProtocol(cls.protImport, wait=True)
        if cls.protImport.isFailed():
            raise Exception("Protocol has failed. Error: ", cls.protImport.getErrorMessage())
        # check that input movies have been imported (a better way to do this?)
        if cls.protImport.outputMovies is None:
            raise Exception('Import of movies: %s, failed. outputMovies is None.' % pattern)
        return cls.protImport
    
    @classmethod
    def runImportMovie1(cls, pattern):
        """ Run an Import movie protocol. """
        return cls.runImportMovie(pattern, samplingRate=1.14, voltage=300, sphericalAberration=2.26, scannedPixelSize=None, magnification=50000)
    
    @classmethod
    def runImportMovie2(cls, pattern):
        """ Run an Import movie protocol. """
        return cls.runImportMovie(pattern, samplingRate=1.4, voltage=300, sphericalAberration=2.7, scannedPixelSize=None, magnification=61000)


class TestXmippAlingMovies(TestMotioncorrBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestMotioncorrBase.setData()
        cls.protImport1 = cls.runImportMovie1(cls.movie1)
        cls.protImport2 = cls.runImportMovie2(cls.movie2)
    
    def testMotioncorr1(self):
        protMotCorr1 = ProtMotionCorr()
        protMotCorr1.inputMovies.set(self.protImport1.outputMovies)
        self.proj.launchProtocol(protMotCorr1, wait=True)
        self.assertIsNotNone(protMotCorr1.outputMicrographs, "SetOfMicrographs has not been created.")
    
    def testMotioncorr2(self):
        # test downsampling a set of micrographs
        protMotioncorr2 = ProtMotionCorr(doSaveMovie=True)
        protMotioncorr2.inputMovies.set(self.protImport2.outputMovies)
        self.proj.launchProtocol(protMotioncorr2, wait=True)
        self.assertIsNotNone(protMotioncorr2.outputMicrographs, "SetOfMicrographs has not been created.")

    def testMotioncorr3(self):
        # test downsampling a set of micrographs
        protMotioncorr3 = ProtMotionCorr(alignFrame0=2, alignFrameN=2, 
                                         cropOffsetX=10, cropOffsetY=10)
        protMotioncorr3.inputMovies.set(self.protImport1.outputMovies)
        self.proj.launchProtocol(protMotioncorr3, wait=True)
        self.assertIsNotNone(protMotioncorr3.outputMicrographs, "SetOfMicrographs has not been created.")
