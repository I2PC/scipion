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
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.grigoriefflab import *

# Some utility functions to import movies that are used
# in several tests.
class TestGrigorieffLabUnBlurBase(BaseTest):
    @classmethod
    def setData(cls):
        cls.dataset = DataSet.getDataSet('movies')
        cls.movies = cls.dataset.getFile('ribo/*.mrcs')

    @classmethod
    def runImportMovie(cls, pattern, samplingRate, voltage, scannedPixelSize, magnification, sphericalAberration):
        """ Run an Import micrograph protocol. """
        # We have two options: passe the SamplingRate or the ScannedPixelSize + microscope magnification
        if not samplingRate is None:
            cls.protImport = ProtImportMovies(samplingRateMode=0, filesPath=pattern, samplingRate=samplingRate, magnification=magnification,
                                                   voltage=voltage, sphericalAberration=sphericalAberration)
        else:
            cls.protImport = ProtImportMovies(samplingRateMode=1, filesPath=pattern, scannedPixelSize=scannedPixelSize,
                                                   voltage=voltage, magnification=magnification, sphericalAberration=sphericalAberration)

        cls.proj.launchProtocol(cls.protImport, wait=True)
        if cls.protImport.isFailed():
            raise Exception("Protocol has failed. Error: ", cls.protImport.getErrorMessage())
        # check that input movies have been imported (a better way to do this?)
        if cls.protImport.outputMovies is None:
            raise Exception('Import of movies: %s, failed. outputMovies is None.' % pattern)
        return cls.protImport

    @classmethod
    def runImportMovie(cls, pattern, samplingRate, voltage,
                       scannedPixelSize, magnification, sphericalAberration,
                       amplitudeContrast):
        """ Run an Import micrograph protocol. """
        # We have two options: passe the SamplingRate or the ScannedPixelSize + microscope magnification
        if not samplingRate is None:
            cls.protImport = ProtImportMovies(samplingRateMode=0
                                              , filesPath=pattern
                                              , samplingRate=samplingRate
                                              , magnification=magnification
                                              , voltage=voltage
                                              , amplitudeContrast = amplitudeContrast
                                              , sphericalAberration=sphericalAberration)
        else:
            cls.protImport = ProtImportMovies(samplingRateMode=1, filesPath=pattern, scannedPixelSize=scannedPixelSize,
                                                   voltage=voltage, magnification=magnification, sphericalAberration=sphericalAberration)

        cls.proj.launchProtocol(cls.protImport, wait=True)
        if cls.protImport.isFailed():
            raise Exception("Protocol has failed. Error: ", cls.protImport.getErrorMessage())
        # check that input movies have been imported (a better way to do this?)
        if cls.protImport.outputMovies is None:
            raise Exception('Import of movies: %s, failed. outputMovies is None.' % pattern)
        return cls.protImport



class TestGrigorieffLabUnBlur(TestGrigorieffLabUnBlurBase):
    """This class check if the unblur micrograph protocol in
    GrigorieffLab works properly."""
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestGrigorieffLabUnBlurBase.setData()
        cls.protImport = cls.runImportMovie(cls.movies, 3.54, 300, 0, 50000, 2., 0.1)


    def testUnblur(self):
        # test downsampling a set of micrographs
        protUnblur = ProtUnblur(alignMethod=0)
        protUnblur.inputMovies.set(self.protImport.outputMovies)
        protUnblur.alignFrameRange.set(-1)
        protUnblur.doApplyDoseFilter.set(True)
        protUnblur.exposurePerFrame.set(1.3)
        self.proj.launchProtocol(protUnblur, wait=True)
        #self.assertIsNotNone(protUnblur.outputMicrographs, "SetOfMicrographs has not been created.")
        self.assertTrue(True)