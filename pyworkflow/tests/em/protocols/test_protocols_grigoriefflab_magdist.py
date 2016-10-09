# **************************************************************************
# *
# * Authors:    Grigory Sharov (sharov@igbmc.fr)
# *
# * L'Institut de genetique et de biologie moleculaire et cellulaire (IGBMC)
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

from pyworkflow.tests import *
from pyworkflow.em.protocol import (ProtImportMovies,
                                    ProtImportMicrographs,
                                    ProtImportCoordinates)
from pyworkflow.em.packages.grigoriefflab import (ProtMagDistEst,
                                                  ProtMagDistCorr,
                                                  ProtMagDistCorrCoord)


# Some utility functions to import files that are used in several tests.
class TestMagDistBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('movies')
        cls.movie1 = cls.dataset.getFile('qbeta/qbeta.mrc')
        cls.movies = cls.dataset.getFile('ribo/*.mrcs')
        cls.dataset2 = DataSet.getDataSet('grigorieff')
        cls.mics = cls.dataset2.getFile('gold/*.mrc')
        cls.coords = cls.dataset2.getFile('gold/')
    
    @classmethod
    def runImportMovie(cls, pattern, objLabel, samplingRate, voltage,
                       scannedPixelSize, magnification, sphericalAberration):
        """ Run an Import movies protocol. """
        # We have two options: pass the SamplingRate or
        # the ScannedPixelSize + microscope magnification
        kwargs = {
            'objLabel': objLabel,
            'filesPath': pattern,
            'magnification': magnification,
            'voltage': voltage,
            'sphericalAberration': sphericalAberration,
            'amplitudeContrast': 0.1
        }

        if samplingRate is not None:
            kwargs.update({'samplingRateMode': 0,
                           'samplingRate': samplingRate})
        else:
            kwargs.update({'samplingRateMode': 1,
                           'scannedPixelSize': scannedPixelSize})

        cls.protImport = cls.newProtocol(ProtImportMovies, **kwargs)
        cls.proj.launchProtocol(cls.protImport, wait=True)

        if cls.protImport.isFailed():
            raise Exception("Protocol has failed. Error: ",
                            cls.protImport.getErrorMessage())

        # Check that input movies have been imported (a better way to do this?)
        if cls.protImport.outputMovies is None:
            raise Exception('Import of movies: %s, failed, '
                            'outputMovies is None.' % pattern)

        return cls.protImport
    
    @classmethod
    def runImportMovie1(cls):
        """ Run an Import movie protocol. """
        return cls.runImportMovie(cls.movie1, 'movie - qbeta.mrc',
                                  samplingRate=1.14, voltage=300,
                                  sphericalAberration=2.26,
                                  scannedPixelSize=None, magnification=50000)
    
    @classmethod
    def runImportMovies(cls):
        """ Run an Import movie protocol. """
        return cls.runImportMovie(cls.movies, 'movies - betagal.mrcs',
                                  samplingRate=3.54, voltage=300,
                                  sphericalAberration=2.,
                                  scannedPixelSize=None, magnification=50000)

    @classmethod
    def runImportMicrograph(cls, pattern, samplingRate, voltage,
                            scannedPixelSize, magnification,
                            sphericalAberration):
        """ Run an Import micrograph protocol. """
        # We have two options: pass the SamplingRate or
        # the ScannedPixelSize + microscope magnification
        kwargs = {
            'filesPath': pattern,
            'magnification': magnification,
            'voltage': voltage,
            'sphericalAberration': sphericalAberration
        }

        if samplingRate is not None:
            kwargs.update({'samplingRateMode': 0,
                           'samplingRate': samplingRate})
        else:
            kwargs.update({'samplingRateMode': 1,
                           'scannedPixelSize': scannedPixelSize})

        cls.protImport2 = ProtImportMicrographs(objLabel='mics - cross grating', **kwargs)
        cls.launchProtocol(cls.protImport2)

        # Check that input micrographs have been imported
        if cls.protImport2.outputMicrographs is None:
            raise Exception('Import of micrograph: %s, failed. '
                            'outputMicrographs is None.' % pattern)

        return cls.protImport2

    @classmethod
    def runImportMics(cls):
        """ Run an Import micrograph protocol. """
        return cls.runImportMicrograph(cls.mics,
                                       samplingRate=1.1,
                                       voltage=300,
                                       sphericalAberration=0.01,
                                       scannedPixelSize=None,
                                       magnification=59000)

    def runBasicTests(self):
        # Expected dimensions of movies and mics
        dims = [(4096, 4096, 7), (1950, 1950, 16), (4096, 4096, 1)]

        # Import 2 different sets of movies and a set of mics
        protImport1 = self.runImportMovie1()
        protImport2 = self.runImportMovies()
        protImport3 = self.runImportMics()

        inputMovies1 = getattr(protImport1, 'outputMovies')
        inputMovies2 = getattr(protImport2, 'outputMovies')
        inputMics = getattr(protImport3, 'outputMicrographs')

        # Check import
        inputSets = [inputMovies1, inputMovies2, inputMics]
        for index, item in enumerate(inputSets):
            self.assertIsNotNone(item)
            self.assertEqual(dims[index], item.getDim())

        # Also import coords
        protImport4 = self.newProtocol(ProtImportCoordinates,
                                       importFrom=ProtImportCoordinates.IMPORT_FROM_XMIPP,
                                       filesPath=self.coords,
                                       filesPattern='*.pos',
                                       boxSize=150)
        protImport4.inputMicrographs.set(inputMics)
        protImport4.setObjLabel('import coords from xmipp')
        self.launchProtocol(protImport4)
        inputCoords = getattr(protImport4, 'outputCoordinates')

        # Run estimation on mics
        protEst = self.newProtocol(ProtMagDistEst,
                                   minAng=0.0,
                                   maxAng=90.0,
                                   angStep=1.0)
        protEst.inputMicrographs.set(inputMics)
        self.launchProtocol(protEst)

        # verify corrected pixSize
        self.assertAlmostEquals(protEst._parseOutputLog()[3], 1.09, delta=0.1)

        # Run correction with input params on movies
        protCorr1 = self.newProtocol(ProtMagDistCorr,
                                     objLabel=ProtMagDistCorr._label + str(1),
                                     scaleMaj=1.02,
                                     scaleMin=0.993,
                                     angDist=30,
                                     newPix=1.13)
        protCorr1.inputMovies.set(inputMovies1)
        self.launchProtocol(protCorr1)

        # Check output movies1
        outputMovies1 = getattr(protCorr1, 'outputMovies', None)
        self.assertIsNotNone(outputMovies1)
        self.assertEqual(protImport1.outputMovies.getSize(),
                         outputMovies1.getSize())

        for movie in outputMovies1:
            movieFn = movie.getFileName()
            self.assertTrue(os.path.exists(self.proj.getPath(movieFn)))

        # Run correction with input = estimation protocol
        protCorr2 = self.newProtocol(ProtMagDistCorr,
                                     objLabel=ProtMagDistCorr._label + str(2),
                                     useEst=True)
        protCorr2.inputMovies.set(inputMovies2)
        protCorr2.inputEst.set(protEst)
        self.launchProtocol(protCorr2)

        # Check output movies2
        outputMovies2 = getattr(protCorr2, 'outputMovies', None)
        self.assertIsNotNone(outputMovies2)
        self.assertEqual(protImport2.outputMovies.getSize(),
                         outputMovies2.getSize())

        for movie in outputMovies2:
            movieFn = movie.getFileName()
            self.assertTrue(os.path.exists(self.proj.getPath(movieFn)))

        # Run correction on input coords
        protCorr3 = self.newProtocol(ProtMagDistCorrCoord,
                                     objLabel=ProtMagDistCorr._label + str(3),
                                     useEst=True)
        protCorr3.inputCoords.set(inputCoords)
        protCorr3.inputEst.set(protEst)
        self.launchProtocol(protCorr3)

        # Check output coords
        outputCoords = getattr(protCorr3, 'outputCoordinates', None)
        self.assertIsNotNone(outputCoords)
        self.assertEqual(protImport4.outputCoordinates.getSize(),
                         outputCoords.getSize())


class TestMagDist(TestMagDistBase):
    def test_movies(self):
        self.runBasicTests()
