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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.grigoriefflab import *


# Some utility functions to import movies that are used in several tests.
class TestMoviesBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('movies')
        cls.movie1 = cls.dataset.getFile('qbeta/qbeta.mrc')
        cls.movie2 = cls.dataset.getFile('cct/cct_1.em')
        cls.movies = cls.dataset.getFile('ribo/*.mrcs')
    
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
            'amplitudeContrast': 0.1,
            'dosePerFrame': 1.3
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
    def runImportMovie2(cls):
        """ Run an Import movie protocol. """
        return cls.runImportMovie(cls.movie2, 'movie - cct.em',
                                  samplingRate=1.4, voltage=300,
                                  sphericalAberration=2.7,
                                  scannedPixelSize=None, magnification=61000)

    @classmethod
    def runImportMovies(cls):
        """ Run an Import movie protocol. """
        return cls.runImportMovie(cls.movies, 'movies - betagal.mrcs',
                                  samplingRate=3.54, voltage=300,
                                  sphericalAberration=2.,
                                  scannedPixelSize=None, magnification=50000)

    def runBasicTests(self, ProtocolClass):
        # Run some basic tests for both Unblur and Summovie

        # Expected dimensions of imported movies
        dims = [(4096, 4096, 7), (4096, 4096, 7), (1950, 1950, 16)]

        # Expected dimensions of imported movies
        dims2 = [(4096, 4096, 5), (4096, 4096, 5), (1950, 1950, 16)]

        # Launch unblur for 3 different set of movies
        for i, protImport in enumerate([self.runImportMovie1(),
                                        self.runImportMovie2(),
                                        self.runImportMovies()]):
            inputMovies = getattr(protImport, 'outputMovies')
            self.assertIsNotNone(inputMovies)
            self.assertEqual(dims[i], inputMovies.getDim())

            prot = self.newProtocol(ProtocolClass,
                                    objLabel=ProtocolClass._label + str(i),
                                    alignFrame0=2,
                                    alignFrameN=6,
                                    sumFrame0=2,
                                    sumFrameN=6,
                                    doApplyDoseFilter=True)
            prot.inputMovies.set(inputMovies)
            self.launchProtocol(prot)
            outputMics = getattr(prot, 'outputMicrographs', None)
            self.assertIsNotNone(outputMics)
            self.assertEqual(protImport.outputMovies.getSize(),
                             outputMics.getSize())

            for mic in outputMics:
                micFn = mic.getFileName()
                self.assertTrue(os.path.exists(self.proj.getPath(micFn)))


class TestUnblur(TestMoviesBase):
    def test_movies(self):
        self.runBasicTests(ProtUnblur)


class TestSummovie(TestMoviesBase):
    def test_movies(self):
        self.runBasicTests(ProtSummovie)

