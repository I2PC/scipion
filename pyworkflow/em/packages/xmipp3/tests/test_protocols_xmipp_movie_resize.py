# **************************************************************************
# *
# * Authors:    Amaya Jimenez (ajimenez@cnb.csic.es)
# *
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


from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3 import ProtImportMovies
from pyworkflow.protocol import getProtocolFromDb
from pyworkflow.em.protocol.protocol_sets import ProtUnionSet
from pyworkflow.em.packages.xmipp3.protocol_preprocess import \
    XmippProtMovieResize
from pyworkflow.em.protocol import ProtCreateStreamData
from pyworkflow.em.protocol.protocol_create_stream_data import \
    SET_OF_MOVIES
from pyworkflow.protocol.constants import STATUS_FINISHED


NUM_MOVIES = 5

# Some utility functions to import movies that are used in several tests.
class TestMovieResize(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('movies')
        cls.movie1 = cls.dataset.getFile('qbeta/qbeta.mrc')
        cls.movie2 = cls.dataset.getFile('cct/cct_1.em')

    def _updateProtocol(self, prot):
        prot2 = getProtocolFromDb(prot.getProject().path,
                                  prot.getDbPath(),
                                  prot.getObjId())
        # Close DB connections
        prot2.getProject().closeMapper()
        prot2.closeMappers()
        return prot2
    
    def runImportMovie(cls, pattern, samplingRate, voltage, scannedPixelSize,
                       magnification, sphericalAberration, dosePerFrame=None):
        """ Run an Import micrograph protocol. """

        kwargs = {
                 'filesPath': pattern,
                 'magnification': magnification,
                 'voltage': voltage,
                 'sphericalAberration': sphericalAberration,
                 'dosePerFrame' : dosePerFrame
                  }

        # We have two options: pass the SamplingRate or
        # the ScannedPixelSize + microscope magnification
        if samplingRate is not None:
            kwargs.update({'samplingRateMode': 0,
                           'samplingRate': samplingRate})
        else:
            kwargs.update({'samplingRateMode': 1,
                           'scannedPixelSize': scannedPixelSize})

        cls.protImport = cls.newProtocol(ProtImportMovies, **kwargs)
        cls.proj.launchProtocol(cls.protImport, wait=False)

        if cls.protImport.isFailed():
            raise Exception("Protocol has failed. Error: ",
                            cls.protImport.getErrorMessage())

        return cls.protImport
    

    def runImportMovie1(cls, pattern):
        """ Run an Import movie protocol. """
        return cls.runImportMovie(pattern, samplingRate=1.14, voltage=300,
                                  sphericalAberration=2.26, dosePerFrame=1.5,
                                  scannedPixelSize=None, magnification=50000)

    def runImportMovie2(cls, pattern):
        """ Run an Import movie protocol. """
        return cls.runImportMovie(pattern, samplingRate=1.4, voltage=300,
                                  sphericalAberration=2.7, dosePerFrame=1.5,
                                  scannedPixelSize=None,
                                  magnification=61000)

    def runImportMoviesRibo(cls):
        args = {'filesPath': cls.dataset.getFile('ribo/'),
                'filesPattern': '*.mrcs',
                'amplitudConstrast': 0.1,
                'sphericalAberration': 2.,
                'voltage': 300,
                'samplingRate': 3.54
                }
        cls.protMovieImport = cls.newProtocol(ProtImportMovies, **args)
        cls.proj.launchProtocol(cls.protMovieImport, wait=False)
        return cls.protMovieImport

    def runUnionSet(cls, input1, input2, input3):
        cls.protUnion = cls.newProtocol(ProtUnionSet)
        cls.protUnion.inputSets.append(input1)
        cls.protUnion.inputSets.append(input2)
        cls.protUnion.inputSets.append(input3)
        cls.proj.launchProtocol(cls.protUnion, wait=False)
        return cls.protUnion


    def importMoviesStr(self, fnMovies):
        kwargs = {'inputMovie': fnMovies,
                  'nDim': NUM_MOVIES,
                  'creationInterval': 30,
                  'delay': 10,
                  'setof': SET_OF_MOVIES  # SetOfMovies
                  }
        protStream = self.newProtocol(ProtCreateStreamData, **kwargs)
        protStream.setObjLabel('create Stream Movies')
        self.proj.launchProtocol(protStream, wait=False)

        return protStream

    def runMovieResize(self, fnMovies, newDim):
        kwargs = {'inputMovies': fnMovies,
                  'resizeDim': newDim
                  }
        protResize = self.newProtocol(XmippProtMovieResize, **kwargs)
        self.proj.launchProtocol(protResize, wait=False)

        return protResize



    def test_pattern(self):

        #protImport1 = self.runImportMovie1(self.movie1)
        #counter = 1
        #while not protImport1.hasAttribute('outputMovies'):
        #    time.sleep(2)
        #    protImport1 = self._updateProtocol(protImport1)
        #    if counter > 100:
        #        break
        #    counter += 1

        #protImport2 = self.runImportMovie2(self.movie2)
        #counter = 1
        #while not protImport2.hasAttribute('outputMovies'):
        #    time.sleep(2)
        #    protImport2 = self._updateProtocol(protImport2)
        #    if counter > 100:
        #        break
        #    counter += 1

        protImport3 = self.runImportMoviesRibo()
        counter = 1
        while not protImport3.hasAttribute('outputMovies'):
            time.sleep(2)
            protImport3 = self._updateProtocol(protImport3)
            if counter > 100:
                break
            counter += 1

        #protUnion = self.runUnionSet(protImport1.outputMovies,
        #                             protImport2.outputMovies,
        #                             protImport3.outputMovies)

        protImportMovsStr = self.importMoviesStr(protImport3.outputMovies)
        counter = 1
        while not protImportMovsStr.hasAttribute('outputMovies'):
            time.sleep(2)
            protImportMovsStr = self._updateProtocol(protImportMovsStr)
            if counter > 100:
                break
            counter += 1
        if protImportMovsStr.isFailed():
            self.assertTrue(False)

        xOrig, yOrig, zOrig = protImportMovsStr.outputMovies.getDim()
        newDim = int(xOrig/2)
        protMovieResize = self.runMovieResize(protImportMovsStr.outputMovies,
                                              newDim)
        counter = 1
        while not protMovieResize.hasAttribute('outputMovies'):
            time.sleep(2)
            protMovieResize = self._updateProtocol(protMovieResize)
            if counter > 100:
                break
            counter += 1

        if not protMovieResize.hasAttribute('outputMovies'):
            self.assertTrue(False)

        while protMovieResize.getStatus() != STATUS_FINISHED:
            protMovieResize = self._updateProtocol(protMovieResize)
            if protMovieResize.isFailed():
                self.assertTrue(False)

        protImportMovsStr = self._updateProtocol(protImportMovsStr)
        protMovieResize = self._updateProtocol(protMovieResize)
        if protMovieResize.outputMovies.getSize() !=\
                protImportMovsStr.outputMovies.getSize():
            self.assertTrue(False)

        xResize, yResize, zResize = protMovieResize.outputMovies.getDim()
        if xResize != newDim or yResize != newDim or zResize != zOrig:
            self.assertTrue(False)


