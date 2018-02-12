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
from pyworkflow.em.packages.xmipp3.protocol_preprocess import \
    XmippProtMovieResize
from pyworkflow.em.protocol import ProtCreateStreamData
from pyworkflow.em.protocol.protocol_create_stream_data import \
    SET_OF_MOVIES
from pyworkflow.protocol.constants import STATUS_FINISHED

RESIZE_SAMPLINGRATE = 0
RESIZE_DIMENSIONS = 1
RESIZE_FACTOR = 2

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



    def importMoviesStr(self, fnMovies):
        kwargs = {'inputMovies': fnMovies,
                  'nDim': NUM_MOVIES,
                  'creationInterval': 30,
                  'delay': 10,
                  'setof': SET_OF_MOVIES  # SetOfMovies
                  }
        protStream = self.newProtocol(ProtCreateStreamData, **kwargs)
        protStream.setObjLabel('create Stream Movies')
        self.proj.launchProtocol(protStream, wait=False)

        return protStream

    def runMovieResize(self, fnMovies, option, paramVal):
        kwargs = {'inputMovies': fnMovies,
                  'resizeOption': option
                  }
        if option==RESIZE_SAMPLINGRATE:
            kwargs['resizeSamplingRate']=paramVal
        if option==RESIZE_DIMENSIONS:
            kwargs['resizeDim']=paramVal
        if option==RESIZE_FACTOR:
            kwargs['resizeFactor']=paramVal

        protResize = self.newProtocol(XmippProtMovieResize, **kwargs)
        self.proj.launchProtocol(protResize, wait=False)

        return protResize



    def test_pattern(self):

        protImport = self.runImportMoviesRibo()
        counter = 1
        while not protImport.hasAttribute('outputMovies'):
            time.sleep(2)
            protImport = self._updateProtocol(protImport)
            if counter > 100:
                break
            counter += 1

        protImportMovsStr = self.importMoviesStr(protImport.outputMovies)
        counter = 1
        while not protImportMovsStr.hasAttribute('outputMovies'):
            time.sleep(2)
            protImportMovsStr = self._updateProtocol(protImportMovsStr)
            if counter > 100:
                break
            counter += 1
        if protImportMovsStr.isFailed():
            self.assertTrue(False)

        newSamplingRate=4.0
        protMovieResize = self.runMovieResize(protImportMovsStr.outputMovies,
                                              RESIZE_SAMPLINGRATE,
                                              newSamplingRate)

        counter = 1
        while not protMovieResize.hasAttribute('outputMovies'):
            time.sleep(2)
            protMovieResize = self._updateProtocol(protMovieResize)
            if counter > 100:
                break
            counter += 1

        newDimensions = 1700
        protMovieResize2 = self.runMovieResize(protMovieResize.outputMovies,
                                              RESIZE_DIMENSIONS,
                                              newDimensions)

        counter = 1
        while not protMovieResize2.hasAttribute('outputMovies'):
            time.sleep(2)
            protMovieResize2 = self._updateProtocol(protMovieResize2)
            if counter > 100:
                break
            counter += 1


        newFactor = 2.15
        protMovieResize3 = self.runMovieResize(protMovieResize2.outputMovies,
                                              RESIZE_FACTOR, newFactor)

        counter = 1
        while not protMovieResize3.hasAttribute('outputMovies'):
            time.sleep(2)
            protMovieResize3 = self._updateProtocol(protMovieResize3)
            if counter > 100:
                break
            counter += 1

        if not protMovieResize3.hasAttribute('outputMovies'):
            self.assertTrue(False)

        while protMovieResize3.getStatus() != STATUS_FINISHED:
            protMovieResize3 = self._updateProtocol(protMovieResize3)
            if protMovieResize3.isFailed():
                self.assertTrue(False)

        protImportMovsStr = self._updateProtocol(protImportMovsStr)
        protMovieResize = self._updateProtocol(protMovieResize)
        protMovieResize2 = self._updateProtocol(protMovieResize2)
        protMovieResize3 = self._updateProtocol(protMovieResize3)
        if protMovieResize3.outputMovies.getSize() != \
                protMovieResize3.outputMovies.getSize():
            self.assertTrue(False)

        xOrig, yOrig, zOrig = protImportMovsStr.outputMovies.getDim()
        samplingRateOrig = protImportMovsStr.outputMovies.getSamplingRate()
        x1, y1, z1 = protMovieResize.outputMovies.getDim()
        factor1 = newSamplingRate / samplingRateOrig
        if x1 != int(xOrig/factor1) or y1 != int(yOrig/factor1) \
                or z1 != zOrig:
            self.assertTrue(False)

        x2, y2, z2 = protMovieResize2.outputMovies.getDim()
        factor2 = float(x1) / float(newDimensions)
        if x2 != int(x1 / factor2) or y2 != int(y1 / factor2) \
                or z2 != z1:
            self.assertTrue(False)

        x3, y3, z3 = protMovieResize3.outputMovies.getDim()
        if x3 != int(x2 / newFactor) or y3 != int(y2 / newFactor) \
                or z3 != z2:
            self.assertTrue(False)

