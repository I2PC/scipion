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

    def runMovieResize(self, fnMovies, resizeVal):
        kwargs = {'inputMovies': fnMovies,
                  'resizeFactor': resizeVal
                  }
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

        resizeFactor=3.5
        protMovieResize = self.runMovieResize(protImportMovsStr.outputMovies,
                                              resizeFactor)
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
        xOrig, yOrig, zOrig = protImportMovsStr.outputMovies.getDim()
        if xResize != int(xOrig/resizeFactor) or \
                yResize != int(yOrig/resizeFactor) or zResize != zOrig:
            self.assertTrue(False)


