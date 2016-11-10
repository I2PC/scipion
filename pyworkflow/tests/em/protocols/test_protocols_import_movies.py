# ***************************************************************************
# * Authors:     Vahid Abrishami (vabrishami@cnb.csic.es)
# *
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
# *  e-mail address 'xmipp@cnb.csic.es'
# ***************************************************************************/

import os
from itertools import izip

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.em.protocol import ProtImportMovies
from pyworkflow.em.data import SetOfMovies


class TestImportBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsMovies = DataSet.getDataSet('movies')


class TestImportMovies(TestImportBase):

    def getArgs(self, filesPath, pattern=''):
        return {'importFrom': ProtImportMovies.IMPORT_FROM_FILES,
                'filesPath': self.dsMovies.getFile(filesPath),
                'filesPattern': pattern,
                'amplitudConstrast': 0.1,
                'sphericalAberration': 2.,
                'voltage': 300,
                'samplingRate': 3.54
                }

    def _checkOutput(self, prot, args, moviesId=[], size=None, dim=None):
        movies = getattr(prot, 'outputMovies', None)
        self.assertIsNotNone(movies)
        self.assertEqual(movies.getSize(), size)

        for i, m in enumerate(movies):
            if moviesId:
                self.assertEqual(m.getObjId(), moviesId[i])
            self.assertAlmostEqual(m.getSamplingRate(),
                                   args['samplingRate'])
            a = m.getAcquisition()
            self.assertAlmostEqual(a.getVoltage(), args['voltage'])

            if dim is not None: # Check if dimensions are the expected ones
                x, y, n = m.getDim()
                self.assertEqual(dim, (x, y, n))

    def test_pattern(self):
        args = self.getArgs('ribo/', pattern='*.mrcs')

        # Id's should be set increasing from 1 if ### is not in the pattern
        protMovieImport = self.newProtocol(ProtImportMovies, **args)
        protMovieImport.setObjLabel('from files')
        self.launchProtocol(protMovieImport)
        self._checkOutput(protMovieImport, args, [1, 2, 3], size=3,
                          dim=(1950, 1950, 16))
        
    def test_em(self):
        args = self.getArgs('cct/cct_1.em')
        args['objLabel'] = 'movie em'

        protMovieImport = self.newProtocol(ProtImportMovies, **args)
        self.launchProtocol(protMovieImport)

        self._checkOutput(protMovieImport, args, size=1, dim=(4096, 4096, 7))

    def test_tif(self):
        args = self.getArgs('c3-adp-se-xyz-0228_200.tif')
        args['objLabel'] = 'movie tif'

        protMovieImport = self.newProtocol(ProtImportMovies, **args)
        self.launchProtocol(protMovieImport)

        self._checkOutput(protMovieImport, args, size=1, dim=(7676, 7420, 38))

    def test_bz2(self):
        args = self.getArgs('Falcon_2012_06_12-14_33_35_0_movie.mrcs.bz2')
        args['objLabel'] = 'movie bz2'

        protMovieImport = self.newProtocol(ProtImportMovies, **args)
        self.launchProtocol(protMovieImport)

        self._checkOutput(protMovieImport, args, size=1)