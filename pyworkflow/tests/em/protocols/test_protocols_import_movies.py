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
    
    def test_pattern(self):
        """ Import several micrographs from a given pattern.
        """
        args = {'importFrom': ProtImportMovies.IMPORT_FROM_FILES,
                'filesPath': self.dsMovies.getFile('ribo/'),
                'filesPattern': '*.mrcs',
                'amplitudConstrast': 0.1,
                'sphericalAberration': 2.,
                'voltage': 300,
                'samplingRate': 3.54
                }
        
        def _checkOutput(prot, moviesId=[], size=None):
            movies = getattr(prot, 'outputMovies', None)
            self.assertIsNotNone(movies)
            self.assertEqual(movies.getSize(), size)
            for i, m in enumerate(movies):
                if moviesId:
                    self.assertEqual(m.getObjId(), moviesId[i])
                self.assertAlmostEqual(m.getSamplingRate(), args['samplingRate'])
                a = m.getAcquisition()
                self.assertAlmostEqual(a.getVoltage(), args['voltage'])

        # Id's should be set increasing from 1 if ### is not in the 
        # pattern
        protMovieImport = self.newProtocol(ProtImportMovies, **args)
        protMovieImport.setObjLabel('from files')
        self.launchProtocol(protMovieImport)
        _checkOutput(protMovieImport, [1, 2, 3], size=3)
        
        # Id's should be taken from filename    
        # args['filesPattern'] = 'BPV_####.mrc'
        # protMovieImport = self.newProtocol(ProtImportMovies, **args)
        # protMovieImport.setObjLabel('from files (with id)')
        # self.launchProtocol(protMovieImport)
        # _checkOutput(protMovieImport, [1386, 1387, 1388], size=3)

