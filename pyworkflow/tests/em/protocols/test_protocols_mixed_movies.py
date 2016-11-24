# ***************************************************************************
# * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# ***************************************************************************/

import os
from itertools import izip

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.em import ImageHandler
from pyworkflow.em.protocol import ProtImportMovies

from pyworkflow.em.packages.xmipp3 import (XmippProtMovieCorr,
                                           XmippProtOFAlignment,
                                           XmippProtMovieAverage)

import xmipp


class TestMixedMovies(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsMovies = DataSet.getDataSet('movies')

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

    def _importMovies(self):
        args = self.getArgs('ribo/', pattern='*.mrcs')

        # Id's should be set increasing from 1 if ### is not in the pattern
        protMovieImport = self.newProtocol(ProtImportMovies, **args)
        protMovieImport.setObjLabel('from files')
        self.launchProtocol(protMovieImport)

        self._checkOutput(protMovieImport, args, [1, 2, 3], size=3,
                          dim=(1950, 1950, 16))

        return protMovieImport

    def _compareMovies(self, micSet1, micSet2):
        print "Comparing micrographs (binary images) from results. "

        ih = ImageHandler()
        img1 = ih.createImage()
        img2 = ih.createImage()

        for mic1, mic2 in izip(micSet1, micSet2):
            img1.read(mic1.getFileName())
            img2.read(mic2.getFileName())
            self.assertTrue(img1.equal(img2, 0.001))

    def _sumShifts(self, movieSet):
        """ Sum all shifts a a movie set """
        xx, yy = movieSet.getFirstItem().getAlignment().getShifts()
        return sum([abs(x) for x in xx]) + sum([abs(y) for y in yy])

    def test_CorrelationOpticalFlow(self):
        protMovieImport = self._importMovies()

        mc1 = self.newProtocol(XmippProtMovieCorr,
                               objLabel='CC (no-write)',
                               alignFrame0=2, alignFrameN=10,
                               useAlignToSum=True,
                               splineOrder=XmippProtMovieCorr.INTERP_LINEAR,
                               numberOfThreads=4)
        mc1.inputMovies.set(protMovieImport.outputMovies)
        self.launchProtocol(mc1)

        f0, fN, ffi = mc1.outputMovies.getFramesRange()
        # Here the first frame index should be 2, since we are not
        # re-writing the movie
        self.assertEqual((2, 10, 2), (f0, fN, ffi))
        # Shifts should be different from zero
        self.assertNotAlmostEqual(0, self._sumShifts(mc1.outputMovies))

        mc2 = self.newProtocol(XmippProtMovieCorr,
                               objLabel='CC (write)',
                               alignFrame0=2, alignFrameN=10,
                               useAlignToSum=True,
                               splineOrder=XmippProtMovieCorr.INTERP_LINEAR,
                               doSaveMovie=True,
                               numberOfThreads=4)
        mc2.inputMovies.set(protMovieImport.outputMovies)
        self.launchProtocol(mc2)

        f0, fN, ffi = mc2.outputMovies.getFramesRange()
        # Here the first frame index should be 1, since we wrote a new movie
        self.assertEqual((2, 10, 1), (f0, fN, ffi))
        # All shifts should be zero, since we have written the move and the
        # shifts have been already applied
        self.assertAlmostEqual(0, self._sumShifts(mc2.outputMovies))

        avg1 = self.newProtocol(XmippProtMovieAverage,
                                objLabel='AVG (1)',
                                sumFrame0=2, sumFrameN=10,
                                splineOrder=XmippProtMovieAverage.INTERP_LINEAR,
                                numberOfThreads=4)
        avg1.inputMovies.set(mc1.outputMovies)
        self.launchProtocol(avg1)

        avg2 = self.newProtocol(XmippProtMovieAverage,
                                objLabel='AVG (2)',
                                sumFrame0=2, sumFrameN=10,
                                splineOrder=XmippProtMovieAverage.INTERP_LINEAR,
                                numberOfThreads=4)

        avg2.inputMovies.set(mc2.outputMovies)
        self.launchProtocol(avg2)

        self._compareMovies(avg1.outputMicrographs, avg2.outputMicrographs)

        of1 = self.newProtocol(XmippProtOFAlignment,
                               objLabel='OF (1)',
                               alignFrame0=2, alignFrameN=10,
                               useAlignment=True,
                               numberOfThreads=4)
        of1.inputMovies.set(mc1.outputMovies)
        self.launchProtocol(of1)

        of2 = self.newProtocol(XmippProtOFAlignment,
                               objLabel='OF (2)',
                               alignFrame0=2, alignFrameN=10,
                               useAlignment=True,
                               useAlignToSum=True,
                               numberOfThreads=4,)
        of2.inputMovies.set(mc2.outputMovies)
        self.launchProtocol(of2)

        self._compareMovies(of1.outputMicrographs, of2.outputMicrographs)