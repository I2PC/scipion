# **************************************************************************
# *
# * Authors:     Grigory Sharov (sharov@igbmc.fr)
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
# *
# **************************************************************************

import os
from os.path import exists

import pyworkflow.protocol.params as params
import pyworkflow.utils.path as pwutils
from pyworkflow.em.protocol import ProtProcessMovies

from grigoriefflab import MAGDISTCORR_PATH


class ProtMagDistCorr(ProtProcessMovies):
    """ This program automatically corrects anisotropic magnification
    distortion using previously estimated parameters
    """
    CONVERT_TO_MRC = 'mrc'
    _label = 'magnification distortion correction'

    # --------------------------- DEFINE params functions --------------------------------------------

    def _defineParams(self, form):
        ProtProcessMovies._defineParams(self, form)

        form.addParam('scaleMaj', params.FloatParam, default=1.0,
                      label='Major scale factor',
                      help='Major scale factor.')
        form.addParam('scaleMin', params.FloatParam, default=1.0,
                      label='Minor scale factor',
                      help='Minor scale factor.')
        form.addParam('angDist', params.FloatParam, default=0.0,
                      label='Distortion angle (deg)',
                      help='Distortion angle, in degrees.')
        form.addParam('doGain', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Do gain correction before undistorting?',
                      help='If Yes, gain reference that you provided during movies import will be chosen.')
        form.addParam('doResample', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Resample images?',
                      help='Resample images after distortion correction and gain correction, '
                           'by cropping their Fourier transforms')
        line = form.addLine('New dimensions (px)',
                            expertLevel=params.LEVEL_ADVANCED,
                            condition='doResample')
        line.addParam('newX', params.IntParam, default=2048,
                      label='X',
                      expertLevel=params.LEVEL_ADVANCED,
                      condition='doResample')
        line.addParam('newY', params.IntParam, default=2048,
                      label='Y',
                      expertLevel=params.LEVEL_ADVANCED,
                      condition='doResample')

        form.addParallelSection(threads=2, mpi=0)

    # --------------------------- STEPS functions --------------------------------------------

    def _processMovie(self, movie):
        inputMovies = self.inputMovies.get()
        outputMovieFn = self._getAbsPath(self._getOutputMovieName(movie))
        logFn = self.getOutputLog()

        self.doSaveMovie = True
        self._createLink(movie)
        self._argsMagDistCor()

        params = {'movieFn': self._getMovieFn(movie),
                  'outputMovieFn': outputMovieFn,
                  'logFn': logFn,
                  'scaleMaj': self.scaleMaj.get(),
                  'scaleMin': self.scaleMin.get(),
                  'angDist': self.angDist.get(),
                  'nthr': self.numberOfThreads.get(),
                  'doGain': 'YES' if self.doGain else 'NO',
                  'doResample': 'YES' if self.doResample else 'NO'
                  }

        if self.doGain:
            params['gainFile'] = inputMovies.getGain()

        if self.doResample:
            params['newX'] = self.newX.get()
            params['newY'] = self.newY.get()

        try:
            self.runJob(self._program % params, self._args % params)
            if self.cleanMovieData:
                pwutils.cleanPath(movie._originalFileName.get())
                print ("Movie %s erased" % movie._originalFileName.get())

        except:
            print("ERROR: Distortion correction for movie %s failed\n" % movie.getFileName())

    def createOutputStep(self):
        # Do nothing now, the output should be ready.
        pass

    # --------------------------- INFO functions ----------------------------------------------------

    def _validate(self):
        validateMsgs = []
        # Check that the program exists
        if not exists(MAGDISTCORR_PATH):
            validateMsgs.append("Binary '%s' does not exits.\n"
                          "Check configuration file: \n"
                          "~/.config/scipion/scipion.conf\n"
                          "and set MAGDIST_HOME variable properly."
                          % MAGDISTCORR_PATH)

        return validateMsgs

    def _citations(self):
        return ["Grant2015"]

    def _summary(self):
        summary = []

        return summary

    def _methods(self):
        txt = []
        txt.append("Anisotropic magnification distortion was corrected using "
                   "Grigorieff's program *mag_distortion_correct*")

        return txt

    # --------------------------- UTILS functions --------------------------------------------

    def getOutputLog(self):
        return self._getExtraPath('mag_dist_correction.log')

    def _argsMagDistCor(self):
        self._program = 'export NCPUS=%(nthr)d ; ' + MAGDISTCORR_PATH

        if self.doGain and self.doResample:
            self._args = """   << eof > %(logFn)s
%(movieFn)s
%(outputMovieFn)s
%(angDist)f
%(scaleMaj)f
%(scaleMin)f
%(doGain)s
%(gainFile)s
%(doResample)s
%(newX)d
%(newY)d
eof
"""
        elif self.doGain and not self.doResample:
            self._args = """   << eof > %(logFn)s
%(movieFn)s
%(outputMovieFn)s
%(angDist)f
%(scaleMaj)f
%(scaleMin)f
%(doGain)s
%(gainFile)s
%(doResample)s
eof
"""
        elif not self.doGain and self.doResample:
            self._args = """   << eof > %(logFn)s
%(movieFn)s
%(outputMovieFn)s
%(angDist)f
%(scaleMaj)f
%(scaleMin)f
%(doGain)s
%(doResample)s
%(newX)d
%(newY)d
eof
"""

        else:
            self._args = """   << eof > %(logFn)s
%(movieFn)s
%(outputMovieFn)s
%(angDist)f
%(scaleMaj)f
%(scaleMin)f
%(doGain)s
%(doResample)s
eof
"""

    def _getMovieFn(self, movie):
        movieFn = movie.getFileName()
        if movieFn.endswith("mrcs"):
            return pwutils.replaceExt(movieFn, self.CONVERT_TO_MRC)
        else:
            return movieFn

    def _createLink(self, movie):
        movieFn = movie.getFileName()
        if movieFn.endswith("mrcs"):
            pwutils.createLink(movieFn, self._getMovieFn(movie))

    def _getAbsPath(self, baseName):
        return os.path.abspath(self._getExtraPath(baseName))

    def _getMovieRoot(self, movie):
        return pwutils.removeBaseExt(movie.getFileName())

    def _getOutputMovieName(self, movie):
        """ Returns the name of the output movie.
        (relative to micFolder)
        """
        return self._getMovieRoot(movie) + '_corrected.mrc'