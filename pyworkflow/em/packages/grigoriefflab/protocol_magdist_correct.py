# **************************************************************************
# *
# * Authors:     Grigory Sharov (sharov@igbmc.fr)
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

import os
from os.path import exists

from pyworkflow import VERSION_1_1
from pyworkflow.object import Set
import pyworkflow.protocol.params as params
import pyworkflow.protocol.constants as cons
import pyworkflow.utils.path as pwutils
from pyworkflow.em.protocol import ProtProcessMovies
from pyworkflow.em.data import SetOfMovies

from grigoriefflab import MAGDISTCORR_PATH
from convert import parseMagCorrInput


class ProtMagDistCorr(ProtProcessMovies):
    """ This program automatically corrects anisotropic magnification
    distortion using previously estimated parameters.
    It works on a set of movies.
    """
    CONVERT_TO_MRC = 'mrc'
    _label = 'mag distortion correction'
    _version = VERSION_1_1
    doSaveAveMic = False
    doSaveMovie = True

    # --------------------------- DEFINE params functions ----------------------

    def _defineParams(self, form):
        ProtProcessMovies._defineParams(self, form)

        form.addParam('useEst', params.BooleanParam, default=False,
                      label='Use previous estimation?',
                      help='Use previously calculated parameters of '
                           'magnification anisotropy.')
        form.addParam('inputEst', params.PointerParam,
                      pointerClass='ProtMagDistEst', condition='useEst',
                      label='Input protocol',
                      help='Select previously executed estimation protocol.')
        form.addParam('scaleMaj', params.FloatParam, default=1.0,
                      condition='not useEst',
                      label='Major scale factor',
                      help='Major scale factor.')
        form.addParam('scaleMin', params.FloatParam, default=1.0,
                      condition='not useEst',
                      label='Minor scale factor',
                      help='Minor scale factor.')
        form.addParam('angDist', params.FloatParam, default=0.0,
                      condition='not useEst',
                      label='Distortion angle (deg)',
                      help='Distortion angle, in degrees.')
        form.addParam('newPix', params.FloatParam,
                      condition='not useEst',
                      label='New pixel size (A)',
                      help='Assign a new corrected pixel size, in Angstrom.')
        form.addParam('doGain', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Do gain correction before undistorting?',
                      help='If Yes, gain reference that you provided during '
                           'movies import will be chosen.')
        form.addParam('doResample', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Resample images?',
                      help='Resample images after distortion correction and gain'
                           ' correction, by cropping their Fourier transforms')
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

        form.addParallelSection(threads=2, mpi=1)

    # --------------------------- STEPS functions ------------------------------

    def _processMovie(self, movie):
        inputMovies = self.inputMovies.get()
        outputMovieFn = self._getAbsPath(self._getOutputMovieName(movie))
        logFn = self._getAbsPath(self.getOutputLog(movie))

        self._createLink(movie)
        self._argsMagDistCor()

        params = {'movieFn': self._getMovieFn(movie),
                  'outputMovieFn': outputMovieFn,
                  'logFn': logFn,
                  'nthr': self.numberOfThreads.get(),
                  'doGain': 'YES' if self.doGain else 'NO',
                  'doResample': 'YES' if self.doResample else 'NO'
                  }

        if self.useEst:
            inputEst = self.inputEst.get().getOutputLog()
            input_params = parseMagCorrInput(inputEst)
            params['angDist'] = input_params[0]
            params['scaleMaj'] = input_params[1]
            params['scaleMin'] = input_params[2]

        else:
            params['angDist'] = self.angDist.get()
            params['scaleMaj'] = self.scaleMaj.get()
            params['scaleMin'] = self.scaleMin.get()

        if self.doGain:
            params['gainFile'] = inputMovies.getGain()

        if self.doResample:
            params['newX'] = self.newX.get()
            params['newY'] = self.newY.get()

        self._storeSummary(movie)

        try:
            self.runJob(self._program % params, self._args % params)
        except:
            print("ERROR: Distortion correction for movie %s failed\n"
                  % movie.getFileName())

    def createOutputStep(self):
        # Do nothing now, the output should be ready.
        pass

    def _loadOutputSet(self, SetClass, baseName, fixSampling=True):
        """
        Load the output set if it exists or create a new one.
        fixSampling: correct the output sampling rate if binning was used,
        except for the case when the original movies are kept.
        """
        setFile = self._getPath(baseName)

        if os.path.exists(setFile):
            outputSet = SetClass(filename=setFile)
            outputSet.loadAllProperties()
            outputSet.enableAppend()
        else:
            outputSet = SetClass(filename=setFile)
            outputSet.setStreamState(outputSet.STREAM_OPEN)

        inputMovies = self.inputMovies.get()
        outputSet.copyInfo(inputMovies)

        if fixSampling:
            outputSet.setSamplingRate(self.calcPixSize())

        return outputSet

    def _checkNewOutput(self):
        if getattr(self, 'finished', False):
            return

        # Load previously done items (from text file)
        doneList = self._readDoneList()
        # Check for newly done items
        newDone = [m for m in self.listOfMovies
                   if m.getObjId() not in doneList and self._isMovieDone(m)]

        # Update the file with the newly done movies
        # or exit from the function if no new done movies
        if newDone:
            self._writeDoneList(newDone)
        else:
            return

        firstTime = len(doneList) == 0
        allDone = len(doneList) + len(newDone)
        # We have finished when there is not more input movies (stream closed)
        # and the number of processed movies is equal to the number of inputs
        self.finished = self.streamClosed and allDone == len(self.listOfMovies)
        streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN

        # FIXME: Even if we save the movie or not, both are aligned
        suffix = '_corrected'
        movieSet = self._loadOutputSet(SetOfMovies,
                                       'movies%s.sqlite' % suffix,
                                       fixSampling=True)

        for movie in newDone:
            newMovie = self._createOutputMovie(movie)
            movieSet.append(newMovie)

        self._updateOutputSet('outputMovies', movieSet, streamMode)

        if firstTime:
            # Probably is a good idea to store a cached summary for the
            # first resulting movie of the processing.
            self._storeSummary(newDone[0])
            self._defineTransformRelation(self.inputMovies, movieSet)

        if self.finished:  # Unlock createOutputStep if finished all jobs
            outputStep = self._getFirstJoinStep()
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(cons.STATUS_NEW)

    # --------------------------- INFO functions -------------------------------

    def _validate(self):
        errors = []
        # Check that the program exists
        if not exists(MAGDISTCORR_PATH):
            errors.append("Binary '%s' does not exits.\n"
                          "Check configuration file: \n"
                          "~/.config/scipion/scipion.conf\n"
                          "and set MAGDIST_HOME variable properly."
                          % MAGDISTCORR_PATH)

        inputMovies = self.inputMovies.get()
        if self.doGain and inputMovies.getGain() is None:
            errors.append('No gain file was provided during movie import!')

        return errors

    def _citations(self):
        return ["Grant2015"]

    def _summary(self):
        return [self.summaryVar.get()]

    def _methods(self):
        txt = []
        txt.append("Anisotropic magnification distortion was corrected using "
                   "Grigorieff's program *mag_distortion_correct*")
        if self.doGain:
            txt.append("Gain reference was applied before undistorting")
        if self.doResample:
            txt.append("Output was resampled to %d x %d pixels" % (self.newX.get(),
                                                                   self.newY.get()))

        return txt

    # --------------------------- UTILS functions ------------------------------

    def getOutputLog(self, movie):
        return 'micrograph_%06d_Log.txt' % movie.getObjId()

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

    def _createOutputMovie(self, movie):
        correctedMovie = movie.clone()
        extraMovieFn = self._getExtraPath(self._getOutputMovieName(movie))
        correctedMovie.setFileName(extraMovieFn)

        return correctedMovie

    def _storeSummary(self, movie):
        """ Implement this method if you want to store the summary. """
        if self.getAttributeValue('useEst', False):
            inputFn = self.getAttributeValue('inputEst', None).getOutputLog()
            input_params = parseMagCorrInput(inputFn)
            self.summaryVar.set("The following magnification distortion parameters "
                                "were used for correction:\n\n"
                                "Distortion Angle: *%0.2f* degrees\n"
                                "Major Scale: *%0.3f*\n"
                                "Minor Scale: *%0.3f*\n"
                                "Corrected pixel size: *%0.3f* A"
                                % (input_params[0], input_params[1],
                                   input_params[2], input_params[3]))
        else:
            self.summaryVar.set("The following magnification distortion parameters "
                                "were used for correction:\n\n"
                                "Distortion Angle: *%0.2f* degrees\n"
                                "Major Scale: *%0.3f*\n"
                                "Minor Scale: *%0.3f*\n"
                                "Corrected pixel size: *%0.3f* A"
                                % (self.getAttributeValue('angDist', 1.0),
                                   self.getAttributeValue('scaleMaj', 1.0),
                                   self.getAttributeValue('scaleMin', 1.0),
                                   self.newPix.get()))

    def calcPixSize(self):
        if self.useEst:
            inputEst = self.inputEst.get().getOutputLog()
            input_params = parseMagCorrInput(inputEst)
            newPix = input_params[3]
            return newPix
        else:
            return self.newPix.get()
