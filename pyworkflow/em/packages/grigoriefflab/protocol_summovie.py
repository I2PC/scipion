# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              Josue Gomez Blanco (jgomez@cnb.csic.es)
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

from os.path import exists, realpath, abspath

import pyworkflow.protocol.params as params
import pyworkflow.protocol.constants as cons
import pyworkflow.utils.path as pwutils
from pyworkflow.em.protocol import ProtAlignMovies, ProtProcessMovies

from grigoriefflab import SUMMOVIE_PATH
from convert import writeShiftsMovieAlignment


class ProtSummovie(ProtAlignMovies):
    """ Summovie generates frame sums that can be used
    in subsequent image processing steps and optionally
    applies an exposure-dependent filter to maximize
    the signal at all resolutions in the frame averages.
    """
    _label = 'summovie'
    CONVERT_TO_MRC = 'mrc'
    doSaveAveMic = True
    
    #--------------------------- DEFINE param functions ------------------------
    def _defineAlignmentParams(self, form):
        form.addHidden('binFactor', params.FloatParam, default=1.)

        group = form.addGroup('Average')
        line = group.addLine('Frames to SUM',
                             help='Frames range to SUM on each movie. The '
                                  'first frame is 1. If you set 0 in the final '
                                  'frame to sum, it means that you will sum '
                                  'until the last frame of the movie.')
        line.addParam('sumFrame0', params.IntParam, default=1,
                      label='from')
        line.addParam('sumFrameN', params.IntParam, default=0,
                      label='to')
        group.addParam('useAlignment', params.BooleanParam, default=True,
                       label="Use movie alignment to Sum frames?")

        form.addParam('doApplyDoseFilter', params.BooleanParam, default=True,
                      label='Apply Dose filter',
                      help='Apply a dose-dependent filter to frames before '
                           'summing them. Pre-exposure and dose per frame were '
                           'specified during movies import.')

        form.addParam('doRestoreNoisePower', params.BooleanParam,
                      default=True,
                      label='Restore Noise Power? ',
                      help='Restore Noise Power? ',
                      expertLevel=cons.LEVEL_ADVANCED)

        form.addParam('cleanInputMovies', params.BooleanParam, default=False,
                      expertLevel=cons.LEVEL_ADVANCED,
                      label='Clean input movies?',
                      help='If set *Yes*, the input movies will be deleted. '
                           'This option is useful if your input movies are '
                           'an intermediates files, obtained with local '
                           'alignment protocols; for instance, *Optical* '
                           '*Flow* protocol or in cases of global alignment '
                           'protocols when you save the movies. If in the '
                           'previous alignment protocol, you save *only* the '
                           'alignment, you *MUST* set *NO*')
        
        form.addParallelSection(threads=1, mpi=1)
    
    #--------------------------- STEPS functions -------------------------------
    def _processMovie(self, movie):
        try:
            self._createLink(movie)
            numberOfFrames = self._getNumberOfFrames(movie)

            if numberOfFrames is None:
                raise Exception("Could not read number of frames.")

            s0, sN = self._getRange(movie, 'sum')

            self._writeMovieAlignment(movie, s0, sN)
            self._argsSummovie()

            params = {'movieFn': self._getMovieFn(movie),
                      'numberOfFrames': numberOfFrames,
                      'micFn': self._getMicFn(movie),
                      'initFrame': s0,
                      'finalFrame': sN,
                      'shiftsFn': self._getShiftsFn(movie),
                      'samplingRate': movie.getSamplingRate(),
                      'voltage': movie.getAcquisition().getVoltage(),
                      'frcFn': self._getFrcFn(movie),
                      'exposurePerFrame': movie.getAcquisition().getDosePerFrame() or 0.0,
                      'doApplyDoseFilter': 'YES' if self.doApplyDoseFilter else 'NO',
                      'doRestoreNoisePower': 'YES' if self.doRestoreNoisePower else 'NO'
                      }

            self.runJob(self._program, self._args % params)
            self._storeSummary(movie)
            if self.cleanInputMovies:
                pwutils.cleanPath(movie._originalFileName.get())
                print ("Movie %s erased" % movie._originalFileName.get())
                
        except Exception as e:
            print("ERROR: Movie %s failed.\n"
                  "       Message: %s\n" % (movie.getFileName(), e))
    
    #--------------------------- INFO functions --------------------------------
    def _citations(self):
        return []
        
    def _methods(self):
        return []

    def _validate(self):

        if self.inputMovies.get() is None:
            return

        errors = []

        inputMovies = self.inputMovies.get()
        firstFrame, lastFrame, _ = inputMovies.getFramesRange()
        frames = lastFrame - firstFrame + 1

        if frames is not None:
            s0, sN = self._getFrameRange(frames, "sum")
            if sN < firstFrame or sN > lastFrame:
                errors.append("Check the selected last frame to *SUM*. "
                              "Last frame (%d) should be in range: %s "
                              % (sN, (firstFrame, lastFrame)))
            if s0 < firstFrame or s0 > lastFrame:
                errors.append("Check the selected first frame to *SUM*. "
                              "First frame (%d) should be in range: %s "
                              % (s0, (firstFrame, lastFrame)))
            if sN < s0:
                errors.append("Check the selected frames range to *SUM*. "
                              "Last frame (%d) should be greater or equal "
                              "than first frame (%d)"
                              % (sN, s0))

        if not exists(SUMMOVIE_PATH):
            errors.append("Cannot find the Summovie program at: %s"
                          % SUMMOVIE_PATH)

        # If this protocol is going to delete the input movies, let's double
        # check that they has been produced by the previous protocol and
        # are not the original input movies data, what can be a serious problem
        if self.cleanInputMovies:
            firstMovie = inputMovies.getFirstItem()
            fullPath = realpath(firstMovie.getFileName())
            prevProt = self.getMapper().getParent(inputMovies)
            protPath = abspath(prevProt.getWorkingDir())
            if not fullPath.startswith(protPath):
                errors.append("You have selected to delete the input movies "
                              "and they are not produced by the previous "
                              "protocol. This is not allowed because you could "
                              "lost your original data.")

        if self.doApplyDoseFilter:
            inputMovies = self.inputMovies.get()
            doseFrame = inputMovies.getAcquisition().getDosePerFrame()

            if doseFrame == 0.0 or doseFrame is None:
                errors.append('Dose per frame for input movies is 0 or not '
                              'set. You cannot apply dose filter.')

        return errors
    
    #--------------------------- UTILS functions -------------------------------
    def _argsSummovie(self):

        # Avoid threads multiplication
        # self._program = 'export OMP_NUM_THREADS=%d; ' % self.numberOfThreads.get()
        self._program = 'export OMP_NUM_THREADS=1; '
        self._program += SUMMOVIE_PATH
        self._args = """ << eof
%(movieFn)s
%(numberOfFrames)s
%(micFn)s
%(shiftsFn)s
%(frcFn)s
%(initFrame)d
%(finalFrame)d
%(samplingRate)f
%(doApplyDoseFilter)s
%(exposurePerFrame)f
%(voltage)f
0
%(doRestoreNoisePower)s
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
        
    def _getMicFn(self, movie):
        return self._getExtraPath(self._getOutputMicName(movie))
    
    def _getShiftsFn(self, movie):
        return self._getExtraPath(self._getMovieRoot(movie) + '_shifts.txt')
    
    def _getFrcFn(self, movie):
        return self._getExtraPath(self._getMovieRoot(movie) + '_frc.txt')
    
    def _writeMovieAlignment(self, movie, s0, sN):
        shiftsFn = self._getShiftsFn(movie)

        if movie.hasAlignment() and self.useAlignment:
            writeShiftsMovieAlignment(movie, shiftsFn, s0, sN)
        else:
            f = open(shiftsFn, 'w')
            frames = sN - s0 + 1
            shift = ("0 " * frames + "\n") * 2
            f.write(shift)
            f.close()
    
    def _createOutputMovies(self):
        """ Returns True if an output set of movies will be generated.
        The most common case is to always generate output movies,
        either with alignment only or the binary aligned movie files.
        Subclasses can override this function to change this behavior.
        """
        return False

    def _storeSummary(self, movie):
        if movie.hasAlignment():
            numOfFrames = self._getNumberOfFrames(movie)
            s0, sN = self._getFrameRange(numOfFrames, 'sum')
            fstFrame, lstFrame, _ = movie.getFramesRange()

            if self.useAlignment and (fstFrame > s0 or lstFrame < sN):
                self.summaryVar.set("Warning!!! You have selected a frame range "
                                    "wider than the range selected to align. All "
                                    "the frames selected without alignment "
                                    "information, will be aligned by setting "
                                    "alignment to 0")

    def _getNumberOfFrames(self, movie):
        _, lstFrame, _ = movie.getFramesRange()

        if movie.hasAlignment():
            _, lastFrmAligned = movie.getAlignment().getRange()
            if lastFrmAligned != lstFrame:
                return lastFrmAligned
            else:
                return movie.getNumberOfFrames()
        else:
            return movie.getNumberOfFrames()

    def _getRange(self, movie, prefix):
        n = self._getNumberOfFrames(movie)
        iniFrame, _, indxFrame = movie.getFramesRange()
        first, last = self._getFrameRange(n, prefix)

        if iniFrame != indxFrame:
            first -= (iniFrame - 1)
            last -= (iniFrame - 1)

        return first, last
