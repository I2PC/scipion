# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
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

import pyworkflow.utils as pwutils
from pyworkflow.em.protocol import ProtAlignMovies
from pyworkflow.protocol.params import  (IntParam,
                                         BooleanParam, FloatParam,
                                         LEVEL_ADVANCED)
from grigoriefflab import UNBLUR_PATH
from convert import writeShiftsMovieAlignment



class ProtUnblur(ProtAlignMovies):
    """ Unblur is used to align the frames of movies recorded on an electron
    microscope to reduce image blurring due to beam-induced motion.
    """
    _label = 'unblur'
    CONVERT_TO_MRC = 'mrc'

    def _defineAlignmentParams(self, form):
        form.addParam('alignFrameRange', IntParam,
                      default=-1,
                      label='Number frames per movie ',
                      help='How many frames per movie. -1 -> all frames ')
        form.addParam('doApplyDoseFilter', BooleanParam, default=True,
                      label='Apply Dose filter',
                      help='Apply a dose-dependent filter to frames before '
                           'summing them')
        form.addParam('exposurePerFrame', FloatParam,
                      label='Exposure per frame (e/A^2)',
                      help='Exposure per frame, in electrons per square '
                           'Angstrom')

        #group = form.addGroup('Expert Options')
        form.addParam('minShiftInitSearch', FloatParam,
                      default=2.,
                      label='Min. Shift Initial search (A)',
                      help='Initial search will be limited to between the '
                           'inner and outer radii',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('OutRadShiftLimit', FloatParam,
                      default=200.,
                      label='Outer radius shift limit (A)',
                      help='The maximum shift of each alignment step will be '
                           'limited to this value',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('Bfactor', FloatParam,
                      default=1500.,
                      label='B-factor (A^2)',
                      help='B-factor to apply to images (A^2)',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('HWVertFourMask', IntParam,
                      default=1,
                      label='Half-width vertical Fourier mask',
                      help='The vertical line mask will be twice this size. '
                           'The central cross mask helps reduce problems by '
                           'line artefacts from the detector',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('HWHoriFourMask', IntParam,
                      default=1,
                      label='Half-width horizontal Fourier mask',
                      help='The horizontal line mask will be twice this size. '
                           'The central cross mask helps reduce problems by '
                           'line artefacts from the detector',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('terminationShiftThreshold', FloatParam,
                      default=0.1,
                      label='Termination shift threshold',
                      help='Alignment will stop at this number, even if the '
                           'threshold shift is not reached',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('maximumNumberIterations', IntParam,
                      default=10,
                      label='Maximum number of iterations',
                      help='Maximum number of iterations',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('doRestoreNoisePower', BooleanParam,
                      default=True,
                      label='Restore Noise Power? ',
                      help='Restore Noise Power? ',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('doVerboseOutput', BooleanParam,
                      default=False,
                      label='Verbose Output?',
                      help='Verbose Output?',
                      expertLevel=LEVEL_ADVANCED)
        form.addParallelSection(threads=1, mpi=1)

    #
    #Input stack filename                [my_movie.mrc] : kk.mrc
    #Number of frames per movie                    [34] :
    #Output aligned sum file       [my_aligned_sum.mrc] :
    #Output shifts file                 [my_shifts.txt] :
    #Pixel size of images (A)                       [1] :
    #Apply Dose filter?                            [NO] : YES
    #Exposure per frame (e/A^2)                   [1.0] :
    #Acceleration voltage (kV)                  [300.0] :
    #Set Expert Options?                           [NO] : YES
    #Output FRC file                       [my_frc.txt] :
    #Minimum shift for initial search (Angstroms)
    #[2.0]                                              :
    #Outer radius shift limit (Angstroms)       [200.0] :
    #B-factor to apply to images (A^2)           [1500] :
    #Half-width of central vertical line of Fourier mask
    #[1]                                                :
    #Half-width of central horizontal line of Fourier mask
    #[1]                                                :
    #Termination shift threshold                  [0.1] :
    #Maximum number of iterations                  [10] :
    #Restore Noise Power?                         [YES] :
    #Verbose Output?                               [NO] : YES

    #--------------------------- STEPS functions ------------------------------

    def _processMovie(self, movie):
        self._createLink(movie)
        numberOfFrames = movie.getNumberOfFrames()
        #FIXME: Figure out how to properly write shifts for unblur
        #self._writeMovieAlignment(movie, numberOfFrames)

        if self.alignFrameRange != -1:

            if self.alignFrameRange > numberOfFrames:
                raise Exception('Frame number (%d) is greater than '
                                'the total frames of the movie (%d)' %
                                (numberOfFrames, self.alignFrameRange))
            numberOfFrames = self.alignFrameRange.get()

        self._argsUnblur(movie, numberOfFrames)

        try:
            self.runJob(self._program, self._args)
        except:
            print("ERROR: Movie %s failed\n" % movie.getName())

    def _argsUnblur(self, movie, numberOfFrames):
        """ Format argument for call unblur program. """
        args = {}

        args['movieName'] = self._getMovieFn(movie)
        args['numberOfFramesPerMovie'] = numberOfFrames
        args['micFnName'] = self._getMicFn(movie)
        args['shiftFnName'] = self._getShiftsFn(movie)
        args['samplingRate'] = self.samplingRate
        args['voltage'] = self.inputMovies.get().getAcquisition().getVoltage()
        args['fscFn'] = self._getFrcFn(movie)
        args['Bfactor'] = self.Bfactor.get()
        args['minShiftInitSearch'] = self.minShiftInitSearch.get()
        args['OutRadShiftLimit'] = self.OutRadShiftLimit.get()
        args['HWVertFourMask'] = self.HWVertFourMask.get()
        args['HWHoriFourMask'] = self.HWHoriFourMask.get()
        args['terminationShiftThreshold'] = self.terminationShiftThreshold.get()
        args['maximumNumberIterations'] = self.maximumNumberIterations.get()
        args['doApplyDoseFilter'] = 'YES' if self.doApplyDoseFilter else 'NO'
        args['doRestoreNoisePower'] = 'YES' if self.doRestoreNoisePower else 'NO'
        args['doVerboseOutput'] = 'YES' if self.doVerboseOutput else 'NO'
        args['exposurePerFrame'] = self.exposurePerFrame.get()

        self._program = 'export OMP_NUM_THREADS=1; ' + UNBLUR_PATH
        self._args = """ << eof
%(movieName)s
%(numberOfFramesPerMovie)s
%(micFnName)s
%(shiftFnName)s
%(samplingRate)f
%(doApplyDoseFilter)s
%(exposurePerFrame)f
%(voltage)f
YES
%(fscFn)s
%(minShiftInitSearch)f
%(OutRadShiftLimit)f
%(Bfactor)f
%(HWVertFourMask)d
%(HWHoriFourMask)d
%(terminationShiftThreshold)f
%(maximumNumberIterations)d
%(doRestoreNoisePower)s
%(doVerboseOutput)s
eof
""" % args

    def _citations(self):
        return []

    def _summary(self):
        return []

    def _methods(self):
        return []

    def _validate(self):
        errors = []
        if not os.path.exists(UNBLUR_PATH):
            errors.append("Cannot find the Unblur program at: " + UNBLUR_PATH)
        return errors

    def _getMicName(self, movieName):
        """ Return the name for the output micrograph given the movie name.
        """
        return self._getExtraPath(self._getNameExt(movieName, '', 'mrc'))

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

    def _writeMovieAlignment(self, movie, numberOfFrames):
        shiftsFn = self._getShiftsFn(movie)
        s0, sN = self._getFrameRange(numberOfFrames)

        if movie.hasAlignment() and self.useAlignment:
            writeShiftsMovieAlignment(movie, shiftsFn, s0, sN)
        else:
            f = open(shiftsFn, 'w')
            frames = sN - s0 + 1
            shift = ("0 " * frames + "\n") * 2
            f.write(shift)
            f.close()

    def _getFrameRange(self, n, prefix=''):
        """
        Params:
        :param n: Number of frames of the movies
        :param prefix: what range we want to consider, either 'align' or 'sum'
        :return: (i, f) initial and last frame range
        """
        last = self.alignFrameRange.get() if self.alignFrameRange > 0 else n

        return 1, min(last, n)