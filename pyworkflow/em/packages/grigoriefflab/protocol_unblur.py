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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
This module contains the protocol for CTF estimation with Unblur3
"""

import os
import sys
import pyworkflow.utils as pwutils
from pyworkflow.em.protocol import ProtProcessMovies
from pyworkflow.em import ImageHandler, DT_FLOAT, Micrograph
from pyworkflow.protocol.params import  (IntParam,
                                         BooleanParam, FloatParam,
                                         LEVEL_ADVANCED)
#from pyworkflow.utils.properties import Message
from grigoriefflab import UNBLUR_PATH
from pyworkflow.utils.path import createLink, relpath


class ProtUnblur(ProtProcessMovies):
    """Unblur is used to align the frames of movies recorded
    on an electron microscope to reduce image blurring due
    to beam-induced motion. It reads stacks of movies that
    are stored in MRC/CCP4 format. Unblur generates frame
    sums that can be used in subsequent image processing
    steps and optionally applies an exposure-dependent
    filter to maximize the signal at all resolutions
    in the frame averages."""
    _label = 'Unblur'

    def _defineParams(self, form):
        ProtProcessMovies._defineParams(self, form)
        form.addParam('alignFrameRange', IntParam,
                      default=-1,
                      label='Number frames per movie ',
                      help='How many frames per monie. -1 -> all frames ')
        form.addParam('doApplyDoseFilter', BooleanParam, default=True,
                      label='Apply Dose filter',
                      help='Apply a dose-dependent filter to frames before summing them')
        form.addParam('exposurePerFrame', FloatParam,
                      label='Exposure per frame (e/A^2)',
                      help='Exposure per frame, in electrons per square Angstrom')

        #group = form.addGroup('Expert Options')
        form.addParam('minShiftInitSearch', FloatParam,
                      default=2.,
                      label='Min. Shift Initial search (A)',
                      help='Initial search will be limited to between the inner and outer radii',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('OutRadShiftLimit', FloatParam,
                      default=200.,
                      label='Outer radius shift limit (A)',
                      help='The maximum shift of each alignment step will be limited to this value',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('Bfactor', FloatParam,
                      default=1500.,
                      label='B-factor (A^2)',
                      help='B-factor to apply to images (A^2)',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('HWVertFourMask', IntParam,
                      default=1,
                      label='Half-width vertical Fourier mask',
                      help='The vertical line mask will be twice this size. The central cross mask helps reduce problems by line artefacts from the detector',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('HWHoriFourMask', IntParam,
                      default=1,
                      label='Half-width horizontal Fourier mask',
                      help='The horizontal line mask will be twice this size. The central cross mask helps reduce problems by line artefacts from the detector',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('terminationShiftThreshold', FloatParam,
                      default=0.1,
                      label='Termination shift threshold',
                      help='Alignment will stop at this number, even if the threshold shift is not reached',
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

    #--------------------------- STEPS functions ----------------------------------------------

    def _getMicFnName(self, movieId,movieFolder):
        return relpath(self._getExtraPath('aligned_sum_%0d5.mrc'%movieId), movieFolder)

    def _getMicFnNameFromRun(self, movieId):
        return self._getExtraPath('aligned_sum_%0d5.mrc'%movieId)

    def _getShiftFnName(self, movieId):
        return 'shifts_%0d5.txt'%movieId

    def _getFSCFnName(self, movieId):
        return 'fsc_%0d5.txt'%movieId

    def _filterMovie(self, movieId, movieFn):
        """I really do not understand how I end writing this function ROB"""
        return True

    def _processMovie(self, movieId, movieName, movieFolder):
        """call program here"""
        # if not mrc convert format to mrc
        # special case is mrc but ends in mrcs
        inMovieName= os.path.join(movieFolder,movieName)
        if movieName.endswith('.mrc'):
            movieNameAux = inMovieName
        elif movieName.endswith('.mrcs'):
            movieNameAux= pwutils.replaceExt(inMovieName, "mrc")
            createLink(inMovieName,movieNameAux)
            movieNameAux = pwutils.replaceExt(movieName, "mrc")
        else:
            micFnMrc = pwutils.replaceExt(inMovieName, "mrc")
            ImageHandler().convert(inMovieName, micFnMrc, DT_FLOAT)
            movieNameAux = pwutils.replaceExt(movieName, "mrc")

        #get number of frames
        if self.alignFrameRange == -1:
            numberOfFramesPerMovie = self.inputMovies.get().getDimensions()[3]
        else:
            numberOfFramesPerMovie = self.alignFrameRange.get()

        doApplyDoseFilter = self.doApplyDoseFilter.get()
        exposurePerFrame = self.exposurePerFrame.get()
        self._argsUnblur(movieNameAux
                        , movieFolder
                        , movieId
                        , numberOfFramesPerMovie
                        , doApplyDoseFilter
                        , exposurePerFrame
                        )
        try:

            self.runJob(self._program, self._args, cwd=movieFolder)
        except:
            print("ERROR: Movie %s failed\n"%movieName)

        logFile = self._getLogFile(movieId)

    def _argsUnblur(self
                    , movieName
                    , movieFolder
                    , movieId
                    , numberOfFramesPerMovie
                    , doApplyDoseFilter
                    , exposurePerFrame):
        """format input as unblur likes it"""
        args = {}

        args['movieName'] = movieName
        args['numberOfFramesPerMovie'] = numberOfFramesPerMovie
        args['micFnName'] = self._getMicFnName(movieId,movieFolder)
        args['shiftFnName'] = self._getShiftFnName(movieId)
        args['samplingRate'] = self.samplingRate
        args['voltage'] = self.inputMovies.get().getAcquisition().getVoltage()
        args['fscFn'] = self._getFSCFnName(movieId)
        args['Bfactor'] = self.Bfactor.get()
        args['minShiftInitSearch'] = self.minShiftInitSearch.get()
        args['OutRadShiftLimit'] = self.OutRadShiftLimit.get()
        args['HWVertFourMask'] = self.HWVertFourMask.get()
        args['HWHoriFourMask'] = self.HWHoriFourMask.get()
        args['terminationShiftThreshold'] = self.terminationShiftThreshold.get()
        args['maximumNumberIterations'] = self.maximumNumberIterations.get()
        doRestoreNoisePower = self.doRestoreNoisePower.get()
        doVerboseOutput = self.doVerboseOutput.get()

        args['HWVertFourMask'] = self.HWVertFourMask.get()

        args['HWVertFourMask'] = self.HWVertFourMask.get()

        if doApplyDoseFilter:
            args['doApplyDoseFilter'] = 'YES'
        else:
            args['doApplyDoseFilter'] = 'NO'
        if doRestoreNoisePower:
            args['doRestoreNoisePower'] = 'YES'
        else:
            args['doRestoreNoisePower'] = 'NO'
        if doVerboseOutput:
            args['doVerboseOutput'] = 'YES'
        else:
            args['doVerboseOutput'] = 'NO'

        args['exposurePerFrame'] = exposurePerFrame

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
"""%args

    def createOutputStep(self):
        inputMovies = self.inputMovies.get()
        micSet = self._createSetOfMicrographs()
        micSet.copyInfo(inputMovies)
        for movie in self.inputMovies.get():
            mic = Micrograph()
            # All micrograph are copied to the 'extra' folder after each step
            mic.setFileName(self._getMicFnNameFromRun(movie.getObjId()))
            micSet.append(mic)
        self._defineOutputs(outputMicrographs=micSet)

        # Also create a Set of Movies with the alignment parameters
        #movieSet = self._createSetOfMovies()
        #movieSet.copyInfo(inputMovies)
        #movieSet.cropOffsetX = Integer(self.cropOffsetX)
        #movieSet.cropOffsetY = Integer(self.cropOffsetY)

    def _citations(self):
        return []

    def _summary(self):
        return []

    def _methods(self):
        return []

    def visualize(self):
        pass
