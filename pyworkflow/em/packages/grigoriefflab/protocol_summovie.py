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
This module contains the protocol for CTF estimation with Summovie3
"""

import os
import sys
import pyworkflow.utils as pwutils
from os.path import exists
from pyworkflow.em.protocol import ProtProcessMovies
from pyworkflow.em import ImageHandler, DT_FLOAT, Micrograph
from pyworkflow.protocol.params import  (IntParam,
                                         BooleanParam, FloatParam,
                                         LEVEL_ADVANCED)
from grigoriefflab import SUMMOVIE_PATH
from pyworkflow.utils.path import createLink, relpath, removeBaseExt


class ProtSummovie(ProtProcessMovies):
    """Summovie is used to align the frames of movies recorded
    on an electron microscope to reduce image blurring due
    to beam-induced motion. It reads stacks of movies that
    are stored in MRC/CCP4 format. Summovie generates frame
    sums that can be used in subsequent image processing
    steps and optionally applies an exposure-dependent
    filter to maximize the signal at all resolutions
    in the frame averages."""
    _label = 'Summovie'

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
        form.addParam('Bfactor', FloatParam,
                      default=1500.,
                      label='B-factor (A^2)',
                      help='B-factor to apply to images (A^2)',
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
    #INPUT_FILENAME kk.mrc
    #number_of_frames_per_movie 34
    #output_filename my_aligned.mrc
    #shifts_filename 0_shifts.txt
    #frc_filename my_frc.txt
    #first_frame 1
    #LAST_FRAME 34
    #PIXEL_SIZE 1.0
    #apply_dose_filter yes
    #dose_per_frame 1.0
    #acceleration_voltage 300.0
    #pre_exposure_amount 0.0
    #restore_power yes

    #--------------------------- STEPS functions ----------------------------------------------

    def _getMicFnName(self, movieId,movieFolder):
        return relpath(self._getExtraPath('aligned_sum_%06d.mrc'%movieId), movieFolder)

    def _getMicFnNameFromRun(self, movieId):
        return self._getExtraPath('aligned_sum_%06d.mrc'%movieId)

    def _getNameExt(self, movieName, postFix, ext):
        if movieName.endswith("bz2"):
            # removeBaseExt function only eliminate the last extension,
            # but if files are compressed, we need to eliminate two extensions:
            # bz2 and its own image extension (e.g: mrcs, em, etc)
            return removeBaseExt(removeBaseExt(movieName)) + postFix + '.' + ext
        else:
            return removeBaseExt(movieName) + postFix + '.' + ext

    def _getShiftFnName(self, movieId):
        return 'shifts_%06d.txt'%movieId

    def _getFSCFnName(self, movieId):
        return 'fsc_%06d.txt'%movieId

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
        #write dummy auxiliary shift file.
        #TODO: this should be done properly when we define how to transfer shift
        #between movies
        shiftFnName= os.path.join(movieFolder,self._getShiftFnName(movieId))
        f=open(shiftFnName,'w')
        shift= ("0 " * numberOfFramesPerMovie + "\n" ) *2
        f.write(shift)
        f.close()

        doApplyDoseFilter = self.doApplyDoseFilter.get()
        exposurePerFrame = self.exposurePerFrame.get()
        self._argsSummovie(movieNameAux
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

    def _argsSummovie(self
                    , movieName
                    , movieFolder
                    , movieId
                    , numberOfFramesPerMovie
                    , doApplyDoseFilter
                    , exposurePerFrame):
        """format input as Summovie likes it"""
        args = {}

        args['movieName'] = movieName
        args['numberOfFramesPerMovie'] = numberOfFramesPerMovie
        ##args['micFnName'] = self._getMicFnName(movieId,movieFolder)
        args['micFnName'] = relpath(self._getExtraPath(self._getNameExt(movieName,'', 'mrc')),movieFolder)
        args['shiftFnName'] = self._getShiftFnName(movieId)
        args['samplingRate'] = self.samplingRate
        args['voltage'] = self.inputMovies.get().getAcquisition().getVoltage()
        args['fscFn'] = self._getFSCFnName(movieId)
        args['Bfactor'] = self.Bfactor.get()
        doRestoreNoisePower = self.doRestoreNoisePower.get()

        if doApplyDoseFilter:
            args['doApplyDoseFilter'] = 'YES'
        else:
            args['doApplyDoseFilter'] = 'NO'
        if doRestoreNoisePower:
            args['doRestoreNoisePower'] = 'YES'
        else:
            args['doRestoreNoisePower'] = 'NO'

        args['exposurePerFrame'] = exposurePerFrame

        self._program = 'export OMP_NUM_THREADS=1; ' + SUMMOVIE_PATH
        self._args = """ << eof
%(movieName)s
%(numberOfFramesPerMovie)s
%(micFnName)s
%(shiftFnName)s
%(fscFn)s
1
%(numberOfFramesPerMovie)s
%(samplingRate)f
%(doApplyDoseFilter)s
%(exposurePerFrame)f
%(voltage)f
0
%(doRestoreNoisePower)s
eof
"""%args

    def createOutputStep(self):
        inputMovies = self.inputMovies.get()
        micSet = self._createSetOfMicrographs()
        micSet.copyInfo(inputMovies)
        for movie in self.inputMovies.get():
            mic = Micrograph()
            # All micrograph are copied to the 'extra' folder after each step
            mic.setFileName(self._getExtraPath(self._getNameExt(movie.getFileName(),'', 'mrc')))
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

    def _validate(self):
        errors=[]
        if not exists(SUMMOVIE_PATH):
            errors.append("Cannot find the Summovie program at: "+SUMMOVIE_PATH)
        return errors

