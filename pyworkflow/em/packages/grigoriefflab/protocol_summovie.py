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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
This module contains the protocol for Summovie
"""

from os.path import exists
import pyworkflow.protocol.params as params
import pyworkflow.protocol.constants as cons
import pyworkflow.utils.path as pwutils
from grigoriefflab import SUMMOVIE_PATH
from pyworkflow.em.protocol import ProtAlignMovies
from convert import writeShiftsMovieAlignment


class ProtSummovie(ProtAlignMovies):
    """Summovie is used to align the frames of movies recorded
    on an electron microscope to reduce image blurring due
    to beam-induced motion. It reads stacks of movies that
    are stored in MRC/CCP4 format. Summovie generates frame
    sums that can be used in subsequent image processing
    steps and optionally applies an exposure-dependent
    filter to maximize the signal at all resolutions
    in the frame averages."""
    
    _label = 'Summovie'
    CONVERT_TO_MRC = 'mrc'
    doSaveAveMic = True
    
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineAlignmentParams(self, form):
        form.addHidden('binFactor', params.FloatParam, default=1.)

        group = form.addGroup('Average')
        line = group.addLine('Remove frames to SUM from',
                             help='How many frames you want remove to sum\n'
                                  'from beginning and/or from the end of each movie.')
        line.addParam('sumFrame0', params.IntParam, default=0, label='beginning')
        line.addParam('sumFrameN', params.IntParam, default=0, label='end')
        
        group.addParam('useAlignment', params.BooleanParam, default=True,
              label="Use movie alignment to Sum frames?")
        group.addParam('doApplyDoseFilter', params.BooleanParam, default=True,
                      label='Apply Dose filter',
                      help='Apply a dose-dependent filter to frames before summing them')
        group.addParam('exposurePerFrame', params.FloatParam,
                      label='Exposure per frame (e/A^2)',
                      help='Exposure per frame, in electrons per square Angstrom')

        group.addParam('doRestoreNoisePower', params.BooleanParam,
                      default=True,
                      label='Restore Noise Power? ',
                      help='Restore Noise Power? ',
                      expertLevel=cons.LEVEL_ADVANCED)
        
        form.addParallelSection(threads=1, mpi=0)
    
    #--------------------------- STEPS functions ----------------------------------------------
    def _processMovie(self, movie):
        self._createLink(movie)
        self._writeMovieAlignment(movie)
        self._argsSummovie()
        s0, sN = self._getFrameRange(self._getNumberOfFrames(movie), 'sum')
        
        params = {'movieFn' : self._getMovieFn(movie),
                  'numberOfFrames' : self._getNumberOfFrames(movie),
                  'micFn' : self._getMicFn(movie),
                  'initFrame' : s0,
                  'finalFrame' : sN,
                  'shiftsFn' : self._getShiftsFn(movie),
                  'samplingRate' : movie.getSamplingRate(),
                  'voltage' : movie.getAcquisition().getVoltage(),
                  'frcFn' : self._getFrcFn(movie),
                  'exposurePerFrame' : self.exposurePerFrame.get()
                  }
        
        if self.doApplyDoseFilter:
            params['doApplyDoseFilter'] = 'YES'
        else:
            params['doApplyDoseFilter'] = 'NO'
        
        if self.doRestoreNoisePower:
            params['doRestoreNoisePower'] = 'YES'
        else:
            params['doRestoreNoisePower'] = 'NO'
        
        try:
            self.runJob(self._program, self._args % params)
        except:
            print("ERROR: Movie %s failed\n" % movie.getFileName())
    
    #--------------------------- INFO functions --------------------------------------------
    def _citations(self):
        return []
    
    def _summary(self):
        summary = []
        movie = self.inputMovies.get().getFirstItem()
        s0, sN = self._getFrameRange(self._getNumberOfFrames(movie), 'sum')
        if movie.hasAlignment():
            fstFrame, lstFrame = movie.getAlignment().getRange()
            if self.useAlignment and (fstFrame > s0 or lstFrame < sN):
                summary.append("Warning!!! You have selected a frame range wider than"
                               " the range selected to align. All the frames selected"
                               " without alignment information, will be aligned by"
                               " setting alignment to 0")
        return summary
    
    def _methods(self):
        return []

    def _validate(self):
        errors=[]
        if not exists(SUMMOVIE_PATH):
            errors.append("Cannot find the Summovie program at: "+SUMMOVIE_PATH)
        return errors
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _argsSummovie(self):
        self._program = 'export OMP_NUM_THREADS=1; ' + SUMMOVIE_PATH
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
    
    def _getNumberOfFrames(self, movie):
        _, _, n = movie.getDim()
        return n
    
    def _writeMovieAlignment(self, movie):
        shiftsFn = self._getShiftsFn(movie)
        s0, sN = self._getFrameRange(self._getNumberOfFrames(movie), 'sum')
        
        if movie.hasAlignment() and self.useAlignment:
            writeShiftsMovieAlignment(movie, shiftsFn, s0, sN)
        else:
            f=open(shiftsFn,'w')
            frames = sN - s0 + 1
            shift= ("0 " * frames + "\n" ) *2
            f.write(shift)
            f.close()
    
    def _doGenerateOutputMovies(self):
        """ Returns True if an output set of movies will be generated.
        The most common case is to always generate output movies,
        either with alignment only or the binary aligned movie files.
        Subclasses can override this function to change this behavior.
        """
        return False
