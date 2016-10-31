# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Vahid Abrishami (vabrishami@cnb.csic.es)
# *              Josue Gomez Blanco (jgomez@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
import pyworkflow.protocol.params as params
import pyworkflow.protocol.constants as cons
from pyworkflow.em.protocol import ProtAlignMovies
import pyworkflow.em.metadata as md
from convert import writeMovieMd



class XmippProtMovieCorr(ProtAlignMovies):
    """
    Wrapper protocol to Xmipp Movie Alignment by cross-correlation
    """
    OUTSIDE_WRAP = 0
    OUTSIDE_AVG = 1
    OUTSIDE_VALUE = 2

    INTERP_LINEAR = 0
    INTERP_CUBIC = 1

    # Map to xmipp interpolation values in command line
    INTERP_MAP = {INTERP_LINEAR: 1, INTERP_CUBIC: 3}

    _label = 'correlation alignment'

    #--------------------------- DEFINE param functions ------------------------

    def _defineAlignmentParams(self, form):
        ProtAlignMovies._defineAlignmentParams(self, form)

        form.addParam('splineOrder', params.EnumParam,
                      default=self.INTERP_CUBIC, choices=['linear', 'cubic'],
                      expertLevel=cons.LEVEL_ADVANCED,
                      label='Interpolation',
                      help="linear (faster but lower quality), "
                           "cubic (slower but more accurate).")

        form.addParam('maxFreq', params.FloatParam, default=4,
                       label='Filter at (A)',
                       help="For the calculation of the shifts with Xmipp, "
                            "micrographs are filtered (and downsized "
                            "accordingly) to this resolution. Then shifts are "
                            "calculated, and they are applied to the original "
                            "frames without any filtering and downsampling.")

        form.addParam('maxShift', params.IntParam, default=30,
                      expertLevel=cons.LEVEL_ADVANCED,
                      label="Maximum shift (pixels)",
                      help='Maximum allowed distance (in pixels) that each '
                           'frame can be shifted with respect to the next.')
        
        form.addParam('outsideMode', params.EnumParam,
                      choices=['Wrapping','Average','Value'],
                      default=self.OUTSIDE_WRAP,
                      expertLevel=cons.LEVEL_ADVANCED,
                      label="How to fill borders",
                      help='How to fill the borders when shifting the frames')

        form.addParam('outsideValue', params.FloatParam, default=0.0,
                       expertLevel=cons.LEVEL_ADVANCED,
                       condition="outsideMode==2",
                       label='Fill value',
                       help="Fixed value for filling borders")

        form.addParallelSection(threads=1, mpi=1)
    
    #--------------------------- STEPS functions -------------------------------

    def _processMovie(self, movie):
        movieFolder = self._getOutputMovieFolder(movie)

        x, y, n = movie.getDim()
        a0, aN = self._getFrameRange(n, 'align')
        s0, sN = self._getFrameRange(n, 'sum')

        inputMd = os.path.join(movieFolder, 'input_movie.xmd')
        writeMovieMd(movie, inputMd, a0, aN, useAlignment=False)

        args  = '-i %s ' % inputMd
        args += '-o %s ' % self._getShiftsFile(movie)
        args += '--sampling %f ' % movie.getSamplingRate()
        args += '--max_freq %f ' % self.maxFreq
        args += '--Bspline %d ' % self.INTERP_MAP[self.splineOrder.get()]

        if self.binFactor > 1:
            args += '--bin %f ' % self.binFactor
        # Assume that if you provide one cropDim, you provide all
        
        offsetX = self.cropOffsetX.get()
        offsetY = self.cropOffsetY.get()
        cropDimX = self.cropDimX.get()
        cropDimY = self.cropDimY.get()
        
        args += '--cropULCorner %d %d ' % (offsetX, offsetY)
        
        if cropDimX <= 0:
            dimX = x - 1
        else:
            dimX = offsetX + cropDimX - 1
        
        if cropDimY <= 0:
            dimY = y - 1
        else:
            dimY = offsetY + cropDimY - 1
        
        args += '--cropDRCorner %d %d ' % (dimX, dimY)
        
        if self.outsideMode == self.OUTSIDE_WRAP:
            args += "--outside wrap"
        elif self.outsideMode == self.OUTSIDE_AVG:
            args += "--outside avg"
        elif self.outsideMode == self.OUTSIDE_AVG:
            args += "--outside value %f" % self.outsideValue
        
        args += ' --frameRange %d %d ' % (0, aN-a0)
        args += ' --frameRangeSum %d %d ' % (s0-a0, sN-s0)
        args += ' --max_shift %d ' % self.maxShift

        if self.doSaveAveMic:
            args += ' --oavg %s' % self._getExtraPath(self._getOutputMicName(movie))

        if self.doSaveMovie:
            args += ' --oaligned %s' % self._getExtraPath(self._getOutputMovieName(movie))

        if self.inputMovies.get().getDark():
            args += ' --dark ' + self.inputMovies.get().getDark()

        if self.inputMovies.get().getGain():
            args += ' --gain ' + self.inputMovies.get().getGain()

        self.runJob('xmipp_movie_alignment_correlation', args, numberOfMpi=1)


    #--------------------------- UTILS functions ------------------------------
    def _getShiftsFile(self, movie):
        return self._getExtraPath(self._getMovieRoot(movie) + '_shifts.xmd')

    def _getMovieShifts(self, movie):
        from convert import readShiftsMovieAlignment
        """ Returns the x and y shifts for the alignment of this movie.
         The shifts should refer to the original micrograph without any binning.
         In case of a bining greater than 1, the shifts should be scaled.
        """
        
        shiftsMd = md.MetaData(self._getShiftsFile(movie))
        return readShiftsMovieAlignment(shiftsMd)

    def _storeSummary(self, movie):
        if self.doSaveAveMic and movie.hasAlignment():
            s0, sN = self._getFrameRange(movie.getNumberOfFrames(), 'sum')
            fstFrame, lstFrame = movie.getAlignment().getRange()
            if fstFrame > s0 or lstFrame < sN:
                self.summaryVar.set("Warning!!! You have selected a frame range "
                                    "wider than the range selected to align. All "
                                    "the frames selected without alignment "
                                    "information, will be aligned by setting "
                                    "alignment to 0")

