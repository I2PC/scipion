# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (jgomez@cnb.csic.es)
# *              J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
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
from pyworkflow.em.protocol import ProtAlignMovies
from convert import writeMovieMd

CROP_NONE = 0
CROP_ALIGNMENT = 1
CROP_NEW = 2


class XmippProtMovieAverage(ProtAlignMovies):
    """
    Protocol to average movies
    """
    _label = 'movie average'
    CONVERT_TO_MRC = 'mrcs'
    doSaveAveMic = True

    INTERP_LINEAR = 0
    INTERP_CUBIC = 1

    # Map to xmipp interpolation values in command line
    INTERP_MAP = {INTERP_LINEAR: 1, INTERP_CUBIC: 3}
    
    #--------------------------- DEFINE param functions ------------------------
    def _defineAlignmentParams(self, form):
        group = form.addGroup('Average')

        line = group.addLine('Frames to SUM',
                             help='Frames range to SUM on each movie. The '
                                  'first frame is 1. If you set 0 in the final '
                                  'frame to sum, it means that you will sum '
                                  'until the last frame of the movie.')
        line.addParam('sumFrame0', params.IntParam, label='from')
        line.addParam('sumFrameN', params.IntParam, label='to')

        group.addParam('binFactor', params.FloatParam, default=1,
                       label='Binning factor',
                       help='Binning factor, it may be any floating number '
                            'Binning in Fourier is the first operation, so '
                            'that crop parameters are referred to the binned '
                            'images. ')
        
        group.addParam('cropRegion', params.EnumParam, 
                       choices=['None', 'From Alignment', 'New'],
                       label="Define crop region", default=CROP_NONE,
                       help="Select if you want to crop the final micrograph. "
                            "If you select: \n"
                            "*None*, the final micrographs will have the same "
                            "dimensions as the input movies. The region of "
                            "interest, if the input movies have alignment, is "
                            "ignored. \n"
                            "*from Alignment*, if the movies has alignment, "
                            "the region of interest is defined and will be "
                            "apply; else, micrographs will have the same "
                            "dimensions as the input movies. \n"
                            "*New*, All crop parameters should be defined below."
                            "The region of interest, if the input movies have "
                            "alignment, is ignored. ")
        
        line = group.addLine('Crop offsets (px)',
                             condition='cropRegion==%d' % CROP_NEW)
        line.addParam('cropOffsetX', params.IntParam, default=0, label='X')
        line.addParam('cropOffsetY', params.IntParam, default=0, label='Y')
        
        line = group.addLine('Crop dimensions (px)',
                             condition='cropRegion==%d' % CROP_NEW,
                             help='How many pixels to crop from offset\n'
                                  'If equal to 0, use maximum size.')
        line.addParam('cropDimX', params.IntParam, default=0, label='X')
        line.addParam('cropDimY', params.IntParam, default=0, label='Y')

        group.addParam('useAlignment', params.BooleanParam, default=True,
                       label="Use previous movie alignment to SUM frames?",
                       help="Input movies could have alignment information from"
                            "a previous protocol. If you select *Yes*, the "
                            "previous alignment will be taken into account.")

        form.addParam('splineOrder', params.EnumParam,
                      default=self.INTERP_CUBIC, choices=['linear', 'cubic'],
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Interpolation',
                      help="linear (faster but lower quality), "
                           "cubic (slower but more accurate).")

        form.addParallelSection(threads=1, mpi=0)
    
    #--------------------------- STEPS functions -------------------------------
    def _processMovie(self, movie):
        movieFolder = self._getOutputMovieFolder(movie)

        x, y, n = movie.getDim()
        s0, sN = self._getFrameRange(n, 'sum')

        inputMd = os.path.join(movieFolder, 'input_movie.xmd')
        writeMovieMd(movie, inputMd, s0, sN,
                     useAlignment=(movie.hasAlignment() and self.useAlignment))

        outputMicFn = self._getExtraPath(self._getOutputMicName(movie))
        
        if self.cropRegion == CROP_ALIGNMENT and movie.hasAlignment():
            roi = movie.getAlignment().getRoi()
        elif self.cropRegion == CROP_NEW:
            roi = [self.cropOffsetX.get(), self.cropOffsetY.get(),
                   self.cropDimX.get(), self.cropDimY.get()]
        else:
            roi = None

        self.averageMovie(movie, inputMd, outputMicFn, self.binFactor.get(),
                          roi, self.inputMovies.get().getDark(),
                          self.inputMovies.get().getGain(),
                          splineOrder=self.INTERP_MAP[self.splineOrder.get()])

        self._storeSummary(movie)
    
    #--------------------------- INFO functions --------------------------------
    # def _validate(self):
    #     errors = []
    #     if (self.cropDimX > 0 and self.cropDimY <= 0 or
    #         self.cropDimY > 0 and self.cropDimX <= 0):
    #         errors.append("If you give cropDimX, you should also give cropDimY "
    #                       "and viceversa")
    #     return errors
    
    #--------------------------- UTILS functions -------------------------------
    def _getShiftsFile(self, movie):
        return self._getExtraPath(self._getMovieRoot(movie) + '_shifts.xmd')
    
    def _doGenerateOutputMovies(self):
        """ Returns True if an output set of movies will be generated.
        The most common case is to always generate output movies,
        either with alignment only or the binary aligned movie files.
        Subclasses can override this function to change this behavior.
        """
        return False
    
    def _storeSummary(self, movie):
        if movie.hasAlignment():
            s0, sN = self._getFrameRange(movie.getNumberOfFrames(), 'sum')
            fstFrame, lstFrame = movie.getAlignment().getRange()
            if self.useAlignment and (fstFrame > s0 or lstFrame < sN):
                self.summaryVar.set("Warning!!! You have selected a frame range "
                                    "wider than the range selected to align. All "
                                    "the frames selected without alignment "
                                    "information, will be aligned by setting "
                                    "alignment to 0")
        else:
            if self.cropRegion == CROP_ALIGNMENT:
                self.summaryVar.set("Warning!!! You select *from Alignment* crop "
                                    "region, but your movies have not alignment."
                                    "Your resulting micrographs were not cropped."
                                    "If you want to crop, please use *New* option.")
                                    


