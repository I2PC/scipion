# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (jgomez@cnb.csic.es)
# *
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
Protocol wrapper around the xmipp correlation alignment only for movie average.
"""

import pyworkflow.protocol.params as params
from pyworkflow.em.protocol import ProtAlignMovies
import pyworkflow.em.metadata as md
from convert import getMovieFileName, writeShiftsMovieAlignment


class XmippProtMovieAverage(ProtAlignMovies):
    """
    Protocol to average movies
    """
    
    _label = 'movie average'
    CONVERT_TO_MRC = 'mrcs'
    doSaveAveMic = True
    
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineAlignmentParams(self, form):
        group = form.addGroup('Average')
        line = group.addLine('Remove frames to SUM from',
                             help='How many frames you want remove to sum\n'
                                  'from beginning and/or from the end of each movie.')
        line.addParam('sumFrame0', params.IntParam, default=0, label='beginning')
        line.addParam('sumFrameN', params.IntParam, default=0, label='end')
        group.addParam('binFactor', params.FloatParam, default=1,
                       label='Binning factor',
                       help='Binning factor, it may be any floating number '
                            'Binning in Fourier is the first operation, so that ' 
                            'crop parameters are referred to the binned images. ')
        
        group.addParam('cropRegion', params.EnumParam, 
                       choices=['None', 'from Alignment', 'New'],
                       label="Define crop region", default=0,
                       help="""Select if you want to crop the final micrograph.
                               If select:
                                     *None*, the final micrographs will have the same dimensions as
                                     the input movies. The region of interest, if the input movies 
                                     have alignment, is ignored.
                                     *from Alignment*, if the movies has alignment, the region of
                                     interest is defined and will be apply; else, micrographs
                                     will have the same dimensions as the input movies.
                                     *New*, All crop parameters should be defined below.
                                     The region of interest, if the input movies have alignment,
                                     is ignored.
                            """)
        
        line = group.addLine('Crop offsets (px)', condition='cropRegion==2')
        line.addParam('cropOffsetX', params.IntParam, default=0, label='X')
        line.addParam('cropOffsetY', params.IntParam, default=0, label='Y')
        
        line = group.addLine('Crop dimensions (px)', condition='cropRegion==2',
                             help='How many pixels to crop from offset\n'
                                  'If equal to 0, use maximum size.')
        line.addParam('cropDimX', params.IntParam, default=0, label='X')
        line.addParam('cropDimY', params.IntParam, default=0, label='Y')
        
        form.addParam('useAlignment', params.BooleanParam, default=True,
                      label="Use movie alignment to Sum frames?",
                      help="If set Yes, the alignment information (if"
                           " it exists) will take into account to align"
                           " your movies.")
        form.addParallelSection(threads=1, mpi=0)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def _processMovie(self, movie):
        inputMd = self._getMovieOrMd(movie)
        x, y, n = movie.getDim()
        args  = '-i %s ' % inputMd
        args += '--sampling %f ' % movie.getSamplingRate()
        args += '--useInputShifts '
        
        if self.binFactor > 1:
            args += '--bin %f ' % self.binFactor
        
        cropRegion = self.cropRegion.get()
        
        if cropRegion == 1 and movie.hasAlignment():
            roi = movie.getAlignment().getRoi()
            args += '--cropULCorner %d %d ' % (roi[0], roi[1])
            
            if roi[2] <= 0:
                dimX = x - 1
            else:
                dimX = roi[0] + roi[2] -1
            
            if roi[3] <= 0:
                dimY = y - 1
            else:
                dimY = roi[1] + roi[3] -1

        elif cropRegion == 2:
            
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
        
        if not (movie.hasAlignment() and self.useAlignment):
            s0, sN = self._getFrameRange(self._getNumberOfFrames(movie), 'sum')
            args += ' --frameRangeSum %d %d ' % (s0-1, sN-1)
        
        args += ' --oavg %s ' % self._getExtraPath(self._getOutputMicName(movie))
        
        if self.inputMovies.get().getDark() is not None:
            args += ' --dark ' + self.inputMovies.get().getDark()
        
        if self.inputMovies.get().getGain() is not None:
            args += ' --gain ' + self.inputMovies.get().getGain()
        
        self.runJob('xmipp_movie_alignment_correlation', args)
    
    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        if (self.cropDimX > 0 and self.cropDimY <= 0 or
            self.cropDimY > 0 and self.cropDimX <= 0):
            errors.append("If you give cropDimX, you should also give cropDimY "
                          "and viceversa")
        return errors
    
    def _summary(self):
        summary = []
        movie = self.inputMovies.get().getFirstItem()
        s0, sN = self._getFrameRange(self._getNumberOfFrames(movie), 'sum')
        
        if self.cropRegion.get() == 1 and not movie.hasAlignment():
            summary.append("Warning!!! You select *from Alignment* crop"
                           " region, but your movies have not alignment."
                           " Your resulting micrographs were not cropped."
                           " If you want to crop, please use *New* option.")
        
        if movie.hasAlignment():
            fstFrame, lstFrame = movie.getAlignment().getRange()
            if self.useAlignment and (fstFrame > s0 or lstFrame < sN):
                summary.append("Warning!!! You have selected a frame range wider than"
                               " the range selected to align. All the frames selected"
                               " without alignment information, will be aligned by"
                               " setting alignment to 0")
        
        return summary
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _getNumberOfFrames(self, movie):
        _, _, n = movie.getDim()
        return n
    
    def _getShiftsFile(self, movie):
        return self._getExtraPath(self._getMovieRoot(movie) + '_shifts.xmd')
    
    def _getMovieOrMd(self, movie):
        if movie.hasAlignment() and self.useAlignment:
            shiftsMd = self._getShiftsFile(movie)
            s0, sN = self._getFrameRange(self._getNumberOfFrames(movie), 'sum')
            writeShiftsMovieAlignment(movie, shiftsMd, s0, sN)
            return shiftsMd
        else:
            return getMovieFileName(movie)
    
    def _doGenerateOutputMovies(self):
        """ Returns True if an output set of movies will be generated.
        The most common case is to always generate output movies,
        either with alignment only or the binary aligned movie files.
        Subclasses can override this function to change this behavior.
        """
        return False
