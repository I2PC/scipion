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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
Protocol wrapper around the xmipp correlation alignment for movie alignment
"""

import sys

import pyworkflow.em as em
import pyworkflow.protocol.params as params
import pyworkflow.protocol.constants as cons
from pyworkflow.em.protocol import ProtProcessMovies

# from convert import parseMovieAlignment

class ProtMovieCorr(ProtProcessMovies):
    """
    Wrapper protocol to Xmipp Movie Alignment by cross-correlation
    """
    
    _label = 'xmipp movie alignment'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        ProtProcessMovies._defineParams(self, form)
        
        # Specific parameters
        form.addParam('splineOrder', params.IntParam, default=3, expertLevel=cons.LEVEL_ADVANCED,
                      label='B-spline order',
                      help="1 for linear interpolation (faster but lower quality), 3 for cubic interpolation (slower but more accurate).")
        
        
        form.addParam('filterFactor', params.FloatParam, default=4,
                       condition="not doGPU", label='Filter at (A)',
                       help="For the calculation of the shifts with Xmipp, micrographs are "
                            "filtered (and downsized accordingly) to this resolution. "
                            "Then shifts are calculated, and they are applied to the "
                            "original frames without any filtering and downsampling.")
        form.addParam('maxShift', params.IntParam, default=30, expertLevel=cons.LEVEL_ADVANCED,
                      label="Maximum shift (pixels)", condition="not doGPU",
                      help='Maximum allowed distance (in pixels) that each frame can be shifted'
                           'with respect to the next.')
        
        form.addParallelSection(threads=1, mpi=1)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def _processMovie(self, movieId, movieName, movieFolder):
        movieSet = self.inputMovies.get()
        micFn = self._getNameExt(movieName, '_aligned', 'mrc')
        alignedMovieFn = self._getCorrMovieName(movieId)
        
        mdCorrelation = self._getNameExt(movieName, '_alignedCorrelation.xmd')
        mdName = self._getNameExt(movieName, '_aligned', 'xmd')
        
        # Some movie have .mrc or .mrcs format but it is recognized as a volume
        if movieName.endswith('.mrcs') or movieName.endswith('.mrc'):
            movieSuffix = ':mrcs'
        else:
            movieSuffix = ''
        
        command  = '-i %s%s ' % (movieName, movieSuffix)
        command += '-o %s ' % mdCorrelation
        command += '--sampling %f ' % self.samplingRate
        command += '--max_freq %f ' % self.maxFreq
        command += '--cropULCorner %d %d ' % (self.cropOffsetX.get(), self.cropOffsetY.get())
        command += '--cropDRCorner %d %d ' % (self.cropOffsetX.get() + self.cropDimX.get() -1
                                             ,self.cropOffsetY.get() + self.cropDimY.get() -1)
        
        command += ' --frameRange %d %d '%(self.alignFrame0.get(), self._getFinalFrame(movieName,"align"))
        command += ' --frameRangeSum %d %d '%(self.sumFrame0.get(), self._getFinalFrame(movieName,"sum"))
        command += ' --max_shift %d ' % self.maxShift
        if self.doSaveAveMic:
            command += ' --oavg %s' % micName
        if self.doSaveMovie:
            command += ' --oaligned %s' % alignedMovieFn
        if self.inputMovies.get().getDark():
            command += ' --dark ' + self.inputMovies.get().getDark()
        if self.inputMovies.get().getGain():
            command += ' --gain ' + self.inputMovies.get().getGain()
        self.runJob('xmipp_movie_alignment_correlation',command, numberOfMpi=1)
        
        
    def createOutputStep(self):
        inputMovies = self.inputMovies.get()
        micSet = self._createSetOfMicrographs()
        micSet.copyInfo(inputMovies)
        # Also create a Set of Movies with the alignment parameters
        movieSet = self._createSetOfMovies()
        movieSet.copyInfo(inputMovies)
        
        for movie in inputMovies:
            movieId = movie.getObjId()
            micName = self._getMicName(movieId)
            movieFolder = self._getMovieFolder(movieId)
          
            mic = micSet.ITEM_TYPE()
            mic.setObjId(movieId)
            mic.setFileName(self._getExtraPath(micName))
            micSet.append(mic)
            
            # Parse the alignment parameters and store the log files
            alignedMovie = movie.clone()
            logFile = self._getExtraPath(self._getLogFile(movieId))
#             alignment = parseMovieAlignment(logFile)
            alignedMovie.setAlignment(alignment)
            movieSet.append(alignedMovie)
            
        self._defineOutputs(outputMicrographs=micSet)
        self._defineTransformRelation(inputMovies, micSet)
        
        self._defineOutputs(outputMovies=movieSet)
        self._defineTransformRelation(inputMovies, movieSet)
    
    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        return summary
    
    def _validate(self):
        errors = []
        if max(self.numberOfThreads, self.numberOfMpi) > 1 and self.doGPU:
            errors.append("GPU and Parallelization can not be used together")
            
        return errors
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _getLogFile(self, movieId):
        return 'micrograph_%06d_log.txt' % movieId
    
    def _getFinalFrame(self, movieName, met="align"):
        ih = em.ImageHandler()
        _, _, z, n = ih.getDimensions(movieName)
        totalFrames = max(z, n) - 1
        if met == "align":
            frameN = self.alignFrameN.get()
        else:
            frameN = self.sumFrameN.get()
        
        if frameN == 0:
            return 0
        else:
            return totalFrames - frameN
    
    
    
    
    
    
    
    
    
    
            
