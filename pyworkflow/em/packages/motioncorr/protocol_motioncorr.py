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
Protocol wrapper around the MotionCorr for movie alignment
"""

import os, sys

import pyworkflow.em as em
import pyworkflow.utils.path as putils
import pyworkflow.protocol.params as params
import pyworkflow.protocol.constants as cons
from pyworkflow.em.protocol import ProtProcessMovies

from convert import parseMovieAlignment


class ProtMotionCorr(ProtProcessMovies):
    """
    Wrapper protocol to Dose Fractionation Tool: Flat fielding and Drift correction
    Wrote by Xueming Li @ Yifan Cheng Lab, UCSF   
    """
    
    _label = 'motioncorr alignment'
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        ProtProcessMovies._defineParams(self, form)
        
        form.addParam('gpuMsg', params.LabelParam, default=True,
                      label='WARNING! You need to have installed CUDA'
                            ' libraries and a nvidia GPU')
        
        group = form.addGroup('Alignment')
        line = group.addLine('Remove frames to ALIGN from',
                            help='How many frames remove'
                                 ' from movie alignment.')
        line.addParam('alignFrame0', params.IntParam, default=0, label='beginning')
        line.addParam('alignFrameN', params.IntParam, default=0, label='end')
        line = group.addLine('Remove frames to SUM from',
                             help='How many frames you want remove to sum\n'
                                  'from beginning and/or from the end of each movie.')
        line.addParam('sumFrame0', params.IntParam, default=0, label='beginning')
        line.addParam('sumFrameN', params.IntParam, default=0, label='end')
        group.addParam('binFactor', params.IntParam, default=1,
                       label='Binning factor',
                       help='1x or 2x. Bin stack before processing.')
        
        line = group.addLine('Crop offsets (px)')
        line.addParam('cropOffsetX', params.IntParam, default=0, label='X')
        line.addParam('cropOffsetY', params.IntParam, default=0, label='Y')
        
        line = group.addLine('Crop dimensions (px)',
                             help='How many pixels to crop from offset\n'
                                  'If equal to 0, use maximum size.')
        line.addParam('cropDimX', params.IntParam, default=0, label='X')
        line.addParam('cropDimY', params.IntParam, default=0, label='Y')
        
        form.addParam('doSaveAveMic', params.BooleanParam, default=True,
                      label="Save aligned micrograph", expertLevel=cons.LEVEL_ADVANCED)
        
        form.addParam('doSaveMovie', params.BooleanParam, default=False,
                      label="Save movie", expertLevel=cons.LEVEL_ADVANCED,
                      help="Save Aligned movie")
        
        group = form.addGroup('GPU', expertLevel=cons.LEVEL_ADVANCED)
        
        group.addParam('GPUCore', params.IntParam, default=0,
                      label="Choose GPU core",
                      help="GPU may have several cores. Set it to zero"
                           " if you do not know what we are talking about."
                           " First core index is 0, second 1 and so on.")
        group.addParam('extraParams', params.StringParam, default='',
                       label='Additional parameters',
                       help="""
-bft       150               BFactor in pix^2.
-pbx       96                Box dimension for searching CC peak.
-fod       2                 Number of frame offset for frame comparison.
-nps       0                 Radius of noise peak.
-sub       0                 1: Save as sub-area corrected sum. 0: Not.
-srs       0                 1: Save uncorrected sum. 0: Not.
-scc       0                 1: Save CC Map. 0: Not.
-slg       1                 1: Save Log. 0: Not.
-atm       1                 1: Align to middle frame. 0: Not.
-dsp       1                 1: Save quick results. 0: Not.
-fsc       0                 1: Calculate and log FSC. 0: Not.
                            """)
        
        form.addParallelSection(threads=0, mpi=0)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def _processMovie(self, movieId, movieName, movieFolder):
        movieSet = self.inputMovies.get()
        self._getFinalFrame(movieSet, movieId)
        
        if not (movieName.endswith("mrc") or movieName.endswith("mrcs")):
            movieEm = os.path.join(movieFolder, movieName)
            movieMrc = os.path.join(movieFolder, self._getNameExt(movieName, '', 'mrc'))
            ih = em.ImageHandler()
            ih.convertStack(movieEm, movieMrc)
            movieName = os.path.basename(movieMrc)
        
        micFn = self._getNameExt(movieName, '_aligned', 'mrc')
        alignedMovieFn = self._getCorrMovieName(movieId)
        
        logFile = self._getLogFile(movieId)
        args = {'-crx': self.cropOffsetX.get(),
                '-cry': self.cropOffsetY.get(),
                '-cdx': self.cropDimX.get(),
                '-cdy': self.cropDimY.get(),
                '-bin': self.binFactor.get(),
                '-nst': self.alignFrame0.get(),
                '-ned': self.frameN,
                '-nss': self.sumFrame0.get(),
                '-nes': self.frameNSum,
                '-gpu': self.GPUCore.get(),
                '-flg': logFile,
                }
        
        command = '%s -fcs %s ' % (movieName, micFn)
        command += ' '.join(['%s %s' % (k, v) for k, v in args.iteritems()])
        
        if movieSet.getGain():
            command += " -fgr " + movieSet.getGain()
            grayCorrected=True
        
        if movieSet.getDark():
            command += " -fdr " + movieSet.getDark()
            grayCorrected=True
        
        if self.doSaveMovie:
            command += " -fct %s -ssc 1" % alignedMovieFn
        command += ' ' + self.extraParams.get()
        program = 'dosefgpu_driftcorr'
        try:
            self.runJob(program, command, cwd=movieFolder)
        except:
            print >> sys.stderr, program, " failed for movie %(movieName)s" % locals()
        
        putils.cleanPattern(os.path.join(movieFolder, movieName))
        putils.moveTree(self._getTmpPath(), self._getExtraPath())
    
    def createOutputStep(self):
        inputMovies = self.inputMovies.get()
        if self.doSaveMovie:
            suffix = "_aligned"
        else:
            suffix = "_original"
        movieSet = self._createSetOfMovies(suffix)
        movieSet.copyInfo(inputMovies)
        if self.binFactor != 1:
            movieSet.setSamplingRate(inputMovies.getSamplingRate()*self.binFactor)
        
        if self.doSaveAveMic:
            micSet = self._createSetOfMicrographs()
            micSet.copyInfo(inputMovies)
            if self.binFactor != 1:
                micSet.setSamplingRate(inputMovies.getSamplingRate()*self.binFactor)
        else:
            micSet = None
        
        for movie in inputMovies:
            self._createOutputMovie(movie, movieSet, micSet)
        
        self._defineOutputs(outputMovies=movieSet)
        self._defineTransformRelation(inputMovies, movieSet)
        
        if self.doSaveAveMic:
            self._defineOutputs(outputMicrographs=micSet)
            self._defineSourceRelation(inputMovies, micSet)
    
    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        return summary
    
    def _validate(self):
        errors = []
        if max(self.numberOfThreads, self.numberOfMpi) > 1:
            errors.append("GPU and Parallelization can not be used together")
        if not (self.binFactor == 1 or self.binFactor == 2):
            errors.append("Binning factor can only be 1 or 2")
            
        return errors
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _getLogFile(self, movieId):
        return 'micrograph_%06d_Log.txt' % movieId
    
    def _getFinalFrame(self, movieSet, movieId):
        if getattr(self, 'frameN', None) is None:
            movie = movieSet[movieId]
            movieName = movie.getFileName()
            
            ih = em.ImageHandler()
            _, _, z, n = ih.getDimensions(movieName)
            self.totalFrames = max(z, n) - 1
            frameN = self.alignFrameN.get()
            
            if frameN == 0:
                self.frameN = 0
            else:
                self.frameN = self.totalFrames - frameN
    
            frameNSum = self.sumFrameN.get()
            if frameNSum == 0:
                self.frameNSum = 0
            else:
                self.frameNSum = self.totalFrames - frameNSum
    
    def _getExtraMovieFolder(self, movieId):
        """ Create a Movie folder where to work with it. """
        return self._getExtraPath('movie_%06d' % movieId)
    
    def _createOutputMovie(self, movie, movieSet, micSet=None):
        movieId = movie.getObjId()
        
        movieFolder = self._getExtraMovieFolder(movieId)
        movieName = os.path.join(movieFolder, self._getCorrMovieName(movieId))
        micFn = os.path.join(movieFolder, self._getNameExt(movie.getFileName(),'_aligned', 'mrc'))
        
        # Parse the alignment parameters and store the log files
        alignedMovie = movie.clone()
        
        if self.doSaveMovie:
            alignedMovie.setFileName(movieName)
            if self.frameN == 0:
                totFrames = self.totalFrames
            else:
                totFrames = self.frameN
            
            diff = totFrames - self.alignFrame0.get() + 1
            
            alignment = em.MovieAlignment(first=self.alignFrame0.get()+1, 
                                          last=totFrames+1,
                                          shifts=[0, 0]*diff)
        else:
            alignment = parseMovieAlignment(os.path.join(movieFolder, 
                                                         self._getLogFile(movieId)))
        alignment.setRoi([self.cropOffsetX.get(), self.cropOffsetY.get(),
                          self.cropDimX.get(), self.cropDimY.get()])
        alignment.setScale(self.binFactor.get())
        
        alignedMovie.setAlignment(alignment)
        movieSet.append(alignedMovie)
        
        if self.doSaveAveMic:
            mic = micSet.ITEM_TYPE()
            mic.setObjId(movieId)
            mic.setFileName(micFn)
            micSet.append(mic)




