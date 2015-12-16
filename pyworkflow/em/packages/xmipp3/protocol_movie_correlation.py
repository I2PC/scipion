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

import sys

import pyworkflow.em as em
import pyworkflow.protocol.params as params
import pyworkflow.protocol.constants as cons
from pyworkflow.em.protocol import ProtProcessMovies

from convert import parseMovieAlignment


# Alignment methods enum
AL_MOTIONCORR = 0
AL_CROSSCORRELATION = 1


class ProtMotionCorr(ProtProcessMovies):
    """
    Wrapper protocol to Dose Fractionation Tool: Flat fielding and Drift correction
    Wrote by Xueming Li @ Yifan Cheng Lab, UCSF   
    """
    
    _label = 'motioncorr alignment'
             
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        ProtProcessMovies._defineParams(self, form)
        
        group = form.addGroup('Common parameters')
        line = group.addLine('How many frames remove to align',
                            help='How many frames remove'
                                 ' from movie alignment.')
        line.addParam('alignFrame0', params.IntParam, default=0, label='from first')
        line.addParam('alignFrameN', params.IntParam, default=0, label='to final')
        
        line = group.addLine('Crop offsets (px)')
        line.addParam('cropOffsetX', params.IntParam, default=0, label='X')
        line.addParam('cropOffsetY', params.IntParam, default=0, label='Y')
        
        line = group.addLine('Crop dimensions (px)',
                      help='How many pixels to crop from offset\n'
                           'If equal to 0, use maximum size.')
        line.addParam('cropDimX', params.IntParam, default=0, label='X')
        line.addParam('cropDimY', params.IntParam, default=0, label='Y')
                
        form.addParam('doGPU', params.BooleanParam, default=False,
                      label="Use GPU (vs CPU)",
                      help="Set to true if you want the GPU implementation of Motioncorr")
        
        # GROUP GPU PARAMETERS
        group = form.addGroup('GPU', condition="doGPU")
        line = group.addLine('How many frames remove to sum',
                             condition="doGPU",
                             help='How many frames you want remove to sum\n'
                                  'from beginning and/or to the end of each movie.')
        line.addParam('sumFrame0', params.IntParam, default=0, label='from first')
        line.addParam('sumFrameN', params.IntParam, default=0, label='to final')

        group.addParam('GPUCore', params.IntParam, default=0, expertLevel=cons.LEVEL_ADVANCED,
                      label="Choose GPU core",
                      condition="doGPU",
                      help="GPU may have several cores. Set it to zero"
                           " if you do not know what we are talking about."
                           " First core index is 0, second 1 and so on.")
        group.addParam('binFactor', params.IntParam, default=1,
                       condition="doGPU", label='Binning factor',
                       help='1x or 2x. Bin stack before processing.')
        group.addParam('extraParams', params.StringParam, default='',
                       expertLevel=cons.LEVEL_ADVANCED, condition="doGPU",
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
        
        # GROUP CPU PARAMETERS
        group = form.addGroup('CPU', condition="not doGPU")
        group.addParam('filterFactor', params.FloatParam, default=4,
                       condition="not doGPU", label='Filter at (A)',
                       help="For the calculation of the shifts with Xmipp, micrographs are "
                            "filtered (and downsized accordingly) to this resolution. "
                            "Then shifts are calculated, and they are applied to the "
                            "original frames without any filtering and downsampling.")
        form.addParam('maxShift', params.IntParam, default=30, expertLevel=cons.LEVEL_ADVANCED,
                      label="Maximum shift (pixels)", condition="not doGPU",
                      help='Maximum allowed distance (in pixels) that each frame can be shifted'
                           'with respect to the next.')
        form.addParam('doSaveAveMic', params.BooleanParam, default=True,
                      label="Save aligned micrograph", expertLevel=cons.LEVEL_ADVANCED)
        
        form.addParam('doSaveMovie', params.BooleanParam, default=False,
                      label="Save movie", expertLevel=cons.LEVEL_ADVANCED,
                      help="Save Aligned movie")
        
        form.addParallelSection(threads=1, mpi=1)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def _processMovie(self, movieId, movieName, movieFolder):
        movieSet = self.inputMovies.get()
        micFn = self._getNameExt(movieName, '_aligned', 'mrc')
        alignedMovieFn = self._getCorrMovieName(movieId)
        
        if self.doGPU:
            logFile = self._getLogFile(movieId)
            args = {'-crx': self.cropOffsetX.get(),
                    '-cry': self.cropOffsetY.get(),
                    '-cdx': self.cropDimX.get(),
                    '-cdy': self.cropDimY.get(),
                    '-bin': self.binFactor.get(),
                    '-nst': self.alignFrame0.get(),
                    '-ned': self._getFinalFrame(movieName),
                    '-nss': self.sumFrame0.get(),
                    '-nes': self._getFinalFrame(movieName, "sum"),
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
                command += " -fct %s -ssc 1" (alignedMovieFn)

            command += ' ' + self.extraParams.get()
            
        else:
            from pyworkflow.em.packages.xmipp3 import getEnviron
            environ = getEnviron()
            
            program = 'xmipp_movie_alignment_correlation'
            
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
            
            command += ' --frameRange %d %d '%(firstFrame, self._getFinalFrame(movieName))
            command += ' --max_shift %d ' % self.maxShift
            if self.doSaveAveMic:
                command += ' --oavg %s' % micName
            if self.doSaveMovie:
                command += ' --oaligned %s' % alignedMovieFn
            if self.inputMovies.get().getDark():
                command += ' --dark ' + self.inputMovies.get().getDark()
                grayCorrected=True
            if self.inputMovies.get().getGain():
                command += ' --gain ' + self.inputMovies.get().getGain()
                grayCorrected=True
        
        try:
            self.runJob('dosefgpu_driftcorr', command, cwd=movieFolder, env=environ)
        except:
            print >> sys.stderr, program, " failed for movie %(movieName)s" % locals()


        
        
        
        
        
#         inputName = movieName
#         micName = self._getMicName(movieId)
#         logFile = self._getLogFile(movieId)
#         gainFile = self.inputMovies.get().getGain()
#         gpuId = self.gpuId.get()
# 
# # TODO Check better way to handle gain correction
# #         if gainFile is not None:
# #             # Apply the gain correction to flat the raw movie
# #             correctedName = movieName.replace('.mrc', '_corrected.mrc')
# #             
# #             self.runJob('dosefgpu_flat', 
# #                         '%(inputName)s %(correctedName)s %(gainFile)s %(gpuId)s' % locals(),
# #                         cwd=movieFolder)
# #            
# #            inputName = correctedName
#         
#         args = {'-crx': self.cropOffsetX.get(),
#                 '-cry': self.cropOffsetY.get(),
#                 '-cdx': self.cropDimX.get(),
#                 '-cdy': self.cropDimY.get(),
#                 '-bin': self.binFactor.get(),
#                 '-nst': self.alignFrame0.get(),
#                 '-ned': self.alignFrameN.get(),
#                 '-nss': self.sumFrame0.get(),
#                 '-nes': self.sumFrameN.get(),
#                 '-gpu': gpuId,
#                 '-flg': logFile,
#                 }
#         
#         #TODO: check the gain can be handle in dosefgpu_driftcoor program
#         #if gainFile is not None:
#         #    args['-fgr'] = gainFile
#         
#         command = '%(inputName)s -fcs %(micName)s ' % locals()
#         command += ' '.join(['%s %s' % (k, v) for k, v in args.iteritems()])
#         command += ' ' + self.extraParams.get()
# 
#         self.runJob('dosefgpu_driftcorr', command, cwd=movieFolder)
#         # Move the micrograph and alignment text file
#         # before clean of movie folder
#         moveFile(join(movieFolder, micName), self._getExtraPath())
#         moveFile(join(movieFolder, logFile), self._getExtraPath())        

        
        
        
        
        
        
        
        
        
        
        
        
        
        
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
            alignment = parseMovieAlignment(logFile)
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
        return 'micrograph_%06d_Log.txt' % movieId
    
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
    
    
    
    
    
    
    
    
    
    
            
