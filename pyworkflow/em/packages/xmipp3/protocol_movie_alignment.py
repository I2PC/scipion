# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Vahid Abrishami (vabrisahmi@cnb.csic.es)
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
In this module are protocol base classes related to EM Micrographs
"""

import sys
from os.path import join

from pyworkflow.object import String, Integer
from pyworkflow.protocol.params import IntParam, FloatParam, StringParam, BooleanParam, LEVEL_ADVANCED, LEVEL_ADVANCED, EnumParam
from pyworkflow.utils.path import moveFile
import pyworkflow.em as em
from pyworkflow.em.protocol import ProtProcessMovies
from pyworkflow.gui.plotter import Plotter
from PIL import Image
import matplotlib.pyplot as plt
import numpy as np

# Alignment methods enum
AL_OPTICAL = 0
AL_DOSEFGPU = 1
AL_DOSEFGPUOPTICAL = 2
AL_AVERAGE = 3
AL_CROSSCORRELATION = 4
AL_CROSSCORRELATIONOPTICAL = 5

class ProtMovieAlignment(ProtProcessMovies):
    """ Aligns movies, from direct detectors cameras, into micrographs.
    """
    _label = 'movie alignment'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        ProtProcessMovies._defineParams(self, form)

        form.addParam('alignMethod', EnumParam, choices=['optical flow'
                                                       , 'dosefgpu'
                                                       , 'dosefgpu + optical flow'
                                                       , 'average'
                                                       , 'croscorrelation'
                                                       , 'croscorrelation + optical flow'
                                                        ],
                      label="Alignment method", default=AL_OPTICAL,
                      display=EnumParam.DISPLAY_COMBO,
                      help='Method to use for movie alignment. '
                           'dosefgpu requires a GPU. Croscorrelation '
                           'and dosefgpu with default parameters are equivalent. '
                           'Croscorrelation is a CPU implementation of dosefgpu  '
                           'with limited functionality (only aligns whole frames)'
                      )

        # GROUP COMMON PARAMETERS
        group = form.addGroup('Common parameters')
        
        line = group.addLine('Used in alignment',
                            help='First and last frames used in alignment.\n'
                                  'The first frame in the stack is *0*.' )
        line.addParam('alignFrame0', IntParam, default=0, label='First')
        line.addParam('alignFrameN', IntParam, default=0, label='Last',
                      help='If *0*, use maximum value')
        
        # GROUP GPU PARAMETERS
        group = form.addGroup('GPU',condition="alignMethod==%d or (alignMethod==%d and expertLevel==%d)"
                                                         " or (alignMethod==%d and expertLevel==%d)"
                                                         % (AL_OPTICAL, AL_DOSEFGPUOPTICAL, LEVEL_ADVANCED, AL_DOSEFGPU, LEVEL_ADVANCED))
        group.addParam('doGPU', BooleanParam, default=False,
                      label="Use GPU (vs CPU)",
                      condition="alignMethod==%d or alignMethod==%d" % (AL_OPTICAL, AL_DOSEFGPUOPTICAL),
                      help="Set to true if you want the GPU implementation of Optical Flow")
        group.addParam('GPUCore', IntParam, default=0, expertLevel=LEVEL_ADVANCED,
                      label="Choose GPU core",
                      condition="doGPU  or alignMethod==%d or alignMethod==%d  " % (AL_DOSEFGPU, AL_DOSEFGPUOPTICAL),
                      help="GPU may have several cores. Set it to zero if you do not know what we are talking about. First core index is 0, second 1 and so on.")
        
        # GROUP OPTICAL FLOW PARAMETERS
        group = form.addGroup('Optical Flow parameters', expertLevel=LEVEL_ADVANCED, condition="alignMethod==%d or alignMethod==%d or alignMethod==%d " % (AL_OPTICAL, AL_DOSEFGPUOPTICAL, AL_CROSSCORRELATIONOPTICAL))
        
        group.addParam('winSize', IntParam, default=150,
                      label="Window size", expertLevel=LEVEL_ADVANCED,
                      help="Window size (shifts are assumed to be constant within this window).")

        group.addParam('groupSize', IntParam, default=1,
                      label="Group Size", expertLevel=LEVEL_ADVANCED,
                      help="In cases with low SNR, the average of a number of frames can be used in alignment")

        group.addParam('doSaveMovie', BooleanParam, default=False,
                      label="Save movie", expertLevel=LEVEL_ADVANCED,
                      help="Save Aligned movie")

        #---------------------------------- DosefGPU Params--------------------------------
        # GROUP DOSEFGPU PARAMETERS
        group = form.addGroup('DosefGPU/Croscorrelation parameters',condition="alignMethod==%d "
                                                              "or alignMethod==%d"
                                                              "or alignMethod==%d"
                                                              "or alignMethod==%d" %
                                                              (AL_DOSEFGPU,\
                                                               AL_DOSEFGPUOPTICAL,\
                                                               AL_CROSSCORRELATION,\
                                                               AL_CROSSCORRELATIONOPTICAL)
                              )
        
        line = group.addLine('Used in final sum',
                             help='First and last frames used in alignment.\n'
                                  'The first frame in the stack is *0*.' )
        line.addParam('sumFrame0', IntParam, default=0, label='First')
        line.addParam('sumFrameN', IntParam, default=0, label='Last',
                      help='If *0*, use maximum value')
        
        line = group.addLine('Crop offsets (px)')
        line.addParam('cropOffsetX', IntParam, default=0, label='X')
        line.addParam('cropOffsetY', IntParam, default=0, label='Y')
        
        line = group.addLine('Crop dimensions (px)',
                      help='How many pixels to crop from offset\n'
                           'If equal to 0, use maximum size.')
        line.addParam('cropDimX', IntParam, default=0, label='X')
        line.addParam('cropDimY', IntParam, default=0, label='Y')

        group.addParam('binFactor', IntParam, default=1,condition="alignMethod==%d "
                                                              "or alignMethod==%d" %
                                                                  (AL_DOSEFGPU,\
                                                                   AL_DOSEFGPUOPTICAL\
                       ),
                       label='Binning factor',
                       help='1x or 2x. Bin stack before processing.')

        group.addParam('filterFactor', FloatParam, default=4,condition="alignMethod==%d "
                                                              "or alignMethod==%d" %
                                                                  (AL_CROSSCORRELATION,\
                                                                   AL_CROSSCORRELATIONOPTICAL\
                       ),
                       label='Filter at (A)',
                       help='1x or 2x. Bin stack before processing.')


        group.addParam('extraParams', StringParam, default='',
                      expertLevel=LEVEL_ADVANCED,
                      label='Additional parameters',
                      help="""
-bft       150               BFactor in pix^2.
-pbx       96                Box dimension for searching CC peak.
-fod       2                 Number of frame offset for frame comparison.
-nps       0                 Radius of noise peak.
-sub       0                 1: Save as sub-area corrected sum. 0: Not.
-srs       0                 1: Save uncorrected sum. 0: Not.
-ssc       0                 1: Save aligned stack. 0: Not.
-scc       0                 1: Save CC Map. 0: Not.
-slg       1                 1: Save Log. 0: Not.
-atm       1                 1: Align to middle frame. 0: Not.
-dsp       1                 1: Save quick results. 0: Not.
-fsc       0                 1: Calculate and log FSC. 0: Not.
                      """)
        form.addParallelSection(threads=1, mpi=1)


    #--------------------------- STEPS functions ---------------------------------------------------
    def createOutputStep(self):
        inputMovies = self.inputMovies.get()
        micSet = self._createSetOfMicrographs()
        micSet.copyInfo(inputMovies)
        # Also create a Set of Movies with the alignment parameters
        if self.doSaveMovie:
            movieSet = self._createSetOfMovies()
            movieSet.copyInfo(inputMovies)
            movieSet.cropOffsetX = Integer(self.cropOffsetX)
            movieSet.cropOffsetY = Integer(self.cropOffsetY)
            movieSet.cropDimX = Integer(self.cropDimX)
            movieSet.cropDimY = Integer(self.cropDimY)
            movieSet.sumFrame0 = Integer(self.sumFrame0)
            movieSet.sumFrameN = Integer(self.sumFrameN)
            
        alMethod = self.alignMethod.get()
        for movie in self.inputMovies.get():
            micName = self._getNameExt(movie.getFileName(),'_aligned', 'mrc')
            metadataName = self._getNameExt(movie.getFileName(), '_aligned', 'xmd')
            plotCartName = self._getNameExt(movie.getFileName(), '_plot_cart', 'png')
            psdCorrName = self._getNameExt(movie.getFileName(),'_aligned_corrected', 'psd')
            # Parse the alignment parameters and store the log files
            alignedMovie = movie.clone()
            
            if self.run:
                alignedMovie.setFileName(self._getExtraPath(self._getNameExt(movie.getFileName(),'_aligned', 'mrcs')))
            ####>>>This is wrong. Save an xmipp metadata
            alignedMovie.alignMetaData = String(self._getExtraPath(metadataName))
            alignedMovie.plotCart = self._getExtraPath(plotCartName)
            alignedMovie.psdCorr = self._getExtraPath(psdCorrName)
            
            '''if (alMethod == AL_OPTICAL or
                alMethod == AL_DOSEFGPUOPTICAL or 
                alMethod == AL_CROSSCORRELATIONOPTICAL):
                movieCreatePlot(alignedMovie, True)'''
            
            if self.doSaveMovie:
                movieSet.append(alignedMovie)

            mic = em.Micrograph()
            # All micrograph are copied to the 'extra' folder after each step
            mic.setFileName(self._getExtraPath(micName))

# The micName of a micrograph MUST be the same as the original movie
#             mic.setMicName(micName)
            mic.setMicName(movie.getMicName())
            
            if (alMethod == AL_OPTICAL or 
                alMethod == AL_DOSEFGPUOPTICAL or 
                alMethod == AL_CROSSCORRELATIONOPTICAL):

                mic.alignMetaData = String(self._getExtraPath(metadataName))
                mic.plotCart = self._getExtraPath(plotCartName)
                movieCreatePlot(mic, True)
                mic.plotCart = em.Image()
                mic.plotCart.setFileName(self._getExtraPath(plotCartName))
            #if alMethod != AL_DOSEFGPU and alMethod != AL_CROSSCORRELATION:
            mic.psdCorr = em.Image()
            mic.psdCorr.setFileName(self._getExtraPath(psdCorrName))
            micSet.append(mic)


            # TODO: Methods for dosefgpu should be transferred to here
            """
            if alMethod == AL_DOSEFGPU:
                # Parse the alignment parameters and store the log files
                alignedMovie = movie.clone()
                logFile = self._getExtraPath(self._getLogFile(movie.getObjId()))
                import pyworkflow.em.packages.dosefgpu as dosefgpu
                alignment = dosefgpu.parseMovieAlignment(logFile)
                alignedMovie.setAlignment(alignment)
                movieSet.append(alignedMovie)
            """
        self._defineOutputs(outputMicrographs=micSet)
        self._defineSourceRelation(self.inputMovies, micSet)
        if self.doSaveMovie:
            self._defineOutputs(outputMovies=movieSet)
        """
        if alMethod == AL_DOSEFGPU:
            self._defineTransformRelation(inputMovies, micSet)
            self._defineOutputs(outputMovies=movieSet)
            self._defineTransformRelation(inputMovies, movieSet)
        """


    #--------------------------- UTILS functions ---------------------------------------------------
    #TODO: In methods calling 2 protocols we should:
    #      1) work with the original movie and not the resampled one
    #      2) output metadata with shifts should be given over
    #         orignal movie not intermediate one
    def _processMovie(self, movieId, movieName, movieFolder):
        """ Process the movie actions, remember to:
        1) Generate all output files inside movieFolder (usually with cwd in runJob)
        2) Copy the important result files after processing (movieFolder will be deleted!!!)
        """

        program = self._getProgram()

        # Read the parameters
        #micName = self._getMicName(movieId)
        micName = self._getNameExt(movieName, '_aligned', 'mrc')
        metadataNameInterMediate = self._getNameExt(movieName, '_alignedIntermediate', 'xmd')
        metadataName = self._getNameExt(movieName, '_aligned', 'xmd')
        psdCorrName = self._getNameExt(movieName,'_aligned_corrected', 'psd')

        firstFrame = self.alignFrame0.get()
        lastFrame = self.alignFrameN.get()
        gpuId = self.GPUCore.get()
        alMethod = self.alignMethod.get()
        doSaveMovie = False
        
        # Some movie have .mrc or .mrcs format but it is recognized as a volume
        if movieName.endswith('.mrcs') or movieName.endswith('.mrc'):
            movieSuffix = ':mrcs'
        else:
            movieSuffix = ''
            
        # For simple average execution
        grayCorrected = False
        if alMethod == AL_AVERAGE:
            command = '-i %(movieName)s%(movieSuffix)s -o %(micName)s' % locals()
            command += ' --nst %d --ned %d --simpleAverage' % (firstFrame, lastFrame)

            if self.inputMovies.get().getDark():
                command += " --dark "+self.inputMovies.get().getDark()
                grayCorrected=True
            if self.inputMovies.get().getGain():
                command += " --gain "+self.inputMovies.get().getGain()
                grayCorrected=True
            try:
                self.runJob(program, command, cwd=movieFolder)
            except:
                print >> sys.stderr, program, " failed for movie %(movieName)s" % locals()

        # For DosefGPU Execution (and combination with optical flow)
        if alMethod == AL_DOSEFGPU or alMethod == AL_DOSEFGPUOPTICAL:
            logFile = self._getLogFile(movieId)
            #gainFile = self.inputMovies.get().getGain()
            args = {'-crx': self.cropOffsetX.get(),
                    '-cry': self.cropOffsetY.get(),
                    '-cdx': self.cropDimX.get(),
                    '-cdy': self.cropDimY.get(),
                    '-bin': self.binFactor.get(),
                    '-nst': self.alignFrame0.get(),
                    '-ned': self.alignFrameN.get(),
                    '-nss': self.sumFrame0.get(),
                    '-nes': self.sumFrameN.get(),
                    '-gpu': gpuId,
                    '-flg': logFile,
                    }
            command = '%(movieName)s -fcs %(micName)s ' % locals()
            command += ' '.join(['%s %s' % (k, v) for k, v in args.iteritems()])
            if alMethod == AL_DOSEFGPUOPTICAL:
                program = 'dosefgpu_driftcorr'
                corrMovieName = self._getCorrMovieName(movieId)
                command += ' ' + '-fct %(corrMovieName)s -ssc 1 ' % locals()
            command += ' ' + self.extraParams.get()
            import pyworkflow.em.packages.dosefgpu as dosefgpu
            try:
                self.runJob(program, command, cwd=movieFolder,
                            env=dosefgpu.getEnviron())
            except:
                print >> sys.stderr, program, " failed for movie %(movieName)s" % locals()
        
        if alMethod == AL_CROSSCORRELATION or alMethod == AL_CROSSCORRELATIONOPTICAL: #not dosefgpu
            program = 'xmipp_movie_alignment_correlation'
            corrMovieName = self._getCorrMovieName(movieId)
            command  = '-i %s%s ' % (movieName, movieSuffix)
            command += '-o %s '% metadataNameInterMediate
            command += '--sampling %f ' % self.samplingRate
            command += '--max_freq %f ' % self.filterFactor
            command += '--cropULCorner %d %d '%(self.cropOffsetX.get(),self.cropOffsetY.get())
            command += '--cropDRCorner %d %d '%(self.cropOffsetX.get() + self.cropDimX.get() -1
                                               ,self.cropOffsetY.get() + self.cropDimY.get() -1)
            _lastFrame = -1
            if lastFrame != 0:
                _lastFrame = lastFrame
            command += '--frameRange %d %d '%(firstFrame,_lastFrame)
            command += '--max_shift %d ' % 16#TODO expert param
            command += '--oavg %s ' % micName
            command += ' --oaligned %s ' % corrMovieName
            if self.inputMovies.get().getDark():
                command += " --dark "+self.inputMovies.get().getDark()
                grayCorrected=True
            if self.inputMovies.get().getGain():
                command += " --gain "+self.inputMovies.get().getGain()
                grayCorrected=True
            
            try:
                self.runJob(program, command, cwd=movieFolder)
            except:
                print >> sys.stderr, program, " failed for movie %(movieName)s" % locals()

        # For Optical Flow execution (and combination with DosefGPU)
        if alMethod == AL_OPTICAL or\
           alMethod == AL_DOSEFGPUOPTICAL or\
           alMethod == AL_CROSSCORRELATIONOPTICAL:
            winSize = self.winSize.get()
            doSaveMovie = self.doSaveMovie.get()
            groupSize = self.groupSize.get()
            if alMethod == AL_DOSEFGPUOPTICAL:
                program = 'xmipp_movie_optical_alignment_gpu'
                corrMovieName = self._getCorrMovieName(movieId)
                command = '-i %(corrMovieName)s%(movieSuffix)s  ' % locals()
                # Set to Zero for Optical Flow (output movie of dosefgpu)
                firstFrame = 0
                lastFrame = 0
            elif alMethod == AL_CROSSCORRELATIONOPTICAL:
                program = 'xmipp_movie_optical_alignment_cpu'
                command = '-i %(corrMovieName)s ' % locals()
            else:
                command = '-i %(movieName)s%(movieSuffix)s ' % locals()
                if self.doGPU:
                    program = 'xmipp_movie_optical_alignment_gpu'
                    command += '--gpu %d ' % gpuId
                else:
                    program = 'xmipp_movie_optical_alignment_cpu'
            # Set to Zero for Optical Flow (output movie of dosefgpu)
            firstFrame = 0
            lastFrame = 0
            if doSaveMovie:
                command += '--ssc '
            command += '-o %(micName)s --winSize %(winSize)d --groupSize %(groupSize)d ' % locals()
            command += '--nst %d --ned %d ' % (firstFrame, lastFrame)
            if self.inputMovies.get().getDark() and not grayCorrected:
                command += " --dark "+self.inputMovies.get().getDark()
                grayCorrected=True
            if self.inputMovies.get().getGain() and not grayCorrected:
                command += " --gain "+self.inputMovies.get().getGain()
                grayCorrected=True
            if doSaveMovie:
                command += '--ssc '
            try:
                self.runJob(program, command, cwd=movieFolder)
            except:
                print >> sys.stderr, program, " failed for movie %(movieName)s" % locals()
            moveFile(join(movieFolder, metadataName), self._getExtraPath())

        # Compute half-half PSD
        ih = em.ImageHandler()
        print join(movieFolder, '%(movieName)s' % locals())
        avg = ih.computeAverage(join(movieFolder, movieName))
        avg.write(join(movieFolder, 'uncorrectedmic.mrc'))
        command = '--micrograph uncorrectedmic.mrc --oroot uncorrectedpsd --dont_estimate_ctf --pieceDim 400 --overlap 0.7'
        program = 'xmipp_ctf_estimate_from_micrograph'
        self.runJob(program, command, cwd=movieFolder)
        command = '--micrograph %(micName)s --oroot correctedpsd --dont_estimate_ctf --pieceDim 400 --overlap 0.7' % locals()
        self.runJob(program, command, cwd=movieFolder)
        correctedPSD = em.ImageHandler().createImage()
        unCorrectedPSD = em.ImageHandler().createImage()
        correctedPSD.read(join(movieFolder, 'correctedpsd.psd'))
        unCorrectedPSD.read(join(movieFolder, 'uncorrectedpsd.psd'))
        x, y, z, n = correctedPSD.getDimensions()
        for i in range(1,y):
            for j in range(1,x//2):
                unCorrectedPSD.setPixel(i, j, correctedPSD.getPixel(i,j))
        unCorrectedPSD.write(join(movieFolder, psdCorrName))

        # Move output micrograph and related information to 'extra' folder
        moveFile(join(movieFolder, micName), self._getExtraPath())
        if doSaveMovie:
            outMovieName = self._getNameExt(movieName,'_aligned', 'mrcs')
            moveFile(join(movieFolder, outMovieName), self._getExtraPath())

        if alMethod == AL_DOSEFGPU:
            # Copy the log file to have shifts information
            moveFile(join(movieFolder, logFile), self._getExtraPath())
        elif alMethod == AL_CROSSCORRELATION:
            # Copy metadatafile otherwise it will be deleted
            #TODO: create a proper scipion object
            moveFile(join(movieFolder, metadataNameInterMediate), self._getExtraPath())
            #moveFile(join(movieFolder, corrMovieName), self._getExtraPath())
        moveFile(join(movieFolder, psdCorrName), self._getExtraPath())


    def _getProgram(self):
        alMethod = self.alignMethod.get()
        if alMethod == AL_AVERAGE:
            return 'xmipp_movie_optical_alignment_cpu'
        if alMethod == AL_OPTICAL:
            if self.doGPU:
		print 'doGPU is set to 1'
                return 'xmipp_movie_optical_alignment_gpu'
            else:
		print 'doGPU is set to 0'
                return 'xmipp_movie_optical_alignment_cpu'
        elif alMethod == AL_DOSEFGPU:
            return 'dosefgpu_driftcorr'

    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        numThreads = self.numberOfThreads;
        alMethod = self.alignMethod.get()
        if numThreads>1:
            if self.doGPU:
                errors.append("GPU and Parallelization can not be used together")

        if self.doGPU and (alMethod == AL_CROSSCORRELATION or \
                           alMethod == AL_CROSSCORRELATIONOPTICAL):
                self.doGPU.set(False)
                #errors.append("Crosscorrelation is not implemente in GPU")

        # FIXME: JMRT: GPU should only be disable for web and temporarly
        #alignMethod = self.alignMethod.get()
        #if alignMethod == 1 or alignMethod == 2:
        #    errors.append("GPU methods are not available at the moment.")
        return errors

    def _citations(self):
        alMethod = self.alignMethod.get()
        if alMethod == AL_OPTICAL:
            return ['Abrishami2015']
        if alMethod == AL_DOSEFGPU:
            return ['Li2013']
        if alMethod == AL_DOSEFGPUOPTICAL:
            return ['Abrishami2015', 'Li2013']

    def _methods(self):
        methods = []
        alMethod = self.alignMethod.get()
        gpuId = self.GPUCore.get()
        if alMethod == AL_AVERAGE:
            methods.append('Aligning method: Simple average')
        if alMethod == AL_DOSEFGPU or alMethod == AL_DOSEFGPUOPTICAL:
            methods.append('Aligning method: DosefGPU')

        if alMethod == AL_OPTICAL or alMethod == AL_DOSEFGPUOPTICAL:
            methods.append('Aligning method: Optical Flow')
            methods.append('- Used a window size of: *%d*' % self.winSize.get())
            methods.append('- Used a pyramid size of: *6*')
            if self.doGPU:
                methods.append('- Used GPU *%d* for processing' % gpuId)

        return methods

    def _summary(self):
        firstFrame = self.alignFrame0.get()
        lastFrame = self.alignFrameN.get()
        summary = []
        if self.inputMovies.get():
            summary.append('Number of input movies: *%d*' % self.inputMovies.get().getSize())
        if lastFrame == 0:
            summary.append('Frames used in alignment: *%d* to *%s* (first frame is 0)' % (firstFrame, 'Last Frame'))
        else:
            summary.append('Frames used in alignment: *%d* to *%d* (first frame is 0)' % (firstFrame, lastFrame))

        return summary

def createPlots(plotType, protocol, micId):
    print "output Micrographs to create Plot %s" % protocol.outputMicrographs
    mic = protocol.outputMicrographs[micId]
    return movieCreatePlot(mic, False).show()

def movieCreatePlot(mic, saveFig):

    import xmipp
    meanX = []
    meanY = []
    figureSize = (8, 6)

    #alignedMovie = mic.alignMetaData
    md = xmipp.MetaData(mic.alignMetaData)
    plotter = Plotter(*figureSize)
    figure = plotter.getFigure()

    preX = 0.0
    preY = 0.0
    meanX.append(0.0)
    meanY.append(0.0)
    ax = figure.add_subplot(111)
    ax.grid()
    ax.set_title('Cartesian representation')
    ax.set_xlabel('Drift x (pixels)')
    ax.set_ylabel('Drift y (pixels)')
    ax.plot(0, 0, 'yo-')
    for objId in md:
        preX += md.getValue(xmipp.MDL_OPTICALFLOW_MEANX, objId)
        preY += md.getValue(xmipp.MDL_OPTICALFLOW_MEANY, objId)
        meanX.append(preX)
        meanY.append(preY)
        ax.plot(preX, preY, 'yo-')
        ax.text(preX-0.02, preY+0.01, str(objId+1))
    ax.plot(np.asarray(meanX), np.asarray(meanY))
    if saveFig:
        plotter.savefig(mic.plotCart)
    return plotter



#class ProtMovieAlignmentWeb(ProtMovieAlignment):
#    """ Aligns a set of volumes using cross correlation.
#    Based on Xmipp protocol for aligning volumes, but
#    the parameters are restricted for ease of use.
#    """
#    _label = 'movie alignment web'
#    
#    def _defineParams(self, form):
#        ProtMovieAlignment._defineParams(self, form)
#        
#        gpuParamsGroup = form.getParam('GPU')
#        gpuParamsGroup.config(condition='False')
#        
#    def getOutputFiles(self):
#        # Redefine the default method to avoid download of movie files
#        return self.outputMicrographs.getFiles()
