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
import numpy as np

from pyworkflow.object import String
import pyworkflow.protocol.params as params
from pyworkflow.utils.path import moveFile
import pyworkflow.em as em
from pyworkflow.gui.plotter import Plotter



class ProtMovieAlignment(em.ProtProcessMovies):
    """ Aligns movies, from direct detectors cameras, into micrographs.
    """
    _label = 'optical alignment'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        em.ProtProcessMovies._defineParams(self, form)

        # Alignment parameters
        group = form.addGroup('Alignment parameters')
        line = group.addLine('Skip frames for alignment',
                            help='Skip frames for alignment.\n'
                                  'The first frame in the stack is *0*.' )
        line.addParam('alignFrame0', params.IntParam, default=0, label='Begin')
        line.addParam('alignFrameN', params.IntParam, default=0, label='End',
                      help='The number of frames to cut from the front and end')

        line = group.addLine('Crop offsets (px)')
        line.addParam('cropOffsetX', params.IntParam, default=0, label='X')
        line.addParam('cropOffsetY', params.IntParam, default=0, label='Y')

        line = group.addLine('Crop dimensions (px)',
                             help='How many pixels to crop from offset\n'
                                  'If equal to 0, use maximum size.')
        line.addParam('cropDimX', params.IntParam, default=0, label='X')
        line.addParam('cropDimY', params.IntParam, default=0, label='Y')

        # GROUP GPU PARAMETERS
        group = form.addGroup('GPU')
        group.addParam('doGPU', params.BooleanParam, default=False,
                       label="Use GPU (vs CPU)",
                       help="Set to true if you want the GPU implementation of "
                           "Optical Flow")
        group.addParam('GPUCore', params.IntParam, default=0,
                       expertLevel=params.LEVEL_ADVANCED,
                       label="Choose GPU core",
                       help="GPU may have several cores. Set it to zero if you "
                            "do not know what we are talking about. First core "
                            "index is 0, second 1 and so on.")
        
        # GROUP OPTICAL FLOW PARAMETERS
        group = form.addGroup('Additional Parameters',
                              expertLevel=params.LEVEL_ADVANCED)
        group.addParam('winSize', params.IntParam, default=150,
                       label="Window size",
                       help="Window size (shifts are assumed to be constant "
                            "within this window).")
        group.addParam('groupSize', params.IntParam, default=1,
                       label="Group Size",
                       help="The number of frames in each group at the last step")
        group.addParam('doSaveMovie', params.BooleanParam, default=False,
                       label="Save movies",
                       help="Save aligned movies")

    #--------------------------- STEPS functions ---------------------------------------------------
    def createOutputStep(self):

        inputMovies = self.inputMovies.get()
        micSet = self._createSetOfMicrographs()
        micSet.copyInfo(inputMovies)
        # Also create a Set of Movies with the alignment parameters
        if self.doSaveMovie:
            movieSet = self._createSetOfMovies()
            movieSet.copyInfo(inputMovies)
        for movie in inputMovies:
            micName = self._getNameExt(movie.getFileName(),'_aligned', 'mrc')
            metadataName = self._getNameExt(movie.getFileName(), '_aligned', 'xmd')
            plotCartName = self._getNameExt(movie.getFileName(), '_plot_cart', 'png')
            psdCorrName = self._getNameExt(movie.getFileName(),'_aligned_corrected', 'psd')
            alignedMovie = movie.clone()
            alignedMovie.alignMetaData = String(self._getExtraPath(metadataName))
            alignedMovie.plotCart = self._getExtraPath(plotCartName)
            alignedMovie.psdCorr = self._getExtraPath(psdCorrName)
            movieCreatePlot(alignedMovie, True)

            if self.doSaveMovie:
                movieSet.append(alignedMovie)
            mic = em.Micrograph()
            # All micrograph are copied to the 'extra' folder after each step
            mic.setFileName(self._getExtraPath(micName))
            # The micName of a micrograph MUST be the same as the original movie
            #mic.setMicName(micName)
            mic.setMicName(movie.getMicName())
            mic.plotCart = em.Image()
            mic.plotCart.setFileName(self._getExtraPath(plotCartName))

            mic.psdCorr = em.Image()
            mic.psdCorr.setFileName(self._getExtraPath(psdCorrName))
            micSet.append(mic)

        self._defineOutputs(outputMicrographs=micSet)
        self._defineSourceRelation(self.inputMovies, micSet)
        if self.doSaveMovie:
            self._defineOutputs(outputMovies=movieSet)

    #--------------------------- UTILS functions ---------------------------------------------------
    def _processMovie(self, movieId, movieName, movieFolder, movieAlignment):
        """ Process the movie actions, remember to:
        1) Generate all output files inside movieFolder (usually with cwd in runJob)
        2) Copy the important result files after processing (movieFolder will be deleted!!!)
        """

        movieSet = self.inputMovies.get()
        # Read the parameters
        #micName = self._getMicName(movieId)
        micName = self._getNameExt(movieName, '_aligned', 'mrc')
        metadataName = self._getNameExt(movieName, '_aligned', 'xmd')
        psdCorrName = self._getNameExt(movieName,'_aligned_corrected', 'psd')
        firstFrame = self.alignFrame0.get()
        lastFrame = self.alignFrameN.get()
        gpuId = self.GPUCore.get()

        # Check if we have global shifts
        if movieAlignment is not None:
            #firstFrame, lastFrame = movieAlignment.getRange()
            #regionInterest = movieAlignment.getRoi()
            shiftListX, shiftListY= movieAlignment.getShifts()
            for shiftX, shiftY in shiftListX, shiftListY:
                print shiftX, shiftY

        # Some movie have .mrc or .mrcs format but it is recognized as a volume
        if movieName.endswith('.mrcs') or movieName.endswith('.mrc'):
            movieSuffix = ':mrcs'
        else:
            movieSuffix = ''
        command = '-i %(movieName)s%(movieSuffix)s -o %(micName)s ' % locals()
        command += '--cutf %d --cute %d ' % (firstFrame, lastFrame)
        if self.inputMovies.get().getDark():
            command += '--dark '+self.inputMovies.get().getDark()
        if self.inputMovies.get().getGain():
            command += '--gain '+self.inputMovies.get().getGain()
        winSize = self.winSize.get()
        doSaveMovie = self.doSaveMovie.get()
        groupSize = self.groupSize.get()
        command += '--winSize %(winSize)d --groupSize %(groupSize)d ' % locals()
        if self.doGPU:
            program = 'xmipp_movie_optical_alignment_gpu'
            command += '--gpu %d ' % gpuId
        else:
            program = 'xmipp_movie_optical_alignment_cpu'
        if doSaveMovie:
            command += '--ssc '
        command += '--crx %d --cry %d --cdx %d --cdy %d' % (self.cropOffsetX,
                                                            self.cropOffsetY,
                                                            self.cropDimX,
                                                            self.cropDimY)
        try:
            self.runJob(program, command, cwd=movieFolder)
        except:
            print >> sys.stderr, program, " failed for movie %(movieName)s" % locals()
        moveFile(join(movieFolder, metadataName), self._getExtraPath())
        if doSaveMovie:
            outMovieName = self._getNameExt(movieName,'_aligned', 'mrcs')
            moveFile(join(movieFolder, outMovieName), self._getExtraPath())

        # Compute half-half PSD
        ih = em.ImageHandler()
        print join(movieFolder, '%(movieName)s' % locals())
        avg = ih.computeAverage(join(movieFolder, movieName))
        avg.write(join(movieFolder, 'uncorrectedmic.mrc'))
        command = '--micrograph uncorrectedmic.mrc --oroot uncorrectedpsd ' \
                  '--dont_estimate_ctf --pieceDim 400 --overlap 0.7'
        program = 'xmipp_ctf_estimate_from_micrograph'
        self.runJob(program, command, cwd=movieFolder)

        command = '--micrograph %(micName)s --oroot correctedpsd ' \
                  '--dont_estimate_ctf --pieceDim 400 --overlap 0.7' % locals()
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
        moveFile(join(movieFolder, psdCorrName), self._getExtraPath())

     #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        numThreads = self.numberOfThreads;
        if numThreads>1:
            if self.doGPU:
                errors.append("GPU and Parallelization can not be used together")
        return errors

    def _citations(self):
        return ['Abrishami2015']

    def _methods(self):
        methods = []
        if self.doGPU:
            gpuId = self.GPUCore.get()
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
        summary.append('The number of frames to cut from the front: '
                       '*%d* to *%s* (first frame is 0)' % (firstFrame, 'Last Frame'))

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


