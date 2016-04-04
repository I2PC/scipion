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

from os.path import join
import sys
import numpy as np

import pyworkflow.em as em
from pyworkflow.em.protocol import ProtAlignMovies
import pyworkflow.protocol.params as params
from pyworkflow.gui.plotter import Plotter
from pyworkflow.utils.path import moveFile
from convert import writeShiftsMovieAlignment, getMovieFileName


class XmippProtOFAlignment(ProtAlignMovies):
    """
    Wrapper protocol to Xmipp Movie Alignment by Optical Flow
    """
    _label = 'optical alignment'
    CONVERT_TO_MRC = 'mrcs'
    

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineAlignmentParams(self, form):
        ProtAlignMovies._defineAlignmentParams(self, form)
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
        
        form.addParam('useAlignment', params.BooleanParam, default=True,
              label="Use movie alignment to Sum frames?",
              help="If set Yes, the alignment information (if"
                   " it exists) will take into account to align"
                   " your movies.")
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def _processMovie(self, movie):
        movieName = getMovieFileName(movie)
        movieFolder = self._getOutputMovieFolder(movie)
        outputMicFn = self._getOutputMicName(movie)

        # Get the number of frames and the range to be used for alignment and sum
        x, y, n = movie.getDim()
        a0, aN = self._getFrameRange(n, 'align')
        s0, sN = self._getFrameRange(n, 'sum')
        
        metadataName = self._getNameExt(movie, '_aligned', 'xmd')
        fnGlobalShifts = self._getNameExt(movie, '_shifts', 'xmd')
        psdCorrName = self._getNameExt(movie,'_aligned_corrected', 'psd')
        gpuId = self.GPUCore.get()

        command = ' -i %s -o %s ' % (movieName, outputMicFn)
        command += '--frameRange %d %d ' % (a0-1, aN-1)
        
        if self.inputMovies.get().getDark():
            command += '--dark %s ' % self.inputMovies.get().getDark()
        
        if self.inputMovies.get().getGain():
            command += '--gain %s ' % self.inputMovies.get().getGain()
        
        winSize = self.winSize.get()
        doSaveMovie = self.doSaveMovie.get()
        groupSize = self.groupSize.get()
        command += ' --winSize %(winSize)d --groupSize %(groupSize)d ' % locals()
        
        # Check if we have global shifts
        if movie.hasAlignment() and self.useAlignment:
            mdFn = join(movieFolder, fnGlobalShifts)
            writeShiftsMovieAlignment(movie, mdFn, s0, sN)
            
            command += '--useInputShifts %(fnGlobalShifts)s ' % locals()
        if self.doGPU:
            program = 'xmipp_movie_optical_alignment_gpu'
            command += '--gpu %d ' % gpuId
        else:
            program = 'xmipp_movie_optical_alignment_cpu'
        if doSaveMovie:
            command += '--ssc '
        command += '--cropULCorner %d %d ' % (self.cropOffsetX, self.cropOffsetY)
        command += '--cropDRCorner %d %d'  % (self.cropOffsetX.get() + self.cropDimX.get() -1,
                                                  self.cropOffsetY.get() + self.cropDimY.get() -1)
        try:
            self.runJob(program, command, cwd=movieFolder)
        except:
            print >> sys.stderr, program, " failed for movie %(movieName)s" % locals()
        moveFile(join(movieFolder, metadataName), self._getExtraPath())
        if doSaveMovie:
            moveFile(join(movieFolder, self._getOutputMovieName(movie)),
                     self._getExtraPath())

        # Compute half-half PSD
        ih = em.ImageHandler()
        avg = ih.computeAverage(join(movieFolder, movieName))
        avg.write(join(movieFolder, 'uncorrectedmic.mrc'))
        command = '--micrograph uncorrectedmic.mrc --oroot uncorrectedpsd ' \
                  '--dont_estimate_ctf --pieceDim 400 --overlap 0.7'
        program = 'xmipp_ctf_estimate_from_micrograph'
        self.runJob(program, command, cwd=movieFolder)

        command = '--micrograph %(outputMicFn)s --oroot correctedpsd ' \
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
        moveFile(join(movieFolder, outputMicFn), self._getExtraPath())
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
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _preprocessOutputMicrograph(self, mic, movie):
        plotCartName = self._getNameExt(movie, '_plot_cart', 'png')
        psdCorrName = self._getNameExt(movie,'_aligned_corrected', 'psd')
        metadataName = self._getNameExt(movie, '_aligned', 'xmd')
        
        mic.alignMetaData = self._getExtraPath(metadataName)
        mic.plotCart = self._getExtraPath(plotCartName)
        # Create plot
        movieCreatePlot(mic, True)
        mic.plotCart = em.Image()
        mic.plotCart.setFileName(self._getExtraPath(plotCartName))
        mic.psdCorr = em.Image()
        mic.psdCorr.setFileName(self._getExtraPath(psdCorrName))
    
    def _getNameExt(self, movie, postFix, ext):
        return self._getMovieRoot(movie) + '_aligned_mic' + '.' + ext
    
    def _doGenerateOutputMovies(self):
        return True


def createPlots(plotType, protocol, micId):
    print "output Micrographs to create Plot %s" % protocol.outputMicrographs
    mic = protocol.outputMicrographs[micId]
    return movieCreatePlot(mic, False).show()


def movieCreatePlot(mic, saveFig):
    import pyworkflow.em.metadata as md

    meanX = []
    meanY = []
    figureSize = (8, 6)

    mdAlign = md.MetaData(mic.alignMetaData)
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
    for objId in mdAlign:
        preX += mdAlign.getValue(md.MDL_OPTICALFLOW_MEANX, objId)
        preY += mdAlign.getValue(md.MDL_OPTICALFLOW_MEANY, objId)
        meanX.append(preX)
        meanY.append(preY)
        ax.plot(preX, preY, 'yo-')
        ax.text(preX-0.02, preY+0.01, str(objId+1))
    ax.plot(np.asarray(meanX), np.asarray(meanY))
    if saveFig:
        plotter.savefig(mic.plotCart)
    return plotter

# Just for backwards compatibility
ProtMovieAlignment = XmippProtOFAlignment
