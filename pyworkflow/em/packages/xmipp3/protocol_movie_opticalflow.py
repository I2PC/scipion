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

from os.path import join, exists
import sys
import numpy as np

import pyworkflow.em as em
from pyworkflow.em.protocol import ProtAlignMovies
from pyworkflow.utils.path import moveFile
import pyworkflow.protocol.params as params
from pyworkflow.gui.plotter import Plotter
from convert import writeShiftsMovieAlignment, getMovieFileName


class XmippProtOFAlignment(ProtAlignMovies):
    """
    Wrapper protocol to Xmipp Movie Alignment by Optical Flow
    """
    _label = 'optical alignment'
    CONVERT_TO_MRC = 'mrcs'
    

    #--------------------------- DEFINE param functions ------------------------
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
        form.addParallelSection(threads=8, mpi=0)
    
    #--------------------------- STEPS functions -------------------------------
    def _processMovie(self, movie):

        inputFn = self._getMovieOrMd(movie)
        outputMicFn = self._getFnInMovieFolder(movie, self._getOutputMicName(movie))
        outMovieName = self._getExtraPath(self._getOutputMovieName(movie))
        #outputMicFn = self._getExtraPath(self._getOutputMicName(movie))
        metadataName = self._getNameExt(movie, '_aligned_mic', 'xmd')
        dark = self.inputMovies.get().getDark()
        gain = self.inputMovies.get().getGain()
        # Get the number of frames and the range to be used for alignment and sum
        x, y, n = movie.getDim()
        a0, aN = self._getFrameRange(n, 'align')
        psdCorrName = self._getExtraPath(self._getNameExt(movie,'_aligned_corrected', 'psd'))
        gpuId = self.GPUCore.get()

        command = '-i %s ' % inputFn
        command += '-o %s ' % self._getExtraPath(metadataName)
        command += '--frameRange %d %d ' % (a0-1, aN-1)

        if dark:
            command += '--dark %s ' % dark
        if gain:
            command += '--gain %s ' % gain
        winSize = self.winSize.get()
        doSaveMovie = self.doSaveMovie.get()
        groupSize = self.groupSize.get()
        command += ' --winSize %(winSize)d --groupSize %(groupSize)d ' % locals()
        if self.doGPU:
            program = 'xmipp_movie_optical_alignment_gpu'
            command += '--gpu %d ' % gpuId
        else:
            program = 'xmipp_movie_optical_alignment_cpu'
        command += '--oavg %s ' % outputMicFn
        if doSaveMovie:
            command += '--outMovie %s ' % outMovieName
        
        roi = [self.cropOffsetX.get(), self.cropOffsetY.get(),
               self.cropDimX.get(), self.cropDimY.get()]
        
        command += '--cropULCorner %d %d ' % (roi[0], roi[1])
        command += '--cropDRCorner %d %d ' % (roi[0] + roi[2] -1,
                                              roi[1] + roi[3] -1)
        try:
            self.runJob(program, command)
            if not exists(outputMicFn):
                raise Exception("Micrograph %s not produced after running %s " %(outputMicFn, program))
            if self.doSaveMovie:
                if not exists(outMovieName):
                    raise Exception("Movie %s not produced after running %s " %(outMovieName, program))
            
            aveMic = self._getFnInMovieFolder(movie, "uncorrected_mic.mrc")
            self.averageMovie(movie, inputFn, aveMic, self.binFactor.get(),
                              roi, dark, gain)
            uncorrectedPSD = self._getFnInMovieFolder(movie, "uncorrected")
            correctedPSD = self._getFnInMovieFolder(movie, "corrected")
            
            self.computePSD(aveMic, uncorrectedPSD)
            self.computePSD(outputMicFn, correctedPSD)
            self.composePSD(uncorrectedPSD + ".psd",
                            correctedPSD + ".psd", psdCorrName)
            if self.doSaveAveMic:
                moveFile(outputMicFn, self._getExtraPath(self._getOutputMicName(movie)))
        except:
            print >> sys.stderr, program, " failed for movie %s" % inputFn

    
    #--------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        numThreads = self.numberOfThreads
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
            summary.append('Number of input movies: '
                           '*%d*' % self.inputMovies.get().getSize())
        summary.append('The number of frames to cut from the front: '
                       '*%d* to *%s* (first frame is 0)' % (firstFrame, 'Last Frame'))

        return summary
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _preprocessOutputMicrograph(self, mic, movie):

        plotCartName = self._getNameExt(movie, '_plot_cart', 'png')
        psdCorrName = self._getExtraPath(self._getNameExt(movie,'_aligned_corrected', 'psd'))
        metadataName = self._getNameExt(movie, '_aligned_mic', 'xmd')

        mic.alignMetaData = self._getExtraPath(metadataName)
        mic.plotCart = self._getExtraPath(plotCartName)
        # Create plot
        movieCreatePlot(mic, True)
        mic.plotCart = em.Image()
        mic.plotCart.setFileName(self._getExtraPath(plotCartName))
        mic.psdCorr = em.Image()
        mic.psdCorr.setFileName(psdCorrName)
    
    def _getNameExt(self, movie, postFix, ext):
        return self._getMovieRoot(movie) + postFix + '.' + ext
    
    def _doGenerateOutputMovies(self):
        """ Returns True if an output set of movies will be generated.
        The most common case is to always generate output movies,
        either with alignment only or the binary aligned movie files.
        Subclasses can override this function to change this behavior.
        """
        if self.doSaveMovie:
            return True
        else:
            return False
    
    def _getFnInMovieFolder(self, movie, filename):
        movieFolder = self._getOutputMovieFolder(movie)
        return join(movieFolder, filename)

    def _getMovieOrMd(self, movie):
        if movie.hasAlignment() and self.useAlignment:
            shiftsMd = self._getShiftsFile(movie)
            numberOfFrames = movie.getNumberOfFrames()
            s0, sN = self._getFrameRange(numberOfFrames, 'sum')
            writeShiftsMovieAlignment(movie, shiftsMd, s0, sN)
            return shiftsMd
        else:
            return getMovieFileName(movie)

    def _getShiftsFile(self, movie):
        movieFolder = self._getOutputMovieFolder(movie)
        return join(movieFolder, self._getMovieRoot(movie) + '_shifts.xmd')

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

