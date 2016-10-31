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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from os.path import join, exists
from glob import glob
from itertools import izip

import pyworkflow.object as pwobj
import pyworkflow.utils as pwutils
import pyworkflow.em as em
from pyworkflow.gui.project import ProjectWindow
from pyworkflow.em.protocol import ProtAlignMovies
import pyworkflow.protocol.params as params
from pyworkflow.gui.plotter import Plotter
import pyworkflow.em.metadata as md
from convert import writeShiftsMovieAlignment, getMovieFileName


PLOT_CART = 0
PLOT_POLAR = 1
OBJCMD_MOVIE_ALIGNCARTESIAN = "Display Cartesian Presentation"


class XmippProtOFAlignment(ProtAlignMovies):
    """
    Wrapper protocol to Xmipp Movie Alignment by Optical Flow
    """
    _label = 'optical alignment'
    CONVERT_TO_MRC = 'mrcs'

    #--------------------------- DEFINE param functions ------------------------
    def _defineAlignmentParams(self, form):
        ProtAlignMovies._defineAlignmentParams(self, form)
        form.addSection("Aditional Parameters")
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
        form.addParam('winSize', params.IntParam, default=150,
                       label="Window size",
                       help="Window size (shifts are assumed to be constant "
                            "within this window).")
        form.addParam('groupSize', params.IntParam, default=1,
                       label="Group Size",
                       help="The number of frames in each group at the "
                            "last step")
        
        form.addParam('useAlignment', params.BooleanParam, default=True,
                      label="Use previous movie alignment to Sum frames?",
                      help="If set Yes, the alignment information (if"
                      " it exists) will take into account to align"
                      " your movies.")

        form.addParam('doComputePSD', params.BooleanParam, default=True,
                      label="Compute PSD (before/after)?",
                      help="If Yes, the protocol will compute for each movie "
                           "the PSD of the average micrograph (without OF "
                           "alignement) and after that, to compare each PSDs")

        form.addParallelSection(threads=8, mpi=0)
    
    #--------------------------- STEPS functions -------------------------------
    def _processMovie(self, movie):

        inputFn = self._getMovieOrMd(movie)
        outMovieFn = self._getExtraPath(self._getOutputMovieName(movie))
        outMicFn = self._getExtraPath(self._getOutputMicName(movie))
        aveMic = self._getFnInMovieFolder(movie, "uncorrected_mic.mrc")
        dark = self.inputMovies.get().getDark()
        gain = self.inputMovies.get().getGain()
        # Get the number of frames and the range
        # to be used for alignment and sum
        x, y, n = movie.getDim()
        a0, aN = self._getFrameRange(n, 'align')
        gpuId = self.GPUCore.get()

        args = '-i %s ' % inputFn
        args += '-o %s ' % self._getOutputShifts(movie)
        args += '--frameRange %d %d ' % (a0-1, aN-1)

        if dark:
            args += '--dark %s ' % dark
        if gain:
            args += '--gain %s ' % gain
        winSize = self.winSize.get()
        doSaveMovie = self.doSaveMovie.get()
        groupSize = self.groupSize.get()
        args += ' --winSize %(winSize)d --groupSize %(groupSize)d ' % locals()
        
        if self.doGPU:
            program = 'xmipp_movie_optical_alignment_gpu'
            args += '--gpu %d ' % gpuId
        else:
            program = 'xmipp_movie_optical_alignment_cpu'

        # We should save the movie either if the user selected it (default)
        # or if the PSD is going to be computed
        if self.doSaveAveMic or self.doComputePSD:
            args += '--oavg %s ' % outMicFn

        if self.doComputePSD:
            args  += '--oavgInitial %s ' % aveMic

        if doSaveMovie:
            args += '--outMovie %s ' % outMovieFn
        
        roi = [self.cropOffsetX.get(), self.cropOffsetY.get(),
               self.cropDimX.get(), self.cropDimY.get()]
        
        args += '--cropULCorner %d %d ' % (roi[0], roi[1])
        args += '--cropDRCorner %d %d ' % (roi[0] + roi[2] -1,
                                           roi[1] + roi[3] -1)
        try:
            self.runJob(program, args)
            
            if self.doSaveAveMic:
                if not exists(outMicFn):
                    raise Exception("Micrograph %s not produced after "
                                    "running %s " % (outMicFn, program))

            if self.doSaveMovie:
                if not exists(outMovieFn):
                    raise Exception("Movie %s not produced after running %s "
                                    % (outMovieFn, program))

            if self.doComputePSD:
                uncorrectedPSD = self._getFnInMovieFolder(movie, "uncorrected")
                correctedPSD = self._getFnInMovieFolder(movie, "corrected")
                # TODO: Compute the PSD inside the OF program?
                self.computePSD(aveMic, uncorrectedPSD)
                self.computePSD(outMicFn, correctedPSD)
                self.composePSD(uncorrectedPSD + ".psd",
                                correctedPSD + ".psd",
                                self._getPsdCorr(movie))
                # If the micrograph was only saved for computing the PSD
                # we can remove it
                if not self.doSaveAveMic:
                    pwutils.cleanPath(outMicFn)
            
            self._saveAlignmentPlots(movie)

        except Exception as e:
            print ("ERROR: %s failed for movie %s.\n  Exception: %s"
                   % (program, inputFn, e))
    
    #--------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = ProtAlignMovies._validate(self)
        if self.numberOfThreads > 1 and self.doGPU:
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
        summary = []
        inputSet = self.inputMovies.get()

        if inputSet:
            summary.append('Number of input movies: *%d*' % inputSet.getSize())
            _, _, n = inputSet.getDim()
            a0, aN = self._getFrameRange(n, 'align')
            summary.append("Frames from *%d* to *%d* were aligned" % (a0, aN))

        return summary
    
    #--------------------------- UTILS functions -------------------------------
    def _setAlignmentInfo(self, movie, obj):
        """ Set alignment info such as plot and psd filename, and
        the cumulative shifts values.
        Params:
            movie: Pass the reference movie
            obj: should pass either the created micrograph or movie
        """
        obj.plotCart = em.Image()
        obj.plotCart.setFileName(self._getPlotCart(movie))
        if self.doComputePSD:
            obj.psdCorr = em.Image()
            obj.psdCorr.setFileName(self._getPsdCorr(movie))

        meanX, meanY = self._loadMeanShifts(movie)
        obj._xmipp_OFMeanX = pwobj.CsvList()
        obj._xmipp_OFMeanX.set(meanX)
        obj._xmipp_OFMeanY = pwobj.CsvList()
        obj._xmipp_OFMeanY.set(meanY)

    def _preprocessOutputMicrograph(self, mic, movie):
        self._setAlignmentInfo(movie, mic)

    def _createOutputMovie(self, movie):
        alignedMovie = ProtAlignMovies._createOutputMovie(self, movie)
        self._setAlignmentInfo(movie, alignedMovie)
        return alignedMovie

    def _getNameExt(self, movie, postFix, ext, extra=False):
        fn = self._getMovieRoot(movie) + postFix + '.' + ext
        if extra:
            return self._getExtraPath(fn)
        else:
            return fn
    
    def _doGenerateOutputMovies(self):
        """ Returns True if an output set of movies will be generated.
        The most common case is to always generate output movies,
        either with alignment only or the binary aligned movie files.
        Subclasses can override this function to change this behavior.
        """
        return self.doSaveMovie.get()

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
        return self._getNameExt(movie, '_shifts', '.xmd', extra=False)

    def _getOutputShifts(self, movie):
        return self._getNameExt(movie, '_aligned_mic', 'xmd', extra=True)

    def _getPlotCart(self, movie):
        return self._getNameExt(movie, '_plot_cart', 'png', extra=True)

    def _getPsdCorr(self, movie):
        return self._getNameExt(movie, '_aligned_corrected', 'psd', extra=True)

    def _loadMeanShifts(self, movie):
        alignMd = md.MetaData(self._getOutputShifts(movie))
        meanX = alignMd.getColumnValues(md.MDL_OPTICALFLOW_MEANX)
        meanY = alignMd.getColumnValues(md.MDL_OPTICALFLOW_MEANY)

        return meanX, meanY

    def _saveAlignmentPlots(self, movie):
        """ Compute alignment shifts plot and save to file as a png image. """
        meanX, meanY = self._loadMeanShifts(movie)
        plotter = createAlignmentPlot(meanX, meanY)
        plotter.savefig(self._getPlotCart(movie))


def showCartesianShiftsPlot(inputSet, itemId):
    item = inputSet[itemId]
    if item.hasAttribute('_xmipp_OFMeanX'):
        meanX = [float(x) for x in item._xmipp_OFMeanX]
        meanY = [float(y) for y in item._xmipp_OFMeanY]
        plotter = createAlignmentPlot(meanX, meanY)
        plotter.show()
    else:
        print "This items does not have OF alignment set. "


ProjectWindow.registerObjectCommand(OBJCMD_MOVIE_ALIGNCARTESIAN,
                                    showCartesianShiftsPlot)

def createAlignmentPlot(meanX, meanY):
    """ Create a plotter with the cumulative shift per frame. """
    sumMeanX = []
    sumMeanY = []
    figureSize = (8, 6)
    plotter = Plotter(*figureSize)
    figure = plotter.getFigure()

    preX = 0.0
    preY = 0.0
    sumMeanX.append(0.0)
    sumMeanY.append(0.0)
    ax = figure.add_subplot(111)
    ax.grid()
    ax.set_title('Cartesian representation')
    ax.set_xlabel('Drift x (pixels)')
    ax.set_ylabel('Drift y (pixels)')
    ax.plot(0, 0, 'yo-')
    i = 1
    for x, y in izip(meanX, meanY):
        preX += x
        preY += y
        sumMeanX.append(preX)
        sumMeanY.append(preY)
        #ax.plot(preX, preY, 'yo-')
        ax.text(preX-0.02, preY+0.02, str(i))
        i += 1

    ax.plot(sumMeanX, sumMeanY, color='b')
    ax.plot(sumMeanX, sumMeanY, 'yo')

    plotter.tightLayout()

    return plotter

# Just for backwards compatibility
ProtMovieAlignment = XmippProtOFAlignment

