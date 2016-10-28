# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Vahid Abrishami (vabrishami@cnb.csic.es)
# *              Josue Gomez Blanco (jgomez@cnb.csic.es)
# *              Grigory Sharov (sharov@igbmc.fr)
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

import os, sys
from pyworkflow.gui.plotter import plt
from itertools import izip

import pyworkflow.protocol.params as params
import pyworkflow.protocol.constants as cons
import pyworkflow.utils as pwutils
import pyworkflow.em as em
from pyworkflow.gui.project import ProjectWindow
from pyworkflow.em.data import MovieAlignment
from pyworkflow.em.packages.xmipp3.convert import writeShiftsMovieAlignment
from pyworkflow.em.protocol import ProtAlignMovies
from pyworkflow.gui.plotter import Plotter
from convert import (MOTIONCORR_PATH, MOTIONCOR2_PATH, getVersion,
                     parseMovieAlignment, parseMovieAlignment2,
                     parseMovieAlignmentLocal, convertShifts)


OBJCMD_MOVIE_ALIGNLOCAL = "Display patch alignment plot"


class ProtMotionCorr(ProtAlignMovies):
    """
    Wrapper protocol to movie alignment programs developed at UCSF:
    motioncorr: Flat fielding and Drift correction
        (written by Xueming Li @ Yifan Cheng Lab)
    motioncor2: anisotropic drift correction and dose weighting
        (written by Shawn Zheng @ David Agard lab)
    """

    _label = 'motioncorr alignment'
    CONVERT_TO_MRC = 'mrc'

    #--------------------------- DEFINE param functions ------------------------
    def _defineAlignmentParams(self, form):
        form.addParam('gpuMsg', params.LabelParam, default=True,
                      label='WARNING! You need to have installed CUDA'
                            ' libraries and a Nvidia GPU')

        form.addParam('GPUIDs', params.StringParam, default='0',
                      expertLevel=cons.LEVEL_ADVANCED,
                      label="Choose GPU IDs",
                      help="GPU may have several cores. Set it to zero"
                           " if you do not know what we are talking about."
                           " First core index is 0, second 1 and so on."
                           " Motioncor2 can use multiple GPUs - in that case"
                           " set to i.e. *0 1 2*.")

        ProtAlignMovies._defineAlignmentParams(self, form)

        form.addParam('doComputePSD', params.BooleanParam, default=False,
                      expertLevel=cons.LEVEL_ADVANCED,
                      label="Compute PSD (before/after)?",
                      help="If Yes, the protocol will compute for each movie "
                           "the average PSD before and after alignment, for comparison")

        form.addParam('extraParams', params.StringParam, default='',
                      expertLevel=cons.LEVEL_ADVANCED,
                      label='Additional parameters',
                      help="""Extra parameters for motioncorr (NOT motioncor2)\n
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

        form.addSection(label="Motioncor2")
        form.addParam('useMotioncor2', params.BooleanParam, default=False,
                      label='Use motioncor2',
                      help='Use new *motioncor2* program with local '
                           'patch-based motion correction and dose weighting.')

        line = form.addLine('Number of patches', condition='useMotioncor2',
                            help='Number of patches to be used for patch based '
                                 'alignment. Set to *0 0* to do only global motion '
                                 'correction.')
        line.addParam('patchX', params.IntParam, default=5, label='X')
        line.addParam('patchY', params.IntParam, default=5, label='Y')

        form.addParam('frameDose', params.FloatParam, default='0.0',
                      label='Frame dose (e/A^2)', condition='useMotioncor2',
                      help='Frame dose in e/A^2. If set to *0.0*, dose '
                            'weighting will be skipped.')

        form.addParam('initDose', params.FloatParam, default='0.0',
                      expertLevel=cons.LEVEL_ADVANCED,
                      label='Pre-exposure (e/A^2)',
                      condition='useMotioncor2 and frameDose and _isNewMotioncor2',
                      help='Initial dose received before stack is acquired, in e/A^2.')

        form.addParam('group', params.IntParam, default='1',
                      label='Group N frames', condition='useMotioncor2',
                      help='Group every specified number of frames by adding '
                           'them together. The alignment is then performed on '
                           'the summed frames. By default, no grouping is '
                           'performed.')

        form.addParam('tol', params.FloatParam, default='0.5',
                      label='Tolerance (px)', condition='useMotioncor2',
                      help='Tolerance for iterative alignment, default *0.5px*.')

        form.addParam('extraParams2', params.StringParam, default='',
                      expertLevel=cons.LEVEL_ADVANCED, condition='useMotioncor2',
                      label='Additional parameters',
                      help="""Extra parameters for motioncor2\n
        -Bft       100        BFactor for alignment, in px^2.
        -Iter      5          Maximum iterations for iterative alignment.
        -MaskCent  0 0        Center of subarea that will be used for alignment,
                              default *0 0* corresponding to the frame center.
        -MaskSize  1.0 1.0    The size of subarea that will be used for alignment,
                              default *1.0 1.0* corresponding full size.
        -Align     1          Generate aligned sum (1) or simple sum (0).
        -FmRef     0          Specify which frame to be the reference to which
                              all other frames are aligned, by default *0* all aligned
                              to the first frame,
                              other value aligns to the central frame.
        -RotGain   0          Rotate gain reference counter-clockwise: 0 - no rotation,
                              1 - 90 degrees, 2 - 180 degrees, 3 - 270 degrees.
        -FlipGain  0          Flip gain reference after gain rotation: 0 - no flipping,
                              1 - flip upside down, 2 - flip left right.
        -Tilt      0 0        Tilt angle range for a dose fractionated tomographic
                              tilt series, e.g. *-60 60*
                              """)

        # Since only runs on GPU, do not allow neither threads nor mpi
        form.addParallelSection(threads=0, mpi=0)

    #--------------------------- STEPS functions -------------------------------
    def _processMovie(self, movie):
        inputMovies = self.inputMovies.get()
        movieFolder = self._getOutputMovieFolder(movie)
        outputMicFn = self._getRelPath(self._getOutputMicName(movie),
                                       movieFolder)
        aveMicFn = pwutils.removeExt(movie.getFileName()) + '_uncorrected_avg.mrc'

        if not self.useMotioncor2:
            outputMovieFn = self._getRelPath(self._getOutputMovieName(movie),
                                             movieFolder)
            logFile = self._getRelPath(self._getMovieLogFile(movie),
                                       movieFolder)

            # Get the number of frames and the range to be used
            # for alignment and sum
            numberOfFrames = movie.getNumberOfFrames()
            a0, aN = self._getFrameRange(numberOfFrames, 'align')
            s0, sN = self._getFrameRange(numberOfFrames, 'sum')

            argsDict = {'-crx': self.cropOffsetX.get(),
                        '-cry': self.cropOffsetY.get(),
                        '-cdx': self.cropDimX.get(),
                        '-cdy': self.cropDimY.get(),
                        '-bin': self.binFactor.get(),
                        '-nst': '%d' % (a0-1),
                        '-ned': '%d' % (aN-1),
                        '-nss': '%d' % (s0-1),
                        '-nes': '%d' % (sN-1),
                        '-gpu': self.GPUIDs.get(),
                        '-flg': logFile,
                        }

            args = '"%s" ' % movie.getBaseName()
            args += ' '.join(['%s %s' % (k, v) for k, v in argsDict.iteritems()])

            if inputMovies.getGain():
                args += ' -fgr "%s"' % inputMovies.getGain()

            if inputMovies.getDark():
                args += ' -fdr "%s"' % inputMovies.getDark()

            if self.doSaveAveMic:
                args += ' -fcs "%s" ' % outputMicFn

            if self.doSaveMovie:
                args += ' -fct "%s" -ssc 1' % outputMovieFn

            args += ' ' + self.extraParams.get()
            program = MOTIONCORR_PATH

        else:
            logFileFn = self._getRelPath(self._getMovieLogFile(movie),
                                         movieFolder)
            logFileBase = logFileFn.replace('0-Full.log', '').replace('0-Patch-Full.log', '')

            # Get the number of frames and the range to be used
            # for alignment and sum
            numberOfFrames = movie.getNumberOfFrames()
            a0, aN = self._getFrameRange(numberOfFrames, 'align')

            # default values for motioncor2 are (1, 1)
            cropDimX = self.cropDimX.get() or 1
            cropDimY = self.cropDimY.get() or 1

            argsDict = {'-OutMrc': '"%s"' % outputMicFn,
                        '-Patch': '%d %d' % (self.patchX, self.patchY),
                        '-MaskCent': '%d %d' % (self.cropOffsetX,
                                                self.cropOffsetY),
                        '-MaskSize': '%d %d' % (cropDimX, cropDimY),
                        '-FtBin': self.binFactor.get(),
                        '-Tol': self.tol.get(),
                        '-Group': self.group.get(),
                        '-FmDose': self.frameDose.get(),
                        '-Throw': '%d' % (a0 - 1),
                        '-Trunc': '%d' % (abs(aN - numberOfFrames)),
                        '-PixSize': inputMovies.getSamplingRate(),
                        '-kV': inputMovies.getAcquisition().getVoltage(),
                        '-Gpu': self.GPUIDs.get(),
                        '-LogFile': logFileBase,
                        }
            if getVersion('MOTIONCOR2') != '03162016':
                argsDict['-InitDose'] = self.initDose.get()

            args = ' -InMrc "%s" ' % movie.getBaseName()
            args += ' '.join(['%s %s' % (k, v) for k, v in argsDict.iteritems()])

            if inputMovies.getGain():
                args += ' -Gain "%s" ' % inputMovies.getGain()

            args += ' ' + self.extraParams2.get()
            program = MOTIONCOR2_PATH

        try:
            self.runJob(program, args, cwd=movieFolder)

            #if self.doComputePSD:
            #    uncorrectedPSD = pwutils.removeExt(movie.getFileName()) + '_uncorrected'
            #    correctedPSD = pwutils.removeExt(movie.getFileName()) + '_corrected'
            #    # Compute uncorrected avg mic
            #    roi = [self.cropOffsetX.get(), self.cropOffsetY.get(),
            #           self.cropDimX.get(), self.cropDimY.get()]
            #    fakeShiftsFn = self.writeZeroShifts(movie)
            #    self.averageMovie(movie, fakeShiftsFn, aveMicFn,
            #                      binFactor=self.binFactor.get(),
            #                      roi=roi, dark=None,
            #                      gain=inputMovies.getGain())
            #    # Compute PSDs
            #    self.computePSD(aveMicFn, uncorrectedPSD)
            #    self.computePSD(outputMicFn, correctedPSD)
            #    self.composePSD(uncorrectedPSD + ".psd",
            #                    correctedPSD + ".psd",
            #                    self._getPsdCorr(movie))
            #    # Remove avg that was used only for computing PSD
            #    pwutils.cleanPath(aveMicFn)

            self._saveAlignmentPlots(movie)

        except:
            print >> sys.stderr, program, " failed for movie %s" % movie.getName()

    #--------------------------- INFO functions --------------------------------
    def _summary(self):
        summary = []
        return summary

    def _validate(self):
        # Check base validation before the specific ones for Motioncorr
        errors = ProtAlignMovies._validate(self)

        program = MOTIONCOR2_PATH if self.useMotioncor2 else MOTIONCORR_PATH

        if not os.path.exists(program):
            errors.append('Missing %s' % program)

        gpu = self.GPUIDs.get()

        if not self.useMotioncor2:
            bin = self.binFactor.get()

            if not (bin == 1.0 or bin == 2.0):
                errors.append("Binning factor can only be 1 or 2")

            if len(gpu) > 1:
                errors.append("Old motioncorr2.1 does not support multiple "
                              "GPUs, use motioncor2.")
        else:
            if not self.doSaveAveMic:
                errors.append('Option not supported. Please select Yes for '
                              'Save aligned micrograph. ')
                errors.append('Optionally you could add -Align 0 to additional'
                              ' parameters so that protocol ')
                errors.append('produces simple movie sum.')

            if self.doSaveMovie:
                errors.append('Saving aligned movies is not supported by '
                              'motioncor2. ')
                errors.append('By default, the protocol will produce '
                              'outputMovies equivalent to the input ')
                errors.append('however containing alignment information.')

            if (self.alignFrame0.get() != self.sumFrame0.get() or 
                self.alignFrameN.get() != self.sumFrameN.get()):
                errors.append('Frame range for ALIGN and SUM must be '
                              'equivalent in case of motioncor2.')

        return errors

    #--------------------------- UTILS functions ------------------------------
    def _getMovieLogFile(self, movie):
        if not self.useMotioncor2:
            return 'micrograph_%06d_Log.txt' % movie.getObjId()
        else:
            if self.patchX == 0 and self.patchY == 0:
                return 'micrograph_%06d_0-Full.log' % movie.getObjId()
            else:
                return 'micrograph_%06d_0-Patch-Full.log' % movie.getObjId()

    def _getAbsPath(self, baseName):
        return os.path.abspath(self._getExtraPath(baseName))

    def _getRelPath(self, baseName, refPath):
        return os.path.relpath(self._getExtraPath(baseName), refPath)

    def _getNameExt(self, movie, postFix, ext, extra=False):
        fn = self._getMovieRoot(movie) + postFix + '.' + ext
        if extra:
            return self._getExtraPath(fn)
        else:
            return fn

    def _preprocessOutputMicrograph(self, mic, movie):
        self._setPlotInfo(movie, mic)

    def _getMovieShifts(self, movie):
        """ Returns the x and y shifts for the alignment of this movie.
        The shifts should refer to the original micrograph without any binning.
        In case of a binning greater than 1, the shifts should be scaled.
        """
        logPath = self._getExtraPath(self._getMovieLogFile(movie))
        binning = self.binFactor.get()
        if not self.useMotioncor2:
            xShifts, yShifts = parseMovieAlignment(logPath)
        else:
            xShifts, yShifts = parseMovieAlignment2(logPath)
        xSfhtsCorr = [x * binning for x in xShifts]
        ySfhtsCorr = [y * binning for y in yShifts]
        return xSfhtsCorr, ySfhtsCorr

    def _getLocalShifts(self, movie):
        """ Returns local shifts per patch for each frame.
        First shift is (0,0). """
        patchFn = 'micrograph_%06d_0-Patch-Patch.log' % movie.getObjId()
        logPath = self._getExtraPath(patchFn)
        binning = self.binFactor.get()
        patchX = self.patchX.get()
        patchY = self.patchY.get()

        # parse log files, multiply by bin factor
        xShifts, yShifts = parseMovieAlignmentLocal(logPath)
        xSfhtsCorr = [x * binning for x in xShifts]
        ySfhtsCorr = [y * binning for y in yShifts]

        # convert shifts for plot format
        xSfhtsFinal, ySfhtsFinal = convertShifts(xSfhtsCorr, ySfhtsCorr, patchX, patchY)
        return xSfhtsFinal, ySfhtsFinal, patchX, patchY

    def _setPlotInfo(self, movie, obj):
        obj.plotGlobal = em.Image()
        obj.plotGlobal.setFileName(self._getPlotGlobal(movie))

        #if self.checkPatchALign():
        #    obj.plotLocal = em.Image(location=None)

    def _getPlotGlobal(self, movie):
        return self._getNameExt(movie, '_global_shifts', 'png', extra=True)

    def _getPsdCorr(self, movie):
        return self._getNameExt(movie, '_aligned_corrected', 'psd')

    def _saveAlignmentPlots(self, movie):
        """ Compute alignment shift plots and save to file as png images. """
        shiftsX, shiftsY = self._getMovieShifts(movie)
        plotter = createGlobalAlignmentPlot(shiftsX, shiftsY)
        plotter.savefig(self._getPlotGlobal(movie))

    def _isNewMotioncor2(self):
        return True if getVersion('MOTIONCOR2') != '03162016' else False

    def writeZeroShifts(self, movie):
        #TODO: find another way to do this
        shiftsMd = self._getTmpPath('zero_shifts.xmd')
        pwutils.cleanPath(shiftsMd)
        xshifts = [0] * movie.getNumberOfFrames()
        yshifts = xshifts
        alignment = MovieAlignment(first=1, last=movie.getNumberOfFrames(),
                                   xshifts=xshifts, yshifts=yshifts)
        roiList = [0, 0, 0, 0]
        alignment.setRoi(roiList)
        movie.setAlignment(alignment)
        writeShiftsMovieAlignment(movie, shiftsMd,
                                  1, movie.getNumberOfFrames())
        return shiftsMd

    def checkPatchALign(self):
        """ Check if we used motioncor2 with patch-based alignment. """
        if self.useMotioncor2 and self.patchX != 0 and self.patchY != 0:
            return True
        else:
            return False


def createGlobalAlignmentPlot(meanX, meanY):
    """ Create a plotter with the cumulative shift per frame. """
    sumMeanX = []
    sumMeanY = []
    preX = 0.0
    preY = 0.0

    figureSize = (6, 4)
    plotter = Plotter(*figureSize)
    figure = plotter.getFigure()
    ax = figure.add_subplot(111)
    ax.grid()
    ax.set_title('Global frame shift (cumulative)')
    ax.set_xlabel('Shift x (pixels)')
    ax.set_ylabel('Shift y (pixels)')
    if meanX[0] != 0 or meanY[0] != 0:
        raise Exception("First frame shift must be (0,0)!")

    i = 1
    for x, y in izip(meanX, meanY):
        preX += x
        preY += y
        sumMeanX.append(preX)
        sumMeanY.append(preY)
        ax.text(preX-0.02, preY+0.02, str(i))
        i += 1

    ax.plot(sumMeanX, sumMeanY, color='b')
    ax.plot(sumMeanX, sumMeanY, 'yo')

    plotter.tightLayout()

    return plotter


def showLocalShiftsPlot(inputSet, itemId):
    item = inputSet[itemId]
    shiftsX, shiftsY, patchX, patchY = ProtMotionCorr._getLocalShifts(item)

    # the values below are for testing
    #meanX, meanY = range(42), range(42)
    #from pyworkflow.object import Integer
    #patchX, patchY = Integer(2), Integer(3)

    plotter = createLocalAlignmentPlot(shiftsX, shiftsY, patchX, patchY)
    plotter.show()


ProjectWindow.registerObjectCommand(OBJCMD_MOVIE_ALIGNLOCAL,
                                    showLocalShiftsPlot)
#FIXME: this opens two figures for some reason... First is empty


def createLocalAlignmentPlot(shiftsX, shiftsY, patchX, patchY):
    figureSize = (8, 6)
    plotter = Plotter(*figureSize)
    cm = plt.get_cmap('winter')
    frames = 0
    fig, aN = plt.subplots(nrows=patchY, ncols=patchX, sharex=True, sharey=True)

    for rows in range(patchY):
        for cols in range(patchX):
            plotdataX = shiftsX[rows][cols].tolist()
            plotdataY = shiftsY[rows][cols].tolist()
            frames = len(plotdataX)
            ax = aN[rows][cols]
            ax.grid()
            ax.tick_params(axis='x', labelsize=8)
            ax.tick_params(axis='y', labelsize=8)
            # Choose a color map, loop through the colors, and assign them to the color
            # cycle. You need len(plotX)-1 colors, because you'll plot that many lines
            # between pairs. In other words, your line is not cyclic, so there's
            # no line from end to beginning
            ax.set_color_cycle([cm(1. * i / (frames - 1))
                                for i in range(frames - 1)])
            for i in range(frames - 1):
                ax.plot(plotdataX[i:i + 2], plotdataY[i:i + 2])

    ax1 = fig.add_subplot(111, frameon=False)
    ax1.set_title('Local patch motion (' + str(patchX) + 'x' + str(patchY) + ' patches)')
    ax1.set_xlabel('Shift x (pixels)')
    ax1.set_ylabel('Shift y (pixels)')
    ax1.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    fig.subplots_adjust(right=0.85)  # FIXME: does not work!!!

    # Add colorbar with frame numbers
    ax99 = fig.add_axes([0.9, 0.15, 0.03, 0.7])
    import matplotlib as mpl
    norm = mpl.colors.Normalize(vmin=1, vmax=frames)
    cb = mpl.colorbar.ColorbarBase(ax99, cmap=cm,
                                   boundaries=range(1, frames + 1),
                                   norm=norm,
                                   orientation='vertical')
    cb.set_label('Frame number')

    plotter.figure = fig
    return plotter
