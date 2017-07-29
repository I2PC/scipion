# ******************************************************************************
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
# ******************************************************************************

import os
from itertools import izip

import pyworkflow.protocol.params as params
import pyworkflow.protocol.constants as cons
import pyworkflow.utils as pwutils
import pyworkflow.em as em
from pyworkflow import VERSION_1_1
from pyworkflow.em.data import MovieAlignment
from pyworkflow.em.packages.xmipp3.convert import writeShiftsMovieAlignment
from pyworkflow.em.packages.grigoriefflab.convert import (parseMagCorrInput,
                                                          parseMagEstOutput)
from pyworkflow.em.protocol import ProtAlignMovies
from pyworkflow.gui.plotter import Plotter
from convert import (MOTIONCORR_PATH, MOTIONCOR2_PATH, getVersion, getEnviron,
                     parseMovieAlignment, parseMovieAlignment2, getCudaLib,
                     MOTIONCORR_CUDA_LIB, MOTIONCOR2_CUDA_LIB, CUDA_LIB)
from pyworkflow.protocol import STEPS_PARALLEL


class ProtMotionCorr(ProtAlignMovies):
    """
    Wrapper protocol to movie alignment programs developed at UCSF:
    motioncorr: Flat fielding and Drift correction
        (written by Xueming Li @ Yifan Cheng Lab)
    motioncor2: anisotropic drift correction and dose weighting
        (written by Shawn Zheng @ David Agard lab)
    """

    _label = 'motioncorr alignment'
    _lastUpdateVersion = VERSION_1_1
    CONVERT_TO_MRC = 'mrc'

    def __init__(self, **args):
        ProtAlignMovies.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL

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
                           "the average PSD before and after alignment, "
                           "for comparison")

        form.addParam('doComputeMicThumbnail', params.BooleanParam,
                      expertLevel=cons.LEVEL_ADVANCED,
                      default=False, condition='doSaveAveMic',
                      label='Compute micrograph thumbnail?',
                      help='When using this option, we will compute a '
                           'micrograph thumbnail and keep it with the '
                           'micrograph object for visualization purposes. ')

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
        form.addParam('doApplyDoseFilter', params.BooleanParam, default=True,
                      condition='useMotioncor2',
                      label='Apply Dose filter',
                      help='Apply a dose-dependent filter to frames before '
                           'summing them. Pre-exposure and dose per frame '
                           'should  be specified during movies import.')

        line = form.addLine('Number of patches', condition='useMotioncor2',
                            help='Number of patches to be used for patch based '
                                 'alignment. Set to *0 0* to do only global motion '
                                 'correction.')
        line.addParam('patchX', params.IntParam, default=5, label='X')
        line.addParam('patchY', params.IntParam, default=5, label='Y')

        form.addParam('group', params.IntParam, default='1',
                      label='Group N frames', condition='useMotioncor2',
                      help='Group every specified number of frames by adding '
                           'them together. The alignment is then performed on '
                           'the summed frames. By default, no grouping is '
                           'performed.')
        form.addParam('tol', params.FloatParam, default='0.5',
                      label='Tolerance (px)', condition='useMotioncor2',
                      help='Tolerance for iterative alignment, default *0.5px*.')
        if self._supportsMagCorrection():
            group = form.addGroup('Magnification correction')
            group.addParam('doMagCor', params.BooleanParam, default=False,
                           label='Correct anisotropic magnification?',
                           condition='useMotioncor2',
                           help='Correct anisotropic magnification by '
                                'stretching image along the major axis, '
                                'the axis where the lower magnification is '
                                'detected.')
            group.addParam('useEst', params.BooleanParam, default=True,
                           label='Use previous estimation?',
                           condition='useMotioncor2 and doMagCor',
                           help='Use previously calculated parameters of '
                                'magnification anisotropy (from magnification '
                                'distortion estimation protocol).')
            group.addParam('inputEst', params.PointerParam,
                           pointerClass='ProtMagDistEst',
                           condition='useEst and useMotioncor2 and doMagCor',
                           label='Input protocol',
                           help='Select previously executed estimation protocol.')
            group.addParam('scaleMaj', params.FloatParam, default=1.0,
                           condition='not useEst and useMotioncor2 and doMagCor',
                           label='Major scale factor',
                           help='Major scale factor.')
            group.addParam('scaleMin', params.FloatParam, default=1.0,
                           condition='not useEst and useMotioncor2 and doMagCor',
                           label='Minor scale factor',
                           help='Minor scale factor.')
            group.addParam('angDist', params.FloatParam, default=0.0,
                           condition='not useEst and useMotioncor2 and doMagCor',
                           label='Distortion angle (deg)',
                           help='Distortion angle, in degrees.')
        else:
            form.addParam('motioncor2Version', params.LabelParam,
                          condition='useMotioncor2',
                          label='Scipion supports some versions of motioncor2 '
                                'that can do magnification correction, '
                                'but they do not seems to be installed. Check '
                                'available versions with: \n'
                                'scipion install --help.\n'
                                'Also, make sure MOTIONCOR2_CUDA_LIB or '
                                'CUDA_LIB point to cuda-8.0/lib path')

        if getVersion('MOTIONCOR2') in ['1.0.0']:
            form.addParam('defectFile', params.FileParam,
                          expertLevel=cons.LEVEL_ADVANCED,
                          condition='useMotioncor2',
                          label='Camera defects file',
                          help='Defect file that stores entries of defects on camera.\n'
                               'Each entry corresponds to a rectangular region in image. '
                               'The pixels in such a region are replaced by '
                               'neighboring good pixel values. Each entry contains '
                               '4 integers x, y, w, h representing the x, y '
                               'coordinates, width, and height, respectively.')

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
        -FmRef     -1         Specify which frame to be the reference to which
                              all other frames are aligned, by default (-1) the
                              the central frame is chosen. The central frame is
                              at N/2 based upon zero indexing where N is the
                              number of frames that will be summed, i.e., not
                              including the frames thrown away.
        -RotGain   0          Rotate gain reference counter-clockwise: 0 - no rotation,
                              1 - 90 degrees, 2 - 180 degrees, 3 - 270 degrees.
        -FlipGain  0          Flip gain reference after gain rotation: 0 - no flipping,
                              1 - flip upside down, 2 - flip left right.
        -Tilt      0 0        Tilt angle range for a dose fractionated tomographic
                              tilt series, e.g. *-60 60*
                              """)
        form.addParam('doSaveUnweightedMic', params.BooleanParam, default=True,
                      condition='doSaveAveMic and useMotioncor2 and doApplyDoseFilter',
                      label="Save unweighted micrographs?",
                      help="Yes by default, if you have selected to apply a "
                           "dose-dependent filter to the frames")

        # Since only runs on GPU, do not allow neither threads nor mpi
        form.addParallelSection(threads=1, mpi=1)

    #--------------------------- STEPS functions -------------------------------
    def _processMovie(self, movie):
        inputMovies = self.inputMovies.get()
        movieFolder = self._getOutputMovieFolder(movie)
        outputMicFn = self._getRelPath(self._getOutputMicName(movie),
                                       movieFolder)
        outputMovieFn = self._getRelPath(self._getOutputMovieName(movie),
                                         movieFolder)
        movieBaseName = pwutils.removeExt(movie.getFileName())
        aveMicFn =  movieBaseName + '_uncorrected_avg.mrc'
        logFile = self._getRelPath(self._getMovieLogFile(movie),
                                   movieFolder)

        a0, aN = self._getRange(movie, 'align')

        if not self.useMotioncor2:
            # Get the number of frames and the range to be used
            # for alignment and sum
            s0, sN = self._getRange(movie, 'sum')

            argsDict = {'-crx': self.cropOffsetX.get(),
                        '-cry': self.cropOffsetY.get(),
                        '-cdx': self.cropDimX.get(),
                        '-cdy': self.cropDimY.get(),
                        '-bin': self.binFactor.get(),
                        '-nst': '%d' % a0,
                        '-ned': '%d' % aN,
                        '-nss': '%d' % s0,
                        '-nes': '%d' % sN,
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
            logFileBase = (logFile.replace('0-Full.log', '').replace(
                '0-Patch-Full.log', ''))
            # default values for motioncor2 are (1, 1)
            cropDimX = self.cropDimX.get() or 1
            cropDimY = self.cropDimY.get() or 1

            numbOfFrames = self._getNumberOfFrames(movie)

            if self.doApplyDoseFilter:
                preExp, dose = self._getCorrectedDose(inputMovies)
            else:
                preExp, dose = 0.0, 0.0

            argsDict = {'-OutMrc': '"%s"' % outputMicFn,
                        '-Patch': '%d %d' % (self.patchX, self.patchY),
                        '-MaskCent': '%d %d' % (self.cropOffsetX,
                                                self.cropOffsetY),
                        '-MaskSize': '%d %d' % (cropDimX, cropDimY),
                        '-FtBin': self.binFactor.get(),
                        '-Tol': self.tol.get(),
                        '-Group': self.group.get(),
                        '-FmDose': dose,
                        '-Throw': '%d' % a0,
                        '-Trunc': '%d' % (abs(aN - numbOfFrames + 1)),
                        '-PixSize': inputMovies.getSamplingRate(),
                        '-kV': inputMovies.getAcquisition().getVoltage(),
                        '-Gpu': self.GPUIDs.get(),
                        '-LogFile': logFileBase,
                        }
            if getVersion('MOTIONCOR2') != '03162016':
                argsDict['-InitDose'] = preExp
                argsDict['-OutStack'] = 1 if self.doSaveMovie else 0

            if self._supportsMagCorrection() and self.doMagCor:
                if self.useEst:
                    inputEst = self.inputEst.get().getOutputLog()
                    if getVersion() == '01302017':
                        input_params = parseMagCorrInput(inputEst)
                        # this version uses stretch parameters as following:
                        # 1/maj, 1/min, -angle
                        argsDict['-Mag'] = '%0.2f %0.2f %0.2f' % (1.0 / input_params[1],
                                                                  1.0 / input_params[2],
                                                                  -1 * input_params[0])
                    else:
                        # While motioncor2 >=1.0.0 uses estimation params  AS IS
                        input_params = parseMagEstOutput(inputEst)
                        argsDict['-Mag'] = '%0.2f %0.2f %0.2f' % (input_params[1],
                                                                  input_params[2],
                                                                  input_params[0])
                else:
                    argsDict['-Mag'] = '%0.2f %0.2f %0.2f' % (self.scaleMaj.get(),
                                                              self.scaleMin.get(),
                                                              self.angDist.get())

            args = ' -InMrc "%s" ' % movie.getBaseName()
            args += ' '.join(['%s %s' % (k, v) for k, v in argsDict.iteritems()])

            if inputMovies.getGain():
                args += ' -Gain "%s" ' % inputMovies.getGain()

            args += ' ' + self.extraParams2.get()
            program = MOTIONCOR2_PATH

        try:
            self.runJob(program, args, cwd=movieFolder,
                        env=getEnviron(self.useMotioncor2))
            self._fixMovie(movie)

            # Compute PSDs
            outMicFn = self._getExtraPath(self._getOutputMicName(movie))
            if not os.path.exists(outMicFn):
                # if only DW mic is saved
                outMicFn = self._getExtraPath(self._getOutputMicWtName(movie))

            if self.doComputePSD:
                uncorrectedPSD = movieBaseName + '_uncorrected'
                correctedPSD = movieBaseName + '_corrected'
                # Compute uncorrected avg mic
                roi = [self.cropOffsetX.get(), self.cropOffsetY.get(),
                       self.cropDimX.get(), self.cropDimY.get()]
                fakeShiftsFn = self.writeZeroShifts(movie)
                self.averageMovie(movie, fakeShiftsFn, aveMicFn,
                                  binFactor=self.binFactor.get(),
                                  roi=roi, dark=None,
                                  gain=inputMovies.getGain())

                self.computePSD(aveMicFn, uncorrectedPSD)
                self.computePSD(outMicFn, correctedPSD)
                self.composePSD(uncorrectedPSD + ".psd",
                                correctedPSD + ".psd",
                                self._getPsdCorr(movie))

            self._saveAlignmentPlots(movie)

            if self._doComputeMicThumbnail():
                self.computeThumbnail(outMicFn,
                                      outputFn=self._getOutputMicThumbnail(movie))
        except:
            print("ERROR: Movie %s failed\n" % movie.getName())

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

        # Check CUDA paths
        cudaLib = getCudaLib(useMC2=self.useMotioncor2)
        cudaConst = (MOTIONCOR2_CUDA_LIB if self.useMotioncor2 else
                     MOTIONCORR_CUDA_LIB)
        
        if cudaLib is None:
            errors.append("Do not know where to find CUDA lib path. "
                          " %s or %s variables have None value or are not"
                          " present in scipion configuration."
                          % (cudaConst, CUDA_LIB))

        elif not pwutils.existsVariablePaths(cudaLib):
            errors.append("Either %s or %s variables points to a non existing "
                          "path (%s). Please, check scipion configuration."
                          % (cudaConst, CUDA_LIB, cudaLib))

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
                              'Save aligned micrograph. '
                              'Optionally you could add -Align 0 to additional '
                              'parameters so that protocol '
                              'produces simple movie sum.')

            if self.doSaveMovie and not self._isOutStackSupport:
                errors.append('Saving aligned movies is not supported by '
                              'this version of motioncor2. '
                              'By default, the protocol will produce '
                              'outputMovies equivalent to the input '
                              'however containing alignment information.')

            if not self.useAlignToSum:
                errors.append('Frame range for ALIGN and SUM must be '
                              'equivalent in case of motioncor2. \n Please, '
                              'set *YES* _Use ALIGN frames range to SUM?_ '
                              'flag or use motioncorr')

            if self.doApplyDoseFilter and self.inputMovies.get():
                inputMovies = self.inputMovies.get()
                doseFrame = inputMovies.getAcquisition().getDosePerFrame()

                if doseFrame == 0.0 or doseFrame is None:
                    errors.append('Dose per frame for input movies is 0 or not '
                                  'set. You cannot apply dose filter.')

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
        return self._getExtraPath(fn) if extra else fn

    def _getPlotGlobal(self, movie):
        return self._getNameExt(movie, '_global_shifts', 'png', extra=True)

    def _getPsdCorr(self, movie):
        return self._getNameExt(movie, '_psd_comparison', 'psd', extra=True)

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

    def _setPlotInfo(self, movie, mic):
        mic.plotGlobal = em.Image(location=self._getPlotGlobal(movie))
        if self.doComputePSD:
            mic.psdCorr = em.Image(location=self._getPsdCorr(movie))
        if self._doComputeMicThumbnail():
            mic.thumbnail = em.Image(location=self._getOutputMicThumbnail(movie))

    def _saveAlignmentPlots(self, movie):
        """ Compute alignment shift plots and save to file as png images. """
        shiftsX, shiftsY = self._getMovieShifts(movie)
        first, _ = self._getFrameRange(movie.getNumberOfFrames(), 'align')
        plotter = createGlobalAlignmentPlot(shiftsX, shiftsY, first)
        plotter.savefig(self._getPlotGlobal(movie))

    def _isOutStackSupport(self):
        # checks if output aligned movies can be saved by motioncor2
        return True if getVersion('MOTIONCOR2') != '03162016' else False

    def _fixMovie(self, movie):
        if self.doSaveMovie and self.useMotioncor2 and self._isOutStackSupport():
            outputMicFn = self._getExtraPath(self._getOutputMicName(movie))
            outputMovieFn = self._getExtraPath(self._getOutputMovieName(movie))
            movieFn = outputMicFn.replace('_aligned_mic.mrc',
                                          '_aligned_mic_Stk.mrc')
            pwutils.moveFile(movieFn, outputMovieFn)

        if self.useMotioncor2 and not self.doSaveUnweightedMic:
            fnToDelete = self._getExtraPath(self._getOutputMicName(movie))
            pwutils.cleanPath(fnToDelete)

    def writeZeroShifts(self, movie):
        # TODO: find another way to do this
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

    def _getRange(self, movie, prefix):

        n = self._getNumberOfFrames(movie)
        iniFrame, _, indxFrame = movie.getFramesRange()
        first, last = self._getFrameRange(n, prefix)

        if iniFrame != indxFrame:
            first -= iniFrame
            last -= iniFrame
        else:
            first -= 1
            last -= 1

        return first, last

    def _getNumberOfFrames(self, movie):
        _, lstFrame, _ = movie.getFramesRange()

        if movie.hasAlignment():
            _, lastFrmAligned = movie.getAlignment().getRange()
            if lastFrmAligned != lstFrame:
                return lastFrmAligned
            else:
                return movie.getNumberOfFrames()
        else:
            return movie.getNumberOfFrames()

    def _createOutputMicrographs(self):
        createWeighted = self._createOutputWeightedMicrographs()
        # To create the unweighted average micrographs
        # we only consider the 'doSaveUnweightedMic' flag if the
        # weighted ones should be created.
        return (self.doSaveAveMic and
                (not createWeighted or self.doSaveUnweightedMic))

    def _createOutputWeightedMicrographs(self):
        return (self.doSaveAveMic and self.useMotioncor2 and
                self.doApplyDoseFilter)

    def _doComputeMicThumbnail(self):
        return (self.doSaveAveMic and self.doComputeMicThumbnail)
    
    def _supportsMagCorrection(self):
        return getVersion('MOTIONCOR2') not in ['03162016', '10192016']


def createGlobalAlignmentPlot(meanX, meanY, first):
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

    i = first
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
