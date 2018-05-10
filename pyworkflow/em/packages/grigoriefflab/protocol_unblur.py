# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              Josue Gomez Blanco (josue.gomez-blanco@mcgill.ca)
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

import os
from itertools import izip
from math import ceil

import pyworkflow.utils as pwutils
import pyworkflow.em as em
from pyworkflow.em.data import MovieAlignment
from pyworkflow.em.packages.xmipp3.convert import writeShiftsMovieAlignment
from pyworkflow.em.protocol import ProtAlignMovies
import pyworkflow.protocol.params as params
from pyworkflow.gui.plotter import Plotter
from grigoriefflab import UNBLUR_PATH, getVersion, UNBLUR_HOME
from convert import readShiftsMovieAlignment



class ProtUnblur(ProtAlignMovies):
    """ Unblur is used to align the frames of movies recorded on an electron
    microscope to reduce image blurring due to beam-induced motion.
    """
    _label = 'unblur'
    CONVERT_TO_MRC = 'mrc'
    CORRECT_GAIN = True

    @classmethod
    def validateInstallation(cls):
        """ Check if the installation of this protocol is correct.
        Can't rely on package function since this is a "multi package" package
        Returning an empty list means that the installation is correct
        and there are not errors. If some errors are found, a list with
        the error messages will be returned.
        """
        missingPaths = []

        if not os.path.exists(UNBLUR_PATH):
            missingPaths.append("%s : %s" % (UNBLUR_HOME, UNBLUR_PATH))
        return missingPaths

    def _defineAlignmentParams(self, form):
        group = form.addGroup('Alignment')
        line = group.addLine('Frames to ALIGN',
                             help='Frames range to ALIGN on each movie. The '
                                  'first frame is 1. If you set 0 in the final '
                                  'frame to align, it means that you will '
                                  'align until the last frame of the movie.')
        line.addParam('alignFrame0', params.IntParam, default=1,
                      label='from')
        line.addParam('alignFrameN', params.IntParam, default=0,
                      label='to')
        form.addParam('doApplyDoseFilter', params.BooleanParam, default=True,
                      label='Apply Dose filter',
                      help='Apply a dose-dependent filter to frames before '
                           'summing them. Pre-exposure and dose per frame were '
                           'specified during movies import.')

        form.addParam('doComputePSD', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Compute PSD (before/after)?",
                      help="If Yes, the protocol will compute for each movie "
                           "the average PSD before and after alignment, "
                           "for comparison")

        form.addParam('doComputeMicThumbnail', params.BooleanParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      default=False,
                      label='Compute micrograph thumbnail?',
                      help='When using this option, we will compute a '
                           'micrograph thumbnail and keep it with the '
                           'micrograph object for visualization purposes. ')

        form.addParam('openmpThreads', params.IntParam,
                      default=1, expertLevel=params.LEVEL_ADVANCED,
                      label='OpenMP threads',
                      help='The number of processors (using OpenMP threads) '
                           'that will work on each movie. Take into account '
                           'that you will be using a total of '
                           'openmpThreads * numberOfThreads (or MPIs) '
                           'processors. ')

        #group = form.addGroup('Expert Options')
        form.addParam('minShiftInitSearch', params.FloatParam,
                      default=2.,
                      label='Min. Shift Initial search (A)',
                      help='Initial search will be limited to between the '
                           'inner and outer radii',
                      expertLevel=params.LEVEL_ADVANCED)
        form.addParam('OutRadShiftLimit', params.FloatParam,
                      default=200.,
                      label='Outer radius shift limit (A)',
                      help='The maximum shift of each alignment step will be '
                           'limited to this value',
                      expertLevel=params.LEVEL_ADVANCED)
        form.addParam('bfactor', params.FloatParam,
                      default=1500.,
                      label='B-factor (A^2)',
                      help='B-factor to apply to images (A^2)',
                      expertLevel=params.LEVEL_ADVANCED)
        form.addParam('HWVertFourMask', params.IntParam,
                      default=1,
                      label='Half-width vertical Fourier mask',
                      help='The vertical line mask will be twice this size. '
                           'The central cross mask helps reduce problems by '
                           'line artefacts from the detector',
                      expertLevel=params.LEVEL_ADVANCED)
        form.addParam('HWHoriFourMask', params.IntParam,
                      default=1,
                      label='Half-width horizontal Fourier mask',
                      help='The horizontal line mask will be twice this size. '
                           'The central cross mask helps reduce problems by '
                           'line artefacts from the detector',
                      expertLevel=params.LEVEL_ADVANCED)
        form.addParam('terminShiftThreshold', params.FloatParam,
                      default=0.1,
                      label='Termination shift threshold',
                      help='Alignment will stop at this number, even if the '
                           'threshold shift is not reached',
                      expertLevel=params.LEVEL_ADVANCED)
        form.addParam('maximumNumberIterations', params.IntParam,
                      default=10,
                      label='Maximum number of iterations',
                      help='Maximum number of iterations',
                      expertLevel=params.LEVEL_ADVANCED)
        form.addParam('doRestoreNoisePwr', params.BooleanParam,
                      default=True,
                      label='Restore Noise Power? ',
                      help='Restore Noise Power? ',
                      expertLevel=params.LEVEL_ADVANCED)
        form.addParam('doVerboseOutput', params.BooleanParam,
                      default=False,
                      label='Verbose Output?',
                      help='Verbose Output?',
                      expertLevel=params.LEVEL_ADVANCED)
        form.addParallelSection(threads=1, mpi=1)

    #
    #Input stack filename                [my_movie.mrc] : kk.mrc
    #Number of frames per movie                    [34] :
    #Output aligned sum file       [my_aligned_sum.mrc] :
    #Output shifts file                 [my_shifts.txt] :
    #Pixel size of images (A)                       [1] :
    #Apply Dose filter?                            [NO] : YES
    #Exposure per frame (e/A^2)                   [1.0] :
    #Acceleration voltage (kV)                  [300.0] :

    #Pre-exposure amount(e/A^2)                   [0.0] :
    #Save Aligned Frames?                          [NO] : NO

    #Set Expert Options?                           [NO] : YES
    #Output FRC file                       [my_frc.txt] :
    #Minimum shift for initial search (Angstroms)
    #[2.0]                                              :
    #Outer radius shift limit (Angstroms)       [200.0] :
    #B-factor to apply to images (A^2)           [1500] :
    #Half-width of central vertical line of Fourier mask
    #[1]                                                :
    #Half-width of central horizontal line of Fourier mask
    #[1]                                                :
    #Termination shift threshold                  [0.1] :
    #Maximum number of iterations                  [10] :
    #Restore Noise Power?                         [YES] :
    #Verbose Output?                               [NO] : YES

    #--------------------------- STEPS functions -------------------------------

    def _processMovie(self, movie):
        numberOfFrames = self._getNumberOfFrames(movie)
        #FIXME: Figure out how to properly write shifts for unblur
        #self._writeMovieAlignment(movie, numberOfFrames)
        
        a0, aN = self._getRange(movie, 'align')
        _, lstFrame, _ = movie.getFramesRange()

        movieBaseName = pwutils.removeExt(movie.getFileName())
        aveMicFn = movieBaseName + '_uncorrected_avg.mrc'
        
        if a0 > 1 or aN < lstFrame:
            from pyworkflow.em import ImageHandler
            ih = ImageHandler()
            movieInputFn = movie.getFileName()
            
            if movieInputFn.endswith("mrc"):
                movieInputFn += ":mrcs"
            
            movieConverted = pwutils.removeExt(movieInputFn) + "_tmp.mrcs"
            ih.convertStack(movieInputFn, movieConverted, a0, aN)
            # Here, only temporal movie file (or link) stored in
            # tmp/movie_?????? is removed before move the converted file. It
            #  is necessary 'cause if it is overwritten you may lost your
            # original data.
            os.remove(movie.getFileName())
            pwutils.moveFile(movieConverted, movie.getFileName())
        
        movieSet = self.inputMovies.get()
        self._createLink(movie)
        range = aN - a0 + 1
        self._argsUnblur(movie, range)
        
        try:
            self.runJob(self._program, self._args)

            outMicFn = self._getExtraPath(self._getOutputMicName(movie))
            if not os.path.exists(outMicFn):
                # if only DW mic is saved
                outMicFn = self._getExtraPath(self._getOutputMicWtName(movie))

            if self.doComputePSD:
                # Compute uncorrected avg mic
                roi = [0, 0, 0, 0]
                fakeShiftsFn = self.writeZeroShifts(movie)
                self.averageMovie(movie, fakeShiftsFn, aveMicFn,
                                  binFactor=1,
                                  roi=roi, dark=None,
                                  gain=movieSet.getGain())

                self.computePSDs(movie, aveMicFn, outMicFn,
                                 outputFnCorrected=self._getPsdJpeg(movie))

            self._saveAlignmentPlots(movie)

            if self._doComputeMicThumbnail():
                self.computeThumbnail(outMicFn,
                                      outputFn=self._getOutputMicThumbnail(
                                          movie))
        except:
            print("ERROR: Movie %s failed\n" % movie.getFileName())

    # --------------------------- INFO functions -------------------------------
    def _citations(self):
        return ['Campbell2012', 'Grant2015b']

    def _summary(self):
        return []

    def _methods(self):
        return []

    def _validate(self):

        errors = []
        if not os.path.exists(UNBLUR_PATH):
            errors.append(
                "Cannot find the Unblur program at: " + UNBLUR_PATH)
        if self.inputMovies.get():

            if self.doApplyDoseFilter:
                inputMovies = self.inputMovies.get()
                doseFrame = inputMovies.getAcquisition().getDosePerFrame()

                if doseFrame == 0.0 or doseFrame is None:
                    errors.append('Dose per frame for input movies is 0 or not '
                                  'set. You cannot apply dose filter.')

        return errors

    #--------------------------- UTILS functions -------------------------------
    def _argsUnblur(self, movie, numberOfFrames):
        """ Format argument for call unblur program. """
        args = {'movieName': self._getMovieFn(movie),
                'numberOfFramesPerMovie': numberOfFrames,
                'micFnName': self._getMicFn(movie),
                'shiftFnName': self._getShiftsFn(movie),
                'samplingRate': self.samplingRate,
                'voltage': movie.getAcquisition().getVoltage(),
                'frcFn': self._getFrcFn(movie),
                'bfactor': self.bfactor.get(),
                'minShiftInitSearch': self.minShiftInitSearch.get(),
                'OutRadShiftLimit': self.OutRadShiftLimit.get(),
                'HWVertFourMask': self.HWVertFourMask.get(),
                'HWHoriFourMask': self.HWHoriFourMask.get(),
                'terminShiftThreshold': self.terminShiftThreshold.get(),
                'maximumNumberIterations': self.maximumNumberIterations.get(),
                'doApplyDoseFilter': 'YES' if self.doApplyDoseFilter else 'NO',
                'doRestoreNoisePwr': 'YES' if self.doRestoreNoisePwr else 'NO',
                'doVerboseOutput': 'YES' if self.doVerboseOutput else 'NO',
                'exposurePerFrame': movie.getAcquisition().getDosePerFrame() or 0.0
                }

        # Avoid threads multiplication
        # self._program = 'export OMP_NUM_THREADS=%d; ' % self.numberOfThreads.get()
        self._program = 'export OMP_NUM_THREADS=%d; ' % self.openmpThreads
        self._program += UNBLUR_PATH

        if getVersion('UNBLUR') != '1.0_150529':
            args['preExposureAmount'] = movie.getAcquisition().getDoseInitial() or 0.0
            self._args = """ << eof
%(movieName)s
%(numberOfFramesPerMovie)s
%(micFnName)s
%(shiftFnName)s
%(samplingRate)f
%(doApplyDoseFilter)s
%(exposurePerFrame)f
%(voltage)f
%(preExposureAmount)f
NO
YES
%(frcFn)s
%(minShiftInitSearch)f
%(OutRadShiftLimit)f
%(bfactor)f
%(HWVertFourMask)d
%(HWHoriFourMask)d
%(terminShiftThreshold)f
%(maximumNumberIterations)d
%(doRestoreNoisePwr)s
%(doVerboseOutput)s
eof
""" % args

        else:
            self._args = """ << eof
%(movieName)s
%(numberOfFramesPerMovie)s
%(micFnName)s
%(shiftFnName)s
%(samplingRate)f
%(doApplyDoseFilter)s
%(exposurePerFrame)f
%(voltage)f
YES
%(frcFn)s
%(minShiftInitSearch)f
%(OutRadShiftLimit)f
%(bfactor)f
%(HWVertFourMask)d
%(HWHoriFourMask)d
%(terminShiftThreshold)f
%(maximumNumberIterations)d
%(doRestoreNoisePwr)s
%(doVerboseOutput)s
eof
""" % args
    
    def _getMicName(self, movieName):
        """ Return the name for the output micrograph given the movie name.
        """
        return self._getExtraPath(self._getNameExt(movieName, '', 'mrc'))

    def _getMovieFn(self, movie):
        movieFn = movie.getFileName()
        if movieFn.endswith("mrcs"):
            return pwutils.replaceExt(movieFn, self.CONVERT_TO_MRC)
        else:
            return movieFn

    def _createLink(self, movie):
        movieFn = movie.getFileName()
        if movieFn.endswith("mrcs"):
            pwutils.createLink(movieFn, self._getMovieFn(movie))

    def _getMicFn(self, movie):
        return self._getExtraPath(self._getOutputMicName(movie))

    def _getShiftsFn(self, movie):
        return self._getExtraPath(self._getMovieRoot(movie) + '_shifts.txt')

    def _getFrcFn(self, movie):
        return self._getExtraPath(self._getMovieRoot(movie) + '_frc.txt')

    def _getMovieShifts(self, movie):
        """ Returns the x and y shifts for the alignment of this movie.
        """
        pixSize = movie.getSamplingRate()
        shiftFn = self._getShiftsFn(movie)
        xShifts, yShifts = readShiftsMovieAlignment(shiftFn)
        # convert shifts from Angstroms to px
        # (e.g. Summovie requires shifts in px)
        xShiftsCorr = [x / pixSize for x in xShifts]
        yShiftsCorr = [y / pixSize for y in yShifts]

        return xShiftsCorr, yShiftsCorr
    
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
            first -= (iniFrame - 1)
            last -= (iniFrame - 1)

        return first, last

    def _isNewUnblur(self):
        return True if getVersion('UNBLUR') != '1.0.150529' else False

    def _doComputeMicThumbnail(self):
        return self.doComputeMicThumbnail.get()

    def _preprocessOutputMicrograph(self, mic, movie):
        if self.doComputePSD:
            mic.psdCorr = em.Image(location=self._getPsdCorr(movie))
            mic.psdJpeg = em.Image(location=self._getPsdJpeg(movie))
        mic.plotGlobal = em.Image(location=self._getPlotGlobal(movie))
        if self._doComputeMicThumbnail():
            mic.thumbnail = em.Image(location=self._getOutputMicThumbnail(movie))

    def _getNameExt(self, movie, postFix, ext, extra=False):
        fn = self._getMovieRoot(movie) + postFix + '.' + ext
        return self._getExtraPath(fn) if extra else fn

    def _getPlotGlobal(self, movie):
        return self._getNameExt(movie, '_global_shifts', 'png', extra=True)

    def _getPsdCorr(self, movie):
        return self._getNameExt(movie, '_psd_comparison', 'psd', extra=True)

    def _getPsdJpeg(self, movie):
        return self._getNameExt(movie, '_psd', 'jpeg', extra=True)

    def _saveAlignmentPlots(self, movie):
        """ Compute alignment shift plots and save to file as png images. """
        shiftsX, shiftsY = self._getMovieShifts(movie)
        first, _ = self._getFrameRange(movie.getNumberOfFrames(), 'align')
        plotter = createGlobalAlignmentPlot(shiftsX, shiftsY, first)
        plotter.savefig(self._getPlotGlobal(movie))

def createGlobalAlignmentPlot(meanX, meanY, first):
    """ Create a plotter with the shift per frame. """
    sumMeanX = []
    sumMeanY = []
    preX = 0.0
    preY = 0.0

    figureSize = (6, 4)
    plotter = Plotter(*figureSize)
    figure = plotter.getFigure()
    ax = figure.add_subplot(111)
    ax.grid()
    ax.set_title('Alignment based upon full frames')
    ax.set_xlabel('Shift x (pixels)')
    ax.set_ylabel('Shift y (pixels)')
    
    i = first
    skipLabels = ceil(len(meanX)/10.0)
    labelTick = 1

    for x, y in izip(meanX, meanY):
        preX = x
        preY = y
        sumMeanX.append(preX)
        sumMeanY.append(preY)
        if labelTick == 1:
            ax.text(preX - 0.02, preY + 0.02, str(i))
            labelTick = skipLabels
        else:
            labelTick -= 1
        i += 1

    ax.plot(sumMeanX, sumMeanY, color='b')
    ax.plot(sumMeanX, sumMeanY, 'yo')

    plotter.tightLayout()

    return plotter

