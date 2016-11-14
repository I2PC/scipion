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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os

from pyworkflow.object import Set
import pyworkflow.utils.path as pwutils
import pyworkflow.protocol.params as params
import pyworkflow.protocol.constants as cons
from pyworkflow.em.convert import ImageHandler
from pyworkflow.em.data import MovieAlignment, SetOfMovies, SetOfMicrographs
from pyworkflow.em.protocol import ProtProcessMovies


class ProtAlignMovies(ProtProcessMovies):
    """
    Base class for movie alignment protocols such as:
    motioncorr, crosscorrelation and optical flow

    Alignment parameters are defined in common. For example,
    the frames range used for alignment and final sum, the binning factor
    or the cropping options (region of interest)
    """

    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        ProtProcessMovies._defineParams(self, form)
        self._defineAlignmentParams(form)

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
        group.addParam('useAlignToSum', params.BooleanParam, default=True,
                       label='Use ALIGN frames range to SUM?',
                       help="If *Yes*, the same frame range will be used to "
                            "ALIGN and to SUM. If *No*, you can selected a "
                            "different range for SUM (must be a subset).")
        line = group.addLine('Frames to SUM', condition="not useAlignToSum",
                             help='Frames range to SUM on each movie. The '
                                  'first frame is 1. If you set 0 in the final '
                                  'frame to sum, it means that you will sum '
                                  'until the last frame of the movie.')
        line.addParam('sumFrame0', params.IntParam, default=1,
                      label='from')
        line.addParam('sumFrameN', params.IntParam, default=0,
                      label='to')
        group.addParam('binFactor', params.FloatParam, default=1.,
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
                      label="Save aligned micrograph",
                      expertLevel=cons.LEVEL_ADVANCED)

        form.addParam('doSaveMovie', params.BooleanParam, default=False,
                      label="Save movie", expertLevel=cons.LEVEL_ADVANCED,
                      help="Save Aligned movie")

    # --------------------------- STEPS functions ----------------------------

    # FIXME: Methods will change when using the streaming for the output
    def createOutputStep(self):
        # Do nothing now, the output should be ready.
        pass

    def _loadOutputSet(self, SetClass, baseName, fixSampling=True):
        """
        Load the output set if it exists or create a new one.
        fixSampling: correct the output sampling rate if binning was used,
        except for the case when the original movies are kept and shifts
        refers to that one.
        """
        setFile = self._getPath(baseName)

        if os.path.exists(setFile):
            outputSet = SetClass(filename=setFile)
            outputSet.loadAllProperties()
            outputSet.enableAppend()
        else:
            outputSet = SetClass(filename=setFile)
            outputSet.setStreamState(outputSet.STREAM_OPEN)

        inputMovies = self.inputMovies.get()
        outputSet.copyInfo(inputMovies)

        if fixSampling:
            newSampling = inputMovies.getSamplingRate() * self._getBinFactor()
            outputSet.setSamplingRate(newSampling)

        return outputSet

    def _checkNewOutput(self):
        if getattr(self, 'finished', False):
            return

        # Load previously done items (from text file)
        doneList = self._readDoneList()
        # Check for newly done items
        newDone = [m for m in self.listOfMovies
                   if m.getObjId() not in doneList and self._isMovieDone(m)]

        # Update the file with the newly done movies
        # or exit from the function if no new done movies
        self.debug('_checkNewOutput: ')
        self.debug('   listOfMovies: %s, doneList: %s, newDone: %s'
                   % (len(self.listOfMovies), len(doneList), len(newDone)))

        firstTime = len(doneList) == 0
        allDone = len(doneList) + len(newDone)
        # We have finished when there is not more input movies (stream closed)
        # and the number of processed movies is equal to the number of inputs
        self.finished = self.streamClosed and allDone == len(self.listOfMovies)
        streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN

        if newDone:
            self._writeDoneList(newDone)
        elif not self.finished:
            # If we are not finished and no new output have been produced
            # it does not make sense to proceed and updated the outputs
            # so we exit from the function here
            return

        self.debug('   finished: %s ' % self.finished)
        self.debug('        self.streamClosed (%s) AND' % self.streamClosed)
        self.debug('        allDone (%s) == len(self.listOfMovies (%s)'
                   % (allDone, len(self.listOfMovies)))
        self.debug('   streamMode: %s' % streamMode)

        if self._createOutputMovies():
            saveMovie = self.getAttributeValue('doSaveMovie', False)
            movieSet = self._loadOutputSet(SetOfMovies, 'movies.sqlite',
                                           fixSampling=saveMovie)

            for movie in newDone:
                newMovie = self._createOutputMovie(movie)
                movieSet.append(newMovie)

            self._updateOutputSet('outputMovies', movieSet, streamMode)

            if firstTime:
                # Probably is a good idea to store a cached summary for the
                # first resulting movie of the processing.
                self._storeSummary(newDone[0])
                self._defineTransformRelation(self.inputMovies, movieSet)

        if self._createOutputMicrographs():
            micSet = self._loadOutputSet(SetOfMicrographs, 'micrographs.sqlite')

            for movie in newDone:
                mic = micSet.ITEM_TYPE()
                mic.copyObjId(movie)
                mic.setMicName(movie.getMicName())
                # The subclass protocol is responsible of generating the output
                # micrograph file in the extra path with the required name
                extraMicFn = self._getExtraPath(self._getOutputMicName(movie))
                mic.setFileName(extraMicFn)
                self._preprocessOutputMicrograph(mic, movie)
                micSet.append(mic)

            self._updateOutputSet('outputMicrographs', micSet, streamMode)

            if firstTime:
                # We consider that Movies are 'transformed' into the Micrographs
                # This will allow to extend the CTF associated to a set of
                # micrographs to another set of micrographs generated from a
                # different movie alignment
                self._defineTransformRelation(self.inputMovies, micSet)

        if self._createOutputWeightedMicrographs():

            micSet2 = self._loadOutputSet(SetOfMicrographs,
                                          'micrographs_dose-weighted.sqlite')

            for movie in newDone:
                mic2 = micSet2.ITEM_TYPE()
                mic2.copyObjId(movie)
                mic2.setMicName(movie.getMicName())
                # The subclass protocol is responsible of generating the output
                # micrograph file in the extra path with the required name
                extraMicFn2 = self._getExtraPath(self._getOutputMicWtName(movie))
                mic2.setFileName(extraMicFn2)
                self._preprocessOutputMicrograph(mic2, movie)
                # FIXME The micSet is not setting properly dimensions (No-Dim)
                micSet2.append(mic2)

            self._updateOutputSet('outputMicrographsDoseWt',
                                  micSet2, streamMode)

            if firstTime:
                self._defineTransformRelation(self.inputMovies, micSet2)

        if self.finished:  # Unlock createOutputStep if finished all jobs
            outputStep = self._getFirstJoinStep()
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(cons.STATUS_NEW)

    # --------------------------- INFO functions --------------------------------

    def _validate(self):
        errors = []

        if (self.cropDimX > 0 and self.cropDimY <= 0 or
                        self.cropDimY > 0 and self.cropDimX <= 0):
            errors.append("If you give cropDimX, you should also give cropDimY"
                          " and vice versa")

        # movie = self.inputMovies.get().getFirstItem()
        # # Close movies db because the getFirstItem open it
        # # we do not want to leave the file open
        # self.inputMovies.get().close()
        # frames = movie.getNumberOfFrames()

        firstFrame, lastFrame, _ = self.inputMovies.get().getFramesRange()
        if lastFrame == 0:
            # Although getFirstItem is not remonended in general, here it is
            # used olny once, for validation purposes, so performance
            # problems not should be apprear.
            frames = self.inputMovies.get().getFirstItem().getNumberOfFrames()
            lastFrame = frames
        else:
            frames = lastFrame - firstFrame + 1

        if frames is not None:
            def _validateRange(prefix):
                # Avoid validation when the range is not defined
                if not hasattr(self, '%sFrame0' % prefix):
                    return

                f0, fN = self._getFrameRange(frames, prefix)
                if fN < firstFrame or fN > lastFrame:
                    errors.append("Check the selected last frame to *%s*. "
                                  "Last frame (%d) should be in range: %s "
                                  % (prefix.upper(), fN, (firstFrame,
                                                          lastFrame)))
                if f0 < firstFrame or f0 > lastFrame:
                    errors.append("Check the selected first frame to *%s*. "
                                  "First frame (%d) should be in range: %s "
                                  % (prefix.upper(), f0, (firstFrame,
                                                          lastFrame)))
                if fN < f0:
                    errors.append("Check the selected frames range to *%s*. "
                                  "Last frame (%d) should be greater or equal "
                                  "than first frame (%d)"
                                  % (prefix.upper(), fN, f0))

            _validateRange("align")
            _validateRange("sum")

        return errors

    # --------------------------- INFO functions -------------------------------

    def _summary(self):
        return [self.summaryVar.get('')]

    # --------------------------- UTILS functions ----------------------------
    def _useAlignToSum(self):
        return self.getAttributeValue('useAlignToSum', False)

    def _getFrameRange(self, n, prefix):
        """
        Params:
        :param n: Number of frames of the movies
        :param prefix: what range we want to consider, either 'align' or 'sum'
        :return: (i, f) initial and last frame range
        """
        # In case that the user select the same range for ALIGN and SUM
        # we also use the 'align' prefix
        if self._useAlignToSum():
            prefix = 'align'

        first = self.getAttributeValue('%sFrame0' % prefix)
        last = self.getAttributeValue('%sFrameN' % prefix)

        if first <= 1:
            first = 1

        if last <= 0:
            last = n

        return first, last

    def _createOutputMovie(self, movie):
        movieId = movie.getObjId()

        # Parse the alignment parameters and store the log files
        alignedMovie = movie.clone()
        n = movie.getNumberOfFrames()
        first, last = self._getFrameRange(n, 'align')
        framesRange = alignedMovie.getFramesRange()
        framesRange.setFirstFrame(first)
        framesRange.setLastFrame(last)
        # Check if user selected to save movie, use the getAttributeValue
        # function for allow the protocol to not define this flag
        # and use False as default
        if self.getAttributeValue('doSaveMovie', False):
            # The subclass protocol is responsible of generating the output
            # movie file in the extra path with the required name
            extraMovieFn = self._getExtraPath(self._getOutputMovieName(movie))
            alignedMovie.setFileName(extraMovieFn)
            # When the output movies are saved, the shifts
            # will be set to zero since they are aligned
            totalFrames = last - first + 1
            xshifts = [0] * totalFrames
            yshifts = xshifts
            # If we save the movies, we need to modify which are the index
            # of the first frame in the stack, now is 1 since the stack is
            # written only with the given frames
            firstFrameIndex = 1
        else:
            xshifts, yshifts = self._getMovieShifts(movie)
            firstFrameIndex = first

        framesRange.setFirstFrameIndex(firstFrameIndex)
        alignment = MovieAlignment(first=first, last=last, xshifts=xshifts,
                                   yshifts=yshifts)

        roiList = [self.getAttributeValue(s, 0) for s in
                   ['cropOffsetX', 'cropOffsetY', 'cropDimX', 'cropDimY']]
        alignment.setRoi(roiList)
        alignedMovie.setAlignment(alignment)

        return alignedMovie

    # ---------- Hook functions that need to be implemented in subclasses ------

    def _getBinFactor(self):
        return self.getAttributeValue('binFactor', 1.0)

    def _getMovieRoot(self, movie):
        return pwutils.removeBaseExt(movie.getFileName())

    def _getOutputMovieName(self, movie):
        """ Returns the name of the output movie.
        (relative to micFolder)
        """
        return self._getMovieRoot(movie) + '_aligned_movie.mrcs'

    def _getOutputMicName(self, movie):
        """ Returns the name of the output micrograph
        (relative to micFolder)
        """
        return self._getMovieRoot(movie) + '_aligned_mic.mrc'

    def _getOutputMicWtName(self, movie):
        """ Returns the name of the output dose-weighted micrograph
        (relative to micFolder)
        """
        return self._getMovieRoot(movie) + '_aligned_mic_DW.mrc'

    def _getMovieShifts(self, movie):
        """ Returns the x and y shifts for the alignment of this movie.
         The shifts should refer to the original micrograph without any binning.
         In case of a binning greater than 1, the shifts should be scaled.
        """
        return [], []

    def _createOutputMovies(self):
        """ Returns True if an output set of movies will be generated.
        The most common case is to always generate output movies,
        either with alignment only or the binary aligned movie files.
        Subclasses can override this function to change this behavior.
        """
        return True

    def _createOutputMicrographs(self):
        """ By default check if the user have selected 'doSaveAveMic'
        property. Subclasses can override this method to implement different
        behaviour.
        """
        return self.getAttributeValue('doSaveAveMic', True)

    def _createOutputWeightedMicrographs(self):
        return False

    def _preprocessOutputMicrograph(self, mic, movie):
        """ Hook function that will be call before adding the micrograph
        to the output set of micrographs.
        """
        pass

    def _storeSummary(self, movie):
        """ Implement this method if you want to store the summary. """
        pass

    def __runXmippProgram(self, program, args):
        """ Internal shortcut function to launch a Xmipp program. """
        import pyworkflow.em.packages.xmipp3 as xmipp3
        xmipp3.runXmippProgram(program, args)

    def averageMovie(self, movie, inputFn, outputMicFn, binFactor=1, roi=None,
                     dark=None, gain=None, splineOrder=None):
        """ Average a movie (using xmipp) taking into account the
         possible shifts and other alignment parameters.
         Params:
            inputFn: input filename, either the movie file or a metadata
                with the shifts and other info.
            dark: dark file
            gain: gain correction file.

         The output will be the averaged micrograph.
        """
        args = '-i %s ' % inputFn
        args += '--sampling %f ' % movie.getSamplingRate()
        args += '--useInputShifts '

        if binFactor > 1:
            args += '--bin %f ' % binFactor

        if roi is not None:
            x, y, _ = movie.getDim()
            offsetX, offsetY, cropDimX, cropDimY = roi
            # cropDim value is <= 0 we should take the whole size
            if cropDimX <= 0:
                dimX = x - 1
            else:
                dimX = offsetX + cropDimX - 1

            if cropDimY <= 0:
                dimY = y - 1
            else:
                dimY = offsetY + cropDimY - 1

            args += '--cropULCorner %d %d ' % (offsetX, offsetY)
            args += '--cropDRCorner %d %d ' % (dimX, dimY)

        args += ' --oavg %s ' % outputMicFn

        if dark is not None:
            args += ' --dark ' + dark

        if gain is not None:
            args += ' --gain ' + gain

        if splineOrder is not None:
            args += '--Bspline %d ' % splineOrder

        self.__runXmippProgram('xmipp_movie_alignment_correlation', args)

    def computePSD(self, inputMic, oroot, dim=400, overlap=0.7):
        args = '--micrograph %s --oroot %s ' % (inputMic, oroot)
        args += '--dont_estimate_ctf --pieceDim %d --overlap %f' % (dim, overlap)

        self.__runXmippProgram('xmipp_ctf_estimate_from_micrograph', args)

    def composePSD(self, psd1, psd2, outputFn):
        """ Compose a single PSD image:
         left part from psd1 (corrected PSD),
         right-part from psd2 (uncorrected PSD)
        """
        ih = ImageHandler()
        psd = ih.read(psd1)
        data1 = psd.getData()
        data2 = ih.read(psd2).getData()
        # Compute middle index
        x, _, _, _ = psd.getDimensions()
        m = int(round(x / 2.))
        data1[:, m:] = data2[:, m:]
        psd.setData(data1)
        psd.write(outputFn)
