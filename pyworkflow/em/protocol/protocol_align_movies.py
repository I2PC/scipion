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
    motioncorr, crosscrorrelation and optical flow

    Alignment parameters are defined in common. For example,
    the frames range used for alignment and final sum, the binning factor
    or the cropping options (region of interest)
    """
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        ProtProcessMovies._defineParams(self, form)
        self._defineAlignmentParams(form)

    def _defineAlignmentParams(self, form):
        group = form.addGroup('Alignment')
        line = group.addLine('Remove frames to ALIGN from',
                            help='How many frames remove'
                                 ' from movie alignment.')
        line.addParam('alignFrame0', params.IntParam, default=0, label='beginning')
        line.addParam('alignFrameN', params.IntParam, default=0, label='end')
        line = group.addLine('Remove frames to SUM from',
                             help='How many frames you want remove to sum\n'
                                  'from beginning and/or from the end of each movie.')
        line.addParam('sumFrame0', params.IntParam, default=0, label='beginning')
        line.addParam('sumFrameN', params.IntParam, default=0, label='end')
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
                      label="Save aligned micrograph", expertLevel=cons.LEVEL_ADVANCED)
        
        form.addParam('doSaveMovie', params.BooleanParam, default=False,
                      label="Save movie", expertLevel=cons.LEVEL_ADVANCED,
                      help="Save Aligned movie")


    #--------------------------- STEPS functions ----------------------------

    #FIXME: Methods will change when using the streaming for the output
    def createOutputStep(self):
        # Do nothing now, the output should be ready.
        pass

    def _loadOutputSet(self, SetClass, baseName, fixSampling=True):
        """
        Load the output set if it exists or create a new one.
        fixSampling: correct the output sampling rate if binning was used,
        except for the case when the original movies are kepts and shifts
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
        if newDone:
            self._writeDoneList(newDone)
        else:
            return

        firstTime = len(doneList) == 0
        allDone = len(doneList) + len(newDone)
        # We have finished when there is not more input movies (stream closed)
        # and the number of processed movies is equal to the number of inputs
        self.finished = self.streamClosed and allDone == len(self.listOfMovies)
        streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN

        if self._doGenerateOutputMovies():
            # FIXME: Even if we save the movie or not, both are aligned
            saveMovie = self.getAttributeValue('doSaveMovie', False)
            suffix = '_aligned' if saveMovie else '_original'
            movieSet = self._loadOutputSet(SetOfMovies,
                                           'movies%s.sqlite' % suffix,
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

        if self.getAttributeValue('doSaveAveMic', True):
            micSet = self._loadOutputSet(SetOfMicrographs,'micrographs.sqlite')

            for movie in newDone:
                mic = micSet.ITEM_TYPE()
                mic.copyObjId(movie)
                mic.setMicName(movie.getMicName())
                # The subclass protocol is responsible of generating the output
                # micrograph file in the extra path with the required name
                extraMicFn = self._getExtraPath(self._getOutputMicName(movie))
                mic.setFileName(extraMicFn)
                self._preprocessOutputMicrograph(mic, movie)
                # FIXME The micSet is not setting properly dimensions (No-Dim)
                micSet.append(mic)

            self._updateOutputSet('outputMicrographs', micSet, streamMode)
            if firstTime:
                self._defineSourceRelation(self.inputMovies, micSet)

        if self.finished: # Unlock createOutputStep if finished all jobs
            outputStep = self._getFirstJoinStep()
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(cons.STATUS_NEW)

    #--------------------------- INFO functions --------------------------------

    def _validate(self):
        errors = []
        if (self.cropDimX > 0 and self.cropDimY <= 0 or
            self.cropDimY > 0 and self.cropDimX <= 0):
            errors.append("If you give cropDimX, you should also give cropDimY"
                          " and viceversa")
        movie = self.inputMovies.get().getFirstItem()
        frames = movie.getNumberOfFrames()
        if frames is not None:
            a0, aN = self._getFrameRange(frames, "align")
            s0, sN = self._getFrameRange(frames, "sum")
            if aN < a0:
                errors.append("The final frame must be greater than the initial"
                              " frame. Please, check the range to remove frames"
                              " to *align* the movies.")
            if sN < s0:
                errors.append("The final frame must be greater than the initial"
                              " frame. Please, check the range to remove frames"
                              " to *sum* the movies.")
            
        return errors

    #--------------------------- INFO functions -------------------------------

    def _summary(self):
        return [self.summaryVar.get('')]

    #--------------------------- UTILS functions ----------------------------
    def _getFrameRange(self, n, prefix):
        """
        Params:
        :param n: Number of frames of the movies
        :param prefix: what range we want to consider, either 'align' or 'sum'
        :return: (i, f) initial and last frame range
        """
        first = 1 + self.getAttributeValue('%sFrame0' % prefix)
        last = n - self.getAttributeValue('%sFrameN' % prefix)

        return first, last
    
    def _createOutputMovie(self, movie):
        movieId = movie.getObjId()

        # Parse the alignment parameters and store the log files
        alignedMovie = movie.clone()
        n = movie.getNumberOfFrames()
        first, last = self._getFrameRange(n, 'align')

        # Check if user selected to save movie, use the getAttriveValue
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
            # If we save the movies, we need to modify which are
            # first and last frames if only a subset of frames have
            # been used to align
            first = 1
            last = totalFrames
        else:
            xshifts, yshifts = self._getMovieShifts(movie)


        alignment = MovieAlignment(first=first, last=last,
                                   xshifts=xshifts, yshifts=yshifts)

        roiList = [self.getAttributeValue(s, 0) for s in
                   ['cropOffsetX', 'cropOffsetY', 'cropDimX', 'cropDimY']]
        alignment.setRoi(roiList)

        alignedMovie.setAlignment(alignment)

        return alignedMovie

    #---------- Hook functions that need to be implemented in subclasses ------

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

    def _getMovieShifts(self, movie):
        """ Returns the x and y shifts for the alignment of this movie.
         The shifts should refer to the original micrograph without any binning.
         In case of a binning greater than 1, the shifts should be scaled.
        """
        return [], []

    def _doGenerateOutputMovies(self):
        """ Returns True if an output set of movies will be generated.
        The most common case is to always generate output movies,
        either with alignment only or the binary aligned movie files.
        Subclasses can override this function to change this behavior.
        """
        return True

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
                     dark=None, gain=None):
        """ Average a movie (using xmipp) taking into account the
         possible shifts and other alignment parameters.
         Params:
            inputFn: input filename, either the movie file or a metadata
                with the shifts and other info.
            dark: dark file
            gain: gain correction file.

         The output will be the averaged micrograph.
        """
        args  = '-i %s ' % inputFn
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
        m = int(round(x/2.))
        data1[:, m:] = data2[:, m:]
        psd.setData(data1)
        psd.write(outputFn)
