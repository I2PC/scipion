# **************************************************************************
# *
# * Authors:     Amaya Jimenez (ajimenez@cnb.csic.es)
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
from pyworkflow import VERSION_1_2
from pyworkflow.utils.properties import Message
from pyworkflow.protocol.params import PointerParam, FloatParam, EnumParam, \
    IntParam
from pyworkflow.em.protocol import EMProtocol, ProtProcessMovies
from pyworkflow.object import Set
import pyworkflow.protocol.constants as cons
from pyworkflow.em import STEPS_PARALLEL
from pyworkflow.em.data import SetOfMovies, Movie
from os.path import exists

RESIZE_SAMPLINGRATE = 0
RESIZE_DIMENSIONS = 1
RESIZE_FACTOR = 2


class XmippProtMovieResize(ProtProcessMovies):
    """
    Resize a set of movies. Only downsampling is allowed.
    """
    _label = 'movie resize'
    _lastUpdateVersion = VERSION_1_2


    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)

        form.addParam('inputMovies', PointerParam, pointerClass='SetOfMovies',
                      label=Message.LABEL_INPUT_MOVS,
                      help='Select a set of movies to be resized.')
        form.addParam('resizeOption', EnumParam,
                      choices=['Sampling Rate', 'Dimensions', 'Factor'],
                      default=RESIZE_SAMPLINGRATE,
                      label="Resize option", display=EnumParam.DISPLAY_COMBO,
                      help='Select an option to resize the images: \n '
                           '_Sampling Rate_: Set the desire sampling rate '
                           'to resize. \n'
                           '_Dimensions_: Set the output dimensions. \n'
                           '_Factor_: Set a resize factor to resize. \n ')
        form.addParam('resizeSamplingRate', FloatParam, default=4.0,
                      condition='resizeOption==%d' % RESIZE_SAMPLINGRATE,
                      label='Resize sampling rate (A/px)',
                      help='Set the new output sampling rate.')
        form.addParam('resizeDim', IntParam, default=1024,
                      condition='resizeOption==%d' % RESIZE_DIMENSIONS,
                      label='New image size (px)',
                      help='Size in pixels of the particle images '
                           '<x> <y=x> <z=x>.')
        form.addParam('resizeFactor', FloatParam, default=2.0,
                      condition='resizeOption==%d' % RESIZE_FACTOR,
                      label='Downsampling factor',
                      help='New size is the old one x resize factor.')
        form.addParallelSection(threads=1, mpi=1)

    # ------------------------ INSERT STEPS functions --------------------------

    def _insertNewMoviesSteps(self, insertedDict, inputMovies):
        """ Insert steps to process new movies (from streaming)
        Params:
            insertedDict: contains already processed movies
            inputMovies: input movies set to be check
        """
        deps = []
        # For each movie insert the step to process it
        for idx, movie in enumerate(self.inputMovies.get()):
            if movie.getObjId() not in insertedDict:
                stepId = self._insertMovieStep(movie)
                deps.append(stepId)
                insertedDict[movie.getObjId()] = stepId
        return deps

    # --------------------------- STEPS functions ------------------------------

    def _processMovie(self, movie):
        movieId = movie.getObjId()
        fnMovie = movie.getFileName()

        dim, _, numMov = self.inputMovies.get().getDim()
        samplingRate = self.inputMovies.get().getSamplingRate()
        if self.resizeOption==RESIZE_FACTOR:
            factor = self.resizeFactor.get()
            newDim = int(dim/factor)
            self.newSamplingRate =  samplingRate*factor
        if self.resizeOption==RESIZE_DIMENSIONS:
            newDim = self.resizeDim.get()
            self.newSamplingRate = float(samplingRate)*(float(dim)/float(
                newDim))
        if self.resizeOption==RESIZE_SAMPLINGRATE:
            self.newSamplingRate = self.resizeSamplingRate.get()
            factor = self.newSamplingRate/samplingRate
            newDim = int(dim/factor)

        args = "-i %s -o %s --fourier %d %d %d" % \
               (fnMovie, self._getPath("movie_%06d_resize.mrcs" % movieId),
                newDim, newDim, numMov)

        self.runJob("xmipp_image_resize", args, numberOfMpi=1)

    # -------------- Methods for Streaming --------------------

    def createOutputStep(self):
        # Do nothing now, the output should be ready.
        pass

    def _checkNewOutput(self):
        if getattr(self, 'finished', False):
            return

        # Load previously done items (from text file)
        doneList = self._readDoneList()
        # Check for newly done items
        newDone = [m.clone() for m in self.listOfMovies
                   if int(m.getObjId()) not in doneList
                   and self._isMovieDone(m)]

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

        imageSet = self._loadOutputSet(SetOfMovies, 'movies.sqlite')

        for movie in newDone:
            imgOut = Movie()
            imgOut.setObjId(movie.getObjId())
            imgOut.setFileName(self._getPath(
                "movie_%06d_resize.mrcs" % movie.getObjId()))
            imgOut.setAcquisition(movie.getAcquisition())
            imgOut.setSamplingRate(self.newSamplingRate)
            imgOut.setFramesRange(self.inputMovies.get().getFramesRange())
            imageSet.append(imgOut)

        self._updateOutputSet('outputMovies', imageSet, streamMode)

        if self.finished:  # Unlock createOutputStep if finished all jobs
            outputStep = self._getFirstJoinStep()
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(cons.STATUS_NEW)

    # --------------------------- UTILS functions ------------------------------

    def _loadOutputSet(self, SetClass, baseName):
        """
        Load the output set if it exists or create a new one.
        """
        setFile = self._getPath(baseName)
        if exists(setFile):
            outputSet = SetClass(filename=setFile)
            outputSet.loadAllProperties()
            outputSet.enableAppend()
        else:
            outputSet = SetClass(filename=setFile)
            outputSet.setStreamState(outputSet.STREAM_OPEN)

        inputMovies = self.inputMovies.get()
        outputSet.copyInfo(inputMovies)

        outputSet.setSamplingRate(self.newSamplingRate)

        return outputSet

    def _updateOutputSet(self, outputName, outputSet, state=Set.STREAM_OPEN):
        outputSet.setStreamState(state)

        if self.hasAttribute(outputName):
            outputSet.write()  # Write to commit changes
            outputAttr = getattr(self, outputName)
            # Copy the properties to the object contained in the protcol
            outputAttr.copy(outputSet, copyId=False)
            # Persist changes
            self._store(outputAttr)
        else:
            # Here the defineOutputs function will call the write() method
            self._defineOutputs(**{outputName: outputSet})
            self._store(outputSet)

        # Close set database to avoid locking it
        outputSet.close()

# --------------------------- INFO functions -------------------------------
    def _summary(self):
        if not hasattr(self, 'outputMovies'):
            summary = ["No summary information yet."]
        else:
            xIn, yIn, zIn = self.inputMovies.get().getDim()
            xOut, yOut, zOut = self.outputMovies.getDim()
            summary = ["%d input movies of size %dx%dx%d resized to %d "
                       "output movies of size %dx%dx%d."
                       %(self.inputMovies.get().getSize(), xIn, yIn, zIn,
                       self.outputMovies.getSize(), xOut, yOut, zOut)]
        return summary

    def validate(self):
        """ Try to find errors on define params. """
        errors = []
        dim, _, numMov = self.inputMovies.get().getDim()
        samplingRate = self.inputMovies.get().getSamplingRate()
        if self.resizeOption == RESIZE_FACTOR:
            if self.resizeFactor.get()<1.0:
                errors.append('Please provide a resizeFactor higher than 1, '
                              'only downsampling is allowed.')
        if self.resizeOption == RESIZE_DIMENSIONS:
            if self.resizeDim.get() > dim:
                errors.append('Please provide a resizeDim higher than the '
                              'size of the input set, only downsampling is '
                              'allowed.')
        if self.resizeOption == RESIZE_SAMPLINGRATE:
            if self.resizeSamplingRate.get() <= samplingRate:
                errors.append('Please provide a resizeSamplingRate higher '
                              'than the original one, only downsampling is '
                              'allowed.')
        return errors