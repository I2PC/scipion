# **************************************************************************
# *
# * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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
from pyworkflow import VERSION_1_1
from pyworkflow.utils.properties import Message
from pyworkflow.utils.path import cleanPath
from pyworkflow.protocol.params import PointerParam, IntParam, BooleanParam, LEVEL_ADVANCED
from pyworkflow.em.protocol import EMProtocol, ProtProcessMovies
from pyworkflow.em.data import SetOfMovies, Movie
from pyworkflow.object import Set
import pyworkflow.protocol.constants as cons
import pyworkflow.em as em
import numpy as np
import os
import math
import xmipp


class XmippProtMovieGain(ProtProcessMovies):
    """
    Estimate the gain image of a camera, directly analyzing one of its movies.
    """
    _label = 'movie gain'
    _lastUpdateVersion = VERSION_1_1

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        self.stepsExecutionMode = em.STEPS_PARALLEL

    #--------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
 
        form.addParam('inputMovies', PointerParam, pointerClass='SetOfMovies, Movie',
                      label=Message.LABEL_INPUT_MOVS,
                      help='Select one or several movies. A gain image will '
                           'be calculated for each one of them.')
        form.addParam('frameStep', IntParam, default=1,
                      label="Frame step", expertLevel=LEVEL_ADVANCED,
                      help='By default, every frame (frameStep=1) is used to '
                           'compute the movie gain. If you set '
                           'this parameter to 2, 3, ..., then only every 2nd, '
                           '3rd, ... frame will be used.')
        form.addParam('movieStep', IntParam, default=1,
                      label="Movie step", expertLevel=LEVEL_ADVANCED,
                      help='By default, every movie (movieStep=1) is used to '
                           'compute the movie gain. If you set '
                           'this parameter to 2, 3, ..., then only every 2nd, '
                           '3rd, ... movie will be used.')
        form.addParam('useExistingGainImage', BooleanParam, default=None,
                      label="Use existing gain image",
                      help='If there is a gain image associated with input '
                           'movies, you can decide to use it instead of '
                           'estimating raw/residual gain image. Location of '
                           'this gain image needs to be indicated in import '
                           'movies protocol.')
        form.addParallelSection(threads=1, mpi=1)


    #--------------------------- STEPS functions -------------------------------

    def createOutputStep(self):
        # Do nothing now, the output should be ready.
        pass


    def _insertNewMoviesSteps(self, insertedDict, inputMovies):
        """ Insert steps to process new movies (from streaming)
        Params:
            insertedDict: contains already processed movies
            inputMovies: input movies set to be check
        """
        deps = []
        if isinstance(self.inputMovies.get(), Movie):
            movie = self.inputMovies.get()
            if movie.getObjId() not in insertedDict:
                stepId = self._insertMovieStep(movie)
                deps.append(stepId)
                insertedDict[movie.getObjId()] = stepId
        else:
            # For each movie insert the step to process it
            for idx, movie in enumerate(self.inputMovies.get()):
                if idx % self.movieStep.get() != 0:
                    continue
                if movie.getObjId() not in insertedDict:
                    stepId = self._insertMovieStep(movie)
                    deps.append(stepId)
                    insertedDict[movie.getObjId()] = stepId
        return deps


    def _processMovie(self, movie):
        movieId = movie.getObjId()
        fnMovie = movie.getFileName()
        gain = self.inputMovies.get().getGain()
        args = "-i %s --oroot %s --iter 1 --singleRef --frameStep %d " \
               "--gainImage %s" % \
               (fnMovie, self._getPath("movie_%06d" % movieId),
                self.frameStep, gain)
        if self.useExistingGainImage.get():
            args += " --applyGain"

        self.runJob("xmipp_movie_estimate_gain", args, numberOfMpi=1)
        cleanPath(self._getPath("movie_%06d_correction.xmp" % movieId))

        fnSummary = self._getPath("summary.txt")
        fnMonitorSummary = self._getPath("summaryForMonitor.txt")
        if not os.path.exists(fnSummary):
            fhSummary = open(fnSummary, "w")
            fnMonitorSummary = open(fnMonitorSummary, "w")
        else:
            fhSummary = open(fnSummary, "a")
            fnMonitorSummary = open(fnMonitorSummary, "a")
        fnGain = self._getPath("movie_%06d_gain.xmp" % movieId)
        if os.path.exists(fnGain):
            G = xmipp.Image()
            G.read(fnGain)
            mean, dev, min, max = G.computeStats()
            Gnp = G.getData()
            p = np.percentile(Gnp, [2.5, 25, 50, 75, 97.5])
            fhSummary.write("movie_%06d: mean=%f std=%f [min=%f,max=%f]\n" % (movieId, mean, dev, min, max))
            fhSummary.write(
                "            2.5%%=%f 25%%=%f 50%%=%f 75%%=%f 97.5%%=%f\n" % (p[0], p[1], p[2], p[3], p[4]))
            fhSummary.close()
            fnMonitorSummary.write("movie_%06d: %f %f %f %f\n" % (movieId, dev, p[0], p[4], max))
            fnMonitorSummary.close()


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

        if isinstance(self.inputMovies.get(), SetOfMovies):
            inputMovies = self.inputMovies.get()
            outputSet.copyInfo(inputMovies)

        if fixSampling:
            newSampling = inputMovies.getSamplingRate() * self._getBinFactor()
            outputSet.setSamplingRate(newSampling)

        return outputSet

    def _checkNewInput(self):
        if isinstance(self.inputMovies.get(), SetOfMovies):
            ProtProcessMovies._checkNewInput(self)

    def _checkNewOutput(self):
        if getattr(self, 'finished', False):
            return
        if isinstance(self.inputMovies.get(), Movie):
            saveMovie = self.getAttributeValue('doSaveMovie', False)
            imageSet = self._loadOutputSet(em.data.SetOfImages, 'movies.sqlite', fixSampling=saveMovie)

            movie = self.inputMovies.get()
            imgOut = em.data.Image()
            imgOut.setObjId(movie.getObjId())
            imgOut.setSamplingRate(movie.getSamplingRate())
            imgOut.setFileName(self._getPath("movie_%06d_gain.xmp" % movie.getObjId()))
            imageSet.setSamplingRate(movie.getSamplingRate())
            imageSet.append(imgOut)

            self._updateOutputSet('outputMovies', imageSet, Set.STREAM_CLOSED)
            outputStep = self._getFirstJoinStep()
            outputStep.setStatus(cons.STATUS_NEW)
            self.finished = True
        else:

            # Load previously done items (from text file)
            doneList = self._readDoneList()
            # Check for newly done items
            newDone = [m.clone() for m in self.listOfMovies
                       if int(m.getObjId()) not in doneList and self._isMovieDone(m)]

            # Update the file with the newly done movies
            # or exit from the function if no new done movies
            self.debug('_checkNewOutput: ')
            self.debug('   listOfMovies: {0}, doneList: {1}, newDone: {2}'
                       .format(int(math.ceil(len(self.listOfMovies)/float(self.movieStep.get()))), len(doneList), len(newDone)))

            allDone = len(doneList) + len(newDone)
            # We have finished when there is not more input movies (stream closed)
            # and the number of processed movies is equal to the number of inputs
            self.finished = self.streamClosed and allDone == int(math.ceil(len(self.listOfMovies)/float(self.movieStep.get())))
            streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN

            if newDone:
                self._writeDoneList(newDone)
            elif not self.finished:
                # If we are not finished and no new output have been produced
                # it does not make sense to proceed and updated the outputs
                # so we exit from the function here
                return

            self.debug('   finished: %s ' % self.finished)
            self.debug('        self.streamClosed ({0}) AND' .format(self.streamClosed))
            self.debug('        allDone ({0}) == len(self.listOfMovies ({1})'
                       .format(allDone, int(math.ceil(len(self.listOfMovies)/float(self.movieStep.get())))))
            self.debug('   streamMode: %s' % streamMode)

            saveMovie = self.getAttributeValue('doSaveMovie', False)
            imageSet = self._loadOutputSet(em.data.SetOfImages, 'movies.sqlite', fixSampling=saveMovie)

            for movie in newDone:
                imgOut = em.data.Image()
                imgOut.setObjId(movie.getObjId())
                imgOut.setSamplingRate(movie.getSamplingRate())
                imgOut.setFileName(self._getPath("movie_%06d_gain.xmp" % movie.getObjId()))
                imageSet.setSamplingRate(movie.getSamplingRate())
                imageSet.append(imgOut)

            self._updateOutputSet('outputMovies', imageSet, streamMode)

            if self.finished:  # Unlock createOutputStep if finished all jobs
                outputStep = self._getFirstJoinStep()
                if outputStep and outputStep.isWaiting():
                    outputStep.setStatus(cons.STATUS_NEW)


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

        # Close set databaset to avoid locking it
        outputSet.close()

    #--------------------------- INFO functions -------------------------------
    def _summary(self):
        fnSummary = self._getPath("summary.txt")
        if not os.path.exists(fnSummary):
            summary = ["No summary information yet."]
        else:
            fhSummary = open(fnSummary,"r")
            summary = []
            for line in fhSummary.readlines():
                summary.append(line.rstrip())
            fhSummary.close()
        return summary