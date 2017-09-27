# **************************************************************************
# *
# * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
# *              David Maluenda    (dmaluenda@cnb.csic.es)
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

from os.path import getmtime, exists
from datetime import datetime
import numpy as np

import pyworkflow.em as em
from pyworkflow.em.protocol.protocol_movies import ProtProcessMovies
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils

from pyworkflow.object import Set
from pyworkflow.protocol.constants import (STATUS_NEW)
from pyworkflow.protocol.params import PointerParam
from pyworkflow.utils.properties import Message
from pyworkflow.em.data import (MovieAlignment, SetOfMovies, SetOfMicrographs,
                                Image)


class XmippProtMovieMaxShift(ProtProcessMovies):
    """
    Protocol to make an automatic selection of those movies whose
    frames move less than a given threshold.
    """
    _label = 'movie maxshift'
    
    def __init__(self, **args):
        ProtProcessMovies.__init__(self, **args)
        self.acceptedIdMoviesList = []
        self.discartedIdMoviesList = []

    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        ProtProcessMovies._defineParams(self, form)

        form.addParam('maxFrameShift', params.FloatParam, default=10, 
                       label='Max. frame shift (A)')
        form.addParam('maxMovieShift', params.FloatParam, default=100,
                       label='Max. movie shift (A)')
        
    #--------------------------- INSERT steps functions ------------------------
    def _processMovie(self, movie):
        """ Create movie only if the alignment is less than the thresshold. """
        movieId = movie.getObjId()
        alignment = movie.getAlignment()

        if alignment is not None:
            shiftListX, shiftListY = alignment.getShifts()
            shiftArrayX = np.asarray(shiftListX)
            shiftArrayY = np.asarray(shiftListY)
            samRate = self.samplingRate
            rejectedByFrame = ( np.max(np.abs(shiftArrayX)) * samRate > 
                                self.maxFrameShift.get()              or
                                np.max(np.abs(shiftArrayY)) * samRate > 
                                self.maxFrameShift.get() )

            cumsumX = np.cumsum(shiftArrayX)
            cumsumY = np.cumsum(shiftArrayY)
            rejectedByMovie = ( np.max(np.abs(cumsumX)) * samRate > 
                                self.maxMovieShift.get()          or
                                np.max(np.abs(cumsumY)) * samRate >
                                self.maxMovieShift.get() )

            if not rejectedByFrame and not rejectedByMovie:
                self._updateLists(movieId,True)
            else:
                self._updateLists(movieId,False)

        else:
            # a no aligned movie is DISCARTED
            # (maybe change this or add a param to control that) 
            self._updateLists(movieId,False)

    def _checkNewOutput(self):
        """ Check for already selected Movies and update the output set. """
        if getattr(self, 'finished', False):
            return

        # Load previously done items (from text file)
        doneList = self._readDoneList()

        # Check for newly done items
        movieListIdAccepted = self._getLists(True)
        movieListIdDiscarted = self._getLists(False)
    
        newDoneAccepted = [movieId for movieId in movieListIdAccepted
                             if movieId not in doneList]
        newDoneDiscarted = [movieId for movieId in movieListIdDiscarted
                             if movieId not in doneList]

        # Update the file with the newly done movies
        # or exit from the function if no new done movies
        self.debug('_checkNewOutput: ')
        self.debug('   listOfMovies: %d,' %len(self.listOfMovies))
        self.debug('   doneList: %d,' %len(doneList))
        self.debug('   newDoneAccepted: %d.' %len(newDoneAccepted))
        self.debug('   newDoneDiscarted: %d,' %len(newDoneDiscarted))

        firstTime = len(doneList) == 0
        # firstTimeDiscarted = len(doneListDiscarted) == 0
        allDone = len(doneList) + len(newDoneAccepted) + len(newDoneDiscarted)
    
        # We have finished when there is not more input movies (stream closed)
        # and the number of processed movie is equal to the number of inputs
        # self.finished = (self.isStreamClosed == Set.STREAM_CLOSED and 
        #                    allDone == len(self.listOfMovies))
        self.finished = self.streamClosed and allDone == len(self.listOfMovies)
        streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN

        # reading the outputs
        if (len(doneList)>0 or len(newDoneAccepted)>0):
            movieSetAccepted = self._loadOutputSet(SetOfMovies, 'movies.sqlite')

        #AJ new subsets with discarted movies
        if (len(doneList)>0 or len(newDoneDiscarted)>0):
            movieSetDiscarted = self._loadOutputSet(SetOfMovies,
                                                       'moviesDiscarted.sqlite')
        if newDoneAccepted:
            inputMovieSet = self._loadInputMovieSet()
            for movieId in newDoneAccepted:
                movieAccepted = inputMovieSet[movieId].clone()
                movieAccepted.setEnabled(True)
                movieSetAccepted.append(movieAccepted)
            inputMovieSet.close()

        if newDoneDiscarted:
            inputMovieSet = self._loadInputMovieSet()
            for movieId in newDoneDiscarted:
                movieDiscarted = inputMovieSet[movieId].clone()
                movieDiscarted.setEnabled(False)
                movieSetDiscarted.append(movieDiscarted)
            inputMovieSet.close()

        if not self.finished and not newDoneDiscarted and not newDoneAccepted:
            # If we are not finished and no new output have been produced
            # it does not make sense to proceed and updated the outputs
            # so we exit from the function here
            return

        if (exists(self._getPath('movies.sqlite'))):
            self._updateOutputSet('outputMovies', movieSetAccepted, streamMode)

        # AJ new subsets with discarted movies
        if (exists(self._getPath('moviesDiscarted.sqlite'))):
            self._updateOutputSet('outputMoviesDiscarted',
                                  movieSetDiscarted, streamMode)
            
        if self.finished:  # Unlock createOutputStep if finished all jobs
            outputStep = self._getFirstJoinStep()
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(STATUS_NEW)

        if (exists(self._getPath('movies.sqlite'))):
            movieSetAccepted.close()

        # AJ new subsets with discarted movies
        if (exists(self._getPath('moviesDiscarted.sqlite'))):
            movieSetDiscarted.close()

    #    # FIXME: Methods will change when using the streaming for the output
    def createOutputStep(self):
        # Do nothing now, the output should be ready.
        pass

    #--------------------------- UTILS functions -------------------------------

    def _createOutputMovies(self):
        """ Returns True if an output set of movies will be generated.
        The most common case is to always generate output movies,
        either with alignment only or the binary aligned movie files.
        Subclasses can override this function to change this behavior.
        """
        return True

    def _loadOutputSet(self, SetClass, baseName, fixSampling=True):
        """
        Load the output set if it exists or create a new one.
        fixSampling: correct the output sampling rate if binning was used,
        except for the case when the original movies are kept and shifts
        refers to that one.
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

        if fixSampling:
            newSampling = inputMovies.getSamplingRate()
            outputSet.setSamplingRate(newSampling)

        return outputSet

    def _loadInputMovieSet(self):
        movieFile = self.inputMovies.get().getFileName()
        self.debug("Loading input db: %s" % movieFile)
        movieSet = SetOfMovies(filename=movieFile)
        movieSet.loadAllProperties()
        return movieSet

    def _getLists(self, accepted):
        if accepted:
            return self.acceptedIdMoviesList
        else:
            return self.discartedIdMoviesList

    def _updateLists(self, movieId, accepted):
        if accepted:
            return self.acceptedIdMoviesList.append(movieId)
        else:
            return self.discartedIdMoviesList.append(movieId)