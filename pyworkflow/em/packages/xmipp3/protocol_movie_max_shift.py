# **************************************************************************
# *
# * Authors:     David Maluenda    (dmaluenda@cnb.csic.es)
# *              Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

from pyworkflow.mapper import Mapper
from pyworkflow.object import Set
from pyworkflow.protocol.constants import (STATUS_NEW)
from pyworkflow.protocol.params import PointerParam
from pyworkflow.utils.properties import Message
from pyworkflow.em.data import (MovieAlignment, SetOfMovies, SetOfMicrographs,
                                Image, Movie)


class XmippProtMovieMaxShift(ProtProcessMovies):
    """
    Protocol to make an automatic selection of those movies whose
    frames move less than a given threshold.
    """
    _label = 'movie maxshift'
    
    REJ_TYPES = ['by frame', 'by whole movie', 'by frame and movie', 
                 'by frame or movie']
    REJ_FRAME = 0
    REJ_MOVIE = 1
    REJ_AND = 2
    REJ_OR = 3

    def __init__(self, **args):
        ProtProcessMovies.__init__(self, **args)
        self.acceptedMoviesList = []
        self.discardedMoviesList = []
        self.acceptedDone = 0
        self.discardedDone = 0

    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        # ProtProcessMovies._defineParams(self, form)
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputMovies', PointerParam, important=True,
                      label=Message.LABEL_INPUT_MOVS,
                      pointerClass='SetOfMovies', 
                      help='Select a set of previously aligned Movies.')

        form.addParam('rejType', params.EnumParam, choices=self.REJ_TYPES,
                      label='Rejection type', default=self.REJ_OR,
                      help='Rejection criteria:\n'
                           ' - *by frame*: Rejects movies with drifts between '
                                    'frames bigger than a certain maximum.\n'
                           ' - *by whole movie*: Rejects movies with a total '
                                    'travel bigger than a certain maximmum.\n'
                           ' - *by frame and movie*: Rejects movies if both '
                                    'conditions above are met.\n'
                           ' - *by frame or movie*: Rejects movies if one of '
                                    'the conditions above are met.')
        form.addParam('maxFrameShift', params.FloatParam, default=1, 
                       label='Max. frame shift (A)',
                       condition='rejType==%s or rejType==%s or rejType==%s'
                                  % (self.REJ_FRAME, self.REJ_AND, self.REJ_OR),
                       help='Maximum drift between consecutive frames '
                            'to evaluate the frame condition.')
        form.addParam('maxMovieShift', params.FloatParam, default=15,
                       label='Max. movie shift (A)',
                       condition='rejType==%s or rejType==%s or rejType==%s'
                                  % (self.REJ_MOVIE, self.REJ_AND, self.REJ_OR),
                       help='Maximum total travel to evaluate the whole movie '
                            'condition.')
        
    #--------------------------- INSERT steps functions ------------------------
    def _insertMovieStep(self, movie):
        """ Insert the processMovieStep for a given movie. """
        # Note1: At this point is safe to pass the movie, since this
        # is not executed in parallel, here we get the params
        # to pass to the actual step that is gone to be executed later on
        # Note2: We are serializing the Movie as a dict that can be passed
        # as parameter for a functionStep
        movieStepId = self._insertFunctionStep('_evaluateMovieAlign',
                                               movie.clone(),
                                               prerequisites=[])
        return movieStepId

    def _evaluateMovieAlign(self, movie):
        """ Fill either the accepted or the rejected list with the movieID """
        alignment = movie.getAlignment()
        sampling = self.samplingRate

        # getShifts() returns the absolute shifts from a certain refference
        shiftListX, shiftListY = alignment.getShifts()

        rejectedByMovie = False
        rejectedByFrame = False
        if any(shiftListX) or any(shiftListY):
            shiftArrayX = np.asarray(shiftListX)
            shiftArrayY = np.asarray(shiftListY)

            evalBoth = self.rejType==self.REJ_AND or self.rejType==self.REJ_OR

            if self.rejType == self.REJ_MOVIE or evalBoth:
                rangeX = np.max(shiftArrayX) - np.min(shiftArrayX)
                rangeY = np.max(shiftArrayY) - np.min(shiftArrayY)
                rejectedByMovie = (rangeX*sampling > self.maxMovieShift.get() or
                                   rangeY*sampling > self.maxMovieShift.get() )

            if self.rejType == self.REJ_FRAME or evalBoth:
                frameShiftX = np.abs(np.diff(shiftArrayX))
                frameShiftY = np.abs(np.diff(shiftArrayY))
                rejectedByFrame = ( np.max(frameShiftX) * sampling > 
                                    self.maxFrameShift.get()          or
                                    np.max(frameShiftY) * sampling >
                                    self.maxFrameShift.get() )

            if self.rejType == self.REJ_AND:
                if rejectedByFrame and rejectedByMovie:
                    self.discardedMoviesList.append(movie)
                else:
                    self.acceptedMoviesList.append(movie)
            else:  # for the OR and also for the individuals evaluations
                if rejectedByFrame or rejectedByMovie:
                    self.discardedMoviesList.append(movie)
                else:
                    self.acceptedMoviesList.append(movie)

    def _checkNewOutput(self):
        """ Check for already selected Movies and update the output set. """
        if getattr(self, 'finished', False):
            return

        # load if first time in order to make dataSets relations
        firstTimeAcc = self.acceptedDone==0
        firstTimeDisc = self.discardedDone==0

        # Load previously done items
        preDone = self.acceptedDone + self.discardedDone

        # Check for newly done items
        newDoneAccepted = self.acceptedMoviesList[self.acceptedDone:]
        newDoneDiscarded = self.discardedMoviesList[self.discardedDone:]

        # Updating the done lists
        self.acceptedDone = len(self.acceptedMoviesList)
        self.discardedDone = len(self.discardedMoviesList)
        
        allDone = preDone + len(newDoneAccepted) + len(newDoneDiscarded)

        # Update the file with the newly done movies
        # or exit from the function if no new done movies
        self.debug('_checkNewOutput: ')
        self.debug('   listOfMovies: %d,' %len(self.listOfMovies))
        self.debug('   doneList: %d,' %preDone)
        self.debug('   newDoneAccepted: %d.' %len(newDoneAccepted))
        self.debug('   newDoneDiscarded: %d,' %len(newDoneDiscarded))
    
        # We have finished when there is not more input movies (stream closed)
        # and the number of processed movie is equal to the number of inputs
        self.finished = self.streamClosed and allDone == len(self.listOfMovies)
        streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN

        if not self.finished and not newDoneDiscarded and not newDoneAccepted:
            # If we are not finished and no new output have been produced
            # it does not make sense to proceed and updated the outputs
            # so we exit from the function here
            return

        if newDoneAccepted:
            movieSetAccepted = self._loadOutputSet(SetOfMovies, 'movies.sqlite')
            for movie in newDoneAccepted:
                movie.setEnabled(True)
                movieSetAccepted.append(movie)
                
            self._updateOutputSet('outputMovies', movieSetAccepted, streamMode)
            if firstTimeAcc:
                # define relation just once
                self._defineTransformRelation(self.inputMovies.get(),
                                           movieSetAccepted)
            movieSetAccepted.close()

        # new subsets with discarded movies
        if newDoneDiscarded:
            movieSetDiscarded = self._loadOutputSet(SetOfMovies,
                                                       'moviesDiscarded.sqlite')
            for movie in newDoneDiscarded:
                movie.setEnabled(False)
                movieSetDiscarded.append(movie)

            self._updateOutputSet('outputMoviesDiscarded',
                                  movieSetDiscarded, streamMode)
            if firstTimeDisc:
                # define relation just once
                self._defineTransformRelation(self.inputMovies.get(),
                                              movieSetDiscarded)
            movieSetDiscarded.close()

        # Unlock createOutputStep if finished all jobs
        if self.finished:  
            outputStep = self._getFirstJoinStep()
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(STATUS_NEW)            


    # FIXME: Methods will change when using the streaming for the output
    def createOutputStep(self):
        # Do nothing now, the output should be ready.
        pass

    #--------------------------- UTILS functions -------------------------------
    def _loadOutputSet(self, SetClass, baseName, fixSampling=True):
        """ Load the output set if it exists or create a new one.
        fixSampling: correct the output sampling rate if binning was used,
        except for the case when the original movies are kept and shifts
        refers to that one. """
        setFile = self._getPath(baseName)

        if exists(setFile):
            outputSet = SetClass(filename=setFile)
            if(outputSet.__len__() is 0):
                pwutils.path.cleanPath(setFile)

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

    # ---------------------- INFO functions ------------------------------------
    def _validate(self):
        errors = []
        if not self.inputMovies.get().getFirstItem().hasAlignment():
            errors.append('The _Input Movies_ must come from an alignment '
                          'protocol.')
        return errors

    def _summary(self):
        moviesAcc = 0 if not self.hasAttribute('outputMovies') else \
                    self.outputMovies.getSize()
        moviesDisc = 0 if not self.hasAttribute('outputMoviesDiscarded') else \
                    self.outputMoviesDiscarded.getSize()

        summary = ['Movies processed: %d'%(moviesAcc+moviesDisc)]
        summary.append('Movies rejected: *%d*' % moviesDisc)

        return summary