# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              Amaya Jimenez    (ajimenez@cnb.csic.es)
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

    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
    #    ProtProcessMovies._defineParams(self, form)
        form.addSection(label='Input')
        form.addParam('inputMovies', PointerParam, pointerClass='SetOfMovies', 
              important=True,
              label=Message.LABEL_INPUT_MOVS,
              help='Select a set of previously imported movies.')
        form.addParam('maxFrameShift', params.FloatParam, default=10, label='Max. frame shift (A)')
        form.addParam('maxMovieShift', params.FloatParam, default=100, label='Max. movie shift (A)')
        
    #--------------------------- INSERT steps functions ------------------------

    def _createOutputMovie(self, movie):
        """ Create movie only if the alignment is less than the thresshold. """
        alignment = movie.getAlignment()

        if alignment is not None:
            shiftListX, shiftListY = alignment.getShifts()
            shiftArrayX = np.asarray(shiftListX)
            shiftArrayY = np.asarray(shiftListY)
            rejectedByFrame = np.max(np.abs(shiftArrayX))*self.samplingRate>self.maxFrameShift.get() or \
                              np.max(np.abs(shiftArrayY))*self.samplingRate>self.maxFrameShift.get()
            cumsumX = np.cumsum(shiftArrayX)
            cumsumY = np.cumsum(shiftArrayY)
            rejectedByMovie = np.max(np.abs(cumsumX))*self.samplingRate>self.maxMovieShift.get() or \
                              np.max(np.abs(cumsumY))*self.samplingRate>self.maxMovieShift.get()

            return not rejectedByFrame and not rejectedByMovie
        else: # If movies are not aligned, this protocal makes not sense
            return False # The protocol returns an empty SetOfMovies to poit out the non-sense.


    #    # FIXME: Methods will change when using the streaming for the output
    def createOutputStep(self):
        # Do nothing now, the output should be ready.
        pass

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

        if os.path.exists(setFile):
            outputSet = SetClass(filename=setFile)
            outputSet.loadAllProperties()
            outputSet.enableAppend()
        else:
            outputSet = SetClass(filename=setFile)
            outputSet.setStreamState(outputSet.STREAM_OPEN)

        inputMovies = self.inputMovies.get()
        outputSet.copyInfo(inputMovies)

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
                if self._createOutputMovie(movie):
                    newMovie = movie.clone()
                    movieSet.append(newMovie)

            self._updateOutputSet('outputMovies', movieSet, streamMode)

            # if firstTime:

            #     # Probably is a good idea to store a cached summary for the
            #     # first resulting movie of the processing.
            #     # self._storeSummary(newDone[0])                 ----------------------

            #     # If the movies are not written out, then dimensions can be
            #     # copied from the input movies

            #     if saveMovie:
            #         movieSet.setDim(self.inputMovies.get().getDim())
                    
            #     self._defineTransformRelation(self.inputMovies, movieSet)

        if self.finished:  # Unlock createOutputStep if finished all jobs
            outputStep = self._getFirstJoinStep()
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(STATUS_NEW)

