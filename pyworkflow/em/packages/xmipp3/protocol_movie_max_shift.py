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


from os.path import getmtime, exists
from datetime import datetime
import numpy as np

import pyworkflow.em as em
from pyworkflow.em.protocol.protocol_movies import ProtProcessMovies
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils

from pyworkflow.object import Set
from pyworkflow.protocol.constants import (STATUS_NEW)


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
        ProtProcessMovies._defineParams(self, form)
        # form.addSection(label='Input')
        form.addParam('maxFrameShift', params.FloatParam, default=10, label='Max. frame shift (A)')
        form.addParam('maxMovieShift', params.FloatParam, default=100, label='Max. movie shift (A)')
        
    #--------------------------- INSERT steps functions ------------------------
    def _filterMovie(self, movie):
        # Try to use the 'original' fileName in case it is present
        # the original could be different from the current filename if
        # we are dealing with compressed movies (e.g., movie.mrc.bz2)
        alignment = movie.getAlignment()
        print(movie.hasAlignment())
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
        # else
        #     print("No alignment")


        print("OK=",not rejectedByFrame and not rejectedByMovie)
        return not rejectedByFrame and not rejectedByMovie

    def createOutputStep(self):
        # Do nothing now, the output should be ready.
        print("Adeu")
        
