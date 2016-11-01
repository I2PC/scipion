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

from pyworkflow.utils.properties import Message
from pyworkflow.utils.path import cleanPath
from pyworkflow.protocol.params import MultiPointerParam
from pyworkflow.em.protocol import ProtPreprocessMicrographs, EMProtocol
import pyworkflow.em as em



class XmippProtMovieGain(ProtPreprocessMicrographs):
    """
    Estimate the gain image of a camera, directly analyzing one of its movies.
    """
    _label = 'movie gain'

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        self.stepsExecutionMode = em.STEPS_PARALLEL

    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
 
        form.addParam('inputMovies', MultiPointerParam, pointerClass='Movie', 
                      label=Message.LABEL_INPUT_MOVS,
                      help='Select one or several movies. A gain image will '
                           'be calculated for each one of them')
        form.addParallelSection(threads=1, mpi=1)


    #--------------------------- STEPS functions -------------------------------
    def _insertAllSteps(self):
        gainSteps = []
        for pointer in self.inputMovies:
            movie = pointer.get()
            movieId = movie.getObjId()
            fnMovie = movie.getFileName()
 
            stepId = self._insertFunctionStep('estimateGain', movieId, fnMovie,
                                              prerequisites=[])
            gainSteps.append(stepId)
             
        self._insertFunctionStep('createOutputStep', prerequisites=gainSteps)

    def createOutputStep(self):
        movie = self.inputMovies[0].get()
        if len(self.inputMovies)==1:
            imgOut = em.data.Image()
            imgOut.setSamplingRate(movie.getSamplingRate())
            imgOut.setFileName(self._getPath("movie_%06d_gain.xmp"
                                             % movie.getObjId()))
    
            self._defineOutputs(outputGainImage=imgOut)
            self._defineSourceRelation(self.inputMovies, imgOut)
        else:
            imgSetOut = self._createSetOfImages()
            imgSetOut.setSamplingRate(movie.getSamplingRate())
            for pointer in self.inputMovies:
                movie = pointer.get()
                imgOut = em.data.Image()
                imgOut.setSamplingRate(movie.getSamplingRate())
                imgOut.setFileName(self._getPath("movie_%06d_gain.xmp"
                                                 % movie.getObjId()))
                imgSetOut.append(imgOut)
    
            self._defineOutputs(outputSetOfGainImages=imgSetOut)
            self._defineSourceRelation(self.inputMovies, imgSetOut)
    
    def estimateGain(self, movieId, fnMovie):
        self.runJob("xmipp_movie_estimate_gain",
                    "-i %s --oroot %s --iter 1 --singleRef"
                    % (fnMovie, self._getPath("movie_%06d" % movieId)),
                    numberOfMpi=1)
        cleanPath(self._getPath("movie_%06d_correction.xmp" % movieId))
    
