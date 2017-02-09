# ******************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
# ******************************************************************************
"""
Protocol wrapper around the ResMap tool for local resolution
"""

from pyworkflow.protocol.params import PointerParam, PathParam
from pyworkflow.em.protocol import ProtProcessMovies, ProtImportFiles

from convert import parseMovieAlignment



class ProtDosefGpuImport(ProtImportFiles, ProtProcessMovies):
    """
    Import alignment from 
    
    Dose Fractionation Tool: Flat fielding and Drift correction
    Wrote by Xueming Li @ Yifan Cheng Lab, UCSF   
    """
    _label = 'import movie alignment'
             
    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputMovies', PointerParam,
                      pointerClass='SetOfMicrographs',
                          label='Input movies',
                          help='Select the movies to which assign the '
                               'alignment parameters')
        form.addParam('pattern', PathParam, label='Alignment pattern',
                      help='Select files containing the dosefgpu alignment '
                           '*_Log.txt files\n'
                           'You should use #### characters in the pattern\n'
                           'to mark where the movie id will be taken from. ')
    
    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('importAlignmentStep', 
                                 self.inputMovies.get().getObjId(),
                                 self.getPattern())
    
    #--------------------------- STEPS functions -------------------------------
    def importAlignmentStep(self, micsId, pattern):
        """ Copy alignment matching the filename pattern
        """
        inputMovies = self.inputMovies.get()
        # Also create a Set of Movies with the alignment parameters
        movieSet = self._createSetOfMovies()
        movieSet.copyInfo(inputMovies)
        
        alignmentFiles = self._getFilePaths(pattern)
        movieDict = {}
        
        for fn in alignmentFiles:
            # Try to match the micrograph id from filename
            # this is set by the user by using #### format in the pattern
            match = self._idRegex.match(fn)
            if match is None:
                raise Exception("File '%s' doesn't match the pattern '%s'" % (fn, self.pattern.get()))
            movieId = int(match.group(1))
            movieDict[movieId] = fn
            
        for movie in inputMovies:
            movieId = movie.getObjId()
            if movieId in movieDict:
                # Parse the alignment parameters and store the log files
                alignedMovie = movie.clone()
                logFileSrc = movieDict[movieId]
                alignment = parseMovieAlignment(logFileSrc)
                alignedMovie.setAlignment(alignment)
                movieSet.append(alignedMovie)
            else:
                self.warning("Alignment for movie with id %d was not found. DISCARDED!!!" % movieId)
            
        self._defineOutputs(outputMovies=movieSet)
        self._defineTransformRelation(inputMovies, movieSet)

    # --------------------------- UTILS functions ------------------------------
    def _summary(self):
        summary = []
        return summary    
    
    def _methods(self):
        return []
