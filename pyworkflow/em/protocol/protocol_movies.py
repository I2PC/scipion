# **************************************************************************
# *
# * Authors:     Airen Zaldivar Peraza (azaldivar@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
# *              J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Vahid Abrishami (vabrishami@cnb.csic.es)
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
"""
In this module are protocol base classes related to EM Micrographs
"""

from os.path import join, basename, exists

from pyworkflow.protocol.params import PointerParam, BooleanParam, LEVEL_ADVANCED
from pyworkflow.protocol.constants import STEPS_PARALLEL
from pyworkflow.utils.path import createLink, removeBaseExt, makePath, cleanPath
from pyworkflow.utils.properties import Message
from pyworkflow.em.convert import ImageHandler

from protocol_micrographs import ProtPreprocessMicrographs
from protocol_particles import ProtExtractParticles
#from pyworkflow.em.packages.dosefgpu.protocol_dosefgpu import ProtDosefGpu
#from pyworkflow.em.protocol import ProtProcessMovies
#from pyworkflow.utils import startDebugger

PLOT_CART = 0
PLOT_POLAR = 1


class ProtProcessMovies(ProtPreprocessMicrographs):
    """
    Protocol base for processing movies from direct detectors cameras.
    This base class will iterate through the movies (extract them if compressed)
    and call a _processMovie method for each one.
    """
    
    def __init__(self, **kwargs):
        ProtPreprocessMicrographs.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        
        form.addParam('inputMovies', PointerParam, pointerClass='SetOfMovies', 
                      important=True,
                      label=Message.LABEL_INPUT_MOVS,
                      help='Select a set of previously imported movies.')
        form.addParam('cleanMovieData', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label='Clean movie data?',
                      help='Movies take a lot of disk space.\n'
                           'So, by default, all protocols that work on\n'
                           'movies will erase the each movie folder after\n'
                           'finishing. Results files should be copied in \n'
                           'the "createOutputStep" function.\n'
                           'If set to *No*, the folder will not be deleted.')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        # Build the list of all processMovieStep ids by 
        # inserting each of the steps for each movie
        self.samplingRate = self.inputMovies.get().getSamplingRate()
        allMovies = [self._insertMovieStep(movie) for movie in self.inputMovies.get()]
        self._insertFunctionStep('createOutputStep', prerequisites=allMovies)

    def _insertMovieStep(self, movie):
        """ Insert the processMovieStep for a given movie. """
        # Note that at this point is safe to pass the movie, since this
        # is not executed in parallel, here we get the params
        # to pass to the actual step that is gone to be executed later on
        movieStepId = self._insertFunctionStep('processMovieStep',
                                               movie.getObjId(), movie.getFileName(),
                                               prerequisites=[])  
        
        return movieStepId      
        
        
    #--------------------------- STEPS functions ---------------------------------------------------
    def processMovieStep(self, movieId, movieFn, *args):
        movieFolder = self._getMovieFolder(movieId)
        movieName = basename(movieFn)
        #export SCIPION_DEBUG=1 # passwd=a
        #startDebugger()

        if self._filterMovie(movieId, movieFn):
            makePath(movieFolder)
            createLink(movieFn, join(movieFolder, movieName))
            toDelete = [movieName]
    
            if movieName.endswith('bz2'):
                movieMrc = movieName.replace('.bz2', '') # we assume that if compressed the name ends with .mrc.bz2
                toDelete.append(movieMrc)
                if not exists(movieMrc):
                    self.runJob('bzip2', '-d -f %s' % movieName, cwd=movieFolder)
            else:
                movieMrc = movieName
            
            self.info("Processing movie: %s" % movieMrc)
            
            if movieMrc.endswith('.em'):
                movieMrc = movieMrc + ":ems"

            self._processMovie(movieId, movieMrc, movieFolder, *args)
            
            if self.cleanMovieData:
                cleanPath(movieFolder)
            else:
                self.info('Clean movie data DISABLED. Movie folder will remain in disk!!!')
        
    #--------------------------- UTILS functions ---------------------------------------------------
    def _filterMovie(self, movieId, movieFn):
        """ Check if process or not this movie.
        """
        return True
    
    def _getMovieFolder(self, movieId):
        """ Create a Movie folder where to work with it. """
        return self._getTmpPath('movie_%06d' % movieId)  
                
    def _getMovieName(self, movieId, ext='.mrc'):
        return self._getExtraPath('movie_%06d%s' % (movieId, ext)) 
    
    def _getNameExt(self, movieName, postFix, ext):
        if movieName.endswith("bz2"):
            # removeBaseExt function only eliminate the last extension,
            # but if files are compressed, we need to eliminate two extensions:
            # bz2 and its own image extension (e.g: mrcs, em, etc)
            return removeBaseExt(removeBaseExt(movieName)) + postFix + '.' + ext
        else:
            return removeBaseExt(movieName) + postFix + '.' + ext
    
    def _getPlotName(self, movieName, plotType):
        if plotType == PLOT_CART:
            return removeBaseExt(movieName) + '_plot_cart.png'
        else:
            return removeBaseExt(movieName) + '_plot_polar.png'

    def _getCorrMovieName(self, movieId, ext='.mrcs'):
        return 'movie_%06d%s' % (movieId, ext)

    def _getLogFile(self, movieId):
        return 'micrograph_%06d_Log.txt' % movieId

    def _processMovie(self, movieId, movieName, movieFolder,shifts):
        """ Process the movie actions, remember to:
        1) Generate all output files inside movieFolder (usually with cwd in runJob)
        2) Copy the important result files after processing (movieFolder will be deleted!!!)
        """
        pass


class ProtExtractMovieParticles(ProtExtractParticles, ProtProcessMovies):
    """ Extract a set of Particles from each frame of a set of Movies.
    """
    pass

         

