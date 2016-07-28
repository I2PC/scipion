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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
from os.path import join, basename, exists
from datetime import datetime

import pyworkflow.object as pwobj
from pyworkflow.protocol.params import PointerParam, BooleanParam, LEVEL_ADVANCED

from pyworkflow.protocol.constants import STEPS_PARALLEL, STATUS_NEW
import pyworkflow.utils as pwutils
from pyworkflow.utils.properties import Message
from pyworkflow.em.data import SetOfMovies, Movie, MovieAlignment, Acquisition
from pyworkflow.em import ImageHandler

from protocol_micrographs import ProtPreprocessMicrographs
from protocol_particles import ProtExtractParticles

PLOT_CART = 0
PLOT_POLAR = 1


class ProtProcessMovies(ProtPreprocessMicrographs):
    """
    Protocol base for processing movies from direct detectors cameras.
    This base class will iterate through the movies (extract them if compressed)
    and call a _processMovie method for each one.
    """
    # Redefine this in subclasses if want to convert the movies to mrc
    # the value should be either 'mrc' or 'mrcs'
    CONVERT_TO_MRC = None

    def __init__(self, **kwargs):
        ProtPreprocessMicrographs.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL
    
    #--------------------------- DEFINE param functions ----------------------
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

    #--------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        # Build the list of all processMovieStep ids by 
        # inserting each of the steps for each movie
        self.insertedDict = {}
        self.samplingRate = self.inputMovies.get().getSamplingRate()
        # FIXME: Not working in scipion-box
        #self.convertStepId = self._insertFunctionStep('convertInputStep')

        movieSteps = self._insertNewMoviesSteps(self.insertedDict,
                                                self.inputMovies.get())
        finalSteps = self._insertFinalSteps(movieSteps)
        self._insertFunctionStep('createOutputStep',
                                 prerequisites=finalSteps, wait=True)

    def _insertFinalSteps(self, deps):
        """ This should be implemented in subclasses"""
        return deps

    def _getFirstJoinStepName(self):
        # This function will be used for streamming, to check which is
        # the first function that need to wait for all micrographs
        # to have completed, this can be overriden in subclasses
        # (ej in Xmipp 'sortPSDStep')
        return 'createOutputStep'

    def _getFirstJoinStep(self):
        for s in self._steps:
            if s.funcName == self._getFirstJoinStepName():
                return s
        return None

    def _loadInputList(self):
        """ Load the input set of movies and create a list. """
        moviesFile = self.inputMovies.get().getFileName()
        self.debug("Loading input db: %s" % moviesFile)
        movieSet = SetOfMovies(filename=moviesFile)
        movieSet.loadAllProperties()
        self.listOfMovies = [m.clone() for m in movieSet]
        self.streamClosed = movieSet.isStreamClosed()
        movieSet.close()
        self.debug("Closed db.")

    def _checkNewInput(self):
        # Check if there are new movies to process from the input set
        localFile = self.inputMovies.get().getFileName()
        now = datetime.now()
        self.lastCheck = getattr(self, 'lastCheck', now)
        mTime = datetime.fromtimestamp(os.path.getmtime(localFile))
        self.debug('Last check: %s, modification: %s'
                  % (pwutils.prettyTime(self.lastCheck),
                     pwutils.prettyTime(mTime)))
        # If the input movies.sqlite have not changed since our last check,
        # it does not make sense to check for new input data
        if self.lastCheck > mTime and hasattr(self, 'listOfMovies'):
            return None

        self.lastCheck = now
        # Open input movies.sqlite and close it as soon as possible
        self._loadInputList()
        newMovies = any(m.getObjId() not in self.insertedDict
                        for m in self.listOfMovies)
        outputStep = self._getFirstJoinStep()

        if newMovies:
            fDeps = self._insertNewMoviesSteps(self.insertedDict,
                                               self.listOfMovies)
            if outputStep is not None:
                outputStep.addPrerequisites(*fDeps)
            self.updateSteps()

    def _checkNewOutput(self):
        pass # To be implemented in sub-classes

    def _stepsCheck(self):
        # Input movie set can be loaded or None when checked for new inputs
        # If None, we load it
        self._checkNewInput()
        self._checkNewOutput()

    def _insertNewMoviesSteps(self, insertedDict, inputMovies):
        """ Insert steps to process new movies (from streaming)
        Params:
            insertedDict: contains already processed movies
            inputMovies: input movies set to be check
        """
        deps = []
        # For each movie insert the step to process it
        for movie in inputMovies:
            if movie.getObjId() not in insertedDict:
                stepId = self._insertMovieStep(movie)
                deps.append(stepId)
                insertedDict[movie.getObjId()] = stepId
        return deps
        
    def _insertMovieStep(self, movie):
        """ Insert the processMovieStep for a given movie. """
        # Note1: At this point is safe to pass the movie, since this
        # is not executed in parallel, here we get the params
        # to pass to the actual step that is gone to be executed later on
        # Note2: We are serializing the Movie as a dict that can be passed
        # as parameter for a functionStep
        movieDict = movie.getObjDict(includeBasic=True)
        movieStepId = self._insertFunctionStep('processMovieStep',
                                               movieDict,
                                               movie.hasAlignment(),
                                               prerequisites=[])
        return movieStepId

    #--------------------------- STEPS functions -----------------------------
    def convertInputStep(self):
        """ Should be implemented in sub-classes if needed. """
        pass

    def processMovieStep(self, movieDict, hasAlignment):
        movie = Movie()
        movie.setAcquisition(Acquisition())

        if hasAlignment:
            movie.setAlignment(MovieAlignment())

        movie.setAttributesFromDict(movieDict, setBasic=True)

        movieFolder = self._getOutputMovieFolder(movie)
        movieFn = movie.getFileName()
        movieName = basename(movieFn)

        # Clean old finished files
        pwutils.cleanPath(self._getMovieDone(movie))

        if self._filterMovie(movie):
            pwutils.makePath(movieFolder)
            pwutils.createLink(movieFn, join(movieFolder, movieName))

            if movieName.endswith('bz2'):
                newMovieName = movieName.replace('.bz2', '')
                # We assume that if compressed the name ends with .mrc.bz2
                if not exists(newMovieName):
                    self.runJob('bzip2', '-d -f %s' % movieName, cwd=movieFolder)

            elif movieName.endswith('tbz'):
                newMovieName = movieName.replace('.tbz', '.mrc')
                # We assume that if compressed the name ends with .tbz
                if not exists(newMovieName):
                    self.runJob('tar', 'jxf %s' % movieName, cwd=movieFolder)

            elif movieName.endswith('.tif'):
                #FIXME: It seems that we have some flip problem with compressed
                # tif files, we need to check that
                newMovieName = movieName.replace('.tif', '.mrc')
                # we assume that if compressed the name ends with .tbz
                if not exists(newMovieName):
                    self.runJob('tif2mrc', '%s %s' % (movieName, newMovieName),
                                                      cwd=movieFolder)
            elif movieName.endswith('.txt'):
                # Support a list of frame as a simple .txt file containing
                # all the frames in a raw list, we could use a xmd as well,
                # but a plain text was choose to simply its generation
                movieTxt = os.path.join(movieFolder, movieName)
                with open(movieTxt) as f:
                    movieOrigin = os.path.basename(os.readlink(movieFn))
                    newMovieName = movieName.replace('.txt', '.mrcs')
                    ih = ImageHandler()
                    for i, line in enumerate(f):
                        if line.strip():
                            inputFrame = os.path.join(movieOrigin, line.strip())
                            ih.convert(inputFrame,
                                       (i+1, os.path.join(movieFolder, newMovieName)))
            else:
                newMovieName = movieName
            
            if (self.CONVERT_TO_MRC and not (newMovieName.endswith("mrc") or
                                             newMovieName.endswith("mrcs"))):
                inputMovieFn = os.path.join(movieFolder, newMovieName)
                if inputMovieFn.endswith('.em'):
                    inputMovieFn += ":ems"
                newMovieName = pwutils.replaceExt(newMovieName,
                                                  self.CONVERT_TO_MRC)
                outputMovieFn = os.path.join(movieFolder, newMovieName)
                self.info("Converting movie '%s' -> '%s'" % (inputMovieFn,
                                                             outputMovieFn))
                ImageHandler().convertStack(inputMovieFn, outputMovieFn)

            # Just store the original name in case it is needed in _processMovie
            movie._originalFileName = pwobj.String(objDoStore=False)
            movie._originalFileName.set(movie.getFileName())
            # Now set the new filename (either linked or converted)
            movie.setFileName(os.path.join(movieFolder, newMovieName))
            self.info("Processing movie: %s" % movie.getFileName())

            self._processMovie(movie)

            if self.cleanMovieData:
                self.info("Erasing.....movieFolder: %s" % movieFolder)
                os.system('rm -rf %s' % movieFolder)
                # cleanPath(movieFolder)
            else:
                self.info('Clean movie data DISABLED. '
                          'Movie folder will remain in disk!!!')

        # Mark this movie as finished
        open(self._getMovieDone(movie), 'w').close()
        
    #--------------------------- UTILS functions ----------------------------
    def _getOutputMovieFolder(self, movie):
        """ Create a Movie folder where to work with it. """
        return self._getTmpPath('movie_%06d' % movie.getObjId())

    def _getMovieName(self, movie, ext='.mrc'):
        return self._getExtraPath('movie_%06d%s' % (movie.getObjId(), ext))

    def _getMovieDone(self, movie):
        return self._getExtraPath('DONE_movie_%06d.TXT' % movie.getObjId())

    def _isMovieDone(self, movie):
        """ A movie is done if the marker file exists. """
        return os.path.exists(self._getMovieDone(movie))

    def _getAllDone(self):
        return self._getExtraPath('DONE_all.TXT')

    def _readDoneList(self):
        """ Read from a text file the id's of the items that have been done. """
        doneFile = self._getAllDone()
        doneList = []
        # Check what items have been previously done
        if os.path.exists(doneFile):
            with open(doneFile) as f:
                doneList += [int(line.strip()) for line in f]

        return doneList

    def _writeDoneList(self, movieList):
        """ Write to a text file the items that have been done. """
        doneFile = self._getAllDone()
        with open(self._getAllDone(), 'a') as f:
            for movie in movieList:
                f.write('%d\n' % movie.getObjId())

    #--------------------------- OVERRIDE functions --------------------------
    def _filterMovie(self, movie):
        """ Check if process or not this movie.
        """
        return True

    def _processMovie(self, movie):
        """ Process the movie actions, remember to:
        1) Generate all output files inside movieFolder
           (usually with cwd in runJob)
        2) Copy the important result files after processing
           (movieFolder will be deleted!!!)
        """
        pass

    # FIXME: check if the following functions could be removed

    def _getNameExt(self, movieName, postFix, ext):
        if movieName.endswith("bz2"):
            # removeBaseExt function only eliminate the last extension,
            # but if files are compressed, we need to eliminate two extensions:
            # bz2 and its own image extension (e.g: mrcs, em, etc)
            return pwutils.removeBaseExt(pwutils.removeBaseExt(movieName)) + postFix + '.' + ext
        else:
            return pwutils.removeBaseExt(movieName) + postFix + '.' + ext
    
    def _getCorrMovieName(self, movieId, ext='.mrcs'):
        return 'movie_%06d%s' % (movieId, ext)

    def _getLogFile(self, movieId):
        return 'micrograph_%06d_Log.txt' % movieId


class ProtExtractMovieParticles(ProtExtractParticles, ProtProcessMovies):
    """ Extract a set of Particles from each frame of a set of Movies.
    """
    pass

    

