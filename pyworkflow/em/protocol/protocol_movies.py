# **************************************************************************
# *
# * Authors:     Airen Zaldivar Peraza (azaldivar@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
# *              J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

from pyworkflow.protocol.params import PointerParam, IntParam, BooleanParam, LEVEL_EXPERT
from pyworkflow.utils.path import copyFile, removeBaseExt, makePath, cleanPath, moveFile
from pyworkflow.utils.properties import Message
from pyworkflow.em.data import Micrograph

from protocol_micrographs import ProtPreprocessMicrographs
from protocol_particles import ProtExtractParticles



class ProtProcessMovies(ProtPreprocessMicrographs):
    """
    Protocol base for processing movies from direct detectors cameras.
    This base class will iterate throught the movies (extract them if compressed)
    and call a _processMovie method for each one.
    """

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        
        form.addParam('inputMovies', PointerParam, pointerClass='SetOfMovies', 
                      important=True,
                      label=Message.LABEL_INPUT_MOVS,
                      help='')
        form.addParam('cleanMovieData', BooleanParam, default=True,
                      expertLevel=LEVEL_EXPERT,
                      label='Clean movie data?',
                      help='Movies take a lot of disk space.\n'
                           'So, by default, all protocols that work on\n'
                           'movies will erase the each movie folder after\n'
                           'finishing. Results files should be copied in \n'
                           'the "createOutputStep" function.\n'
                           'If set to *No*, the folder will not be deleted.')
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        
        for movie in self.inputMovies.get():
            self._insertFunctionStep('processMovieStep', 
                                     movie.getObjId(), movie.getFileName())
        self._insertFunctionStep('createOutputStep')
        self._insertFunctionStep('cleanMovieDataStep')

    #--------------------------- STEPS functions ---------------------------------------------------
    def processMovieStep(self, movieId, movieFn):
        movieFolder = self._getMovieFolder(movieId)
        movieName = basename(movieFn)
        
        makePath(movieFolder)
        copyFile(movieFn, join(movieFolder, movieName))
        
        self._enterDir(movieFolder)

        if movieName.endswith('bz2'):
            movieMrc = movieName.replace('.bz2', '') # we assume that if compressed the name ends with .mrc.bz2
            if not exists(movieMrc):
                self.runJob('bzip2', '-d %s' % movieName)
        else:
            movieMrc = movieName

        self.info("Processing movie: %s" % movieMrc)            
        self._processMovie(movieId, movieMrc)
        
        self._leaveDir()
        
    def cleanMovieDataStep(self):
        if self.cleanMovieData:
            for movie in self.inputMovies.get():
                movieFolder = self._getMovieFolder(movie.getObjId())
                self.info('Cleanning folder: %s' % movieFolder)
                cleanPath(movieFolder)
        else:
            self.info('Clean movie data DISABLED. Movie folder will remain in disk!!!')

    #--------------------------- UTILS functions ---------------------------------------------------
    def _getMovieFolder(self, movieId):
        """ Create a Movie folder where to work with it. """
        return self._getExtraPath('movie_%06d' % movieId)  
                
    def _getMovieName(self, movieId, ext='mrc'):
        return 'movie_%06d%s' % (movieId, ext) 
    
    def _getMicName(self, movieId):
        return 'micrograph_%06d.mrc' % movieId
    
    
class ProtAverageMovies(ProtProcessMovies):
    """ Simple protocol to just average all frames
    in a movie and produce a micrograph per movie.
    """
    _label = 'movie average'
    
    def _processMovie(self, movieId, movieName):
        import xmipp
        movieImg = xmipp.Image()
        movieImg.read(movieName, xmipp.HEADER)
    
        x, y, z, n = movieImg.getDimensions()
    
        # Store all the particles
        img = xmipp.Image()
        sumImg = xmipp.Image()
        
        for i in range(n):
            img.read('%d@%s' % (i+1, movieName))
            if i:
                sumImg += img
            else:
                sumImg = img
        #sumImg *= (1./n)
        sumImg.write(self._getMicName(movieId))
        
    def createOutputStep(self):
        inputMovies = self.inputMovies.get()
        micSet = self._createSetOfMicrographs()
        micSet.copyInfo(inputMovies)
        
        # Create a folder to store the resulting micrographs
        micsFolder = self._getPath('micrographs')
        makePath(micsFolder)
        
        for movie in self.inputMovies.get():
            movieId = movie.getObjId()
            micName = self._getMicName(movieId)
            # Move the resulting micrograph before delete of movies folder
            micNameSrc = join(self._getMovieFolder(movieId), micName)
            micNameDst = join(micsFolder, micName)
            moveFile(micNameSrc, micNameDst)
            
            mic = micSet.ITEM_TYPE()
            mic.setFileName(micNameDst)
            micSet.append(mic)
            
        self._defineOutputs(outputMicrographs=micSet)
        self._defineTransformRelation(inputMovies, micSet)
            
    
class ProtOpticalAlignment(ProtProcessMovies):
    """ Aligns movies, from direct detectors cameras, into micrographs.
    """
    _label = 'movie alignment' 
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        
        form.addParam('inputMovies', PointerParam, important=True,
                      label=Message.LABEL_INPUT_MOVS, pointerClass='SetOfMovies')
        form.addParam('doGPU', BooleanParam, default=False,
                      label="Use GPU (vs CPU)", help="Set to true if you want the GPU implementation")
        form.addParam('GPUCore', IntParam, default=0,
                      label="Choose GPU core",
                      condition="doGPU",
                      help="GPU may have several cores. Set it to zero if you do not know what we are talking about")
        form.addParam('winSize', IntParam, default=150,
                      label="Window size", expertLevel=LEVEL_EXPERT,
                      help="Window size (shifts are assumed to be constant within this window).")
        line = form.addLine('Frames rangeDrop Frames (NOT IMPLEMENTED):',
                      help='Drop first and last frames. set to 0 in orser to keep all\n\n'
                           'NOT IMPLEMENTED YET.')
        line.addParam('firstFrame', IntParam, default='0',
                      label='First')
        line.addParam('lastFrame',  IntParam, default='0',
                      label='Last')
        form.addParallelSection(threads=1, mpi=1)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        movSet = self.inputMovies.get()
        for mov in movSet:
            #TODO What milks is this?
            mov.load()
            movFn = mov.getFirstItem().getFileName()
            self._insertFunctionStep('processMoviesStep', movFn)
        self._insertFunctionStep('createOutputStep')
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def processMoviesStep(self, movFn):
        movName = removeBaseExt(movFn)
        self._createMovWorkingDir(movName)
        self._defineProgram()
        self._micList = []
        micFn = self._getExtraPath( movName+ "_aligned.spi")
        args = '-i %s -o %s --winSize %d'%(movFn, micFn, self.winSize.get())
        if False:
            args += ' --dropFirst %d --dropLast %d' % (self.firstFrames.get(),self.lastFrames.get())
        if self.doGPU:
            args += ' --gpu %d' % self.GPUCore.get()
        self.runJob(self._program, args)

        # micJob = join(movDir, "justtest.mrc")
        # micFn = self._getExtraPath(movName + ".mrc")
        # moveFile(micJob, micFn)
        self._micList.append(micFn)
    
    def createOutputStep(self):
        micSet = self._createSetOfMicrographs()
        movSet = self.inputMovies.get()
        micSet.setAcquisition(movSet.getAcquisition())
        micSet.setSamplingRate(movSet.getSamplingRate())
        
        for m in self._micList:
            mic = Micrograph()
            mic.setFileName(m)
            micSet.append(mic)
        self._defineOutputs(outputMicrographs=micSet)
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _createMovWorkingDir(self, movFn):
        """create a new directory for the movie and change to this directory.
        """
        workDir = self._movWorkingDir(movFn)
        makePath(workDir)   # Create a directory for a current iteration
    
    def _movWorkingDir(self, movFn):
        """ Define which is the directory for the current movie"""
        workDir = self._getTmpPath(movFn)
        return workDir
    
    def _defineProgram(self):
        if self.doGPU:
            self._program = 'xmipp_optical_alignment_gpu'
        else:
            self._program = 'xmipp_optical_alignment_cpu'
    
    
class ProtExtractMovieParticles(ProtExtractParticles, ProtProcessMovies):
    """ Extract a set of Particles from each frame of a set of Movies.
    """
    _label = 'movie extract particles' 
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        
        form.addParam('inputMovies', PointerParam, pointerClass='SetOfMovies',
                      important=True,
                      label=Message.LABEL_INPUT_MOVS)
        form.addParam('inputCoordinates', PointerParam, pointerClass='SetOfCoordinates',
                      important=True,
                      label='Input coordinates')
        
        line = form.addLine('Frames range',
                      help='')
        line.addParam('firstFrame', IntParam, default=0,
                      label='First')
        line.addParam('lastFrame',  IntParam, default=0,
                      label='Last')
        form.addParallelSection(threads=1, mpi=1)    
