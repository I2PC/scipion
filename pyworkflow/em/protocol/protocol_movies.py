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

from pyworkflow.protocol.params import PointerParam, IntParam, BooleanParam, LEVEL_EXPERT
from pyworkflow.protocol.constants import STEPS_PARALLEL
from pyworkflow.utils.path import copyFile, removeBaseExt, makePath, cleanPath, moveFile, getExt
from pyworkflow.utils.properties import Message
from pyworkflow.em.data import Micrograph

from protocol_micrographs import ProtPreprocessMicrographs
from protocol_particles import ProtExtractParticles
#from pyworkflow.em.packages.dosefgpu.protocol_dosefgpu import ProtDosefGpu
#from pyworkflow.em.protocol import ProtProcessMovies



class ProtProcessMovies(ProtPreprocessMicrographs):
    """
    Protocol base for processing movies from direct detectors cameras.
    This base class will iterate throught the movies (extract them if compressed)
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
        allMovies = []
        for movie in self.inputMovies.get():
            movieStepId = self._insertFunctionStep('processMovieStep', 
                                                   movie.getObjId(), movie.getFileName(),
                                                   prerequisites=[])
            allMovies.append(movieStepId)
        self._insertFunctionStep('createOutputStep', prerequisites=allMovies)

    #--------------------------- STEPS functions ---------------------------------------------------
    def _filterMovie(self, movieId, movieFn):
        """ Check if process or not this movie.
        """
        return True
    
    def processMovieStep(self, movieId, movieFn):
        movieFolder = self._getMovieFolder(movieId)
        movieName = basename(movieFn)
        
        if self._filterMovie(movieId, movieFn):
            makePath(movieFolder)
            copyFile(movieFn, join(movieFolder, movieName))
            toDelete = [movieName]
    
            if movieName.endswith('bz2'):
                movieMrc = movieName.replace('.bz2', '') # we assume that if compressed the name ends with .mrc.bz2
                toDelete.append(movieMrc)
                if not exists(movieMrc):
                    self.runJob('bzip2', '-d %s' % movieName)
            else:
                movieMrc = movieName
    
            self.info("Processing movie: %s" % movieMrc) 
            
            self._processMovie(movieId, movieName, movieFolder)
            
            if self.cleanMovieData:
                cleanPath(movieFolder)
            else:
                self.info('Clean movie data DISABLED. Movie folder will remain in disk!!!')

        
    #--------------------------- UTILS functions ---------------------------------------------------
    def _getMovieFolder(self, movieId):
        """ Create a Movie folder where to work with it. """
        return self._getTmpPath('movie_%06d' % movieId)  
                
    def _getMovieName(self, movieId, ext='.mrc'):
        return 'movie_%06d%s' % (movieId, ext) 
    
    def _getMicName(self, movieId):
        return 'micrograph_%06d.mrc' % movieId

    def _getCorrMovieName(self, movieId, ext='.mrcs'):
        return 'movie_%06d%s' % (movieId, ext)

    def _getLogFile(self, movieId):
        return 'micrograph_%06d_Log.txt' % movieId

    def _processMovie(self, movieId, movieName, movieFolder):
        """ Process the movie actions, remember to:
        1) Generate all output files inside movieFolder (usually with cwd in runJob)
        2) Copy the important result files after processing (movieFolder will be deleted!!!)
        """
        pass
    
    
class ProtAverageMovies(ProtProcessMovies):
    """ Simple protocol to just average all frames
    in a movie and produce a micrograph per movie.
    """
    _label = 'movie average'
    
    def _processMovie(self, movieId, movieName, movieFolder):
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
            
    
class ProtExtractMovieParticles(ProtExtractParticles, ProtProcessMovies):
    """ Extract a set of Particles from each frame of a set of Movies.
    """
    _label = 'movie extract particles' 
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        ProtProcessMovies._defineParams(self, form)
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
                     
    def _processMovie(self, movieId, movieName):
        import xmipp
        from pyworkflow.em.packages import xmipp3
        
        boxSize = self.inputCoordinates.get().getBoxSize()
        # Aux image to read and write images        
        img = xmipp.Image()
        
        # Read movie dimensions to iterate through each frame
        movieImg = xmipp.Image()
        movieImg.read(movieName, xmipp.HEADER)
        x, y, z, n = movieImg.getDimensions()
        
        first = self.firstFrame.get()
        last = self.lastFrame.get()
        if last == 0:
            last = n-1
        # Get the shift if the is alignment
        movie = self.inputMovies.get()[movieId]
        if movie.hasAlignment():
            shifts = movie.getAlignment().getShifts()
        else:
            shifts = [0] * (2*n)
            
        stkIndex = 0
        movieStk = self._getMovieName(movieId, '.stk')
        movieMdFile = self._getMovieName(movieId, '.xmd')
        movieMd = xmipp.MetaData()
        frameMd = xmipp.MetaData()
        frameMdImages = xmipp.MetaData()
        frameRow = xmipp3.XmippMdRow()
        
        for frame in range(first, last+1):
            i = frame + 1
            # Get the frame shifts
            index = frame - first
            shiftX = shifts[2*index]
            shiftY = shifts[2*index+1]
             
            frameRoot = 'frame_%02d' % frame
            frameName = frameRoot + '.mrc'
            frameMdFile = frameRoot + '.xmd'
            framePosFile = frameRoot + '_coordinates.xmd'
            coordinatesName = frameRoot + '_coordinates.xmd'
            
            
            hasCoordinates = self._writeXmippPosFile(movieId, movieName, coordinatesName, 
                                                     shiftX, shiftY)
            
            if hasCoordinates:            
                self.info("Writing frame: %s" % frameName)
                img.read('%d@%s' % (i, movieName))
                img.write(frameName)
                
                self.info("Extracting particles")
                frameImages = frameRoot + '_images'
                args = '-i %(frameName)s --pos %(coordinatesName)s -o %(frameRoot)s --Xdim %(boxSize)d' % locals()
                self.runJob('xmipp_micrograph_scissor', args, env=xmipp3._environ)
                cleanPath(frameName)
                
                frameStk = frameRoot + '.stk'
                #img.read(frameStk, xmipp.HEADER)
                
                #x2, y2, z2, n2 = img.getDimensions()
                
                self.info("Combining particles into one stack.")
                
                frameMdImages.read(frameMdFile)
                frameMd.read('particles@%s' % framePosFile)
                frameMd.merge(frameMdImages)
                #frameMd.removeDisabled()
                
                for objId in frameMd:
                    stkIndex += 1
                    frameRow.readFromMd(frameMd, objId)
                    img.read(frameRow.getValue(xmipp.MDL_IMAGE))
                    img.write('%d@%s' % (stkIndex, movieStk))
                    # Fix the name to be accesible from the Project directory
                    # so we know that the movie stack file will be moved
                    # to final particles folder
                    newImageName = '%d@%s' % (stkIndex, self._getPath('particles', movieStk))
                    frameRow.setValue(xmipp.MDL_IMAGE, newImageName)
                    frameRow.setValue(xmipp.MDL_MICROGRAPH_ID, long(movieId))
                    frameRow.setValue(xmipp.MDL_MICROGRAPH, str(movieId))
                    frameRow.setValue(xmipp.MDL_FRAME_ID, long(frame))
                    frameRow.setValue(xmipp.MDL_PARTICLE_ID, frameRow.getValue(xmipp.MDL_ITEM_ID))
                    frameRow.writeToMd(movieMd, movieMd.addObject())
                movieMd.addItemId()
                movieMd.write(movieMdFile)
                cleanPath(frameStk)
            
    def createOutputStep(self):
        import xmipp
        from pyworkflow.em.packages import xmipp3
        
        inputMovies = self.inputMovies.get()
        particleSet = self._createSetOfMovieParticles()        
        particleSet.copyInfo(inputMovies)
        
        # Create a folder to store the resulting micrographs
        particleFolder = self._getPath('particles')
        makePath(particleFolder)
        md = xmipp.MetaData()
        mdAll = xmipp.MetaData()
        
        for movie in inputMovies:
            movieId = movie.getObjId()
            movieName = self._getMovieName(movieId)
            movieFolder = self._getMovieFolder(movieId)
            movieStk = movieName.replace('.mrc', '.stk')
            movieMdFile = movieName.replace('.mrc', '.xmd') 
            # Store particle stack and metadata files in final particles folder
            movieStkSrc = join(movieFolder, movieStk)
            if exists(movieStkSrc):
                moveFile(movieStkSrc, particleFolder)
                moveFile(join(movieFolder, movieMdFile), particleFolder)
                # Read particles from each movie and join in a full set
                md.read(join(particleFolder, movieMdFile))
                mdAll.unionAll(md)
            
        particleMd = self._getPath('movie_particles.xmd')
        mdAll.addItemId()
        mdAll.write(particleMd)
        xmipp3.readSetOfParticles(particleMd, particleSet, removeDisabled=False)
        
        self._defineOutputs(outputParticles=particleSet)
        self._defineSourceRelation(inputMovies, particleSet)
    
    def _filterMovie(self, movieId, movieFn):
        """ Check if process or not this movie.
        In this case if there are not coordinates, do not process it.
        """
        coordinates = self.inputCoordinates.get()
        micrograph = coordinates.getMicrographs()[movieId]
        
        if micrograph is None:
            self.warning("Micrograph with id %d was not found in SetOfCoordinates!!!" % movieId)
            return False
        
        n = len([coord for coord in coordinates.iterCoordinates(micrograph)])
        
        if n == 0:
            self.warning("Skipping micrograph with id %d, not coordinates found!!!" % movieId)
            return False
        
        return True
        
        
    def _writeXmippPosFile(self, movieId, movieName, coordinatesName, shiftX, shiftY):
        """ Create an Xmipp coordinates files to be extracted
        from the frames of the movie.
        """
        import xmipp
        from pyworkflow.em.packages.xmipp3 import coordinateToRow, XmippMdRow
        
        coordinates = self.inputCoordinates.get()
        micrograph = coordinates.getMicrographs()[movieId]
        
        if micrograph is None:
            raise Exception("Micrograph with id %d was not found in SetOfCoordinates!!!" % movieId)
        
        md = xmipp.MetaData()
        coordRow = XmippMdRow()
        
        for coord in coordinates.iterCoordinates(micrograph):
            coord.shiftX(1 * int(round(float(shiftX))))
            coord.shiftY(1 * int(round(float(shiftY))))
            coordinateToRow(coord, coordRow)
            coordRow.writeToMd(md, md.addObject())
            
        if md.isEmpty():
            return False
        else:
            self.info("Writing coordinates metadata: %s, with shifts: %s %s" % (coordinatesName, shiftX, shiftY))
            md.write('particles@' + coordinatesName)
            return True
        
        
         
