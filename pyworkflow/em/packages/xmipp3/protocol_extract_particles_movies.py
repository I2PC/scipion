# *****************************************************************************
# *
# * Authors:     Josue Gomez Blanco (jgomez@cnb.csic.es)
# *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# *
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
# *****************************************************************************

import os

import pyworkflow.em.metadata as md

from pyworkflow.em import SetOfCoordinates
from pyworkflow.em.packages.xmipp3.convert import (readSetOfMovieParticles,
                                                   xmippToLocation)
from pyworkflow.em.convert import ImageHandler
from pyworkflow.em.protocol import ProtExtractMovieParticles
from pyworkflow.protocol.constants import LEVEL_ADVANCED, STEPS_PARALLEL
from pyworkflow.protocol.params import (PointerParam, IntParam, BooleanParam,
                                        Positive, FloatParam, EnumParam)

from pyworkflow.utils.path import cleanPath

from pyworkflow.em.packages.xmipp3 import coordinateToRow, XmippMdRow



class XmippProtExtractMovieParticles(ProtExtractMovieParticles):
    """ Extract a set of Particles from each frame of a set of Movies.
    """
    _label = 'movie extract particles'

    def __init__(self, **kwargs):
        ProtExtractMovieParticles.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL
    
    
    
    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        self.movieFolder = self._getTmpPath('movie_%(movieId)06d/')
        self.frameRoot = self.movieFolder + 'frame_%(frame)02d'
        myDict = {
                  'frameImages' : self.frameRoot + '_images',
                  'frameMic' : self.frameRoot + '.mrc',
                  'frameMdFile' : self.frameRoot + '_images.xmd',
                  'frameCoords' : self.frameRoot + '_coordinates.xmd',
                  'frameStk' : self.frameRoot + '_images.stk'
                 }

        self._updateFilenamesDict(myDict)

    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        ProtExtractMovieParticles._defineParams(self, form)
        form.addParam('inputCoordinates', PointerParam,
                      pointerClass='SetOfCoordinates',
                      important=True,
                      label='Input coordinates')
        form.addParam('boxSize', IntParam, default=0,
                      label='Particle box size', validators=[Positive],
                      help='In pixels. The box size is the size of the boxed '
                           'particles, actual particles may be smaller than '
                           'this.')
        form.addParam('applyAlignment', BooleanParam, default=False,
                      label='Apply movie alignments to extract?',
                      help='If the input movies contains frames alignment, '
                           'you decide whether to use that information '
                           'for extracting the particles taking into account '
                           'the shifts between frames.')
        line = form.addLine('Frames range',
                      help='Specify the frames range to extract particles. '
                           'The first frame is 1. If you set 0 in the '
                           'last frame, it means that you will use until the '
                           'last frame of the movie. If you apply the '
                           'previous alignment of the movies, you only can use '
                           'a frame range equal or less as used to alignment.')
        line.addParam('frame0', IntParam, label='First')
        line.addParam('frameN',  IntParam, label='Last')
        
        form.addSection(label='Preprocess')
        form.addParam('doRemoveDust', BooleanParam, default=True,
                      important=True,
                      label='Dust removal (Recommended)', 
                      help='Sets pixels with unusually large values to random '
                           'values from a Gaussian with zero-mean and '
                           'unity-standard deviation.')
        form.addParam('thresholdDust', FloatParam, default=3.5,
                      condition='doRemoveDust', expertLevel=LEVEL_ADVANCED,
                      label='Threshold for dust removal',
                      help='Pixels with a signal higher or lower than this '
                           'value times the standard deviation of the image '
                           'will be affected. For cryo, 3.5 is a good value. '
                           'For high-contrast negative stain, the signal '
                           'itself may be affected so that a higher value may '
                           'be preferable.')
        form.addParam('doInvert', BooleanParam, default=False,
                      label='Invert contrast', 
                      help='Invert the contrast if your particles are black '
                           'over a white background.')
        form.addParam('doNormalize', BooleanParam, default=True,
                      label='Normalize (Recommended)', 
                      help='It subtract a ramp in the gray values and '
                           'normalizes so that in the  background there is 0 '
                           'mean and standard deviation 1.')
        form.addParam('normType', EnumParam,
                      choices=['OldXmipp','NewXmipp','Ramp'], default=2,
                      condition='doNormalize', expertLevel=LEVEL_ADVANCED,
                      display=EnumParam.DISPLAY_COMBO,
                      label='Normalization type', 
                      help='OldXmipp (mean(Image)=0, stddev(Image)=1). \n'
                           'NewXmipp (mean(background)=0, '
                           'stddev(background)=1)\n'
                           'Ramp (subtract background+NewXmipp).')
        form.addParam('backRadius', IntParam, default=-1,
                      condition='doNormalize',
                      label='Background radius',
                      help='Pixels outside this circle are assumed to be '
                           'noise and their stddev is set to 1. Radius for '
                           'background circle definition (in pix.). If this '
                           'value is 0, then half the box size is used.',
                      expertLevel=LEVEL_ADVANCED)

        form.addParallelSection(threads=3, mpi=1)
    #--------------------------- INSERT steps functions ------------------------

    def _insertAllSteps(self):
        #ROB: deal with the case in which sampling rate use for picking and
        #  movies is different
        
        inputCoords = self.inputCoordinates.get()
        #coordinates sampling mic used for picking
        samplingCoords = inputCoords.getMicrographs().getSamplingRate()
        #coordinates sampling input mic
        samplingMic = self.inputMovies.get().getSamplingRate()
        self.factor = samplingMic / samplingCoords
        
        self._createFilenameTemplates()

        # Build the list of all processMovieStep ids by
        # inserting each of the steps for each movie
        self.insertedDict = {}
        self.samplingRate = self.inputMovies.get().getSamplingRate()
        # FIXME: Not working in scipion-box
        # self.convertStepId = self._insertFunctionStep('convertInputStep')

        movieSteps = self._insertNewMoviesSteps(self.insertedDict,
                                                self.inputMovies.get())
        finalSteps = self._insertFinalSteps(movieSteps)
        self._insertFunctionStep('createOutputStep',
                                 prerequisites=finalSteps, wait=False)

    def _insertMovieStep(self, movie):
        # Redefine this function to add the shifts and factor to the
        # processMovieStep function and run properly in parallel with threads

        #retrive shifts here so there is no conflict
        #if the object is accessed inside at the same time by multiple threads
        movieDict = movie.getObjDict(includeBasic=True)
        movieStepId = self._insertFunctionStep('processMovieStep',
                                               movieDict,
                                               movie.hasAlignment(),
                                               prerequisites=[])
        
        return movieStepId 
    
    #--------------------------- STEPS functions -------------------------------
    def _processMovie(self, movie):
        movId = movie.getObjId()
        x, y, n = movie.getDim()
        iniFrame, lastFrame, _ = movie.getFramesRange()
        frame0, frameN = self._getRange(movie)
        boxSize = self.boxSize.get()
        
        if movie.hasAlignment() and self.applyAlignment:
            shiftX, shiftY = movie.getAlignment().getShifts()  # lists.
        else:
            shiftX = [0] * (lastFrame-iniFrame+1)
            shiftY = shiftX
        
        stkIndex = 0
        movieStk = self._getMovieName(movie, '.stk')
        movieMdFile = self._getMovieName(movie, '.xmd')
        movieMd = md.MetaData()
        frameMd = md.MetaData()
        frameMdImages = md.MetaData()
        frameRow = md.Row()
        
        imgh = ImageHandler()
        for frame in range(frame0, frameN+1):
            indx = frame-iniFrame
            print  "Index: ", indx, shiftX, shiftY
            frameName = self._getFnRelated('frameMic',movId, frame)
            frameMdFile = self._getFnRelated('frameMdFile',movId, frame)
            coordinatesName = self._getFnRelated('frameCoords',movId, frame)
            frameImages = self._getFnRelated('frameImages',movId, frame)
            frameStk = self._getFnRelated('frameStk', movId, frame)
            
            hasCoordinates = self._writeXmippPosFile(movie, coordinatesName,
                                                     shiftX[indx], shiftY[indx])
            
            if hasCoordinates:
                self.info("Writing frame: %s" % frameName)
                #TODO: there is no need to write the frame and then operate
                #the input of the first operation should be the movie
                movieName = imgh.fixXmippVolumeFileName(movie)
                print "Passing Ih: ", frame, movieName, frameName
                imgh.convert(tuple([frame-1, movieName]), frameName)
                
                if self.doRemoveDust:
                    self.info("Removing Dust")
                    self._runNoDust(frameName)

                self.info("Extracting particles")
                args = '-i %(frameName)s --pos %(coordinatesName)s ' \
                       '-o %(frameImages)s --Xdim %(boxSize)d' % locals()
                
                if self.doInvert:
                    args += " --invert"
                
                args += " --downsampling %f " % self.factor
                self.runJob('xmipp_micrograph_scissor', args)
                cleanPath(frameName)
                
                self.info("Combining particles into one stack.")
                
                frameMdImages.read(frameMdFile)
                frameMd.read('particles@%s' % coordinatesName)
                frameMd.merge(frameMdImages)
                 
                for objId in frameMd:
                    stkIndex += 1
                    frameRow.readFromMd(frameMd, objId)
                    location = xmippToLocation(frameRow.getValue(md.MDL_IMAGE))
                    newLocation = (stkIndex, movieStk)
                    imgh.convert(location, newLocation)
                    
                    # Fix the name to be accesible from the Project directory
                    # so we know that the movie stack file will be moved
                    # to final particles folder
                    newImageName = '%d@%s' % newLocation
                    frameRow.setValue(md.MDL_IMAGE, newImageName)
                    frameRow.setValue(md.MDL_MICROGRAPH_ID, long(movId))
                    frameRow.setValue(md.MDL_MICROGRAPH, str(movId))
                    frameRow.setValue(md.MDL_FRAME_ID, long(frame))
                    frameRow.setValue(md.MDL_PARTICLE_ID,
                                      frameRow.getValue(md.MDL_ITEM_ID))
                    frameRow.writeToMd(movieMd, movieMd.addObject())
                movieMd.addItemId()
                movieMd.write(movieMdFile)
                cleanPath(frameStk)
        
        if self.doNormalize:
            numberOfFrames = frameN - frame0 + 1
            self._runNormalize(movieStk, numberOfFrames)
    
    def _runNoDust(self, frameName):
        args = (" -i %s -o %s --bad_pixels outliers %f"
                % (frameName, frameName, self.thresholdDust.get()))
        self.runJob("xmipp_transform_filter", args)

    def _runNormalize(self, stack, frames):
        import math
        program = "xmipp_transform_normalize"
        args = "-i %(stack)s "
        normType = self.normType.get()
        bgRadius = self.backRadius.get()
        
        size, _, _ , _ = ImageHandler().getDimensions(stack)
        if bgRadius <= 0:
            bgRadius = int(size/2)
        
        if normType=="OldXmipp":
            args += "--method OldXmipp"
        elif normType=="NewXmipp":
            args += "--method NewXmipp --background circle %(bgRadius)d"
        else:
            args += "--method Ramp --background circle %(bgRadius)d"
        self.runJob(program, args % locals())
        
        # The background's stddev of the movie
        # particles must be equal to sqrt(NoOfFrames)
        # for the background's stddev of the
        # average particle be equal to 1.
        val = math.sqrt(frames)
        opArgs = "-i %s --mult %f" % (stack, val)
        self.runJob("xmipp_image_operate", opArgs)
    
    def createOutputStep(self):
        inputMovies = self.inputMovies.get()
        particleSet = self._createSetOfMovieParticles()        
        particleSet.copyInfo(inputMovies)
          
        # Create a folder to store the resulting micrographs
         
#         particleFolder = self._getPath('particles')
#         makePath(particleFolder)
        mData = md.MetaData()
        mdAll = md.MetaData()
          
        for movie in inputMovies:
            movieName = self._getMovieName(movie)
            movieStk = movieName.replace('.mrc', '.stk')
            movieMdFile = movieName.replace('.mrc', '.xmd')
            
            # Store particle stack and metadata files in final particles folder
            if os.path.exists(movieStk):
                mData.read(movieMdFile)
                mdAll.unionAll(mData)
              
        particleMd = self._getPath('movie_particles.xmd')
        mdAll.addItemId()
        mdAll.write(particleMd)
        readSetOfMovieParticles(particleMd, particleSet, removeDisabled=False,
                                postprocessImageRow=self._postprocessImageRow)
        
        self._defineOutputs(outputParticles=particleSet)
        self._defineSourceRelation(self.inputMovies, particleSet)
    
    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        
        inputSet = self.inputMovies.get()
        
        # Although getFirstItem is not remonended in general, here it is
        # used olny once, for validation purposes, so performance
        # problems not should be apprear.
        movie = inputSet.getFirstItem()
        if (not movie.hasAlignment()) and self.applyAlignment:
            errors.append("Your movies has not alignment. Please, set *No* "
                          "the parameter _Apply movie alignments to extract?_")

        firstFrame, lastFrame, _ = inputSet.getFramesRange()
        if lastFrame == 0:
            # Although getFirstItem is not remonended in general, here it is
            # used olny once, for validation purposes, so performance
            # problems not should be apprear.
            frames = inputSet.getFirstItem().getNumberOfFrames()
            lastFrame = frames
        else:
            frames = lastFrame - firstFrame + 1
    
        if frames is not None:
            # Avoid validation when the range is not defined
            if not (hasattr(self, 'frame0') or hasattr(self, 'frameN')):
                return
        
            f0, fN = self._getRange(movie)
            if fN < firstFrame or fN > lastFrame:
                errors.append("Check the selected last frame. "
                              "Last frame (%d) should be in range: %s "
                              % (fN, (firstFrame, lastFrame)))
            
            if f0 < firstFrame or f0 > lastFrame:
                errors.append("Check the selected first frame. "
                              "First frame (%d) should be in range: %s "
                              % (f0, (firstFrame, lastFrame)))
            if fN < f0:
                errors.append("Check the selected frames range. Last frame "
                              "(%d) should be greater than first frame (%d)"
                              % (fN, f0))
        return errors

    def _methods(self):
        methods = []
        return methods

    def _summary(self):
        summary = []
        return summary
    
    #--------------------------- UTILS functions -------------------------------
    def _getFnRelated(self, keyFile, movId, frameIndex):
        return self._getFileName(keyFile, movieId=movId, frame=frameIndex)
    
    def _writeXmippPosFile(self, movie, coordinatesName, shiftX, shiftY):
        """ Create an Xmipp coordinates files to be extracted
        from the frames of the movie.
        """
        coordSet = self.inputCoordinates.get()
        
        # to support multiple access to db
        coordToIter = SetOfCoordinates()
        coordToIter.copy(coordSet)
        coordSet.close()

        mData = md.MetaData()
        coordRow = XmippMdRow()

        for coord in coordToIter.iterCoordinates(movie.getObjId()):
            coord.shiftX( int(round(float(shiftX))))
            coord.shiftY( int(round(float(shiftY))))
            coordinateToRow(coord, coordRow)
            coordRow.writeToMd(mData, mData.addObject())
        coordToIter.close()

        if mData.isEmpty():
            return False
        else:
            self.info("Writing coordinates metadata: %s, with shifts: %s %s"
                      % (coordinatesName, shiftX, shiftY))
            mData.write('particles@' + coordinatesName)
            return True
    
    def _filterMovie(self, movie):
        """ Check if process or not this movie. If there are coordinates or
        not to this movie, is cheked in _processMovie step.
        """
        
        movieId = movie.getObjId()
        coordinates = self.inputCoordinates.get()
        micrograph = coordinates.getMicrographs()[movieId]
        
        if micrograph is not None:
            return True
        else:
            self.warning("Micrograph with id %d was not found in "
                         "SetOfCoordinates!!!" % movieId)
            return False
    
    def _postprocessImageRow(self, img, imgRow):
        img.setFrameId(imgRow.getValue(md.MDL_FRAME_ID))
        img.setParticleId(imgRow.getValue(md.MDL_PARTICLE_ID))


    def _useAlignToSum(self):
        return self.getAttributeValue('useAlignToSum', False)

    def _getRange(self, movie):
        n = self._getNumberOfFrames(movie)
        iniFrame, _, indxFrame = movie.getFramesRange()

        first = self.getAttributeValue('frame0')
        last = self.getAttributeValue('frameN')

        if first <= 1:
            first = 1

        if last <= 0:
            last = n

        if iniFrame != indxFrame:
            first -= (iniFrame - 1)
            last -= (iniFrame - 1)

        return first, last

    def _getNumberOfFrames(self, movie):
        _, lstFrame, _ = movie.getFramesRange()

        if movie.hasAlignment():
            _, lastFrmAligned = movie.getAlignment().getRange()
            if lastFrmAligned != lstFrame:
                return lastFrmAligned
            else:
                return movie.getNumberOfFrames()
        else:
            return movie.getNumberOfFrames()
