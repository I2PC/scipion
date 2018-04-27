# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
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
from os.path import relpath

from pyworkflow.protocol.constants import STATUS_FINISHED, STEPS_SERIAL
from pyworkflow.em.protocol import ProtExtractMovieParticles
from pyworkflow.em.constants import ALIGN_PROJ
from pyworkflow.utils.properties import Message
from pyworkflow import VERSION_1_2
import pyworkflow.protocol.params as params
import pyworkflow.em as em
import pyworkflow.em.metadata as md
import pyworkflow.utils as pwutils

from convert import (writeSetOfMovies, isVersion2,
                     writeSetOfParticles, locationToRelion)
from protocol_base import ProtRelionBase


class ProtRelionExtractMovieParticles(ProtExtractMovieParticles,
                                      ProtRelionBase):
    """Protocol to extract particles from a set of coordinates"""
    _label = 'movie particles extraction'
    _lastUpdateVersion = VERSION_1_2

    @classmethod
    def isDisabled(cls):
        return not isVersion2()

    def __init__(self, **kwargs):
        ProtExtractMovieParticles.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_SERIAL
        self.movieDict = {}
        self.ptclDict = {}

    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        myDict = {'inputMovies': 'input_movies.star',
                  'inputParts': 'input_particles.star',
                  'outputDir': pwutils.join('extra', 'output'),
                  'outputParts': self._getExtraPath('output/movie_particles.star')
                  }

        self._updateFilenamesDict(myDict)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputMovies', params.PointerParam,
                      pointerClass='SetOfMovies',
                      important=True,
                      label='Input aligned movies',
                      help='Select a set of ALIGNED movies.')

        form.addParam('inputParticles', params.PointerParam,
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      important=True,
                      label='Input aligned particles',
                      help='Relion requires input particles set, which it will '
                           'use to get the coordinates and X/Y-shifts for '
                           'movie particle extraction.')

        form.addParam('boxSize', params.IntParam,
                      label='Particle box size (px)',
                      validators=[params.Positive],
                      help='This is size of the boxed particles (in pixels).')

        line = form.addLine('Frames range',
                            help='Specify the frames range to extract '
                                 'particles. The first frame is 1. If '
                                 'you set 0 in the last frame, it '
                                 'means that you will use until the '
                                 'last frame of the movie.')
        line.addParam('frame0', params.IntParam, label='First')
        line.addParam('frameN',  params.IntParam, label='Last')

        form.addParam('avgFrames', params.IntParam, default=1,
                      label='Average every so many frames',
                      validators=[params.Positive],
                      help='Average over this number of individual movie frames. '
                           'For Relion movie refinement it is recommended to '
                           'adjust this value to have a dose of at least '
                           'approximately 1 e/A^2 in the averaged frames, '
                           'so use a value higher than 1 if you have a '
                           'dose of less than 0.5-1 e/A^2 in each '
                           'individual movie frame.')

        form.addParam('doRescale', params.BooleanParam, default=False,
                      label='Rescale particles?',
                      help='If set to Yes, particles will be re-scaled. '
                           'Note that the re-scaled size below will be in '
                           'the down-scaled images.')

        form.addParam('rescaledSize', params.IntParam,
                      validators=[params.Positive],
                      condition='doRescale',
                      label='Re-scaled size (px)',
                      help='Final size in pixels of the extracted particles. '
                           'The provided value should be an even number. ')

        form.addSection(label='Preprocess')

        form.addParam('doInvert', params.BooleanParam, default=None,
                      label='Invert contrast?',
                      help='Invert the contrast if your particles are black '
                           'over a white background.  Xmipp, Spider, Relion '
                           'and Eman require white particles over a black '
                           'background. Frealign (up to v9.07) requires black '
                           'particles over a white background')

        form.addParam('doNormalize', params.BooleanParam, default=True,
                      label='Normalize particles?',
                      help='If set to Yes, particles will be normalized in '
                           'the way RELION prefers it.')

        form.addParam('backDiameter', params.IntParam, default=-1,
                      condition='doNormalize',
                      label='Diameter background circle (px)',
                      help='Particles will be normalized to a mean value of '
                           'zero and a standard-deviation of one for all '
                           'pixels in the background area. The background area '
                           'is defined as all pixels outside a circle with '
                           'this given diameter in pixels (before rescaling). '
                           'When specifying a negative value, a default value '
                           'of 75% of the Particle box size will be used.')

        form.addParam('stddevWhiteDust', params.FloatParam, default=-1,
                      condition='doNormalize',
                      label='Stddev for white dust removal: ',
                      help='Remove very white pixels from the extracted '
                           'particles. Pixels values higher than this many '
                           'times the image stddev will be replaced with '
                           'values from a Gaussian distribution. \n'
                           'Use negative value to switch off dust removal.')

        form.addParam('stddevBlackDust', params.FloatParam, default=-1,
                      condition='doNormalize',
                      label='Stddev for black dust removal: ',
                      help='Remove very black pixels from the extracted '
                           'particles. Pixels values higher than this many '
                           'times the image stddev will be replaced with '
                           'values from a Gaussian distribution. \n'
                           'Use negative value to switch off dust removal.')

        form.addParallelSection(threads=0, mpi=4)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._insertFunctionStep('convertInputStep',
                                 self.getInputMovies().getObjId(),
                                 self.getParticles().getObjId())
        args = self._getExtractArgs()
        self._insertFunctionStep('extractNoStreamParticlesStep', args)
        self._insertFunctionStep('createOutputStep')

        # Disable streaming functions:
        self._insertFinalSteps = self._doNothing
        self._stepsCheck = self._doNothing

    def _doNothing(self, *args):
        pass  # used to avoid some streaming functions

    # -------------------------- STEPS functions ------------------------------
    def convertInputStep(self, moviesId, partsId):
        self._ih = em.ImageHandler()  # used to convert movies
        inputMovies = self.getInputMovies()
        firstMovie = inputMovies.getFirstItem()
        self._linkOrConvertMovie = self._getConvertFunc(firstMovie)

        # convert movies to mrcs and write a movie star file
        writeSetOfMovies(inputMovies,
                         self._getPath(self._getFileName('inputMovies')),
                         preprocessImageRow=self._preprocessMovieRow)

        inputParts = self.getParticles()
        firstPart = inputParts.getFirstItem()
        self._linkOrConvertPart = self._getConvertFunc(firstPart)
        # convert particle set to star file
        # we avoid to use the defaul convertBinaryFiles because we want
        # to set a different naming convention to map between movie mrcs
        # and particles stack mrcs
        writeSetOfParticles(inputParts,
                            self._getPath(self._getFileName('inputParts')),
                            None,  # do not use convertBinaryFiles
                            fillMagnification=True,
                            alignType=inputParts.getAlignment(),
                            postprocessImageRow=self._postprocessParticleRow)

    def extractNoStreamParticlesStep(self, args):
        """ Run relion extract particles job. """
        self.runJob(self._getProgram('relion_preprocess'), args,
                    cwd=self.getWorkingDir())

    def createOutputStep(self):
        inputMovies = self.getInputMovies()
        nFrames = inputMovies.getFirstItem().getNumberOfFrames()

        inputParts = self.getParticles()
        movieParticles = self._createSetOfMovieParticles()
        movieParticles.copyInfo(inputParts)
        movieParticles.setSamplingRate(self._getNewSampling())

        self.lastMicName = None
        self.partList = []

        def _addPartsFromMic():
            # To avoid parsing the Relion star files...we are assuming here
            # the order in which Relion is generating the movie-particles per stack
            # it start part 1, 2, N of frame 1, then 1, 2..N of frame 2 and so on.
            # If this way changes in the future, the following code could break.
            # For the sake of performace, I will take the risk now.
            count = 0
            avgFrames = self.avgFrames.get()

            for frame in range(0, nFrames, avgFrames):
                frameId = min(frame + avgFrames, nFrames)

                for mPart in self.partList:
                    mPart.setObjId(None)  # clear objId to insert a new one
                    mPart.setFrameId(frameId)
                    count += 1
                    mPart.setIndex(count)
                    mPart._rlnAverageNrOfFrames = em.Integer(avgFrames)
                    movieParticles.append(mPart)

            del self.partList  # free unnecessary particle list memory
            self.partList = []

        for part in inputParts.iterItems(orderBy='_micId'):
            micName = part.getCoordinate().getMicName()

            if micName != self.lastMicName:
                if self.lastMicName is not None:
                    _addPartsFromMic()
                self.lastMicName = micName
                movieBase = '%s_movie.mrcs' % pwutils.removeBaseExt(micName)

                def _replaceSuffix(suffix):
                    return movieBase.replace('_movie.mrcs', suffix)

                # Move the resulting stack of movie-particles to extra directly
                movieStack = self._getExtraPath('output', 'extra', movieBase)
                self.newMovieStack = self._getExtraPath(_replaceSuffix('_ptcls.mrcs'))
                pwutils.moveFile(movieStack, self.newMovieStack)

                # Clean up intermediate files (either links or converted)
                # plus generated files not needed anymore
                if not pwutils.envVarOn("SCIPION_DEBUG_NOCLEAN"):
                    pwutils.cleanPath(self._getExtraPath(movieBase),
                                      self._getExtraPath(_replaceSuffix('.mrcs')))

            # Create a movie particles based on that one and
            # store in the list of this movie
            mPart = em.MovieParticle()
            mPart.copy(part)  # copy all information from part
            mPart.setParticleId(part.getObjId())

            mPart.setFileName(self.newMovieStack)
            self.partList.append(mPart)

        pwutils.cleanPath(self._getExtraPath('output'))

        _addPartsFromMic()

        self._defineOutputs(outputParticles=movieParticles)
        self._defineSourceRelation(self.inputMovies, movieParticles)
        self._defineSourceRelation(self.inputParticles, movieParticles)

            # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        
        self.validatePackageVersion('RELION_HOME', errors)
        if self.doNormalize and self.backDiameter > self.boxSize:
            errors.append("Background diameter for normalization should "
                          "be equal or less than the box size.")

        # check if movie frame shifts are all zeros
        movie = self.getInputMovies().getFirstItem()
        iniFrame, lastFrame, _ = movie.getFramesRange()
        shiftsAligned = [0] * (lastFrame-iniFrame+1)

        if movie.hasAlignment():
            shiftX, shiftY = movie.getAlignment().getShifts()
            if not (shiftX == shiftY == shiftsAligned):
                errors.append('Input movies must have zero shifts, meaning they '
                              'should be already aligned.')

        return errors
    
    def _citations(self):
        return ['Scheres2012b']
        
    def _summary(self):
        summary = []
        summary.append("Particle box size: %d" % self.boxSize)
        
        if not hasattr(self, 'outputParticles'):
            summary.append("Output images not ready yet.") 
        else:
            summary.append("Movie particles extracted: %d" %
                           self.outputParticles.getSize())
            if self.avgFrames > 1:
                summary.append("Averaged every %d frames." %
                               self.avgFrames.get())
            
        return summary
    
    def _methods(self):
        methodsMsgs = []

        if self.getStatus() == STATUS_FINISHED:
            msg = ("A total of %d movie particles of size %d were extracted"
                   % (self.getOutput().getSize(), self.boxSize))

            msg += self.methodsVar.get('')
            methodsMsgs.append(msg)

            if self.doInvert:
                methodsMsgs.append("Inverted contrast on images.")
            if self._doDownsample():
                methodsMsgs.append("Particles downsampled by a factor of %0.2f."
                                   % self._getDownFactor())

        return methodsMsgs

    # -------------------------- UTILS functions ------------------------------
    def _getExtractArgs(self):
        params = ' --i %s' % self._getFileName('inputMovies')
        params += ' --reextract_data_star %s' % self._getFileName('inputParts')
        params += ' --part_dir %s' % self._getFileName('outputDir')
        params += ' --extract --extract_movies --movie_name movie'
        params += ' --avg_movie_frames %d' % self.avgFrames.get()
        params += ' --first_movie_frame %d --last_movie_frame %d' %\
                  (self.frame0.get(), self.frameN.get())
        params += ' --extract_size %d' % self.boxSize

        if self.backDiameter <= 0:
            diameter = self.boxSize.get() * 0.75 / self._getDownFactor()
        else:
            diameter = self.backDiameter.get()

        params += ' --bg_radius %d' % int(diameter/2)

        if self.doInvert:
            params += ' --invert_contrast'

        if self.doNormalize:
            params += ' --norm'

        if self._doDownsample():
            params += ' --scale %d' % self.rescaledSize

        if self.stddevWhiteDust > 0:
            params += ' --white_dust %0.3f' % self.stddevWhiteDust

        if self.stddevBlackDust > 0:
            params += ' --black_dust %0.3f' % self.stddevBlackDust

        return params

    def _getDownFactor(self):
        if self.doRescale:
            return self.boxSize.get() / self.rescaledSize.get()
        return 1

    def _getNewSampling(self):
        newSampling = self.getInputMovies().getSamplingRate()

        if self._doDownsample():
            # Set new sampling, it should be the input sampling of the used
            # movies multiplied by the downFactor
            newSampling *= self._getDownFactor()

        return newSampling

    def _doDownsample(self):
        return self.doRescale and self.rescaledSize != self.boxSize

    def getInputMovies(self):
        return self.inputMovies.get()

    def getParticles(self):
        return self.inputParticles.get()

    def getOutput(self):
        if (self.hasAttribute('outputParticles') and
                self.outputParticles.hasValue()):
            return self.outputParticles
        else:
            return None

    def _getConvertFunc(self, img):
        """ Return a function that will link or convert the input
        file depending if the input img has a format valid for Relion.
        """
        srcFile = img.getFileName()

        if srcFile.endswith('mrcs'):
            return lambda s, d: pwutils.createLink(s, d)
        else:
            return lambda s, d: self._ih.convert(s, d)

    def _preprocessMovieRow(self, movie, movieRow):
        # 'movie' suffix in newName is required for Relion
        micBase = pwutils.removeBaseExt(movie.getMicName())
        newName = self._getExtraPath('%s_movie.mrcs' % micBase)
        self._linkOrConvertMovie(movie.getFileName(), newName)
        # The command will be launched from the working dir
        # so, let's make the movie path relative to that
        movie.setFileName(relpath(newName, self._getPath()))

    def _postprocessParticleRow(self, part, partRow):
        # Keep a map of the already converted files
        self._stackDict = getattr(self, '_stackDict', {})
        micBase = pwutils.removeBaseExt(part.getCoordinate().getMicName())

        if micBase not in self._stackDict:
            # Same convention without '_movie' suffix
            newName = self._getExtraPath('%s.mrcs' % micBase)
            self._linkOrConvertPart(part.getFileName(), newName)
            self._stackDict[micBase] = newName
        else:
            newName = self._stackDict[micBase]

        # The command will be launched from the working dir
        # so, let's make the stack path relative to that
        newPath = locationToRelion(part.getIndex(),
                                   relpath(newName, self._getPath()))
        partRow.setValue(md.RLN_IMAGE_NAME, newPath)

        if part.hasAttribute('_rlnGroupName'):
            partRow.setValue(md.RLN_MLMODEL_GROUP_NAME,
                             '%s' % part.getAttributeValue('_rlnGroupName'))
        else:
            partRow.setValue(md.RLN_MLMODEL_GROUP_NAME,
                             '%s' % part.getMicId())

        ctf = part.getCTF()

        if ctf is not None and ctf.getPhaseShift():
            partRow.setValue(md.RLN_CTF_PHASESHIFT, ctf.getPhaseShift())

    def _postprocessImageRow(self, img, imgRow):
        img.setFrameId(imgRow.getValue(md.MDL_FRAME_ID))
        img.setParticleId(imgRow.getValue(md.MDL_PARTICLE_ID))
