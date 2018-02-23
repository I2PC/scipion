# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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

import pyworkflow.utils as pwutils
from pyworkflow.protocol.constants import STATUS_FINISHED, STEPS_SERIAL
from pyworkflow.em.protocol import ProtExtractMovieParticles
from pyworkflow.em.constants import ALIGN_PROJ
from pyworkflow import VERSION_1_2
import pyworkflow.protocol.params as params
import pyworkflow.em as em
import pyworkflow.em.metadata as md
from pyworkflow.utils.properties import Message

from convert import (writeSetOfMovies, isVersion2,
                     writeSetOfParticles, readSetOfMovieParticles,
                     MOVIE_EXTRA_LABELS)
from protocol_base import ProtRelionBase



class ProtRelionExtractMovieParticles(ProtExtractMovieParticles, ProtRelionBase):
    """Protocol to extract particles from a set of coordinates"""
    _label = 'movie particles extraction'
    _lastUpdateVersion = VERSION_1_2

    @classmethod
    def isDisabled(cls):
        return not isVersion2()

    def __init__(self, **kwargs):
        ProtExtractMovieParticles.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_SERIAL

    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        myDict = {'inputMovies': 'input_movies.star',
                  'inputParts': 'input_particles.star',
                  'outputDir': pwutils.join('extra', 'output'),
                  'outputParts': self._getExtraPath('output/movie_particles.star')
                  }

        self._updateFilenamesDict(myDict)

    #--------------------------- DEFINE param functions ------------------------
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

    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        inputMovies = self.getInputMovies()
        firstStepId = self._insertFunctionStep('convertInputStep',
                                               inputMovies.getObjId())
        args = self._getExtractArgs()
        self._insertFunctionStep('extractNoStreamParticlesStep', args,
                                 prerequisites=[firstStepId])
        self._insertFunctionStep('createOutputStep')

        # Disable streaming functions:
        self._insertFinalSteps = self._doNothing
        self._stepsCheck = self._doNothing

    def _doNothing(self, *args):
        pass  # used to avoid some streaming functions

    # -------------------------- STEPS functions -------------------------------
    def convertInputStep(self, moviesId):
        self._ih = em.ImageHandler()  # used to convert movies
        self.inputMovieSet = self.getInputMovies()
        self.inputParts = self.getParticles()
        # convert movies to mrcs and write a movie star file
        writeSetOfMovies(self.inputMovieSet,
                         self._getPath(self._getFileName('inputMovies')),
                         preprocessImageRow=self._preprocessMovieRow)

        # convert particle set to star file
        writeSetOfParticles(self.inputParts,
                            self._getPath(self._getFileName('inputParts')),
                            self._getExtraPath(),
                            fillMagnification=True,
                            alignType=self.inputParts.getAlignment(),
                            postprocessImageRow=self._postprocessParticleRow)

    def extractNoStreamParticlesStep(self, args):
        """ Run relion extract particles job. """
        self.runJob(self._getProgram('relion_preprocess'), args,
                    cwd=self.getWorkingDir())

    def createOutputStep(self):
        movieParticles = self._createSetOfMovieParticles()
        movieParticles.copyInfo(self.inputParts)
        movieParticles.setSamplingRate(self._getNewSampling())

        mData = md.MetaData()
        mdAll = md.MetaData()

        for movie in self.inputMovieSet:
            movieName = movie.getMicName()
            # 'aligned_movie' suffix comes from protocol_align_movies
            movieParts = movieName.replace('.mrcs',
                                           '_aligned_movie_extract.star')
            moviePartsStar = pwutils.join(self._getExtraPath('output/extra'),
                                          movieParts)
            # Join output star files with movie particles
            if pwutils.exists(moviePartsStar):
                mData.read(moviePartsStar)
                mData.removeLabel(md.RLN_IMAGE_ID)
                mdAll.unionAll(mData)

        # move output movie particle stacks and change rlnImageName in star file
        from convert import relionToLocation, locationToRelion
        for objId in mdAll:
            index, imgPath = relionToLocation(mdAll.getValue(md.RLN_IMAGE_NAME, objId))
            baseFn = os.path.basename(imgPath)
            newPath = pwutils.join(self._getExtraPath('output'), baseFn)
            newLoc = locationToRelion(index, newPath)
            if not pwutils.exists(newPath):
                pwutils.moveFile(self._getPath(imgPath), newPath)
            mdAll.setValue(md.RLN_IMAGE_NAME, newLoc, objId)

        # delete old imgPath
        pwutils.cleanPath(self._getExtraPath('output/extra'))

        particleMd = self._getFileName('outputParts')
        mdAll.write(particleMd)
        readSetOfMovieParticles(particleMd, movieParticles,
                                removeDisabled=False,
                                alignType=ALIGN_PROJ,
                                extraLabels=MOVIE_EXTRA_LABELS,
                                postprocessImageRow=self._postprocessImageRow)

        self._defineOutputs(outputParticles=movieParticles)
        self._defineSourceRelation(self.inputMovies, movieParticles)

    #--------------------------- INFO functions --------------------------------
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
        shiftX, shiftY = movie.getAlignment().getShifts()

        if movie.hasAlignment() and not (shiftX == shiftY == shiftsAligned):
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

    #--------------------------- UTILS functions -------------------------------
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

    def _preprocessMovieRow(self, img, imgRow):
        newName = pwutils.replaceBaseExt(img.getFileName(), 'mrcs')
        newPath = self._getExtraPath(newName)

        # If the movies are in 'mrcs' format just create a link
        # if not, convert to 'mrcs'
        if img.getFileName().endswith('mrcs'):
            pwutils.createLink(img.getFileName(), newPath)
        else:
            self._ih.convert(img, newPath)
        # The command will be launched from the working dir
        # so, let's make the movie path relative to that
        img.setFileName(pwutils.join('extra', newName))

    def _postprocessParticleRow(self, part, partRow):
        # we need to remove the suffix (from protocol_align_movies) so that
        # Relion can match the rlnImageName (from particle set) and
        # the rlnMicrographMovieName (from movie set)
        imgName = partRow.getValue(md.RLN_IMAGE_NAME)
        newImgName = imgName.replace('_aligned_mic_DW.', '_aligned.')
        newImgName = newImgName.replace('_aligned_mic.', '_aligned.')
        partRow.setValue(md.RLN_IMAGE_NAME, newImgName)

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
        img.setFrameId(imgRow.getValue(md.RLN_IMAGE_FRAME_NR))
        #img.setParticleId(imgRow.getValue(md.RLN_PARTICLE_ID))  # this is set in rowToParticle
