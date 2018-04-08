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
                            None,  # do not convertBinaryFiles here, do it later
                            fillMagnification=True,
                            alignType=self.inputParts.getAlignment(),
                            postprocessImageRow=self._postprocessParticleRow)

        # We have to change rlnImageName in inputParts star file to match
        # rlnMicrographMovieName from inputMovies star file.

        # Create a dict with old and new image names
        imgNameDict = {}

        for movie in self.movieDict:
            for img in self.ptclDict:
                micName = self.ptclDict[img]
                if micName in movie:
                    newImgName = pwutils.removeBaseExt(self.movieDict[movie])
                    newImgName = newImgName.replace('_movie', '.mrcs')
                    if img not in imgNameDict:
                        imgNameDict[img] = newImgName

        mData = md.MetaData(self._getPath(self._getFileName('inputParts')))
        imgNames = mData.getColumnValues(md.RLN_IMAGE_NAME)

        for index, img in enumerate(imgNames):
            loc, name = img.split('@')
            imgNames[index] = "%s@%s" % (loc, self._getExtraPath(imgNameDict[name]))

        mData.setColumnValues(md.RLN_IMAGE_NAME, imgNames)
        mData.write(self._getPath(self._getFileName('inputParts')))

        # link new particle stacks
        for oldName in imgNameDict:
            newName = self._getExtraPath(imgNameDict[oldName])
            if oldName.endswith('mrcs'):
                pwutils.createLink(oldName, newName)
            else:
                img = tuple(oldName.split('@'))
                self._ih.convert(img, newName)


    def extractNoStreamParticlesStep(self, args):
        """ Run relion extract particles job. """
        self.runJob(self._getProgram('relion_preprocess'), args,
                    cwd=self.getWorkingDir())

    def createOutputStep(self):
        movieParticles = self._createSetOfMovieParticles()
        movieParticles.copyInfo(self.inputParts)
        movieParticles.setSamplingRate(self._getNewSampling())

        print "movieDict = ", self.movieDict, "\n\n"
        print "ptclsDict = ", self.ptclDict, "\n\n"

        mData = md.MetaData()
        mdAll = md.MetaData()

        # Move output movie particle stacks and join output star files
        for movieKey in self.movieDict:
            movie = self.movieDict[movieKey]
            outStar = os.path.basename(movie).replace('.mrcs', '_extract.star')
            outStar = self._getExtraPath('output/extra/%s' % outStar)

            if pwutils.exists(outStar):
                mData.read(outStar)
                mData.removeLabel(md.RLN_IMAGE_ID)
                mdAll.unionAll(mData)
                mdAll.addItemId()

            outStack = self._getExtraPath('output/extra/%s' % os.path.basename(movie))
            newOutStack = self._getExtraPath('output/' +
                                             pwutils.removeBaseExt(movieKey) +
                                             '_ptcls.mrcs')
            pwutils.moveFile(outStack, newOutStack)

        # Modify rlnMicrographName, rlnOriginalParticleName back to original values
        # Change rlnImageName to point to newOutStack

        # reverse movieDict and ptclsDict
        revMovieDict = dict((v, k) for k, v in self.movieDict.iteritems())
        revPtclDict = dict((v, k) for k, v in self.ptclDict.iteritems())

        for row in md.iterRows(mdAll):
            micInd = row.getValue(md.RLN_MICROGRAPH_NAME).split('@')[0]
            imgStr = row.getValue(md.RLN_IMAGE_NAME).split('@')
            ptclInd = row.getValue(md.RLN_PARTICLE_ORI_NAME).split('@')[0]
            name = os.path.basename(imgStr[1])
            fnPath = self._getExtraPath(name)
            origMicName = revMovieDict[fnPath]
            row.setValue(md.RLN_MICROGRAPH_NAME, "%s@%s" % (micInd, origMicName))  # restore orig movie name
            row.setValue(md.RLN_PARTICLE_ORI_NAME, "%s@%s" % (ptclInd, revPtclDict[pwutils.removeBaseExt(origMicName)]))  # restore origPtclname
            row.setValue(md.RLN_IMAGE_NAME, "%s@extra/output/%s_ptcls.mrcs" % (imgStr[0], pwutils.removeBaseExt(origMicName)))
            row.setValue(md.MDL_FRAME_ID, long(micInd))
            row.setValue(md.MDL_PARTICLE_ID, long(ptclInd))


        mdAll.write(self._getExtraPath('test.star'))
        raise Exception('Debugging in progress!')

        # delete old imgPath
        pwutils.cleanPath(self._getExtraPath('output/extra'))

        #particleMd = self._getFileName('outputParts')
        #
        #mdAll.write(particleMd)
        #readSetOfMovieParticles(particleMd, movieParticles,
        #                        removeDisabled=False,
        #                        alignType=ALIGN_PROJ,
        #                        extraLabels=MOVIE_EXTRA_LABELS,
        #                        postprocessImageRow=self._postprocessImageRow)
        #
        #self._defineOutputs(outputParticles=movieParticles)
        #self._defineSourceRelation(self.inputMovies, movieParticles)

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
        oldName = img.getFileName()
        movId = img.getObjId()
        # 'movie' suffix in newName is required for Relion
        newName = self._getExtraPath('file_%06d_movie.mrcs' % movId)
        if oldName not in self.movieDict:
            self.movieDict[oldName] = newName
        # If the movies are in 'mrcs' format just create a link
        # if not, convert to 'mrcs'
        if oldName.endswith('mrcs'):
            pwutils.createLink(oldName, newName)
        else:
            self._ih.convert(img, newName)
        # The command will be launched from the working dir
        # so, let's make the movie path relative to that
        img.setFileName(pwutils.join('extra', 'file_%06d_movie.mrcs' % movId))

    def _postprocessParticleRow(self, part, partRow):
        oldImgName = partRow.getValue(md.RLN_IMAGE_NAME)
        oldImgName = str(oldImgName).split('@')[1]
        micName = partRow.getValue(md.RLN_MICROGRAPH_NAME)
        micBase = pwutils.removeBaseExt(micName)

        # fill in ptclDict[rlnImageName] = rlnMicrographName
        if oldImgName not in self.ptclDict:
            self.ptclDict[oldImgName] = micBase

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
