# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (jgomez@cnb.csic.es)
# *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
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
from os.path import exists
from glob import glob
from pyworkflow.protocol.params import (PointerParam, FloatParam, StringParam,
                                        IntParam, BooleanParam, LEVEL_ADVANCED)
from pyworkflow.em.data import Volume
from pyworkflow.em import ALIGN_PROJ
import pyworkflow.em.metadata as md
import pyworkflow.utils as pwutils
from pyworkflow.em.protocol import ProtProcessParticles

from protocol_base import ProtRelionBase
from convert import (createItemMatrix, isVersion2, getVersion, V1_3,
                    readSetOfParticles, writeSetOfParticles, MOVIE_EXTRA_LABELS)


class ProtRelionPolish(ProtProcessParticles, ProtRelionBase):
    """
    This Relion protocol tracks particle movement in movie frames
    (from previous movie refinement run), applies a resolution
    and dose-dependent weighting scheme for each frame
    and finally sums them together, producing so-called
    shiny, or polished particles.
    """
    _label = 'particle polishing'
    
    PREFIXES = ['half1_', 'half2_']

    def _initialize(self):
        """ This function is meant to be called after the
        working dir for the protocol have been set.
        (maybe after recovery from mapper)
        """
        self._createFilenameTemplates()
        self._createIterTemplates()

    
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('inputMovieParticles', PointerParam,
                      pointerClass='SetOfMovieParticles',
                      label='Input refined movie particles',
                      help='Select the refined movie particles '
                           'from Relion 3D auto-refine run.')
        form.addParam('maskRadius', IntParam, default=100,
                      label='Particle mask RADIUS (px)',
                      help='Radius of the circular mask that will be used '
                           'to define the background area.')
        form.addParam('extraParams', StringParam, default='',
                      expertLevel=LEVEL_ADVANCED,
                      label='Additional parameters',
                      help='')

        form.addSection('Movement')
        form.addParam('linerFitParticleMovements', BooleanParam, default=True,
                      label='Linear fit particle movements?',
                      help='If set to Yes, then the program will fit linear '
                           'tracks (in X-Y and in time) through the estimated '
                           'movie tracks in the input STAR file. For small '
                           'particles (e.g. < 1MDa) this will lead to more '
                           'robust beam-induced movement modelling. Because '
                           'particles that are close to each other on a '
                           'micrograph often move in similar directions, the '
                           'estimated tracks from neighbouring particles may '
                           'also be included in the fitting of each particle. '
                           'Again, in particular for smaller particles this '
                           'may improve the robustness of the fits.')
        if not isVersion2():
            form.addParam('movieAvgWindow', FloatParam, default=5,
                          label='Running average window',
                          help='The individual movie frames will be '
                               'averaged using a running average window '
                               'with the specified width. Use an odd '
                               'number. The optimal value will depend on '
                               'the SNR in the individual movie frames. '
                               'For ribosomes, we used a value of 5, '
                               'where each movie frame integrated '
                               'approximately 1 electron per squared '
                               'Angstrom.')
        form.addParam('stddevParticleDistance', IntParam, default=100,
                      condition='linerFitParticleMovements',
                      label='Stddev on particle distance (px)',
                      help='This value determines how much neighbouring '
                           'particles contribute to the fit of the movements '
                           'of each particle. This value is the standard '
                           'deviation of a Gaussian on the inter-particle '
                           'distance. Larger values mean that particles that '
                           'are further away still contribute more. Particles '
                           'beyond 3 standard deviations are excluded from '
                           'the fit. Very large values will lead to all '
                           'fitted tracks pointing in the same direction. '
                           'A value of zero means that each particle is '
                           'fitted independently.')
    
        form.addSection('Damage')
        form.addParam('performBfactorWeighting', BooleanParam, default=True,
                      label='Perform B-factor weighting?',
                      help='If set to Yes, then the program will estimate '
                           'a resolution-dependent weight for each movie '
                           'frames by calculating independent '
                           'half-reconstructions for each movie frame '
                           'separately. Gold-standard FSCs between these '
                           'are then converted into relative Guinier plots, '
                           'through which straight lines are fitted. Linear '
                           'fits are often suitable beyond 20A resolution. '
                           'Small particles may not yield such high '
                           'resolutions for the individual-frame '
                           'reconstructions. Therefore, in some cases '
                           'it may be better to skip this step. It is '
                           'however recommended to always try and perform '
                           'B-factor weighting, and to inspect the output '
                           'bfactors.star and guinier.star files, as '
                           'adequate weighting may significantly improve '
                           'resolution in the final map.')
        form.addParam('highresLimitPerFrameMaps', FloatParam, default=5,
                      condition='performBfactorWeighting',
                      label='Highres-limit per-frame maps (A)',
                      help='To estimate the resolution and dose dependency of '
                           'the radiation damage, the program will calculate '
                           'reconstructions from all first movie frames, '
                           'second movie frames, etc. These per-frame '
                           'reconstructions will have lower resolution than '
                           'the reconstruction from all-frames. To speed up '
                           'the calculations (and reduce memory '
                           'requirements), the per-frame reconstructions may '
                           'be limited in resolution using this parameter. '
                           'One can inspect the output STAR files of the '
                           'per-frame reconstructions to check afterwards '
                           'that this value was not chosen lower than the '
                           'actual resolution of these reconstruction')
        form.addParam('lowresLimitBfactorEstimation', FloatParam, default=20,
                      condition='performBfactorWeighting',
                      label='Lowres-limit B-factor estimation (A)',
                      help='This value describes the lowest resolution that '
                           'is included in the B-factor estimation of the '
                           'per-frame reconstructions. Because the power '
                           'spectrum of per-frame reconstructions is compared '
                           'to the power spectrum of the reconstruction from '
                           'all frames, a much lower value than the 10A '
                           'described in the Rosenthal and Henderson (2003) '
                           'paper in JMB can be used. Probably a value around '
                           '20A is still OK.')
        if isVersion2():
            form.addParam('avgFramesBfac', IntParam, default=1,
                          condition='performBfactorWeighting',
                          label='Average frames B-factor estimation',
                          help='B-factors for each movie frame will be estimated '
                               'from reconstructions of all particles for that '
                               'movie frame. Single-frame reconstructions '
                               'sometimes give not enough signal to estimate '
                               'reliable B-factors. This option allows one to '
                               'calculate the B-factors from running averages of '
                               'movie frames. The value specified should be an '
                               'odd number. Calculating B-factors from multiple '
                               'movie frames improves the SNR in the reconstructions, '
                               'but limits the estimation of sudden changes in '
                               'B-factors throughout the movie, for example in '
                               'the first few frames when beam-induced movement '
                               'is very rapid. Therefore, one should not use '
                               'higher values than strictly necessary.')
           
        form.addParam('maskForReconstructions', PointerParam,
                      pointerClass='VolumeMask',
                      label='Mask for the reconstructions', allowsNull=True,
                      help='A continuous mask with values between 0 (solvent) '
                           'and 1 (protein). You may provide the same mask '
                           'that was used in the post-processing of the '
                           'corresponding 3D auto-refine jobs before the movie '
                           'processing.')
        self.addSymmetry(form)

        form.addParallelSection(threads=0, mpi=3)
    
    #--------------------------- INSERT steps functions --------------------------------------------  

    def _insertAllSteps(self): 
        self._initialize()
        self._insertFunctionStep('convertInputStep',
                                 self._getInputParticles().getObjId())
        self._insertPolishStep()
        self._insertFunctionStep('organizeDataStep')
        self._insertFunctionStep('createOutputStep')
        
    def _insertPolishStep(self):
        imgSet = self._getInputParticles()
        imgStar = self._getFileName('movie_particles')
        
        params = ' --i %s' % imgStar
        if isVersion2():
            params += ' --o %s/shiny' % self._getExtraPath()
        else:
            params += ' --o shiny'
        params += ' --angpix %0.3f' % imgSet.getSamplingRate()
        params += ' --sym %s' % self.symmetryGroup.get()
        params += ' --sigma_nb %d' % self.stddevParticleDistance.get()
        params += ' --perframe_highres %0.3f' % self.highresLimitPerFrameMaps.get()
        params += ' --autob_lowres %0.3f' % self.lowresLimitBfactorEstimation.get()

        if not isVersion2():
            params += ' --movie_frames_running_avg %d' % self.movieAvgWindow.get()
            params += ' --dont_read_old_files'
        else:
            params += ' --bfactor_running_avg %d' % self.avgFramesBfac.get()

        x, _, _ = imgSet.getDimensions()
        if self.maskRadius.get() >= x/2 or self.maskRadius.get() < 0:
            params += ' --bg_radius %d' % (x/2)
        else:
            params += ' --bg_radius %d' % self.maskRadius.get()
        
        if self.performBfactorWeighting:
            params += ' --mask %s' % self.maskForReconstructions.get().getFileName()
        else:
            params += ' --skip_bfactor_weighting'
        params += ' ' + self.extraParams.get()

        self._insertFunctionStep('polishStep', params)
        
    #--------------------------- STEPS functions --------------------------------------------
    def convertInputStep(self, particlesId):
        """ Create the input file in STAR format as expected by Relion.
        If the input particles comes from Relion, just link the file.
        Params:
            particlesId: use this parameter just to force redo of convert if
                the input particles are changed.
        """
        imgSet = self._getInputParticles()
        imgStar = self._getFileName('movie_particles')
        imgStarTmp = self._getTmpPath('movie_particles.star')
        self.info("Converting set from '%s' into '%s'" %
                  (imgSet.getFileName(), imgStarTmp))

        writeSetOfParticles(imgSet, imgStarTmp, self._getExtraPath(),
                            alignType=imgSet.getAlignment(),
                            extraLabels=MOVIE_EXTRA_LABELS)
        mdImg = md.MetaData(imgStarTmp)

        # replace mdColumn from *.stk to *.mrcs as Relion2 requires
        if getVersion() == V1_3:
            mdColumn = md.RLN_PARTICLE_NAME
        else:
            mdColumn = md.RLN_PARTICLE_ORI_NAME

        from convert import relionToLocation, locationToRelion
        for objId in mdImg:
            index, imgPath = relionToLocation(mdImg.getValue(mdColumn, objId))
            if not imgPath.endswith('mrcs'):
                newName = pwutils.replaceBaseExt(os.path.basename(imgPath), 'mrcs')
                newPath = self._getTmpPath(newName)
                newLoc = locationToRelion(index, newPath)
                if not exists(newPath):
                    pwutils.createLink(imgPath, newPath)
                mdImg.setValue(mdColumn, newLoc, objId)
        mdImg.write(imgStar)

    def polishStep(self, params):
        self.runJob(self._getProgram('relion_particle_polish'), params)
    
    def organizeDataStep(self):
        from convert import relionToLocation, locationToRelion
        if getVersion() == V1_3:
            mdColumn = md.RLN_PARTICLE_NAME
        else:
            mdColumn = md.RLN_PARTICLE_ORI_NAME

        shinyStar = self._getFileName('shiny')
        newDir = self._getExtraPath('polished_particles')
        pwutils.makePath(newDir)

        if not isVersion2():
            pwutils.makePath(self._getExtraPath('shiny'))
            shinyOld = "shiny.star"
            inputFit = "movie_particles_shiny.star"
            try:
                pwutils.moveFile(shinyOld, shinyStar)
                pwutils.moveFile(self._getPath(inputFit), self._getExtraPath("shiny/all_movies_input_fit.star"))
                for half in self.PREFIXES:
                    pwutils.moveFile(self._getPath('movie_particles_shiny_%sclass001_unfil.mrc' % half),
                                     self._getExtraPath('shiny/shiny_%sclass001_unfil.mrc' % half))
                self._renameFiles('movie_particles_shiny_post*', 'movie_particles_')
                self._renameFiles('movie_particles_shiny*', 'movie_particles_shiny_')
            except:
                raise Exception('ERROR: some file(s) were not found!')

        # move polished particles from Tmp to Extra path
        # and restore previous mdColumn
        mdShiny = md.MetaData(shinyStar)
        oldPath = ""

        for objId in mdShiny:
            index, imgPath = relionToLocation(mdShiny.getValue(md.RLN_IMAGE_NAME, objId))
            newPath = pwutils.join(newDir, str(imgPath).split('/')[-1])
            newLoc = locationToRelion(index, newPath)
            mdShiny.setValue(md.RLN_IMAGE_NAME, newLoc, objId)
            if oldPath != imgPath and exists(imgPath):
                pwutils.moveFile(imgPath, newPath)
                oldPath = imgPath

            index2, imgPath2 = relionToLocation(mdShiny.getValue(mdColumn, objId))
            absPath = os.path.realpath(imgPath2)
            newPath2 = 'Runs' + str(absPath).split('Runs')[1]
            newLoc2 = locationToRelion(index2, newPath2)
            mdShiny.setValue(mdColumn, newLoc2, objId)

        mdShiny.write(shinyStar, md.MD_OVERWRITE)
        pwutils.cleanPath(self._getExtraPath('shiny/Runs'))

    def createOutputStep(self):
        imgSet = self._getInputParticles()
        vol = Volume()
        vol.setFileName(self._getFileName('volume_shiny'))
        vol.setSamplingRate(imgSet.getSamplingRate())
        shinyPartSet = self._createSetOfParticles()
        shinyPartSet.copyInfo(imgSet)
        shinyPartSet.setAlignmentProj()
        readSetOfParticles(self._getFileName('shiny'), shinyPartSet,
                           alignType=ALIGN_PROJ)

        self._defineOutputs(outputParticles=shinyPartSet)
        self._defineOutputs(outputVolume=vol)
        
        self._defineSourceRelation(imgSet, shinyPartSet)
        self._defineSourceRelation(imgSet, vol)

    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        """ Should be overwritten in subclasses to
        return summary message for NORMAL EXECUTION. 
        """
        errors = []
        self.validatePackageVersion('RELION_HOME', errors)

        if self.performBfactorWeighting:
            if self.maskForReconstructions.get() is None:
                errors.append('You should provide a *mask* when performing B-factor weighting.')

        return errors
    
    def _summary(self):
        """ Should be overwritten in subclasses to
        return summary message for NORMAL EXECUTION. 
        """
        summary = []
        if not hasattr(self, 'outputVolume'):
            summary.append("Output is not ready yet.")
        else:
            if self.performBfactorWeighting:
                summary.append('Performed B-factor weighting')
            row = md.getFirstRow('general@' + self._getExtraPath('shiny/shiny_post.star'))
            res = row.getValue('rlnFinalResolution')
            summary.append('Resolution of reconstructions from shiny particles: *%0.2f A*' % res)

        return summary
    
    #--------------------------- UTILS functions --------------------------------------------
    def _getInputParticles(self):
        return self.inputMovieParticles.get()

    def _lastFrame(self):
        mdFn = md.MetaData(self._getFileName('shiny'))
        nrOfFrames = mdFn.getValue(md.RLN_PARTICLE_NR_FRAMES, mdFn.firstObject())

        return nrOfFrames

    def _createItemMatrix(self, item, row):
        from pyworkflow.em import ALIGN_PROJ
        createItemMatrix(item, row, align=ALIGN_PROJ)

    def _renameFiles(self, pattern1, pattern2):
        # find files by pattern1, move and rename them by replacing pattern2
        filesList = sorted(glob(self._getPath(pattern1)))
        for fn in filesList:
            oldFn = os.path.basename(fn)
            newFn = pwutils.join(self._getExtraPath('shiny'), oldFn.replace(pattern2, ''))
            pwutils.moveFile(fn, newFn)
