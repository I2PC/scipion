# **************************************************************************
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

import os
from os.path import exists
from pyworkflow.protocol.params import (PointerParam, FloatParam, StringParam,
                                        IntParam, BooleanParam, LEVEL_ADVANCED)
from pyworkflow.em.data import Volume 
from pyworkflow.em.protocol import ProtProcessParticles

from protocol_base import ProtRelionBase


class ProtRelionPolish(ProtProcessParticles, ProtRelionBase):
    """    
    Reconstruct a volume using Relion from a given set of particles.
    The alignment parameters will be converted to a Relion star file
    and used as direction projections to reconstruct.
    """
    _label = 'polish particles'
    
    def _initialize(self):
        """ This function is mean to be called after the 
        working dir for the protocol have been set. (maybe after recovery from mapper)
        """
        self._createFilenameTemplates()
        self._createIterTemplates()
        
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('refineRun', PointerParam, pointerClass='ProtRelionRefine3D',              
              label='Movie particles refine run',
              help='Select the Relion 3D auto-refine run where \n'
                   'the movie particles where refined.')
        form.addParam('maskRadius', IntParam, default=100,
                      label='Particle mask RADIUS (px)',
                      help='Radius of the circular mask that will be used '
                           'to define the background area.') 

        group = form.addGroup('Movement')
        group.addParam('linerFitParticleMovements', BooleanParam, default=True,
                      label='Linear fit particle movements?',
                      help='If set to Yes, then the program will fit linear tracks (in X-Y and in time) through the estimated movie tracks in the input STAR file. For small particles (e.g. < 1MDa) this will lead to more robust beam-induced movement modelling. Because particles that are close to each other on a micrograph often move in similar directions, the estimated tracks from neighbouring particles may also be included in the fitting of each particle. Again, in particular for smaller particles this may improve the robustness of the fits.')
        #NOTE: we do not need to ask for running average windows as in relion gui
        # we can take it from the previous run
        group.addParam('stddevParticleDistance', IntParam, default=100,
                       condition='linerFitParticleMovements',
                       label='Stddev on particle distance (px)',
                       help='This value determines how much neighbouring particles contribute to the fit of the movements of each particle. This value is the standard deviation of a Gaussian on the inter-particle distance. Larger values mean that particles that are further away still contribute more. Particles beyond 3 standard deviations are excluded from the fit. Very large values will lead to all fitted tracks pointing in the same direction. A value of zero means that each particle is fitted independently.')
        
    
        group = form.addGroup('Damage')
        group.addParam('performBfactorWeighting', BooleanParam, default=True,
                      label='Perform B-factor weighting?',
                      help='If set to Yes, then running averages of the individual frames of recorded movies will also be aligned rotationally. \n'
                           'If one wants to perform particle polishing, then rotational alignments of the movie frames is NOT necessary and will only take more computing time.')
        #NOTE: we do not need to ask for running average windows as in relion gui
        # we can take it from the previous run
        group.addParam('highresLimitPerFrameMaps', FloatParam, default=5,
                       condition='performBfactorWeighting',
                       label='Highres-limit per-frame maps (A)',
                       help='To estimate the resolution and dose dependency of the radiation damage, the program will calculate reconstructions from all first movie frames, second movie frames, etc. These per-frame reconstructions will have lower resolution than the reconstruction from all-frames. To speed up the calculations (and reduce memory requirements), the per-frame reconstructions may be limited in resolution using this parameter. One can inspect the output STAR files of the per-frame reconstructions to check afterwards that this value was not chosen lower than the actual resolution of these reconstruction')
        group.addParam('lowresLimitBfactorEstimation', FloatParam, default=20,
                       condition='performBfactorWeighting',
                       label='Lowres-limit B-factor estimation (A)',
                       help='This value describes the lowest resolution that is included in the B-factor estimation of the per-frame reconstructions. Because the power spectrum of per-frame reconstructions is compared to the power spectrum of the reconstruction from all frames, a much lower value than the 10A described in the Rosenthal and Henderson (2003) paper in JMB can be used. Probably a value around 20A is still OK.')
           
        group.addParam('maskForReconstructions', PointerParam, pointerClass='VolumeMask',
                      label='Mask for the reconstructions', allowsNull=True,
                      help='A continuous mask with values between 0 (solvent) and 1 (protein). You may provide the same map that was obtained in the post-processing of the corresponding auto-refine jobs before the movie processing.') 

        self.addSymmetry(group)
        
        form.addParam('extraParams', StringParam, default='',
                      expertLevel=LEVEL_ADVANCED,
                      label='Additional parameters',
                      help='')
        
        form.addParallelSection(threads=1, mpi=1) 
        #TODO: Add an option to allow the user to decide if copy binary files or not        
            
    #--------------------------- INSERT steps functions --------------------------------------------  

    def _insertAllSteps(self): 
        self._initialize()
        self._insertPolishStep()
        self._insertFunctionStep('organizeDataStep')
        self._insertFunctionStep('createOutputStep')
        
    def _insertPolishStep(self):
        refineRun = self.refineRun.get()
        imgSet = refineRun._getInputParticles()
        refineRun._initialize() # load filenames stuff
        
        imgStar = self._getFileName('data', iter=refineRun._lastIter())
        
        # `which relion_particle_polish_mpi` --i Refine3D/run2_ct21_it021_data.star 
        # --o folder/shiny  --angpix 3.54 --movie_frames_running_avg 5 --dont_read_old_files  
        # --sigma_nb 100 --perframe_highres 6 --autob_lowres 20 --sym C1 --bg_radius 28 
        #--white_dust -1 --black_dust -1
        
        params = ' --i %s' % imgStar
        params += ' --o shiny'
        params += ' --angpix %0.3f' % imgSet.getSamplingRate()
        params += ' --movie_frames_running_avg %d' % refineRun.movieAvgWindow.get()
        params += ' --dont_read_old_files'
        params += ' --sym %s' % self.symmetryGroup.get()
        params += ' --sigma_nb %d' % self.stddevParticleDistance.get()
        params += ' --perframe_highres %0.3f' % self.highresLimitPerFrameMaps.get()
        params += ' --autob_lowres %0.3f' % self.lowresLimitBfactorEstimation.get()
        
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
    def polishStep(self, params):
        """ Create the input file in STAR format as expected by Relion.
        If the input particles comes from Relion, just link the file. 
        """
        from pyworkflow.utils.path import copyTree
        copyTree(self.refineRun.get()._getExtraPath(), self._getExtraPath())
        self.runJob(self._getProgram('relion_particle_polish'), params)
    
    def organizeDataStep(self):
        from pyworkflow.utils import moveFile
        import pyworkflow.em.metadata as md
        from convert import relionToLocation, locationToRelion
        
        # moving shiny.star form project base path to the current protocol extra path.
        shinyStar = "shiny.star"
        pathFixedShiny = self._getExtraPath(shinyStar)
        
        if exists(shinyStar):
            moveFile(shinyStar, pathFixedShiny)
        mdShiny = md.MetaData(pathFixedShiny)
        
        oldImgPath = ""
        for objId in mdShiny:
            index, imgPath = relionToLocation(mdShiny.getValue(md.RLN_IMAGE_NAME, objId))
            newPath = self._getExtraPath(os.path.basename(imgPath))
            newLoc = locationToRelion(index, newPath)
            print "newLoc ", newLoc
            mdShiny.setValue(md.RLN_IMAGE_NAME, newLoc, objId)
            if oldImgPath != imgPath and exists(imgPath):
                moveFile(imgPath, newPath)
                oldImgPath = imgPath
        mdShiny.write(pathFixedShiny)
    
    def createOutputStep(self):
        from pyworkflow.em.packages.relion.convert import readSetOfParticles
        imgSet = self.refineRun.get()._getInputParticles()
        vol = Volume()
        vol.setFileName(self._getFileName('volume_shiny', iter=self._lastIter()))
        vol.setSamplingRate(imgSet.getSamplingRate())
        
        shinyPartSet = self._createSetOfParticles()
        
        shinyPartSet.copyInfo(imgSet)
        shinyPartSet.setAlignmentProj()
        readSetOfParticles(self._getExtraPath("shiny.star"), shinyPartSet)
        
        self._defineOutputs(outputParticles=shinyPartSet)
        self._defineSourceRelation(imgSet, shinyPartSet)
        self._defineOutputs(outputVolume=vol)
        self._defineSourceRelation(imgSet, vol)

    
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        errors = []
        if self.performBfactorWeighting:
            if self.maskForReconstructions.get() is None:
                errors.append('You should select a *mask* if performing b-factor weighting.')
        return errors
    
    def _summary(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        return []
    
    #--------------------------- UTILS functions --------------------------------------------
