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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
This module contains the protocol for 3d refinement with Relion.
"""
import pyworkflow.em.metadata as md
from pyworkflow.em.data import(Volume,FSC)
from pyworkflow.em.protocol import ProtRefine3D

from pyworkflow.em.packages.relion.protocol_base import ProtRelionBase


class ProtRelionRefine3D(ProtRefine3D, ProtRelionBase):
    """Protocol to refine a 3D map using Relion. Relion employs an empirical
Bayesian approach to refinement of (multiple) 3D reconstructions
or 2D class averages in electron cryo-microscopy (cryo-EM). Many
parameters of a statistical model are learned from the data,which
leads to objective and high-quality results.
    """    
    _label = '3D auto-refine'
    IS_CLASSIFY = False
    CHANGE_LABELS = [md.RLN_OPTIMISER_CHANGES_OPTIMAL_ORIENTS, 
                     md.RLN_OPTIMISER_CHANGES_OPTIMAL_OFFSETS]

    PREFIXES = ['half1_', 'half2_']
    
    def __init__(self, **args):        
        ProtRelionBase.__init__(self, **args)
        
    def _initialize(self):
        """ This function is mean to be called after the 
        working dir for the protocol have been set. (maybe after recovery from mapper)
        """
        ProtRelionBase._initialize(self)
        self.ClassFnTemplate = '%(ref)03d@%(rootDir)s/relion_it%(iter)03d_classes.mrcs'
        self.numberOfClasses.set(1) # For refinement we only need just one "class"
    
    #--------------------------- INSERT steps functions --------------------------------------------  
    def _setSamplingArgs(self, args):
        """ Set sampling related params"""
        # Sampling stuff
        args['--auto_local_healpix_order'] = self.localSearchAutoSamplingDeg.get()
        
        if not self.doContinue:
            args['--healpix_order'] = self.angularSamplingDeg.get()
            args['--offset_range'] = self.offsetSearchRangePix.get()
            args['--offset_step']  = self.offsetSearchStepPix.get() * 2
            args['--auto_refine'] = ''
            args['--split_random_halves'] = ''
            
            joinHalves = "--low_resol_join_halves"
            if not joinHalves in self.extraParams.get():
                args['--low_resol_join_halves'] = 40
        
        # Set movie refinement arguments
        if self.realignMovieFrames:
            args['--realign_movie_frames'] = self._getFileName('movie_particles')
            args['--movie_frames_running_avg'] = self.movieAvgWindow.get()
            args['--sigma_off'] = self.movieStdTrans.get()
            if not self.movieIncludeRotSearch:
                args['--skip_rotate'] = ''
                args['--skip_maximize'] = ''
            else:
                args['--sigma_ang'] = self.movieStdRot.get()
        
    #--------------------------- STEPS functions --------------------------------------------     
    def createOutputStep(self):
        
        if not self.realignMovieFrames:
            imgSet = self._getInputParticles()
            vol = Volume()
            vol.setFileName(self._getExtraPath('relion_class001.mrc'))
            vol.setSamplingRate(imgSet.getSamplingRate())
            half1 = self._getFileName("final_half1_volume", ref3d=1)
            half2 = self._getFileName("final_half2_volume", ref3d=1)
            vol.setHalfMaps([half1, half2])
            
            outImgSet = self._createSetOfParticles()
            outImgSet.copyInfo(imgSet)
            self._fillDataFromIter(outImgSet, self._lastIter())

            self._defineOutputs(outputVolume=vol)
            self._defineSourceRelation(self.inputParticles, vol)
            self._defineOutputs(outputParticles=outImgSet)
            self._defineTransformRelation(self.inputParticles, outImgSet)

            fsc = FSC(objLabel=self.getRunName())
            blockName = 'model_class_%d@'%1
            fn = blockName + self._getExtraPath("relion_model.star")
            mData = md.MetaData(fn)
            fsc.loadFromMd(mData,
                       md.RLN_RESOLUTION,
                       md.RLN_MLMODEL_FSC_HALVES_REF)
            self._defineOutputs(outputFSC=fsc)
            self._defineSourceRelation(vol,fsc)


        else:
            pass
    
    #--------------------------- INFO functions -------------------------------------------- 
    def _validateNormal(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        errors = []
        self._validateDim(self._getInputParticles(), self.referenceVolume.get(),
                          errors, 'Input particles', 'Reference volume')
        return errors
    
    def _validateContinue(self):
        """ Should be overriden in subclasses to
        return summary messages for CONTINUE EXECUTION.
        """
        errors = []
        continueRun = self.continueRun.get()
        continueRun._initialize()
        lastIter = continueRun._lastIter()
        
        if self.continueIter.get() == 'last':
            continueIter = lastIter
        else:
            continueIter = int(self.continueIter.get())
        
        if continueIter > lastIter:
            errors += ["The iteration from you want to continue must be %01d or less" % lastIter]
        
        return errors
    
    def _summaryNormal(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        summary = []
        it = self._lastIter()
        if it >= 1:
            row = md.getFirstRow('model_general@' + self._getFileName('half1_model', iter=it))
            resol = row.getValue("rlnCurrentResolution")
            summary.append("Current resolution: *%0.2f*" % resol)
        return summary
    
    def _summaryContinue(self):
        """ Should be overriden in subclasses to
        return summary messages for CONTINUE EXECUTION.
        """
        summary = []
        summary.append("Continue from iteration %01d" % self._getContinueIter())
        return summary

    def _summary(self):
        if not hasattr(self, 'outputVolume'):
            return ["Output volume not ready yet."]
        else:
            return ProtRelionBase._summary(self)

    #--------------------------- UTILS functions --------------------------------------------
    def _fillDataFromIter(self, imgSet, iteration):
        outImgsFn = self._getFileName('data', iter=iteration)
        imgSet.setAlignmentProj()
        imgSet.copyItems(self._getInputParticles(),
                         updateItemCallback=self._createItemMatrix,
                         itemDataIterator=md.iterRows(outImgsFn, sortByLabel=md.RLN_IMAGE_ID))
    
    def _createItemMatrix(self, item, row):
        from pyworkflow.em.packages.relion.convert import createItemMatrix
        from pyworkflow.em import ALIGN_PROJ
        
        createItemMatrix(item, row, align=ALIGN_PROJ)
        