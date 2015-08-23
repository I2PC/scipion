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
"""
This module contains the protocol for 3d classification with relion.
"""

import pyworkflow.em.metadata as md
import pyworkflow.em as em
from pyworkflow.em.protocol import ProtClassify3D
from pyworkflow.em.data import Volume

from pyworkflow.em.packages.relion.protocol_base import ProtRelionBase
from pyworkflow.em.packages.relion.convert import relionToLocation, rowToAlignment


class ProtRelionClassify3D(ProtClassify3D, ProtRelionBase):
    """    
    Protocol to classify 3D using Relion. Relion employs an empirical
    Bayesian approach to refinement of (multiple) 3D reconstructions
    or 2D class averages in electron cryo-microscopy (cryo-EM). Many
    parameters of a statistical model are learned from the data,which
    leads to objective and high-quality results.
    """
    _label = '3D classification'
    CHANGE_LABELS = [md.RLN_OPTIMISER_CHANGES_OPTIMAL_ORIENTS, 
                     md.RLN_OPTIMISER_CHANGES_OPTIMAL_OFFSETS]
    
    def __init__(self, **args):        
        ProtRelionBase.__init__(self, **args)
        
    def _initialize(self):
        """ This function is mean to be called after the 
        working dir for the protocol have been set. (maybe after recovery from mapper)
        """
        ProtRelionBase._initialize(self)
    
    #--------------------------- INSERT steps functions --------------------------------------------  
    def _setSamplingArgs(self, args):
        """ Set sampling related params. """
        if self.doImageAlignment:
            args['--healpix_order'] = self.angularSamplingDeg.get()
            args['--offset_range'] = self.offsetSearchRangePix.get()
            args['--offset_step']  = self.offsetSearchStepPix.get() * 2
            if self.localAngularSearch:
                args['--sigma_ang'] = self.localAngularSearchRange.get() / 3.
        else:
            args['--skip_align'] = ''
    
    #--------------------------- STEPS functions --------------------------------------------
    def createOutputStep(self):
        partSet = self.inputParticles.get()
        classes3D = self._createSetOfClasses3D(partSet)
        self._fillClassesFromIter(classes3D, self._lastIter())
        
        self._defineOutputs(outputClasses=classes3D)
        self._defineSourceRelation(self.inputParticles, classes3D)


        # create a SetOfVolumes and define its relations
        volumes = self._createSetOfVolumes()
        volumes.setSamplingRate(partSet.getSamplingRate())
        
        for class3D in classes3D:
            vol = class3D.getRepresentative()
            vol.setObjId(class3D.getObjId())
            volumes.append(vol)
        
        self._defineOutputs(outputVolumes=volumes)
        self._defineSourceRelation(self.inputParticles, volumes)
        
        if not self.doContinue:
            self._defineSourceRelation(self.referenceVolume, classes3D)
            self._defineSourceRelation(self.referenceVolume, volumes)
    
    #--------------------------- INFO functions -------------------------------------------- 
    def _validateNormal(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        errors = []
        partSizeX, _, _ = self._getInputParticles().getDim()
        volSizeX, _, _ = self.referenceVolume.get().getDim()
        if partSizeX != volSizeX:
            errors.append('Volume and particles dimensions must be equal!!!')

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
            row = md.getFirstRow('model_general@' + self._getFileName('model', iter=it))
            resol = row.getValue("rlnCurrentResolution")
            summary.append("Current resolution: *%0.2f*" % resol)
        
        summary.append("Input Particles: *%d*\nClassified into *%d* 3D classes\n" % (self.inputParticles.get().getSize(),
                                                                              self.numberOfClasses.get()))
        
        return summary
    
    def _summaryContinue(self):
        """ Should be overriden in subclasses to
        return summary messages for CONTINUE EXECUTION.
        """
        summary = []
        summary.append("Continue from iteration %01d" % self._getContinueIter())
        return summary
    
    def _methods(self):
        strline=''
        if hasattr(self, 'outputClasses'):
            strline += 'We classified %d particles into %d 3D classes using Relion Classify3d. '%\
                           (self.inputParticles.get().getSize(), self.numberOfClasses.get())
        return [strline]
    
    #--------------------------- UTILS functions --------------------------------------------
    def _loadClassesInfo(self, iteration):
        """ Read some information about the produced Relion 3D classes
        from the *model.star file.
        """
        self._classesInfo = {} # store classes info, indexed by class id
         
        modelStar = md.MetaData('model_classes@' + self._getFileName('model', iter=iteration))
        
        for classNumber, row in enumerate(md.iterRows(modelStar)):
            index, fn = relionToLocation(row.getValue('rlnReferenceImage'))
            # Store info indexed by id, we need to store the row.clone() since
            # the same reference is used for iteration            
            self._classesInfo[classNumber+1] = (index, fn, row.clone())
    
    def _fillClassesFromIter(self, clsSet, iteration):
        """ Create the SetOfClasses3D from a given iteration. """
        self._loadClassesInfo(iteration)
        dataStar = self._getFileName('data', iter=iteration)
        clsSet.classifyItems(updateItemCallback=self._updateParticle,
                             updateClassCallback=self._updateClass,
                             itemDataIterator=md.iterRows(dataStar))
    
    def _updateParticle(self, item, row):
        item.setClassId(row.getValue(md.RLN_PARTICLE_CLASS))
        item.setTransform(rowToAlignment(row, em.ALIGN_3D))
        
        item._rlnLogLikeliContribution = em.Float(row.getValue('rlnLogLikeliContribution'))
        item._rlnMaxValueProbDistribution = em.Float(row.getValue('rlnMaxValueProbDistribution'))
        
    def _updateClass(self, item):
        classId = item.getObjId()
        if  classId in self._classesInfo:
            index, fn, row = self._classesInfo[classId]
            fn = fn + ":mrc"
            item.setAlignment3D()
            item.getRepresentative().setLocation(index, fn)
            item._rlnclassDistribution = em.Float(row.getValue('rlnClassDistribution'))
            item._rlnAccuracyRotations = em.Float(row.getValue('rlnAccuracyRotations'))
            item._rlnAccuracyTranslations = em.Float(row.getValue('rlnAccuracyTranslations'))
