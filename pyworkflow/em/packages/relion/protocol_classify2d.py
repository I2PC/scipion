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

import pyworkflow.em as em
import pyworkflow.em.metadata as md
from pyworkflow.em.data import SetOfClasses2D
from pyworkflow.em.protocol import ProtClassify2D

from pyworkflow.em.packages.relion.protocol_base import ProtRelionBase
from pyworkflow.em.packages.relion.convert import rowToAlignment, relionToLocation


class ProtRelionClassify2D(ProtRelionBase, ProtClassify2D):
    """ Wrapper around the Relion 2D classification program."""
    _label = '2D classification'
    IS_2D = True
    OUTPUT_TYPE = SetOfClasses2D
    
    def __init__(self, **args):        
        ProtRelionBase.__init__(self, **args)
        
    def _initialize(self):
        """ This function is mean to be called after the 
        working dir for the protocol have been set.
        (maybe after recovery from mapper)
        """
        ProtRelionBase._initialize(self)
        self.ClassFnTemplate = '%(ref)03d@%(rootDir)s/relion_it%(iter)03d_classes.mrcs'
        

    #--------------------------- INSERT steps functions ------------------------
    def _setSamplingArgs(self, args):
        """ Set sampling related params. """
        # Sampling stuff
        if self.doImageAlignment:
            args['--offset_range'] = self.offsetSearchRangePix.get()
            args['--offset_step']  = self.offsetSearchStepPix.get() * 2
            args['--psi_step'] = self.inplaneAngularSamplingDeg.get() * 2
        else:
            args['--skip_align'] = ''

    #--------------------------- STEPS functions -------------------------------
    def _loadClassesInfo(self, iteration):
        """ Read some information about the produced Relion 2D classes
        from the *model.star file.
        """
        self._classesInfo = {} # store classes info, indexed by class id
         
        modelStar = md.MetaData('model_classes@%s' %
                                self._getFileName('model', iter=iteration))
        
        for classNumber, row in enumerate(md.iterRows(modelStar)):
            index, fn = relionToLocation(row.getValue('rlnReferenceImage'))
            # Store info indexed by id, we need to store the row.clone() since
            # the same reference is used for iteration            
            self._classesInfo[classNumber+1] = (index, fn, row.clone())
    
    def _fillClassesFromIter(self, clsSet, iteration):
        """ Create the SetOfClasses2D from a given iteration. """
        self._loadClassesInfo(iteration)
        dataStar = self._getFileName('data', iter=iteration)
        clsSet.classifyItems(updateItemCallback=self._updateParticle,
                             updateClassCallback=self._updateClass,
                             itemDataIterator=md.iterRows(dataStar,
                                                          sortByLabel=md.RLN_IMAGE_ID))
        
    def createOutputStep(self):
        partSet = self.inputParticles.get()       
        
        classes2D = self._createSetOfClasses2D(partSet)
        self._fillClassesFromIter(classes2D, self._lastIter())
        
        self._defineOutputs(outputClasses=classes2D)
        self._defineSourceRelation(self.inputParticles, classes2D)
        
    #--------------------------- INFO functions --------------------------------
    def _validateNormal(self):
        """ Should be overwritten in subclasses to
        return summary message for NORMAL EXECUTION. 
        """
        return []
    
    def _validateContinue(self):
        """ Should be overwritten in subclasses to
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
        """ Should be overwritten in subclasses to
        return summary message for NORMAL EXECUTION. 
        """
        summary = []
        summary.append("Input Particles: %s" % self.getObjectTag('inputParticles'))
        summary.append("Classified into *%d* classes." % self.numberOfClasses)
        summary.append("Output set: %s" % self.getObjectTag('outputClasses'))
        
        return summary
    
    def _summaryContinue(self):
        """ Should be overwritten in subclasses to
        return summary messages for CONTINUE EXECUTION.
        """
        summary = []
        summary.append("Continue from iteration %01d" % self._getContinueIter())
        
        return summary
    
    def _methods(self):
        methods = ''
        
        if hasattr(self, 'outputClasses'):
            
            methods += "We classified input particles %s (%d items) " % (self.getObjectTag('inputParticles'),
                                                                         self.inputParticles.get().getSize())
            methods += "into %d classes using Relion Classify2d. " % self.numberOfClasses
            methods += 'Output classes: %s' % self.getObjectTag('outputClasses')
        return [methods]
    
    #--------------------------- UTILS functions --------------------------------------------
    def _updateParticle(self, item, row):
        #print "_updateParticle"
        #row.printDict()
        item.setClassId(row.getValue(md.RLN_PARTICLE_CLASS))
        item.setTransform(rowToAlignment(row, em.ALIGN_2D))
        
        item._rlnNormCorrection = em.Float(row.getValue('rlnNormCorrection'))
        item._rlnLogLikeliContribution = em.Float(row.getValue('rlnLogLikeliContribution'))
        item._rlnMaxValueProbDistribution = em.Float(row.getValue('rlnMaxValueProbDistribution'))
        
    def _updateClass(self, item):
        classId = item.getObjId()
        if  classId in self._classesInfo:
            index, fn, row = self._classesInfo[classId]
            item.setAlignment2D()
            item.getRepresentative().setLocation(index, fn)
            item._rlnclassDistribution = em.Float(row.getValue('rlnClassDistribution'))
            item._rlnAccuracyRotations = em.Float(row.getValue('rlnAccuracyRotations'))
            item._rlnAccuracyTranslations = em.Float(row.getValue('rlnAccuracyTranslations'))

