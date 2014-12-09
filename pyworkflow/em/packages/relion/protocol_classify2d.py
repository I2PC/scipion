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
This module contains the protocol base class for Relion protocols
"""
import pyworkflow.em as em
import pyworkflow.em.metadata as md
from pyworkflow.em.data import SetOfClasses2D
from pyworkflow.em.protocol import ProtClassify2D

from pyworkflow.em.packages.relion.protocol_base import ProtRelionBase
from pyworkflow.em.packages.relion.convert import fillClasses, rowToAlignment, relionToLocation


class ProtRelionClassify2D(ProtRelionBase, ProtClassify2D):
    """ Wrapper to Relion 2D - class averages program.
    """
    _label = '2d classify'
    IS_2D = True
    OUTPUT_TYPE = SetOfClasses2D
    
    def __init__(self, **args):        
        ProtRelionBase.__init__(self, **args)
        
    def _initialize(self):
        """ This function is mean to be called after the 
        working dir for the protocol have been set. (maybe after recovery from mapper)
        """
        ProtRelionBase._initialize(self)
        self.ClassFnTemplate = '%(ref)03d@%(rootDir)s/relion_it%(iter)03d_classes.mrcs'
        

    #--------------------------- INSERT steps functions --------------------------------------------  
    def _setSamplingArgs(self, args):
        """ Set sampling related params. """
        # Sampling stuff            
        args['--psi_step'] = self.inplaneAngularSamplingDeg.get() * 2

    #--------------------------- STEPS functions --------------------------------------------       
    def createOutputStep(self):
        from pyworkflow.utils.path import findRootFrom
        
        outImgSet = self._createSetOfParticles()
        imgSet = self.inputParticles.get()
        
        data = md.MetaData(self._getFileName('data', iter=self._lastIter()))
        outImgSet.copyInfo(imgSet)
        outImgSet.copyItems(imgSet,
                    updateItemCallback=self.updatePartProperties,
                    itemDataIterator=md.iterRows(data))
        
        self._classesDict = {} # store classes info, indexed by class id
         
        modelMd = md.MetaData('model_classes@' + self._getFileName('model', iter=self._lastIter()))
        for classNumber, objId in enumerate(modelMd):
            row = md.Row()
            row.readFromMd(modelMd, objId)
            index, fn = relionToLocation(row.getValue('rlnReferenceImage'))
            
            self._classesDict[classNumber+1] = (index, fn, row)
        
        clsSet = self._createSetOfClasses2D(outImgSet)
        fillClasses(clsSet, updateClassCallback=self._updateClass)
        
        self._defineOutputs(outputParticles=outImgSet)
        self._defineTransformRelation(imgSet, outImgSet)
        self._defineOutputs(outputClasses=clsSet)
        self._defineSourceRelation(imgSet, clsSet)
        
    #--------------------------- INFO functions -------------------------------------------- 
    def _validateNormal(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        return []
    
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
        summary.append("Input Particles: *%d*\nClassified into *%d* classes\n" % (self.inputParticles.get().getSize(),
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
            strline += 'We classified %d particles into %d classes using Relion Classify2d. '%\
                           (self.inputParticles.get().getSize(), self.numberOfClasses.get())
        return [strline]
    
    #--------------------------- UTILS functions --------------------------------------------
    def updatePartProperties(self, item, row):
        
        item.setClassId(row.getValue(md.RLN_PARTICLE_CLASS))
        item.setTransform(rowToAlignment(row, em.ALIGN_2D))
        
        item._rlnNormCorrection = em.Float(row.getValue('rlnNormCorrection'))
        item._rlnLogLikeliContribution = em.Float(row.getValue('rlnLogLikeliContribution'))
        item._rlnMaxValueProbDistribution = em.Float(row.getValue('rlnMaxValueProbDistribution'))
        
    def _updateClass(self, item):
        classId = item.getObjId()
        if  classId in self._classesDict:
            index, fn, row = self._classesDict[classId]
            item.getRepresentative().setLocation(index, fn)
            item._rlnclassDistribution = em.Float(row.getValue('rlnClassDistribution'))
            item._rlnAccuracyRotations = em.Float(row.getValue('rlnAccuracyRotations'))
            item._rlnAccuracyTranslations = em.Float(row.getValue('rlnAccuracyTranslations'))

