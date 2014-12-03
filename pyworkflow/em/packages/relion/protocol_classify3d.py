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

import xmipp
from pyworkflow.em.protocol import ProtClassify3D
from pyworkflow.em.data import Volume

from pyworkflow.em.packages.relion.protocol_base import ProtRelionBase


class ProtRelionClassify3D(ProtClassify3D, ProtRelionBase):
    """    
    Protocol to classify 3D using Relion. Relion employs an empirical
    Bayesian approach to refinement of (multiple) 3D reconstructions
    or 2D class averages in electron cryo-microscopy (cryo-EM). Many
    parameters of a statistical model are learned from the data,which
    leads to objective and high-quality results.
    """
    _label = '3D classify'
    CHANGE_LABELS = [xmipp.MDL_AVG_CHANGES_ORIENTATIONS, 
                     xmipp.MDL_AVG_CHANGES_OFFSETS]
    
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
        args['--healpix_order'] = self.angularSamplingDeg.get()
        if self.localAngularSearch:
            args['--sigma_ang'] = self.localAngularSearchRange.get() / 3.
    
    #--------------------------- STEPS functions --------------------------------------------
    def createOutputStep(self):
        from pyworkflow.em.packages.relion.convert import readSetOfClasses3D
        
        imgSet = self.inputParticles.get()
        
        # create a SetOfClasses3D and define its relations
        classesSqlite = self._getIterClasses(self._lastIter(), clean=True)
        classes = self._createSetOfClasses3D(imgSet)
        readSetOfClasses3D(classes, classesSqlite)
        
        self._defineOutputs(outputClasses=classes)
        self._defineSourceRelation(imgSet, classes)
        self._defineSourceRelation(self.referenceVolume.get(), classes)
        
        # create a SetOfVolumes and define its relations
        volumes = self._createSetOfVolumes()
        volumes.setSamplingRate(imgSet.getSamplingRate())
        
        for ref3d in range(1, self.numberOfClasses.get()+1):
            vol = Volume()
            vol.setFileName(self._getFileName('volume', iter=self._lastIter(), ref3d=ref3d))
            volumes.append(vol)
        
        self._defineOutputs(outputVolumes=volumes)
        self._defineSourceRelation(imgSet, volumes)
    
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
        from pyworkflow.em.packages.xmipp3 import getMdFirstRow
        summary = []
        it = self._lastIter()
        if it >= 1:
            row = getMdFirstRow('model_general@' + self._getFileName('model', iter=it))
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
