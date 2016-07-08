# **************************************************************************
# *
# * Authors:     Mohsen Kazemi  (mkazemi@cnb.csic.es)
# *              Joaquin Oton   (joton@cnb.csic.es)
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
import glob
import pyworkflow.protocol.params as params
from pyworkflow.em import Protocol

#from pyworkflow.em.data import Volume
#from pyworkflow.em.protocol import ProtReconstruct3D
#from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles
#from pyworkflow.utils import getFloatListFromValues
#from pyworkflow.utils.path import cleanPattern, cleanPath, copyFile
#import xmipp
#from pyworkflow.object import Float, String
#from math import sqrt
#from plotter import XmippPlotter


class XmippProtImport(Protocol):
    """    
    ??????
    """
    _label = 'import tilt-series'    
    #--------------------------- DEFINE param functions --------------------------------------------   
   
    def _defineParams(self, form):
        form.addSection(label='Input')

          
                      
        form.addParallelSection(threads=0, mpi=0)
    #--------------------------- INSERT steps functions --------------------------------------------
    
    def _insertAllSteps(self):
        
        
        
        
        
        
        #self._insertFunctionStep('reconstructionStep', ....)
                 
        #self._insertFunctionStep('gatherResultsStep', ....)        
    #--------------------------- STEPS functions --------------------------------------------
    
    #def reconstructionStep(self, ...):
        
            
    
    #def gatherResultsStep(self, ...):
         
    #--------------------------- INFO functions -------------------------------------------- 
    
#    def _summary(self):
#        """ Should be overriden in subclasses to 
#        return summary message for NORMAL EXECUTION. 
#        """
#              
#        msg = []
#        msg.append()        
#        msg.append()
#        return msg
#    
#    def _methods(self):
#        messages = []
#        messages.append('Joton')
#        return messages
#    
#    def _citations(self):
#        return ['Joton']
#    
#    def _validate(self):
#        errors=[]
#        if :
#            errors.append() 
#        if :
#            errors.append()
#        return errors 
             
    #--------------------------- UTILS functions --------------------------------------------
    
