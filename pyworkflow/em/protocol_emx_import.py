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
In this module are protocol base classes related to EM.
Them should be sub-classes in the different sub-packages from
each EM-software packages.
"""

import os
import shutil
from glob import glob
from protocol import *
from constants import *
from data import * 
import emx


class ProtEmxImport(ProtClassify):
    """
    Import micrographs, coordinates or particles from EMX file.
    
    EMX is a joint initiative for data exchange format between different 
    EM software packages. See more about [[http://i2pc.cnb.csic.es/emx][EMX format]]
    """
        
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputEMX', FileParam, 
                      label="Input EMX file",
                      help='Provide the path to a valid EMX file.')
                 
    #--------------------------- INSERT steps functions --------------------------------------------  
    def _insertAllSteps(self):
        self._insertFunctionStep('importDataStep', self.inputSet.get().getObjId())       

    #--------------------------- STEPS functions --------------------------------------------       
    def importDataStep(self, micsId):
        """ Export micrographs to EMX file.
        micsId is only passed to force redone of this step if micrographs change.
        """
        emxFn = self._getPath('emxData')
        emx.exportData(emxDir, self.inputSet.get(), ctfSet=None)
        emx._iterXml(emxFn, {'micrographs': self._handleMicrograph})
    
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        errors = []
        return errors
    
    def _citations(self):
        cites = []
        return cites
    
    def _summary(self):
        summary = []
        return summary
    
    def _methods(self):
        return self._summary()  # summary is quite explicit and serve as methods
            
        
