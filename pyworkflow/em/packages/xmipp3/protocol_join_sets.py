# **************************************************************************
# *
# * Authors:     Adrian Quintana (aquintana@cnb.csic.es)
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
This sub-package contains wrapper around ML2D Xmipp program
"""

from pyworkflow.em import *  
from pyworkflow.utils import *  
import xmipp

class XmippProtJoinSets(ProtPreprocessMicrographs):
    """ Protocol to obtain a set of initial volumes. """
    _label = 'join sets'
    
    def __init__(self, **args):
        ProtPreprocessMicrographs.__init__(self, **args)
        
    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('inputMicrographs', MultiPointerParam, label="Input set of micrographs", important=True, 
                      pointerClass='SetOfMicrographs', minNumObjects=1, maxNumObjects=0,
                      help='Select the input set of micrographs from the project.'
                           'They should 2 or more SetOfMicrographs classes')
        
#        form.addParallelSection(mpi=2)
         
        
    def _defineSteps(self):
        
        self._insertFunctionStep('createOutput')
  
    def createOutput(self):
       
        micsSet = self._createSetOfMicrographs()
        
        # Get pointer to input micrographs 
        for setOfMicrograph in self.inputMicrographs:
            setOfMicrograph.printAll()
            for mic in setOfMicrograph.get():
                micsSet.append(mic)
                
        micsSet.write()
        
        self._defineOutputs(outputMicrographs=micsSet)

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputMicrographs'):
            summary.append("Output micrographs not ready yet.")
        else:
            summary.append("RANSAC iterations: ")

            return summary
            