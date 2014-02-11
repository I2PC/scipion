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
        
        form.addParam('inputSet', MultiPointerParam, label="Input set of images", important=True, 
                      pointerClass='SetOfImages', minNumObjects=2, maxNumObjects=0,
                      help='Select the set of images (it can be a set of micrographs, particles o volumes) from the project.'
                           'They should 2 or more object classes')
        
    def _insertAllSteps(self):
        self._insertFunctionStep('createOutput')
  
    def createOutput(self):
        self.inputType = str(self.inputSet[0].get().getClassName())
       
        func ="_create%s" % self.inputType
        
        print "func", func
       
        outputSetFunction = getattr(self, func)
        
        outputSet = outputSetFunction()
        outputSet.printAll()
        outputSet.copyInfo(self.inputSet[0].get())
       
        for itemSet in self.inputSet:
            for itemObj in itemSet.get():
                itemObj.cleanObjId()
                outputSet.append(itemObj)
                
        outputSet.write()
       
       
#        micsSet = self._createSetOfMicrographs()
#        micsSet.copyInfo(self.inputMicrographs[0].get())
#        
#        # Get pointer to input micrographs 
#        for setOfMicrograph in self.inputMicrographs:
#            for mic in setOfMicrograph.get():
#                mic.cleanObjId()
#                micsSet.append(mic)
#                
#        micsSet.write()
        
        self._defineOutputs(outputSet=outputSet)

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputSet'):
            summary.append("Output set not ready yet.")
        else:
            summary.append("Input sets of type [%s]:" % self.inputType)
            for itemSet in self.inputSet:
                summary.append("[%s]" % itemSet.get().getNameId())
        return summary
        
    def _methods(self):
        methods = []
        if not hasattr(self, 'outputSet'):
            methods.append("Protocol has not finished yet.")
        else:
            methods.append("We have joint the following sets ")
        
        return methods
            