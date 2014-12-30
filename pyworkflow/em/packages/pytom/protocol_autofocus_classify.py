# **************************************************************************
# *
# * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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
# *  e-mail address 'jgomez@cnb.csic.es'
# *
# **************************************************************************
"""
This sub-package contains wrapper around auto_focus_classify algorithm in Pytom
"""

import os
from glob import glob

import pyworkflow.protocol.params as params

import pyworkflow.em as em  
import pytom



class ProtAutofocusClassify(em.ProtClassify3D):
    """ Subtomogram averaging using pytom autofocus """
    _label = 'autofocus'

    #--------------------------- DEFINE param functions --------------------------------------------
    
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputVolumes', params.PointerParam, label="Set of volumes", important=True, 
                      pointerClass='SetOfVolumes',
                      help='Subtomograms to average')
        form.addParam('numberOfReferences',params.IntParam,label='Number of references', default=4,
                      help="How many references are computed at the end of the process")

        # form.addParallelSection(threads=8, mpi=0)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    
    def _insertAllSteps(self):
        """ Mainly prepare the command line for calling reconstruct_significant program"""
        
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('runAutofocus')
        self._insertFunctionStep('createOutputStep')        

    #--------------------------- STEPS functions --------------------------------------------        
    def convertInputStep(self):
        return
        self.inputClasses.get().writeStack(self._getExtraPath("classes.spi:stk"))
            
    def runAutofocus(self):
        args = os.path.join(os.environ['PYTOM_HOME'],"classification","auto_focus_classify.py")
        self.runJob("python", args, cwd=self._getExtraPath())

    def createOutputStep(self):
        pass

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        summary.append("Input classes: %s" % self.inputClasses.get().getNameId())
        summary.append("Starting from: %d random volumes"%self.Nvolumes.get())
        return summary
    
    def _citations(self):
        return ['Chen2014']
    
    def _methods(self):
        return []
