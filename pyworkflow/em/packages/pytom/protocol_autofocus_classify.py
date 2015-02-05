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
from pyworkflow.utils.path import createLink
from pyworkflow.em.constants import NO_INDEX
from pyworkflow.em.convert import ImageHandler

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
        ih = ImageHandler()
        for vol in self.inputVolumes.get():
            newFn = self._getVolTmp(vol)
            index, fn = vol.getLocation()
            if index == NO_INDEX and fn.endswith('.mrc'):
                createLink(fn, newFn)
            else:
                ih.convert(vol, newFn)
        args = os.path.join(os.environ['PYTOM_HOME'],"bin","createParticleListFromDir.py")+\
           " -d tmp -w 0.0 -o volumes.xml"
        self.runJob("python", args, cwd=self._getPath())
            
    def runAutofocus(self):
        args = os.path.join(os.environ['PYTOM_HOME'],"classification","auto_focus_classify.py")
        self.runJob("python", args, cwd=self._getExtraPath())
        # scipion run python /home/coss/usb_linux/scipion/scipion/software/em/pytom/classification/auto_focus_classify.py -p volumes.xml -k 1 -f 2 -m ../phantom.mrc -s 5 -i 10 -n 0.2 -g -2 -t 0.4

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
    
    #--------------------------- UTILS functions --------------------------------------------   
    
    def _getVolTmp(self, vol):
        return self._getTmpPath('volume_%06d.mrc' % vol.getObjId())
