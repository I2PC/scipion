# **************************************************************************
# *
# * Authors:     Airen Zaldivar (azaldivar@cnb.csic.es)
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

import os

from pyworkflow.protocol.params import IntParam, FloatParam, BooleanParam, StringParam, LEVEL_ADVANCED
from pyworkflow.em.protocol import ProtParticlePicking


from convert import readSetOfCoordinates



class SparxGaussianProtPicking(ProtParticlePicking):
    """Protocol to pick particles automatically in a set of micrographs using sparx gaussian picker"""
    _label = 'sparx gaussian picker'
        
    def __init__(self, **args):     
        ProtParticlePicking.__init__(self, **args)
        self.extraParams = 'pixel_input=1:pixel_output=1:invert_contrast=True:use_variance=True'


    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):

        ProtParticlePicking._defineParams(self, form)
        form.addParam('boxSize', IntParam, default=100,
                   label='Box Size')
        line = form.addLine('Picker range',
                            help='')
        line.addParam('lowerThreshold', FloatParam, default='1',
                      label='Lower')
        line.addParam('higherThreshold', FloatParam, default='30',
                      label='Higher')

        form.addParam('gaussWidth', FloatParam, default='1',
              label='Gauss Width')
        
    

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        args = {"lowerThreshold": self.lowerThreshold,
                "higherThreshold": self.higherThreshold,
                "boxSize": self.boxSize,
                "gaussWidth": self.gaussWidth,
                'extraParams': self.extraParams}
        params = 'demoparms --makedb=thr_low=%(lowerThreshold)s:thr_hi=%(higherThreshold)s:boxsize=%(boxSize)s:gauss_width=%(gaussWidth)s:%(extraParams)s'%args
        self.runJob('sxprocess.py', params, cwd=self.getWorkingDir()) 
        deps = [] # Store all steps ids, final step createOutput depends on all of them
        for mic in self.inputMicrographs.get():
            micFile = os.path.relpath(mic.getFileName(), self.workingDir.get())
            stepId = self._insertFunctionStep('executeSparxGaussianPickerStep', micFile)
            deps.append(stepId)


        self._insertFunctionStep('createOutputStep', prerequisites=deps)

    
    #--------------------------- STEPS functions ---------------------------------------------------
    def executeSparxGaussianPickerStep(self, micFile):
        print micFile
        params = '--gauss_autoboxer=demoparms --write_dbbox %s'%micFile
        self.runJob('e2boxer.py', params, cwd=self.getWorkingDir()) 
        


    def createOutputStep(self):
        coordSet = self._createSetOfCoordinates(self.getInputMicrographs())
        self.readSetOfCoordinates(self.workingDir.get(), coordSet)
        coordSet.setBoxSize(self.boxSize.get())
        self._defineOutputs(outputCoordinates=coordSet)
        self._defineSourceRelation(self.inputMicrographs, coordSet)

    
    #--------------------------- UTILS functions --------------------------------------------------
    def getFiles(self):
        filePaths = self.inputMicrographs.get().getFiles() | ProtParticlePicking.getFiles(self)
        return filePaths


    def readSetOfCoordinates(self, workingDir, coordSet):
        readSetOfCoordinates(workingDir, self.inputMicrographs.get(), coordSet)

  

   