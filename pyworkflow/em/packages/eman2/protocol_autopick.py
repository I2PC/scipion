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

        form.addParam('extraParams', StringParam, expertLevel=LEVEL_ADVANCED,
              label='Additional parameters',
              help='Additional parameters for sparx guassian picker: \n  invert_contrast, use_variance,...')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):

#         self._params = {}
#         # diameter must be passed in Armstrongs and therefore should be converted
#         self._params['diam'] = self.diameter.get() * self.getInputMicrographs().getSamplingRate()
#         # self._params['num-slices'] = self.numberSizes.get()
#         # self._params['size-range'] = self.sizeRange.get()
#         self._params['apix'] = self.inputMicrographs.get().getSamplingRate()
#         self._params['thresh'] = self.threshold.get()
#         # self._params['max-thresh'] = self.maxThreshold.get()
#         # self._params['max-area'] = self.maxArea.get()
#         # self._params['max-peaks'] = self.maxPeaks.get()
# 
#         args = ""
#         for par, val in self._params.iteritems():
#             args += " --%s=%s" % (par, str(val))
# 
#         if self.invert:
#             args += " --invert"
# 
#         args += " " + self.extraParams.get('')
# 
#         deps = [] # Store all steps ids, final step createOutput depends on all of them
# 
#         ih = ImageHandler()
# 
#         for mic in self.inputMicrographs.get():
#             # Create micrograph folder
#             micName = mic.getFileName()
#             micDir = self._getTmpPath(removeBaseExt(micName))
#             makePath(micDir)
# 
#             # If needed convert micrograph to mrc format, otherwise link it
#             if getExt(micName) != ".mrc":
#                 fnMicBase = replaceBaseExt(micName, 'mrc')
#                 inputMic = join(micDir, fnMicBase)
#                 ih.convert(mic.getLocation(), inputMic)
#             else:
#                 inputMic = join(micDir, basename(micName))
#                 createLink(micName, inputMic)
# 
#             # Insert step to execute program
#             stepId = self._insertFunctionStep('executeDogpickerStep', inputMic, args)
#             deps.append(stepId)
# 
# 
#         self._insertFinalSteps(deps)
        pass


    def _insertFinalSteps(self, deps):
        # Insert step to create output objects
        self._insertFunctionStep('createOutputStep', prerequisites=deps)


    #--------------------------- STEPS functions ---------------------------------------------------
    def executeSparxGaussianPickerStep(self, inputMic, args):

        # Program to execute and it arguments
#         program = "ApDogPicker.py"
#         outputFile = self._getExtraPath(replaceBaseExt(inputMic, "txt"))
# 
#         args += " --image=%s --outfile=%s" % (inputMic, outputFile)
# 
#         # Run the command with formatted parameters
# 
#         self._log.info('Launching: ' + program + ' ' + args)
#         self.runJob(program, args)
        pass


    def createOutputStep(self):
        coordSet = self._createSetOfCoordinates(self.getInputMicrographs())
        self.readSetOfCoordinates(self._getExtraPath(), coordSet)
        coordSet.setBoxSize(self.diameter.get())
        self._defineOutputs(outputCoordinates=coordSet)
        self._defineSourceRelation(self.inputMicrographs, coordSet)

    
    #--------------------------- UTILS functions --------------------------------------------------
    def getFiles(self):
        filePaths = self.inputMicrographs.get().getFiles() | ProtParticlePicking.getFiles(self)
        return filePaths


    def readSetOfCoordinates(self, workingDir, coordSet):
        readSetOfCoordinates(workingDir, self.inputMicrographs.get(), coordSet)

  

   