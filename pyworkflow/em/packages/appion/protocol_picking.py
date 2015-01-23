# **************************************************************************
# *
# * Authors:     Laura del Cano (ldelcano@cnb.csic.es)
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

from pyworkflow.protocol.params import IntParam, FloatParam, BooleanParam, StringParam, LEVEL_EXPERT
from pyworkflow.em.protocol import ProtParticlePicking

from convert import readSetOfCoordinates

from pyworkflow.utils.path import replaceBaseExt, join


class DogPickerProtPicking(ProtParticlePicking):
    """Protocol to pick particles in a set of micrographs using appion dogpicker"""
    _label = 'dogpicker'
        
    def __init__(self, **args):     
        ProtParticlePicking.__init__(self, **args)


    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):

        ProtParticlePicking._defineParams(self, form)
        form.addParam('diameter', IntParam, default=100,
                   label='Diameter of particle')
        # form.addParam('numberSizes', IntParam, expertLevel=LEVEL_EXPERT,
        #            label='Number of sizes', help = "Number of different sizes to try.")
        # form.addParam('sizeRange', IntParam, expertLevel=LEVEL_EXPERT,
        #            label='Size range', help = "Size range in pixels about diam to search.")
        # form.addParam('pixelSize', IntParam, expertLevel=LEVEL_EXPERT,
        #            label='Pixel size (A)', help = "Pixel size of images in Angstroms.")
        # form.addParam('threshold', FloatParam, expertLevel=LEVEL_EXPERT,
        #            label='Threshold', help = "Threshold in standard deviations above the mean, e.g. --thresh=0.7")
        # form.addParam('maxThreshold', FloatParam, expertLevel=LEVEL_EXPERT,
        #            label='Max threshold', help = "Threshold in standard deviations above the mean, e.g. --thresh=0.7")
        # form.addParam('maxArea', FloatParam, expertLevel=LEVEL_EXPERT,
        #            label='Max area', help = "When thresholded the peak must be less than maxarea*pi*r^2.")
        # form.addParam('maxPeaks', FloatParam, expertLevel=LEVEL_EXPERT,
        #            label='Max peaks', help = "Maximum number of allowed peaks.")
        form.addParam('invert', BooleanParam, default=False,
                   label='Invert', help = "Invert image before picking, DoG normally picks white particles.")
        form.addParam('extraParams', StringParam, expertLevel=LEVEL_EXPERT,
              label='Additional parameters',
              help='Additional parameters for dogpicker: \n  --numberSizes, --sizeRange, --threshold,...')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):

        self._params = {}
        self._params['diam'] = self.diameter.get()
        # self._params['num-slices'] = self.numberSizes.get()
        # self._params['size-range'] = self.sizeRange.get()
        # self._params['apix'] = self.pixelSize.get()
        # self._params['thresh'] = self.threshold.get()
        # self._params['max-thresh'] = self.maxThreshold.get()
        # self._params['max-area'] = self.maxArea.get()
        # self._params['max-peaks'] = self.maxPeaks.get()

        args = ""
        for par, val in self._params.iteritems():
            args += " --%s=%s" % (par, str(val))

        if self.invert:
            args += " --invert"

        args += " " + self.extraParams.get('')

        deps = [] # Store all steps ids, final step createOutput depends on all of them

        for mic in self.inputMicrographs.get():
            stepId = self._insertFunctionStep('executeDogpickerStep', mic.getFileName(), args)
            deps.append(stepId)


        self._insertFinalSteps(deps)


    def _insertFinalSteps(self, deps):
        # Insert step to create output objects
        self._insertFunctionStep('createOutputStep', prerequisites=deps)


    #--------------------------- STEPS functions ---------------------------------------------------
    def executeDogpickerStep(self, inputMic, args):
        # Program to execute and it arguments
        program = "ApDogPicker.py"
        outputFile = join(self.getWorkingDir(), replaceBaseExt(inputMic, "txt"))

        args += " --image=%s --outfile=%s" % (inputMic, outputFile)

        # Run the command with formatted parameters

        self._log.info('Launching: ' + program + ' ' + args)
        self.runJob(program, args)



    def createOutputStep(self):
        micSet = self.getInputMicrographs()
        coordSet = self._createSetOfCoordinates(micSet)
        self.readSetOfCoordinates(self.getWorkingDir(), coordSet)
        coordSet.setBoxSize(self.diameter.get())
        self._defineOutputs(outputCoordinates=coordSet)
        self._defineSourceRelation(micSet, coordSet)

    
    #--------------------------- UTILS functions --------------------------------------------------
    def getFiles(self):
        filePaths = self.inputMicrographs.get().getFiles() | ProtParticlePicking.getFiles(self)
        return filePaths


    def readSetOfCoordinates(self, workingDir, coordSet):
        readSetOfCoordinates(workingDir, self.inputMicrographs.get(), coordSet)