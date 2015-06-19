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

from pyworkflow.protocol.params import IntParam, FloatParam, BooleanParam, StringParam, LEVEL_ADVANCED
from pyworkflow.em.protocol import ProtParticlePicking
from pyworkflow.em.convert import ImageHandler
from pyworkflow.utils.properties import Message

from convert import readSetOfCoordinates

from pyworkflow.utils.path import replaceBaseExt, join, makePath, removeBaseExt, createLink, basename, getExt


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
        form.addParam('invert', BooleanParam, default=False,
                   label='Invert', help = "Invert image before picking, DoG normally picks white particles.")
        form.addParam('threshold', FloatParam, default=0.5,
                    label='Threshold', help = "Threshold in standard deviations above the mean, e.g. --thresh=0.7")
        # form.addParam('numberSizes', IntParam, expertLevel=LEVEL_ADVANCED,
        #            label='Number of sizes', help = "Number of different sizes to try.")
        # form.addParam('sizeRange', IntParam, expertLevel=LEVEL_ADVANCED,
        #            label='Size range', help = "Size range in pixels about diam to search.")
        # form.addParam('maxThreshold', FloatParam, expertLevel=LEVEL_ADVANCED,
        #            label='Max threshold', help = "Threshold in standard deviations above the mean, e.g. --thresh=0.7")
        # form.addParam('maxArea', FloatParam, expertLevel=LEVEL_ADVANCED,
        #            label='Max area', help = "When thresholded the peak must be less than maxarea*pi*r^2.")
        # form.addParam('maxPeaks', FloatParam, expertLevel=LEVEL_ADVANCED,
        #            label='Max peaks', help = "Maximum number of allowed peaks.")
        form.addParam('extraParams', StringParam, expertLevel=LEVEL_ADVANCED,
              label='Additional parameters',
              help='Additional parameters for dogpicker: \n  --numberSizes, --sizeRange, --threshold,...')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):

        self._params = {}
        # diameter must be passed in Armstrongs and therefore should be converted
        self._params['diam'] = self.diameter.get() * self.getInputMicrographs().getSamplingRate()
        # self._params['num-slices'] = self.numberSizes.get()
        # self._params['size-range'] = self.sizeRange.get()
        self._params['apix'] = self.inputMicrographs.get().getSamplingRate()
        self._params['thresh'] = self.threshold.get()
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

        ih = ImageHandler()

        for mic in self.inputMicrographs.get():
            # Create micrograph folder
            micName = mic.getFileName()
            micDir = self._getTmpPath(removeBaseExt(micName))
            makePath(micDir)

            # If needed convert micrograph to mrc format, otherwise link it
            if getExt(micName) != ".mrc":
                fnMicBase = replaceBaseExt(micName, 'mrc')
                inputMic = join(micDir, fnMicBase)
                ih.convert(mic.getLocation(), inputMic)
            else:
                inputMic = join(micDir, basename(micName))
                createLink(micName, inputMic)

            # Insert step to execute program
            stepId = self._insertFunctionStep('executeDogpickerStep', inputMic, args)
            deps.append(stepId)


        self._insertFinalSteps(deps)


    def _insertFinalSteps(self, deps):
        # Insert step to create output objects
        self._insertFunctionStep('createOutputStep', prerequisites=deps)


    #--------------------------- STEPS functions ---------------------------------------------------
    def executeDogpickerStep(self, inputMic, args):

        # Program to execute and it arguments
        program = "ApDogPicker.py"
        outputFile = self._getExtraPath(replaceBaseExt(inputMic, "txt"))

        args += " --image=%s --outfile=%s" % (inputMic, outputFile)

        # Run the command with formatted parameters

        self._log.info('Launching: ' + program + ' ' + args)
        self.runJob(program, args)


    def createOutputStep(self):
        micSet = self.getInputMicrographs()
        coordSet = self._createSetOfCoordinates(micSet)
        self.readSetOfCoordinates(self._getExtraPath(), coordSet)
        coordSet.setBoxSize(self.diameter.get())
        self._defineOutputs(outputCoordinates=coordSet)
        self._defineSourceRelation(micSet, coordSet)

    
    #--------------------------- UTILS functions --------------------------------------------------
    def getFiles(self):
        filePaths = self.inputMicrographs.get().getFiles() | ProtParticlePicking.getFiles(self)
        return filePaths


    def readSetOfCoordinates(self, workingDir, coordSet):
        readSetOfCoordinates(workingDir, self.inputMicrographs.get(), coordSet)

    def _summary(self):
        summary = []
        summary.append("Number of input micrographs: %d" % self.getInputMicrographs().getSize())
        if(self.getOutputsSize() > 0):
            summary.append("Number of particles picked: %d" % self.getCoords().getSize())
            summary.append("Particle size: %d" % self.getCoords().getBoxSize())
            summary.append("Threshold: %0.2f" % self.threshold)
            if self.extraParams.hasValue():
                summary.append("And other parameters: %s" % self.extraParams)
        else:
            summary.append(Message.TEXT_NO_OUTPUT_CO)
        return summary

    def _citations(self):
        return ['Voss2009']

    def _methods(self):
        methodsMsgs = []
        if self.getInputMicrographs() is None:
            return ['Input micrographs not available yet.']
        methodsMsgs.append("Input micrographs %s of size %d." % (self.getObjectTag(self.getInputMicrographs()), self.getInputMicrographs().getSize()))

        if self.getOutputsSize() > 0:
            output = self.getCoords()
            methodsMsgs.append('%s: User picked %d particles with a particle size of %d and threshold %0.2f.'
                               % (self.getObjectTag(output), output.getSize(), output.getBoxSize(), self.threshold.get()))
        else:
            methodsMsgs.append(Message.TEXT_NO_OUTPUT_CO)

        return methodsMsgs