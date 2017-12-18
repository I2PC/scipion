# **************************************************************************
# *
# * Authors:     Laura del Cano (ldelcano@cnb.csic.es) [1]
# *              J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [2]
# *
# * [1] BCU, Centro Nacional de Biotecnologia, CSIC
# * [2] SciLifeLab, Stockholm University
# *
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os

import pyworkflow.protocol.params as params
from pyworkflow.em.protocol import ProtParticlePickingAuto
from pyworkflow.em.convert import ImageHandler
import pyworkflow.utils as pwutils
from pyworkflow.utils.properties import Message

from convert import readSetOfCoordinates



class DogPickerProtPicking(ProtParticlePickingAuto):
    """ Protocol to pick particles in a set of micrographs using appion
    dogpicker.
    """
    _label = 'dogpicker'
        
    def __init__(self, **args):
        ProtParticlePickingAuto.__init__(self, **args)


    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):

        ProtParticlePickingAuto._defineParams(self, form)
        form.addParam('diameter', params.IntParam, default=100,
                   label='Diameter of particle')
        form.addParam('invert', params.BooleanParam, default=False,
                      label='Invert',
                      help = "Invert image before picking, DoG normally picks "
                             "white particles.")
        form.addParam('threshold', params.FloatParam, default=0.5,
                      label='Threshold',
                      help = "Threshold in standard deviations above the mean, "
                             "e.g. --thresh=0.7")
        form.addParam('extraParams', params.StringParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Additional parameters',
                      help='Additional parameters for dogpicker: \n  '
                           '--numberSizes, --sizeRange, --threshold,...')

    #--------------------------- STEPS functions -------------------------------
    def _pickMicrograph(self, mic, args):
        # Prepare mic folder and convert if needed
        micName = mic.getFileName()
        micDir = self._getTmpPath(pwutils.removeBaseExt(micName))
        pwutils.makePath(micDir)

        ih = ImageHandler()
        # If needed convert micrograph to mrc format, otherwise link it
        if pwutils.getExt(micName) != ".mrc":
            fnMicBase = pwutils.replaceBaseExt(micName, 'mrc')
            inputMic = os.path.join(micDir, fnMicBase)
            ih.convert(mic.getLocation(), inputMic)
        else:
            inputMic = os.path.join(micDir, os.path.basename(micName))
            pwutils.createLink(micName, inputMic)

        # Program to execute and it arguments
        program = "ApDogPicker.py"
        outputFile = self._getExtraPath(pwutils.replaceBaseExt(inputMic, "txt"))

        args += " --image=%s --outfile=%s" % (inputMic, outputFile)

        self.runJob(program, args)
    
    def createOutputStep(self):
        pass

    #--------------------------- INFO functions --------------------------------
    def _summary(self):
        summary = []
        summary.append("Number of input micrographs: %d"
                       % self.getInputMicrographs().getSize())
        if self.getOutputsSize() > 0:
            summary.append("Number of particles picked: %d"
                           % self.getCoords().getSize())
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
        methodsMsgs.append("Input micrographs %s of size %d."
                           % (self.getObjectTag(self.getInputMicrographs()),
                              self.getInputMicrographs().getSize()))

        if self.getOutputsSize() > 0:
            output = self.getCoords()
            methodsMsgs.append('%s: User picked %d particles with a particle '
                               'size of %d and threshold %0.2f.'
                               % (self.getObjectTag(output), output.getSize(),
                                  output.getBoxSize(), self.threshold))
        else:
            methodsMsgs.append(Message.TEXT_NO_OUTPUT_CO)

        return methodsMsgs

    #--------------------------- UTILS functions -------------------------------
    def _getPickArgs(self):
        inputSampling = self.getInputMicrographs().getSamplingRate()
        args = "--diam=%0.3f " % (inputSampling * self.diameter.get())
        args += "--apix=%0.3f " % inputSampling
        args += "--thresh=%f" % self.threshold

        if self.invert:
            args += " --invert"

        args += " " + self.extraParams.get('')

        return [args]

    def readCoordsFromMics(self, workingDir, micList, coordSet):
        coordSet.setBoxSize(self.diameter.get())
        readSetOfCoordinates(workingDir, micList, coordSet)
    
    def getCoordsDir(self):
        return self._getExtraPath()
 