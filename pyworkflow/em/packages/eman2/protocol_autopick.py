# **************************************************************************
# *
# * Authors:     Airen Zaldivar (azaldivar@cnb.csic.es) [1]
# *              J.M. de la Rosa Trevin (delarosatrevin@scilifelab.se) [2]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# * [2] Science for Life Laboratory, Stockholm University
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

from pyworkflow import VERSION_1_1
from pyworkflow.protocol.params import IntParam, FloatParam
from pyworkflow.em.protocol import ProtParticlePickingAuto

import eman2
from convert import readSetOfCoordinates



class SparxGaussianProtPicking(ProtParticlePickingAuto):
    """
    Protocol to pick particles automatically in a set of micrographs
    using sparx gaussian picker.
    For more information see http://sparx-em.org/sparxwiki/e2boxer
    """
    _label = 'sparx gaussian picker'
    _version = VERSION_1_1
        
    def __init__(self, **kwargs):
        ProtParticlePickingAuto.__init__(self, **kwargs)
        self.extraParams = 'pixel_input=1:pixel_output=1:invert_contrast=True:use_variance=True'

    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        ProtParticlePickingAuto._defineParams(self, form)
        form.addParam('boxSize', IntParam, default=100,
                      label='Box Size', help='Box size in pixels')
        line = form.addLine('Picker range',
                            help='CCF threshold range for automatic picking')
        line.addParam('lowerThreshold', FloatParam, default='1',
                      label='Lower')
        line.addParam('higherThreshold', FloatParam, default='30',
                      label='Higher')

        form.addParam('gaussWidth', FloatParam, default='1',
                      label='Gauss Width',
                      help='Width of the Gaussian kernel used')

    #--------------------------- INSERT steps functions ------------------------
    def _insertInitialSteps(self):
        initId = self._insertFunctionStep('initSparxDb',
                                 self.lowerThreshold.get(),
                                 self.higherThreshold.get(),
                                 self.boxSize.get(), self.gaussWidth.get())
        return [initId]

    #--------------------------- STEPS functions -------------------------------
    def initSparxDb(self, lowerThreshold, higherThreshold, boxSize, gaussWidth):
        args = {"lowerThreshold": lowerThreshold,
                "higherThreshold": higherThreshold,
                "boxSize": boxSize,
                "gaussWidth": gaussWidth,
                "extraParams": self.extraParams}
        params = 'demoparms --makedb=thr_low=%(lowerThreshold)s:'
        params += 'thr_hi=%(higherThreshold)s:boxsize=%(boxSize)s:'
        params += 'gauss_width=%(gaussWidth)s:%(extraParams)s'

        self.runJob('sxprocess.py', params % args, cwd=self.getWorkingDir())

    def _pickMicrograph(self, mic, *args):
        micFile = os.path.relpath(mic.getFileName(), self.workingDir.get())
        params = ('--gauss_autoboxer=demoparms --write_dbbox --boxsize=%d %s'
                  % (self.boxSize, micFile))
        self.runJob('e2boxer.py', params, cwd=self.getWorkingDir()) 
        
    def createOutputStep(self):
        coordSet = self._createSetOfCoordinates(self.getInputMicrographs())
        self.readSetOfCoordinates(self.workingDir.get(), coordSet)
        coordSet.setBoxSize(self.boxSize.get())
        self._defineOutputs(outputCoordinates=coordSet)
        self._defineSourceRelation(self.inputMicrographs, coordSet)
    
    #--------------------------- INFO functions --------------------------------
    def _validate(self):
        errors = []
        eman2.validateVersion(self, errors)
        return errors

    #--------------------------- UTILS functions -------------------------------
    def getFiles(self):
        return (self.inputMicrographs.get().getFiles() |
                ProtParticlePickingAuto.getFiles(self))

    def readSetOfCoordinates(self, workingDir, coordSet):
        readSetOfCoordinates(workingDir, self.inputMicrographs.get(), coordSet)

