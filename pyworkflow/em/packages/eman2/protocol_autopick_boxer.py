# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
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

from pyworkflow import VERSION_1_2
from pyworkflow.protocol.params import (IntParam, FloatParam,
                                        EnumParam, PointerParam)
from pyworkflow.em.protocol import ProtParticlePickingAuto
from pyworkflow.utils import makePath

import eman2
from convert import readSetOfCoordinates, convertReferences
from constants import *


class EmanProtAutopick(ProtParticlePickingAuto):
    """ Protocol to pick particles automatically in a set of micrographs
    using EMAN2 boxer.
    """
    _label = 'boxer auto'
    _lastUpdateVersion = VERSION_1_2

    @classmethod
    def isDisabled(cls):
        return not eman2.isNewVersion()

    def _createFilenameTemplates(self):
        """ Centralize the names of the files. """

        myDict = {'goodRefsFn': self._getExtraPath('info/boxrefs.hdf'),
                  'badRefsFn': self._getExtraPath('info/boxrefsbad.hdf'),
                  'bgRefsFn': self._getExtraPath('info/bgrefsbad.hdf')
                  }
        self._updateFilenamesDict(myDict)

    def __init__(self, **kwargs):
        ProtParticlePickingAuto.__init__(self, **kwargs)

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        ProtParticlePickingAuto._defineParams(self, form)
        form.addParam('boxSize', IntParam, default=128,
                      label='Box size (px)',
                      help="Box size in pixels.")
        form.addParam('particleSize', IntParam, default=100,
                      label='Particle size (px)',
                      help="Longest axis of particle in pixels (diameter, "
                           "not radius).")
        form.addParam('boxerMode', EnumParam,
                      choices=['local search', 'by ref', 'gauss', 'neural net'],
                      label="Autopicker mode:", default=AUTO_LOCAL,
                      display=EnumParam.DISPLAY_COMBO,
                      help="Choose autopicker mode:\n\n"
                           " _local search_ - Reference based search by "
                           "downsampling and 2-D alignment to references.\n"
                           " _by ref_ - simple reference-based "
                           "cross-correlation picker with exhaustive "
                           "rotational search.\n"
                           " _gauss_ - Gaussian sparx picker.\n"
                           " _neural net_ - convolutional neural network "
                           "boxer.")
        form.addParam('threshold', FloatParam, default='5.0',
                      label='Threshold')

        form.addSection('References')
        form.addParam('goodRefs', PointerParam,
                      pointerClass='SetOfAverages',
                      important=True,
                      label="Good references",
                      help="Good particle references.")
        form.addParam('badRefs', PointerParam,
                      pointerClass='SetOfAverages',
                      allowsNull=True,
                      label="Bad references",
                      help="Bad particle references like ice contamination "
                           "or large aggregation.")
        form.addParam('bgRefs', PointerParam,
                      pointerClass='SetOfAverages',
                      allowsNull=True,
                      label="Background references",
                      help="Pure noise regions in micrograph.")

        form.addParallelSection(threads=1, mpi=0)

    # --------------------------- INSERT steps functions ------------------------
    def _insertInitialSteps(self):
        self._createFilenameTemplates()
        initId = self._insertFunctionStep('convertInputStep')
        return [initId]

    # --------------------------- STEPS functions -------------------------------
    def convertInputStep(self):
        goodRefs = self.goodRefs.get() if self.goodRefs.hasValue() else None
        badRefs = self.badRefs.get() if self.badRefs.hasValue() else None
        bgRefs = self.bgRefs.get() if self.bgRefs.hasValue() else None
        storePath = self._getExtraPath("info")
        makePath(storePath)
        output = [self._getFileName('goodRefsFn'),
                  self._getFileName('badRefsFn'),
                  self._getFileName('bgRefsFn')]

        for i, refs in enumerate([goodRefs, badRefs, bgRefs]):
            if refs is not None:
                convertReferences(refs, output[i])

    def _pickMicrograph(self, mic, *args):
        micFile = os.path.relpath(mic.getFileName(), self.getCoordsDir())
        params = " --apix=%f --no_ctf" % self.inputMicrographs.get().getSamplingRate()
        params += " --boxsize=%d" % self.boxSize.get()
        params += " --ptclsize=%d" % self.particleSize.get()
        params += " --threads=%d" % self.numberOfThreads.get()

        modes = ['auto_local', 'auto_ref', 'auto_gauss', 'auto_convnet']
        params += " --autopick=%s:threshold=%0.2f" % (
            modes[self.boxerMode.get()], self.threshold.get())

        params += ' %s' % micFile
        program = eman2.getBoxerCommand(eman2.getVersion())

        self.runJob(program, params, cwd=self.getCoordsDir())

    def createOutputStep(self):
        pass

    # --------------------------- INFO functions --------------------------------
    def _validate(self):
        errors = []
        eman2.validateVersion(self, errors)
        if self.boxerMode.get() == AUTO_GAUSS:
            errors.append('Gauss mode is not implemented for new e2boxer yet.')
        if self.boxerMode.get() == AUTO_CONVNET:
            if not (self.badRefs.hasValue() and self.bgRefs.hasValue()):
                errors.append('Neural net picker requires all three types of references.')

        return errors

    # --------------------------- UTILS functions -------------------------------
    def getCoordsDir(self):
        return self._getExtraPath()

    def getFiles(self):
        return (self.inputMicrographs.get().getFiles() |
                ProtParticlePickingAuto.getFiles(self))

    def readCoordsFromMics(self, workingDir, micList, coordSet):
        coordSet.setBoxSize(self.boxSize.get())
        readSetOfCoordinates(workingDir, micList, coordSet, newBoxer=True)
