# **************************************************************************
# *
# * Authors:     Grigory Sharov (sharov@igbmc.fr)
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

from os.path import exists

import pyworkflow.protocol.params as params
import pyworkflow.protocol.constants as cons
import pyworkflow.utils.path as pwutils
from pyworkflow.em.protocol import ProtPreprocessMicrographs

from grigoriefflab import MAGDISTEST_PATH
from convert import parseMagEstOutput


class ProtMagDistEst(ProtPreprocessMicrographs):
    """ This program automatically estimates anisotropic magnification
    distortion from a set of images of a standard gold shadowed diffraction
    grating
    """    
    _label = 'magnification distortion estimation'

    def __init__(self, **args):
        ProtPreprocessMicrographs.__init__(self, **args)
        self.stepsExecutionMode = cons.STEPS_SERIAL

    # --------------------------- DEFINE params functions ----------------------

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputMicrographs', params.PointerParam,
                      pointerClass='SetOfMicrographs',
                      label="Input micrographs", important=True,
                      help='Select the SetOfMicrograph containing ~20 images '
                           'of different areas of polycrystalline gold')

        form.addParam('scaleFactor', params.FloatParam, default=0.03,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Scale factor',
                      help='Maximum allowed scale factor.')
        form.addParam('scaleStep', params.FloatParam, default=0.0005,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Scale step',
                      help='Step size for the scale search.')

        line = form.addLine('Resolution limit',
                            expertLevel=params.LEVEL_ADVANCED,
                            help='Resolution limits for the search. '
                                 'Default values are optimized '
                                 'for the gold diffraction ring')
        line.addParam('lowRes', params.FloatParam, default=2.5, label='Low')
        line.addParam('highRes', params.FloatParam, default=2.1, label='High')

        line = form.addLine('Angle range (deg)',
                            expertLevel=params.LEVEL_ADVANCED,
                            help='Allowed angle range for the search.')
        line.addParam('minAng', params.FloatParam, default=0.0, label='Min')
        line.addParam('maxAng', params.FloatParam, default=180.0, label='Max')

        form.addParam('angStep', params.FloatParam, default=0.1,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Angular step (deg)',
                      help='Step size for the angle search.')

        line = form.addLine('Filter radius (freq.)',
                            expertLevel=params.LEVEL_ADVANCED,
                            help='Filter radius for the amplitude bandpass '
                                 'filter.')

        line.addParam('lowp', params.FloatParam, default=0.2, label='Low-pass')
        line.addParam('highp', params.FloatParam, default=0.01,
                      label='High-pass')

        form.addParam('box', params.IntParam, default=512,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Amplitude box size',
                      help='Box size for the calculated amplitudes.')


        form.addParallelSection(threads=2, mpi=0)

    def _defineInputs(self):
        """ Store some of the input parameters in a dictionary for
        an easy replacement in the program command line.
        """
        pixSize = self.inputMicrographs.get().getSamplingRate()
        self.params = {'scaleFactor': self.scaleFactor.get(),
                       'scaleStep': self.scaleStep.get(),
                       'lowRes': self.lowRes.get(),
                       'highRes': self.highRes.get(),
                       'minAng': self.minAng.get(),
                       'maxAng': self.maxAng.get(),
                       'angStep': self.angStep.get(),
                       'lowp': self.lowp.get(),
                       'highp': self.highp.get(),
                       'box': self.box.get(),
                       'pixSize': pixSize,
                       'nthr': self.numberOfThreads.get()}

    # --------------------------- INSERT steps functions -----------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        parameters = self.runMagDistEst()
        self._argsMagDistEst()
        self._insertRunJobStep(self._program % parameters,
                               self._args % parameters)
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions ------------------------------

    def convertInputStep(self):
        """ Convert input micrographs into a single mrcs stack """
        inputMics = self.inputMicrographs.get()
        stackFn = self._getTmpPath('input_stack.mrcs')
        stackFnMrc = self._getTmpPath('input_stack.mrc')

        inputMics.writeStack(stackFn, applyTransform=False)
        # Grigorieff's program recognizes only mrc extension
        pwutils.moveFile(stackFn, stackFnMrc)

    def runMagDistEst(self):
        self._defineInputs()

        stackFnMrc = self._getTmpPath('input_stack.mrc')
        spectraFn = self.getOutputAmplitudes()
        rotAvgFn = self.getOutputAmplitudesRot()
        spectraCorrFn = self.getOutputAmplitudesCorr()
        logFn = self.getOutputLog()

        self.params['stackFnMrc'] = stackFnMrc
        self.params['spectraFn'] = spectraFn
        self.params['rotAvgFn'] = rotAvgFn
        self.params['spectraCorrFn'] = spectraCorrFn
        self.params['logFn'] = logFn

        return self.params

    def createOutputStep(self):
        pass

    # --------------------------- INFO functions -------------------------------

    def _validate(self):
        errors = []
        # Check that the program exists
        if not exists(MAGDISTEST_PATH):
            errors.append("Binary '%s' does not exits.\n"
                          "Check configuration file: \n"
                          "~/.config/scipion/scipion.conf\n"
                          "and set MAGDIST_HOME variable properly."
                          % MAGDISTEST_PATH)

        return errors

    def _citations(self):
        return ["Grant2015"]

    def _summary(self):
        summary = []
        result = self._parseOutputLog()

        if result is None or len(result) < 5:
            summary.append('Output is not ready yet.')
        else:
            result = tuple(result)
            distAngle, majorAxis, minorAxis, pixCorr, totDist = result
            summary.append('This protocol does not generate any output. It only'
                           ' estimates magnification distortion parameters.\n\n'
                           'Total amount of distortion: *%0.2f%%* ' % totDist)
            summary.append('Distortion Angle: *%0.2f* degrees' % distAngle)
            summary.append('Major Scale: *%0.3f*' % majorAxis)
            summary.append('Minor Scale: *%0.3f*' % minorAxis)
            summary.append('Corrected pixel size: *%0.3f* A' % pixCorr)

        return summary

    def _methods(self):
        txt = []
        txt.append("Anisotropic magnification distortion was estimated using "
                   "Grigorieff's program *mag_distortion_estimate*")

        return txt

    # --------------------------- UTILS functions ------------------------------

    def getOutputLog(self):
        return self._getExtraPath('mag_dist_estimation.log')

    def _parseOutputLog(self):
        """ Return the distortion amount, angle and two scale params. """
        fnOut = self.getOutputLog()

        return parseMagEstOutput(fnOut)

    def _argsMagDistEst(self):
        self._program = 'export NCPUS=%(nthr)d ; ' + MAGDISTEST_PATH
        self._args = """   << eof > %(logFn)s
%(stackFnMrc)s
%(spectraFn)s
%(rotAvgFn)s
%(spectraCorrFn)s
%(pixSize)f
YES
%(lowRes)f
%(highRes)f
%(scaleFactor)f
%(scaleStep)f
%(minAng)f
%(maxAng)f
%(angStep)f
%(lowp)f
%(highp)f
%(box)d
eof
"""

    def getOutputAmplitudes(self):
        return self._getExtraPath('output_amp.mrc')

    def getOutputAmplitudesCorr(self):
        return self._getExtraPath('output_amp_corrected.mrc')

    def getOutputAmplitudesRot(self):
        return self._getExtraPath('output_amp_rot.mrc')
