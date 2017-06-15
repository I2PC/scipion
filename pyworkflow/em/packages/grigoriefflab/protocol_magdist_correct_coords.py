# **************************************************************************
# *
# * Authors:     Grigory Sharov (sharov@igbmc.fr)
# *
# * L'Institut de genetique et de biologie moleculaire et cellulaire (IGBMC)
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

import pyworkflow.protocol.params as params
from pyworkflow import VERSION_1_1
from pyworkflow.utils.properties import Message
from pyworkflow.em.protocol import EMProtocol
from convert import parseMagCorrInput, unDistortCoord


class ProtMagDistCorrCoord(EMProtocol):
    """ This program automatically corrects anisotropic magnification
    distortion using previously estimated parameters.
    It works on a set of coordinates.
    """
    _label = 'mag distortion correction (coords)'
    _lastUpdateVersion = VERSION_1_1
    # --------------------------- DEFINE params functions ----------------------

    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)

        form.addParam('inputCoords', params.PointerParam,
                      pointerClass='SetOfCoordinates',
                      important=True,
                      label='Input coordinates',
                      help='Select a set of coordinates.')
        form.addParam('useEst', params.BooleanParam, default=False,
                      label='Use previous estimation?',
                      help='Use previously calculated parameters of '
                           'magnification anisotropy.')
        form.addParam('inputEst', params.PointerParam,
                      pointerClass='ProtMagDistEst', condition='useEst',
                      label='Input protocol',
                      help='Select previously executed estimation protocol.')
        form.addParam('scaleMaj', params.FloatParam, default=1.0,
                      condition='not useEst',
                      label='Major scale factor',
                      help='Major scale factor.')
        form.addParam('scaleMin', params.FloatParam, default=1.0,
                      condition='not useEst',
                      label='Minor scale factor',
                      help='Minor scale factor.')
        form.addParam('angDist', params.FloatParam, default=0.0,
                      condition='not useEst',
                      label='Distortion angle (deg)',
                      help='Distortion angle, in degrees.')

        form.addParallelSection(threads=0, mpi=0)

    # --------------------------- STEPS functions ------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep')

    def createOutputStep(self):
        inputCoords = self.inputCoords.get()
        inputMics = inputCoords.getMicrographs()
        coordSet = self._createSetOfCoordinates(inputMics)

        coordSet.copyInfo(inputCoords)
        coordSet.copyItems(inputCoords, updateItemCallback=self._updateCoordinates)
        coordSet.setBoxSize(inputCoords.getBoxSize())

        self._storeSummary()

        self._defineOutputs(outputCoordinates=coordSet)
        self._defineSourceRelation(inputCoords, coordSet)

    # --------------------------- INFO functions -------------------------------

    def _validate(self):
        errors = []
        return errors

    def _citations(self):
        return ["Grant2015"]

    def _summary(self):
        return [self.summaryVar.get()]

    def _methods(self):
        txt = []
        txt.append("Anisotropic magnification distortion was corrected using "
                   "Grigorieff's program *mag_distortion_correct*")
        return txt

    # --------------------------- UTILS functions ------------------------------

    def _storeSummary(self):
        if self.getAttributeValue('useEst', False):
            inputFn = self.getAttributeValue('inputEst', None).getOutputLog()
            input_params = parseMagCorrInput(inputFn)
            self.summaryVar.set("The following magnification distortion parameters "
                                "were used for correction:\n\n"
                                "Distortion Angle: *%0.2f* degrees\n"
                                "Major Scale: *%0.3f*\n"
                                "Minor Scale: *%0.3f*\n"
                                % (input_params[0],
                                   input_params[1],
                                   input_params[2]))
        else:
            self.summaryVar.set("The following magnification distortion parameters "
                                "were used for correction:\n\n"
                                "Distortion Angle: *%0.2f* degrees\n"
                                "Major Scale: *%0.3f*\n"
                                "Minor Scale: *%0.3f*\n"
                                % (self.getAttributeValue('angDist', 1.0),
                                   self.getAttributeValue('scaleMaj', 1.0),
                                   self.getAttributeValue('scaleMin', 1.0)))

    def _updateCoordinates(self, coord, row):
        params = [coord.getX(), coord.getY()] + self._getParams()
        newX, newY = unDistortCoord(params)
        coord.setPosition(newX, newY)

    def _getParams(self):
        inputMics = self.inputCoords.get().getMicrographs()
        mic_x, mic_y, _ = inputMics.getFirstItem().getDim()

        if self.useEst:
            inputEst = self.inputEst.get().getOutputLog()
            input_params = parseMagCorrInput(inputEst)
            ang = input_params[0]
            major_scale = input_params[1]
            minor_scale = input_params[2]

        else:
            ang = self.angDist.get()
            major_scale = self.scaleMaj.get()
            minor_scale = self.scaleMin.get()

        return [mic_x, mic_y, ang, major_scale, minor_scale]
