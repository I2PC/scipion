# **************************************************************************
# *
# * Authors:     Laura del Cano (ldelcano@cnb.csic.es)
# *              Ignacio Foche (ifoche@cnb.csic.es)
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
In this module are protocol related to EM imports of Masks...
"""

from os.path import exists, basename

from pyworkflow.utils.properties import Message

from pyworkflow.protocol.params import (PathParam, FloatParam,
                                        EnumParam, StringParam, LEVEL_ADVANCED)
from pyworkflow.utils.path import copyFile
from pyworkflow.em.protocol.protocol_import.base import ProtImport
from pyworkflow.em.data import Mask, VolumeMask
from pyworkflow.em.convert import ImageHandler


class ProtImportMask(ProtImport):
    """ Class for import masks from existing files. """
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        
        form.addSection(label='Import')

        form.addParam('maskPath', PathParam, 
                      label="Mask path",
                      help="Select the file path of the mask\n")
        form.addParam('samplingRate', FloatParam, default=1.,
                   label=Message.LABEL_SAMP_RATE)

    def _insertAllSteps(self):
        self._insertFunctionStep('importMaskStep', self.maskPath.get(), self.samplingRate.get())

    #--------------------------- STEPS functions ---------------------------------------------------

    def importMaskStep(self, path, samplingRate):
        """ Copy mask from maskPath.
        """
        self.info("Using mask path: '%s'" % path)

        # Copy the image file into the project
        dst = self._getExtraPath(basename(path))
        copyFile(path, dst)

        # Retrive image dimensions
        imgh = ImageHandler()
        x, y, z, n = imgh.getDimensions(dst)

        # Create a 2D or 3D Mask
        if z == 1:
            mask = Mask()
        else:
            mask = VolumeMask()

        mask.setFileName(dst)
        mask.setSamplingRate(samplingRate)

        self._defineOutputs(outputVolume=mask)

    #--------------------------- INFO functions ----------------------------------------------------
    def _validate(self):
        return []
    
    #--------------------------- UTILS functions ---------------------------------------------------


