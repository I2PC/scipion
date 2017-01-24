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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
In this module are protocol related to EM imports of Masks...
"""

from os.path import exists, basename

from pyworkflow.protocol.params import PathParam, FloatParam
import pyworkflow.utils as pwutils

from base import ProtImport
from pyworkflow.em.data import Mask, VolumeMask
from pyworkflow.em.convert import ImageHandler



class ProtImportMask(ProtImport):
    """ Class for import masks from existing files. """
    _label = 'import mask'
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        
        form.addSection(label='Import')

        form.addParam('maskPath', PathParam, 
                      label="Mask path",
                      help="Select the file path of the mask\n")
        form.addParam('samplingRate', FloatParam, default=1.,
                   label=pwutils.properties.Message.LABEL_SAMP_RATE)

    def _insertAllSteps(self):
        self._insertFunctionStep('importMaskStep', self.maskPath.get(), self.samplingRate.get())

    #--------------------------- STEPS functions ---------------------------------------------------

    def importMaskStep(self, path, samplingRate):
        """ Copy mask from maskPath.
        """
        self.info("Using mask path: '%s'" % path)

        # Copy the image file into the project
        dst = self._getExtraPath(basename(path))
        pwutils.copyFile(path, dst)

        # Retrive image dimensions
        imgh = ImageHandler()
        _, _, z, n = imgh.getDimensions(dst)

        # Create a 2D or 3D Mask, consider the case of n>1
        # as the case of some volume maps in mrc format
        if z > 1 or n > 1:
            mask = VolumeMask()
        else:
            mask = Mask()

        mask.setFileName(dst)
        mask.setSamplingRate(samplingRate)

        self._defineOutputs(outputMask=mask)

    #--------------------------- INFO functions ----------------------------------------------------
    def _validate(self):
        errors = []
        if not self.maskPath.hasValue():
            errors.append("Mask path cannot be empty.")
        elif not exists(self.maskPath.get()):
            errors.append("Mask not found at *%s*" % self.maskPath.get())
        if not self.samplingRate.hasValue():
            errors.append("Sampling rate cannot be empty.")
        return errors

    def _summary(self):
        summary = ['Mask file imported from *%s*' % self.maskPath.get()]
        return summary
