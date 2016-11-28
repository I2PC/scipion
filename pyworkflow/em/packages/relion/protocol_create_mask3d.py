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

from pyworkflow.em import ProtCreateMask3D, VolumeMask
from convert import convertBinaryVol
import pyworkflow.protocol.params as params


AND = 0
OR = 1
AND_NOT = 2
OR_NOT = 3


class ProtRelionCreateMask3D(ProtCreateMask3D):
    """ Create a 3D mask.
    The mask is created from a 3d volume or by comparing two input volumes.
    """
    _label = 'create 3d mask'
    
    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Mask generation')
        form.addParam('inputVolume', params.PointerParam, pointerClass="Volume",
                      label="Input volume",
                      help="Select the volume that will be used to create the mask")
        # TODO: add wizard
        form.addParam('threshold', params.FloatParam, default=0.01,
                      label='Threshold',
                      help='Select a threshold value for input volume binarization')
        form.addParam('doCompare', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Compare with another volume to produce a mask?',
                      help='Logical comparison of two input volumes to produce a mask')
        form.addParam('inputVolume2', params.PointerParam, pointerClass="Volume",
                      condition='doCompare',
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Input volume (second)",
                      help="Select the volume that will be compared to the first one")
        form.addParam('operation', params.EnumParam, default=AND,
                      condition='doCompare',
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Operation',
                      choices=['AND', 'OR', 'AND_NOT', 'OR_NOT'],
                      help='*AND*: Pixels in the initial mask will be one if the input '
                           'AND the second volume are above the threshold value.\n'
                           '*OR*: Pixels in the initial mask will be one if the input '
                           'OR the second volume are above the threshold value.\n'
                           '*AND_NOT*: pixels in the initial mask will be one if '
                           'the input is above the threshold AND the second volume '
                           'is below it.\n*OR_NOT*: pixels in the initial mask will be '
                           'one if the input is above the threshold OR the second '
                           'volume is below it.')
        form.addParam('extend', params.IntParam, default=0,
                      label='Extend the mask by (px)',
                      help='Extend initial binary mask this number of pixels')
        form.addParam('edge', params.IntParam, default=0,
                      label='Soft edge (px)',
                      help='Width (in pixels) of the additional soft edge on the binary mask')
        form.addParam('doInvert', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Invert final mask',
                      help='Invert the final mask')

    # --------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self.maskFile = self._getExtraPath('mask.mrc')
        self._insertFunctionStep('convertInputStep', self.inputVolume.get().getObjId())
        self._insertFunctionStep('createMaskStep')
        self._insertFunctionStep('createOutputStep')
    
    # --------------------------- STEPS functions -------------------------------
    def convertInputStep(self, volId):
        self.inputVolFn = convertBinaryVol(self.inputVolume.get(), self._getTmpPath())

        if self.doCompare:
            self.inputVol2Fn = convertBinaryVol(self.inputVolume2.get(), self._getTmpPath())

    def createMaskStep(self):
        argsDict = {'--i ': self.inputVolFn,
                    '--ini_threshold ': self.threshold.get(),
                    '--extend_inimask ': self.extend.get(),
                    '--width_soft_edge ': self.edge.get()
                    }
        args = ' --o %s ' % self.maskFile
        args += ' '.join(['%s %s' % (k, v) for k, v in argsDict.iteritems()])

        if self.doCompare:
            if self.operation.get() == AND:
                args += ' --and %s' % self.inputVol2Fn
            elif self.operation.get() == OR:
                args += ' --or %s' % self.inputVol2Fn
            elif self.operation.get() == AND_NOT:
                args += ' --and_not %s' % self.inputVol2Fn
            elif self.operation.get() == OR_NOT:
                args += ' --or_not %s' % self.inputVol2Fn

        if self.doInvert:
            args += '--invert'

        self.runJob("relion_mask_create", args)

        return [self.maskFile]

    def createOutputStep(self):
        volMask = VolumeMask()
        volMask.setFileName(self.maskFile)
        volMask.setSamplingRate(self.inputVolume.get().getSamplingRate())

        self._defineOutputs(outputMask=volMask)
        self._defineSourceRelation(self.inputVolume, self.outputMask)
        
    # --------------------------- INFO functions --------------------------------
    def _validate(self):
        errors = []
        self.validatePackageVersion('RELION_HOME', errors)

        if self.doCompare:
            pixVol1 = self.inputVolume.get().getSamplingRate()
            pixVol2 = self.inputVolume2.get().getSamplingRate()

            if pixVol1 != pixVol2:
                errors.append('Pixel size should be the same for both input volumes.')

        return errors

    def _summary(self):
        messages = []
        messages.append("Created a mask from input volume using threshold %0.3f" % self.threshold.get())
        messages.append("*Mask processing*")
        messages.append("   Extend by %d pixels" % self.extend.get())
        messages.append("   Apply soft edge of %d pixels" % self.edge.get())
        if self.doCompare:
            messages.append("   Logical operation: %s" % self.getEnumText('operation'))
        if self.doInvert:
            messages.append("   Inverted")

        return messages

    def _citations(self):
        return []

    def _methods(self):
        messages = []

        messages.append("*Mask creation*")
        messages.append('We processed the volume %s.' % self.inputVolume.get().getNameId())
        messages.append("We binarized it at threshold of %0.3f. " % self.threshold.get())
        messages.append("We extended binary mask by %d voxels." % self.extend.get())
        messages.append("And, we smoothed it by applying a soft edge of %d voxels." % self.edge.get())

        if self.doInvert:
            messages.append("We inverted the mask. ")
            messages.append("And, we smoothed it by applying a soft edge of %d voxels." % self.edge.get())
        if self.hasAttribute('outputMask'):
            messages.append('We refer to the output mask as %s.' % self.outputMask.getNameId())

        return messages
