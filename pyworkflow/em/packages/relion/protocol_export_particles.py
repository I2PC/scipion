# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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

import pyworkflow.em as em
import pyworkflow.utils as pwutils

from pyworkflow.em.protocol.protocol_particles import ProtProcessParticles
import pyworkflow.protocol.params as params
from convert import writeSetOfParticles, locationToRelion
from protocol_base import ProtRelionBase

STACK_NONE = 0
STACK_MULT = 1
STACK_ONE = 2


class ProtRelionExportParticles(ProtProcessParticles, ProtRelionBase):
    """ Export particles from Relion to be used outside Scipion. """
    _label = 'export particles'
    
    # --------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):

        form.addSection(label='Input')

        form.addParam('inputParticles', params.PointerParam,
                      pointerClass='SetOfParticles',
                      label="Input particles", important=True,
                      help='Select the input images from the project.')

        form.addParam('useAlignment', params.BooleanParam, default=True,
                      label='Write alignment information?',
                      help='If *Yes* the alignment information (2D or 3D) '
                           'will be written to the resulting .star file if'
                           'the particles contains such information.')

        form.addParam('stackType', params.EnumParam,
                      choices=["Don't write stacks",
                               "Write multiple stacks",
                               "Write a single stack"], default=STACK_MULT,
                      display=params.EnumParam.DISPLAY_LIST,
                      label="Binary stack files",
                      help="If *Don't write stacks* is choose, only the star "
                           "files will be written out. Alternatively, you can "
                           "select to write images into a single stack file or"
                           " several stacks (one per micrograph). ")

    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        objId = self.inputParticles.get().getObjId()
        self._insertFunctionStep("exportParticlesStep", objId)

    # --------------------------- STEPS functions ------------------------------
    def exportParticlesStep(self, particlesId):
        """ Create the input file in STAR format as expected by Relion.
        If the input particles comes from Relion, just link the file. 
        """
        imgSet = self.inputParticles.get()
        self._stackType = self.stackType.get()
        self._ih = em.ImageHandler()
        self._stackDict = {}
        particlesPath = self._getPath('Particles')
        pwutils.cleanPath(particlesPath)
        pwutils.makePath(particlesPath)

        alignType = imgSet.getAlignment() if self.useAlignment else em.ALIGN_NONE
        # Create links to binary files and write the relion .star file
        writeSetOfParticles(imgSet, self._getPath("particles.star"),
                            outputDir=self._getExtraPath(),
                            alignType=alignType,
                            postprocessImageRow=self._postprocessImageRow,
                            fillMagnification=True)

        pwutils.prettyDict(self._stackDict)
    
    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        validateMsgs = []
        return validateMsgs
    
    def _summary(self):
        summary = []
        return summary
    
    # --------------------------- UTILS functions ------------------------------
    def _postprocessImageRow(self, img, row):
        """ Write the binary image to the final stack
        and update the row imageName. """

        if self._stackType > STACK_NONE:
            rlnImageName = row.getValue('rlnImageName')
            # backup the original name
            row.setValue('rlnOriginalParticleName', rlnImageName)

            if self._stackType == STACK_ONE:
                self._count = getattr(self, '_count', 1)
                index, stackName = (self._count, 'particles.mrcs')
                self._count += 1
            else: # STACK_MULT
                baseName = pwutils.removeBaseExt(img.getFileName())
                if baseName not in self._stackDict:
                    self._stackDict[baseName] = 0
                index = self._stackDict[baseName] + 1
                stackName = baseName + '.mrcs'
                self._stackDict[baseName] = index

            stackFn = self._getPath('Particles', stackName)

            self._ih.convert(img, (index, stackFn))
            # Store relative path in the star file
            relStackFn = os.path.relpath(stackFn, self._getPath())
            row.setValue('rlnImageName', locationToRelion(index, relStackFn))
