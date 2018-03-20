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

from glob import glob

import pyworkflow.em as em
import pyworkflow.em.metadata as md
from pyworkflow.utils.path import moveFile

from pyworkflow.em.protocol.protocol_particles import ProtProcessParticles
from pyworkflow.protocol.params import PointerParam, EnumParam
from pyworkflow.em.packages.relion.convert import (writeSetOfParticles,
                                                   locationToRelion, relionToLocation)
from pyworkflow.em.packages.relion.protocol_base import ProtRelionBase

STACK_NONE = 0
STACK_MULT = 1
STACK_ONE = 2


class ProtRelionExportParticles(ProtProcessParticles, ProtRelionBase):
    """ Export particles from Relion to be used outside Scipion. """
    _label = 'export particles'
    
    # --------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      label="Input particles", important=True,
                      help='Select the input images from the project.')

        form.addParam('stackType', EnumParam,
                      choices=["Don't write stacks",
                               "Write multiple stacks",
                               "Write a single stack"], default=STACK_MULT,
                      display=EnumParam.DISPLAY_LIST,
                      label="Iteration to visualize",
                      help="")

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

        # Create links to binary files and write the relion .star file
        writeSetOfParticles(imgSet, self._getPath("particles.star"),
                            outputDir=self._getExtraPath(),
                            writeAlignment=False,
                            postprocessImageRow=self._postprocessImageRow)
    
    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        validateMsgs = []
        return validateMsgs
    
    def _summary(self):
        summary = []
        return summary
    
    # --------------------------- UTILS functions ------------------------------
    def _postprocessImageRow(self, img, row):
        self._count = getattr(self, '_count', 1)

        if self._count == 1:
            row.printDict()

        if self._stackType > STACK_NONE:
            row.setValue('rlnOriginalParticleName',
                         row.getValue('rlnImageName'))
            newLoc = locationToRelion(self._count, 'particles.mrcs')
            row.setValue('rlnImageName', newLoc)

        self._count += 1