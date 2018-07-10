# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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
import glob
import json

from pyworkflow.em.protocol import EMProtocol
from pyworkflow.protocol.params import PointerParam
from convert import runPhenixProgram, getProgram
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.em import PdbFile


class PhenixProtRunSuperposePDBs(EMProtocol):
    """Superpose two PDBs so that they optimally match """
    _label = 'superpose pdbs'
    _program = ""
    # _version = VERSION_1_2

    # --------------------------- DEFINE param functions -------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputStructureFixed', PointerParam, pointerClass="PdbFile",
                      label='Fixed atomic structure', allowsNull=True,
                      help="The moving PDB will be aligned to the fixed one")
        form.addParam('inputStructureMoving', PointerParam,
                      pointerClass="PdbFile",
                      label='Moving atomic structure',
                      help="PDBx/mmCIF to be aligned")
        
    # --------------------------- INSERT steps functions ---------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runSuperposePDBsStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------

    def runSuperposePDBsStep(self):
        args = os.path.abspath(self.inputStructureFixed.get().getFileName())
        args += " "+os.path.abspath(self.inputStructureMoving.get().getFileName())
        self.runJob(getProgram("phenix.superpose_pdbs"), args, cwd=self._getExtraPath())

    def createOutputStep(self):
        fnPdb = os.path.basename(self.inputStructureMoving.get().getFileName())
        pdb = PdbFile()
        pdb.setFileName(self._getExtraPath(fnPdb+"_fitted.pdb"))
        pdb.setVolume(self.inputStructureFixed.get().getVolume())
        self._defineOutputs(outputPdb=pdb)
        self._defineSourceRelation(self.inputStructureFixed.get(), pdb)
        self._defineSourceRelation(self.inputStructureMoving.get(), pdb)

    # --------------------------- INFO functions ---------------------------

    def _validate(self):
        errors = []
        # Check that the program exists
        program = getProgram("phenix.superpose_pdbs")
        if not os.path.exists(program):
            errors.append("Cannot find phenix.superpose_pdbs")

        # If there is any error at this point it is related to config variables
        if errors:
            errors.append("Check configuration file: "
                          "~/.config/scipion/scipion.conf")
            errors.append("and set PHENIX_HOME variables properly.")
            if program is not None:
                errors.append("Current values:")
                errors.append("PHENIX_HOME = %s" % os.environ['PHENIX_HOME'])

        return errors
