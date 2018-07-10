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
from pyworkflow.protocol.params import BooleanParam, PointerParam, FloatParam, IntParam
from convert import runPhenixProgram, getProgram
from pyworkflow.em.convert_header.CCP4.convert import adaptFileToCCP4, START
from pyworkflow.object import String
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.em import PdbFile


class PhenixProtRunRSRefine(EMProtocol):
    """Autorefine refines a PDB comparing it to a volume.
"""
    _label = 'real space refine'
    _program = ""
    # _version = VERSION_1_2

    # --------------------------- DEFINE param functions -------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVolume', PointerParam, pointerClass="Volume",
                      label='Input Volume', allowsNull=True,
                      help="Set the starting volume")
        form.addParam('inputStructure', PointerParam,
                      pointerClass="PdbFile",
                      label='Input atomic structure',
                      help="PDBx/mmCIF to be refined")
        form.addParam("doSecondary", BooleanParam, label="Secondary structure", default=False, expertLevel=LEVEL_ADVANCED)
        form.addParam("resolution", FloatParam, label="Resolution", default=3.0, help="Current resolution of the map")
        form.addParam("macroCycles", IntParam, label="Macro cycles", default=20, expertLevel=LEVEL_ADVANCED, help="Number of iterations of refinement")
        
    # --------------------------- INSERT steps functions ---------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('runRSrefineStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------

    def convertInputStep(self):
        """ convert 3D maps to MRC '.mrc' format
        """
        vol = self._getInputVolume()
        inVolName = vol.getFileName()
        newFn = self._getExtraPath("volume.mrc")
        origin = vol.getOrigin(force=True).getShifts()
        sampling = vol.getSamplingRate()
        adaptFileToCCP4(inVolName, newFn, origin, sampling, START)  # ORIGIN

    def runRSrefineStep(self):
        args = os.path.abspath(self.inputStructure.get().getFileName())
        args += " "+os.path.abspath(self._getExtraPath("volume.mrc"))
        args += " secondary_structure.enabled=%s"%self.doSecondary
        args += " resolution=%f"%self.resolution
        args += " run=minimization_global+local_grid_search+morphing+simulated_annealing"
        args += " macro_cycles=%d"%self.macroCycles
        self.runJob(getProgram("phenix.real_space_refine"), args, cwd=self._getExtraPath())

    def createOutputStep(self):
        fnPdb = os.path.basename(self.inputStructure.get().getFileName())
        pdb = PdbFile()
        pdb.setFileName(self._getExtraPath(fnPdb.replace(".pdb","_real_space_refined.pdb")))
        if self.inputVolume.get() is not None:
            pdb.setVolume(self.inputVolume.get())
        else:
            pdb.setVolume(self.inputStructure.get().getVolume())
        self._defineOutputs(outputPdb=pdb)
        self._defineSourceRelation(self.inputStructure.get(), pdb)
        if self.inputVolume.get() is not None:
            self._defineSourceRelation(self.inputVolume.get(), pdb)

    # --------------------------- INFO functions ---------------------------

    def _validate(self):
        errors = []
        # Check that the program exists
        program = getProgram("phenix.real_space_refine")
        if not os.path.exists(program):
            errors.append("Cannot find phenix.real_space_refine")

        # If there is any error at this point it is related to config variables
        if errors:
            errors.append("Check configuration file: "
                          "~/.config/scipion/scipion.conf")
            errors.append("and set PHENIX_HOME variables properly.")
            if program is not None:
                errors.append("Current values:")
                errors.append("PHENIX_HOME = %s" % os.environ['PHENIX_HOME'])

        # Check that the input volume exist
        if self._getInputVolume() is None:
            errors.append("Error: You should provide a volume.\n")

        return errors

    def _citations(self):
        return ['Barad_2015']

    # --------------------------- UTILS functions --------------------------
    def _getInputVolume(self):
        if self.inputVolume.get() is None:
            fnVol = self.inputStructure.get().getVolume()
        else:
            fnVol = self.inputVolume.get()
        return fnVol
