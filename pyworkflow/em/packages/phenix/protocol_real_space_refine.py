# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
# *              Marta Martinez (mmmtnez@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
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
from pyworkflow.protocol.params import BooleanParam,  IntParam
from convert import runPhenixProgram, getProgram, REALSPACEREFINE, MOLPROBITY
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.em import PdbFile
from protocol_refinement_base import PhenixProtRunRefinementBase


class PhenixProtRunRSRefine(PhenixProtRunRefinementBase):
    """Tool for extensive real-space refinement of an atomic structure
    against the map provided. The map can be derived from X-ray or neutron
    crystallography, or cryoEM. The program obtains a model that fits the map
    as well as possible having appropriate geometry. The model should not show
    validation outliers, such as Ramachandran plot or rotamer outliers.
    """
    _label = 'real space refine'
    _program = ""
    # _version = VERSION_1_2
    REALSPACEFILE = 'real_space.mrc'

    # --------------------------- DEFINE param functions -------------------
    def _defineParams(self, form):
        super(PhenixProtRunRSRefine, self)._defineParams(form)
        form.addParam("doSecondary", BooleanParam, label="Secondary structure",
                      default=False, expertLevel=LEVEL_ADVANCED,
                      help="Set to TRUE to use secondary structure "
                           "restraints.\n")
        form.addParam("macroCycles", IntParam, label="Macro cycles",
                      default=5, expertLevel=LEVEL_ADVANCED,
                      help="Number of iterations of refinement.\nAlthough 5 "
                           "macro-cycles is usually sufficient, in cases in "
                           "which model geometry or/and model-to-map fit is "
                           "poor the use of more macro-cycles could be "
                           "helpful.\n")

    # --------------------------- INSERT steps functions ---------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep', self.REALSPACEFILE)
        self._insertFunctionStep('runRSrefineStep', self.REALSPACEFILE)
        self._insertFunctionStep('runMolprobityStep', self.REALSPACEFILE)
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------
    def runRSrefineStep(self, tmpMapFile):
        pdb = os.path.abspath(self.inputStructure.get().getFileName())
        args = "" + pdb
        vol = os.path.abspath(self._getExtraPath(tmpMapFile))
        args += " " + vol
        args += " resolution=%f" % self.resolution
        args += " secondary_structure.enabled=%s" % self.doSecondary
        # args += " run=minimization_global+local_grid_search+morphing+" \
        #         "simulated_annealing"
        args += " macro_cycles=%d" % self.macroCycles
        args += " write_pkl_stats=True"

        runPhenixProgram(getProgram(REALSPACEREFINE), args,
                         cwd=self._getExtraPath())

    def runMolprobityStep(self, tmpMapFile):
        # PDBx/mmCIF
        outPdbName = self._getPdbRSRefineOutput()
        args = ""
        args += os.path.abspath(outPdbName)
        # starting volume (.mrc)
        vol = os.path.abspath(self._getExtraPath(tmpMapFile))
        args += " "
        args += "map_file_name=%s" % vol
        args += " "
        args += "d_min=%f" % self.resolution.get()
        args += " "
        args += "pickle=True"
        # args += " wxplots=True" # TODO: Avoid the direct opening of plots
        # script with auxiliary files
        runPhenixProgram(getProgram(MOLPROBITY), args,
                         cwd=self._getExtraPath())

    def createOutputStep(self):
        outPdbName = self._getPdbRSRefineOutput()
        pdb = PdbFile()
        pdb.setFileName(outPdbName)

        if self.inputVolume.get() is not None:
            pdb.setVolume(self.inputVolume.get())
        else:
            self.pdb.setVolume(self.inputStructure.get().getVolume())
        self._defineOutputs(outputPdb=pdb)
        self._defineSourceRelation(self.inputStructure.get(), pdb)
        if self.inputVolume.get() is not None:
            self._defineSourceRelation(self.inputVolume.get(), pdb)

        # fileName = glob.glob(self._getExtraPath("*.log"))[0]
        # self._parseFile(fileName)

        MOLPROBITYOUTFILENAME = self._getExtraPath(
            self.MOLPROBITYOUTFILENAME)
        self._parseFile(MOLPROBITYOUTFILENAME)
        self._store()

    # --------------------------- INFO functions ---------------------------

    def _validate(self):
        errors = self.validateBase(REALSPACEREFINE, 'REALSPACEREFINE')

        # Check that the input volume exist
        if self._getInputVolume() is None:
            errors.append("Error: You should provide a volume.\n")
        return errors

    def _citations(self):
        return ['Barad_2015']

    def _summary(self):
        summary = PhenixProtRunRefinementBase._summary(self)
        summary.append(
            "https://www.phenix-online.org/documentation/reference/"
            "real_space_refine.html")
        return summary

    # --------------------------- UTILS functions --------------------------

    def _getPdbRSRefineOutput(self):
        inPdbName = os.path.basename(self.inputStructure.get().getFileName())
        outPdbName = self._getExtraPath(
            inPdbName.replace(".pdb", "_real_space_refined.pdb"))
        return outPdbName