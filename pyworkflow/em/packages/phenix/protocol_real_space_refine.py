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
from pyworkflow.protocol.params import BooleanParam,  IntParam, EnumParam
from convert import runPhenixProgram, getProgram, REALSPACEREFINE, MOLPROBITY
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.em import PdbFile
from protocol_refinement_base import PhenixProtRunRefinementBase

PDB = 0
mmCIF = 1
OUTPUT_FORMAT = ['pdb', 'mmcif']


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
        form.addParam('outputFormat', EnumParam, choices=OUTPUT_FORMAT,
                      default=PDB, label="Select output format",
                      help="Refined atomic structure is the protocol output. "
                           "You can choose PDB or mmCIF as output format. "
                           "PDB format has been selected by default.")
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
        group = form.addGroup('Optimization strategy options')
        group.addParam('minimizationGlobal', BooleanParam,
                       label="Global minimization: ", default=True,
                       expertLevel=LEVEL_ADVANCED,
                       help="Phenix default parameter to look the global "
                            "minimum of the model.\nGenerally, refinement "
                            "with all defaults is "
                            "sufficient.\nOther options "
                            "of use: run=minimization_global+local_grid_search"
                            "+morphing+simulated_annealing\n")
        group.addParam('rigidBody', BooleanParam,
                       label="Rigid body: ", default=False,
                       expertLevel=LEVEL_ADVANCED,
                       help="Refinement strategy that considers groups of "
                            "atoms that move (rotate and translate) as a "
                            "single body.\n")
        group.addParam('localGridSearch', BooleanParam,
                       label="Local grid search: ", default=False,
                       expertLevel=LEVEL_ADVANCED,
                       help="Refinement strategy that considers "
                            "local rotamer fitting.\n\n Generally, refinement "
                            "with all defaults is sufficient.\n Including "
                            "local fitting, morphing, "
                            "or simulated annealing "
                            "( local_grid_search+morphing+simulated_annealing) "
                            "into refinement may significantly increase "
                            "runtime.\nOther options "
                            "of use: run=minimization_global+local_grid_search"
                            "+morphing+simulated_annealing\n")
        group.addParam('morphing', BooleanParam,
                       label="Morphing ", default=False,
                       expertLevel=LEVEL_ADVANCED,
                       help="Morphing procedure distorts a model to match an "
                            "electron density map.\n\nGenerally, refinement "
                            "with all defaults is "
                            "sufficient.\n Including local fitting, morphing, "
                            "or simulated annealing "
                            "( local_grid_search+morphing+simulated_annealing) "
                            "into refinement may significantly increase "
                            "runtime.\nOther options "
                            "of use: run=minimization_global+local_grid_search"
                            "+morphing+simulated_annealing\n")
        group.addParam('simulatedAnnealing', BooleanParam,
                       label="Simulated annealing ", default=False,
                       expertLevel=LEVEL_ADVANCED,
                       help="Optimization technique known as molecular "
                            "dynamics refinement; it minimizes the energy of "
                            "the model.\n"
                            "Generally, refinement with all defaults is "
                            "sufficient.\n Including local fitting, morphing, "
                            "or simulated annealing "
                            "( local_grid_search+morphing+simulated_annealing) "
                            "into refinement may significantly increase "
                            "runtime.\nOther options "
                            "of use: run=minimization_global+local_grid_search"
                            "+morphing+simulated_annealing\n")
        group.addParam('adp', BooleanParam,
                       label="Atomic Displacement Parameters (ADPs) ",
                       default=True,
                       expertLevel=LEVEL_ADVANCED,
                       help="Phenix default parameter.\nGenerally, refinement "
                            "with all defaults is sufficient.\n\nADP ("
                            "B-factors) refinement against the map is "
                            "performed at the last macro-cycle only. ")

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
        args += " run="
        if self.minimizationGlobal == True:
            args += "minimization_global+"
        if self.rigidBody == True:
            args += "rigid_body+"
        if self.localGridSearch == True:
            args += "local_grid_search+"
        if self.morphing == True:
            args += "morphing+"
        if self.simulatedAnnealing == True:
            args += "simulated_annealing+"
        if self.adp == True:
            args += "adp+"
        args = args[:-1]
        args += " macro_cycles=%d" % self.macroCycles
        args += " model_format=%s" % OUTPUT_FORMAT[int(self.outputFormat)]
        args += " write_pkl_stats=True"
        runPhenixProgram(getProgram(REALSPACEREFINE), args,
                         cwd=self._getExtraPath())

    def runMolprobityStep(self, tmpMapFile):
        # PDBx/mmCIF
        self._getRSRefineOutput()
        args = ""
        args += os.path.abspath(self.outPdbName)
        # starting volume (.mrc)
        vol = os.path.abspath(self._getExtraPath(tmpMapFile))
        args += " "
        args += "map_file_name=%s" % vol
        args += " "
        args += "d_min=%f" % self.resolution.get()
        args += " "
        args += "pickle=True"
        # args += " wxplots=True" # Direct opening of the three wx plots (
                                  # Ramachandran, Chi1-Chi2 and
                                  # Multi-criterion plots)
                                  # script with auxiliary files
        runPhenixProgram(getProgram(MOLPROBITY), args,
                         cwd=self._getExtraPath())

    def createOutputStep(self):
        self._getRSRefineOutput()
        pdb = PdbFile()
        pdb.setFileName(self.outPdbName)

        if self.inputVolume.get() is not None:
            pdb.setVolume(self.inputVolume.get())
        else:
            pdb.setVolume(self.inputStructure.get().getVolume())
        self._defineOutputs(outputPdb=pdb)
        self._defineSourceRelation(self.inputStructure.get(), pdb)
        if self.inputVolume.get() is not None:
            self._defineSourceRelation(self.inputVolume.get(), pdb)

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

    def _getRSRefineOutput(self):
        inPdbName = os.path.basename(self.inputStructure.get().getFileName())
        if self.outputFormat == PDB:
            self.outPdbName = self._getExtraPath(
                inPdbName.replace("." + inPdbName.split(".")[1],
                                  "_real_space_refined.pdb"))
        elif self.outputFormat == mmCIF:
            self.outPdbName = self._getExtraPath(
                inPdbName.replace("." + inPdbName.split(".")[1],
                                  "_real_space_refined.cif"))

