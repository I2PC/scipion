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
import glob
import math
from pyworkflow.object import String, Float, Integer
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.protocol.params import BooleanParam, PointerParam, FloatParam, IntParam
from convert import runPhenixProgram, getProgram, REALSPACEREFINE
from pyworkflow.em.headers import adaptFileToCCP4, START
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
        PhenixProtRunRefinementBase._defineParams(self, form)
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
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------
    def runRSrefineStep(self, tmpMapFile):
        args = os.path.abspath(self.inputStructure.get().getFileName())
        vol = os.path.abspath(self._getExtraPath(tmpMapFile))
        args += " " + vol
        args += " resolution=%f" % self.resolution
        args += " secondary_structure.enabled=%s"%self.doSecondary
        args += " run=minimization_global+local_grid_search+morphing+" \
                "simulated_annealing"
        args += " macro_cycles=%d"%self.macroCycles
        runPhenixProgram(getProgram(REALSPACEREFINE), args,
                         cwd=self._getExtraPath())

    def createOutputStep(self):
        fnPdb = os.path.basename(self.inputStructure.get().getFileName())
        pdb = PdbFile()
        pdb.setFileName(self._getExtraPath(
            fnPdb.replace(".pdb", "_real_space_refined.pdb")))
        if self.inputVolume.get() is not None:
            pdb.setVolume(self.inputVolume.get())
        else:
            pdb.setVolume(self.inputStructure.get().getVolume())
        self._defineOutputs(outputPdb=pdb)
        self._defineSourceRelation(self.inputStructure.get(), pdb)
        if self.inputVolume.get() is not None:
            self._defineSourceRelation(self.inputVolume.get(), pdb)

        fileName = glob.glob(self._getExtraPath("*.log"))[0]
        self._parseFile(fileName)
        self._store()

    # --------------------------- INFO functions ---------------------------

    def _validate(self):
        errors = self.validateBase(REALSPACEREFINE, 'REALSPACEREFINE')
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
    def _parseFile(self, fileName):
        with open(fileName) as f:
            line = f.readline()
            while line:
                words = line.strip().split()
                if len(words) > 1:
                    if (words[0] == 'Final' and words[1] == 'info'
                        and words[2] == '(overall)'):
                        line = f.readline()
                        while line:
                            words = line.strip().split()
                            if len(words) > 1:
                                if (words[0] == 'all-atom' and
                                    words[1] =='clashscore'):
                                    self.clashscore = Float(words[3])
                                elif (words[0] == 'outliers'and
                                      words[1] == ':'):
                                    self.ramachandranOutliers = Float(words[2])
                                elif (words[0] == 'favored'and
                                      words[1] == ':'):
                                    self.ramachandranFavored = Float(words[2])
                                elif (words[0] == 'rotamer' and
                                      words[1] == 'outliers'):
                                    self.rotamerOutliers = Float(words[3])
                                elif (words[0] == 'cbeta' and
                                      words[1] == 'deviations'):
                                    self.cbetaOutliers = Integer(words[3])

                            line = f.readline()
                line = f.readline()
        self.molprobity_score(self.clashscore.get(),
                              self.rotamerOutliers.get(),
                              self.ramachandranFavored.get()
                              )

    def molprobity_score(self, clashscore, rota_out, rama_fav):
        """
        Calculate the overall Molprobity score, as described here:
          http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2877634/?tool=pubmed
          http://kinemage.biochem.duke.edu/suppinfo/CASP8/methods.html
        """
        if (clashscore >= 0) and (rota_out >= 0) and (rama_fav >= 0):
            rama_iffy = 100. - rama_fav
            self.overallScore = Float(((0.426 * math.log(1 + clashscore)) +
                       (0.33 * math.log(1 + max(0, rota_out - 1))) +
                       (0.25 * math.log(1 + max(0, rama_iffy - 2)))) + 0.5)
        else:
            return -1  # FIXME prevents crashing on RNA
