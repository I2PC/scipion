# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              Marta Martinez (mmmtnez@cnb.csic.es)
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
from pyworkflow.object import Float, Integer
from convert import runPhenixProgram, getProgram, MOLPROBITY
from protocol_refinement_base import PhenixProtRunRefinementBase

class PhenixProtRunMolprobity(PhenixProtRunRefinementBase):
    """MolProbity is a Phenix application to validate the geometry of an
atomic structure derived from a cryo-EM density map.
"""
    _label = 'molprobity'
    _program = ""
    #_version = VERSION_1_2
    MOLPROBITYFILE = 'molprobity.mrc'

    # --------------------------- DEFINE param functions -------------------
    def _defineParams(self, form):
        super(PhenixProtRunMolprobity, self)._defineParams(form)

    # --------------------------- INSERT steps functions --------------------

    def _insertAllSteps(self):
        if (self.inputVolume.get() or self.inputStructure.get().getVolume()) \
                is not None:
            self._insertFunctionStep('convertInputStep', self.MOLPROBITYFILE)
        self._insertFunctionStep('runMolprobityStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------

    def runMolprobityStep(self):
        # PDBx/mmCIF
        pdb = os.path.abspath(self.inputStructure.get().getFileName())
        args = ""
        args += pdb
        # starting volume (.mrc)
        tmpMapFile = self.MOLPROBITYFILE
        vol = os.path.abspath(self._getExtraPath(tmpMapFile))
        if vol is not None:
            args += " "
            volume = os.path.abspath(self._getExtraPath(self.MOLPROBITYFILE))
            args += "map_file_name=%s" % volume
            args += " "
            args += "d_min=%f" % self.resolution.get()
        args += " "
        args += "pickle=True"
            # args += " wxplots=True" # TODO: Avoid the direct opening of plots
        # script with auxiliary files
        try:
            runPhenixProgram(getProgram(MOLPROBITY), args,
                         cwd=self._getExtraPath())
        except:
            print "WARNING!!!\nPHENIX error:\n pdb_interpretation.clash_guard" \
                  " failure: Number of nonbonded interaction distances < 0.5: " \
                  "109 This error has been disable by running the same " \
                  "command with the same following additional " \
                  "argument:\npdb_interpretation.clash_guard." \
                  "nonbonded_distance_threshold=None "
            args += " "
            args += "pdb_interpretation.clash_guard." \
                    "nonbonded_distance_threshold=None"
            runPhenixProgram(getProgram(MOLPROBITY), args,
                             cwd=self._getExtraPath())


    def createOutputStep(self):
        MOLPROBITYOUTFILENAME = self._getExtraPath(
            self.MOLPROBITYOUTFILENAME)
        self._parseFile(MOLPROBITYOUTFILENAME)
        self._store()

    # --------------------------- INFO functions ---------------------------
    def _validate(self):
        errors = self.validateBase(MOLPROBITY,'MOLPROBITY')
        return errors

    def _summary(self):
        summary = PhenixProtRunRefinementBase._summary(self)
        summary.append("MolProbity: http://molprobity.biochem.duke.edu/")
        return summary

    def _methods(self):
        methodsMsgs = []
        methodsMsgs.append("TODO")
        return methodsMsgs

    def _citations(self):
        return ['Chen_2010']
