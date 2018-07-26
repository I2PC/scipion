# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
# *              Marta Martinez (mmmtnez@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
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
from pyworkflow.object import String, Float, Integer


from pyworkflow.em.protocol import EMProtocol
from pyworkflow.protocol.params import PointerParam
from convert import runPhenixProgram, getProgram, SUPERPOSE
from pyworkflow.em import PdbFile
import collections


class PhenixProtRunSuperposePDBs(EMProtocol):
    """Superpose two PDBs so that they optimally match """
    _label = 'superpose pdbs'
    _program = ""
    # _version = VERSION_1_2

    # --------------------------- DEFINE param functions -------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputStructureFixed', PointerParam,
                      pointerClass="PdbFile",
                      label='Fixed atomic structure', important=True,
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
        args += " "
        args += os.path.abspath(self.inputStructureMoving.get().getFileName())

        runPhenixProgram(getProgram(SUPERPOSE), args,
                         cwd=self._getExtraPath())

    def createOutputStep(self):
        fnPdb = os.path.basename(self.inputStructureMoving.get().getFileName())
        pdb = PdbFile()
        pdb.setFileName(self._getExtraPath(fnPdb+"_fitted.pdb"))
        if self.inputStructureFixed.get().getVolume() is not None:
            pdb.setVolume(self.inputStructureFixed.get().getVolume())
        self._defineOutputs(outputPdb=pdb)
        self._defineSourceRelation(self.inputStructureFixed.get(), pdb)
        self._defineSourceRelation(self.inputStructureMoving.get(), pdb)

        logFile = os.path.abspath(self._getLogsPath()) + "/run.stdout"
        self._parseLogFile(logFile)
        self._store()

    # --------------------------- INFO functions ---------------------------

    def _validate(self):
        errors = []
        # Check that the program exists
        program = getProgram(SUPERPOSE)
        if not os.path.exists(program):
            errors.append("Cannot find " + program)

        # If there is any error at this point it is related to config variables
            errors.append("Check configuration file: "
                          "~/.config/scipion/scipion.conf")
            errors.append("and set PHENIX_HOME variables properly.")
            if program is not None:
                errors.append("Current values:")
                errors.append("PHENIX_HOME = %s" % os.environ['PHENIX_HOME'])
                errors.append("SUPERPOSE = %s" % SUPERPOSE)

        return errors

    def _summary(self):
        summary = []
        try:
            summary.append("RMSD between fixed and moving atoms (start): " +
                       str(self.startRMSD))
            summary.append("RMSD between fixed and moving atoms (final): " +
                       str(self.finalRMSD))
        except:
            summary.append("RMSD not yet computed" )
        summary.append(
            "http://www.phenix-online.org/documentation/superpose_pdbs.htm")
        summary.append("Peter Zwart, Pavel Afonine, Ralf W. Grosse-Kunstleve")

        return summary

    # --------------------------- UTILS functions --------------------------

    def _parseLogFile(self, logFile):
        with open(logFile) as f:
            line = f.readline()
            while line:
                words = line.strip().split()
                if len(words) > 1:
                    if (words[0] == 'RMSD' and words[1] == 'between'and
                        words[6] == '(start):'):
                        self.startRMSD = Float(words[7])
                    elif (words[0] == 'RMSD' and words[1] == 'between' and
                          words[6] == '(final):'):
                        self.finalRMSD = Float(words[7])
                line = f.readline()
