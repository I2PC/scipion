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
from pyworkflow.object import String, Float, Integer
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.protocol.params import PointerParam, FloatParam
from pyworkflow.em.convert_header.CCP4.convert import adaptFileToCCP4, START
from convert import runPhenixProgram, getProgram

class PhenixProtRunMolprobity(EMProtocol):
    """MolProbity is a Phenix application to validate the geometry of an
atomic structure derived from a cryo-EM density map.
"""
    _label = 'molprobity: model validation'
    _program = ""
    #_version = VERSION_1_2
    MOLPROBITY = 'molprobity.py'
    MOLPROBITYFILE = 'molprobity.mrc'
    MOLPROBITYOUTFILENAME = 'molprobity.out'
    MOLPROBITYCOOTFILENAME = 'molprobity_coot.py'

    # --------------------------- DEFINE param functions -------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVolume', PointerParam, pointerClass="Volume",
                      label='Input Volume', allowsNull=True,
                      help="Set the starting volume (optional).")
        form.addParam('resolution', FloatParam,
                      label='Resolution (A):', allowsNull=True,
                      help='Set the resolution of the input volume.')
        form.addParam('inputStructure', PointerParam,
                      pointerClass="PdbFile", allowsNull=False,
                      label='Input atomic structure',
                      help="Set the PDBx/mmCIF to be validated.")

    # --------------------------- INSERT steps functions --------------------

    def _insertAllSteps(self):
        if self.inputVolume.get() is not None:
            self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('runMolprobityStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------

    def convertInputStep(self):
        """ convert 3D maps to MRC '.mrc' format
        """
        vol = self._getInputVolume()
        inVolName = vol.getFileName()
        newFn = self._getTmpPath(self.MOLPROBITYFILE)
        origin = vol.getOrigin(force=True).getShifts()
        sampling = vol.getSamplingRate()
        adaptFileToCCP4(inVolName, newFn, origin, sampling, START)  # ORIGIN

    def runMolprobityStep(self):

        # PDBx/mmCIF
        pdb = os.path.abspath(self.inputStructure.get().getFileName())
        args = ""
        args += pdb
        # starting volume (.mrc)
        if self.inputVolume.get() is not None:
            args += " "
            volume = os.path.abspath(self._getTmpPath(self.MOLPROBITYFILE))
            args += "map_file_name=%s" % volume
            args += " "
            args += "d_min=%f" % self.resolution.get()

        # script with auxiliary files

        runPhenixProgram(getProgram(self.MOLPROBITY), args,
                         cwd=self._getExtraPath())


    def createOutputStep(self):
        MOLPROBITYOUTFILENAME = self._getExtraPath(
            self.MOLPROBITYOUTFILENAME)
        self._parseFile(MOLPROBITYOUTFILENAME)
        self._store()

    # --------------------------- INFO functions ---------------------------
    def _validate(self):
        errors = []
        # Check that the program exists
        program = getProgram(self.MOLPROBITY)
        if program is None:
            errors.append("Missing variables MOLPROBITY and/or PHENIX_HOME")
        elif not os.path.exists(program):
            errors.append("Binary '%s' does not exists.\n" % program)

        # If there is any error at this point it is related to config variables
        if errors:
            errors.append("Check configuration file: "
                          "~/.config/scipion/scipion.conf")
            errors.append("and set MOLPROBITY and PHENIX_HOME variables "
                          "properly.")
            if program is not None:
                errors.append("Current values:")
                errors.append("PHENIX_HOME = %s" % os.environ['PHENIX_HOME'])
                errors.append("MOLPROBITY = %s" % self.MOLPROBITY)

        # Check that the input volume exist
        # if self._getInputVolume() is None:
        #     errors.append("Error: You should provide a volume.\n")

        return errors

    def _summary(self):
        #  Think on how to update this summary with created PDB
        summary = []
        try:
            summary.append("Ramachandran outliers: %0.2f %%   (Goal: < 0.2%%)  "\
                           "Ramachandran favored: %0.2f %%   (Goal: > 98%%) " %\
                           (self.ramachandranOutliers.get(),
                            self.ramachandranFavored.get())
                           )
            summary.append("Rotamer outliers:           %0.2f"\
                           " %%   (Goal: < 1%%)    "\
                           "C-beta outliers:              %d"\
                           " %%         (Goal: 0%%) " % \
                            (self.rotamerOutliers.get(),
                             self.cbetaOutliers.get()
                            ))
            summary.append("Clashscore:                   %0.2f"\
                           "                             Overall "\
                           "score:                %0.2f" %\
                           (self.clashscore.get(), self.overallScore.get())
                           )
        except:
            summary = ["Overall score not yet computed"]

        return summary

    def _methods(self):
        methodsMsgs = []
        methodsMsgs.append("TODO")
        return methodsMsgs

    def _citations(self):
        return ['Chen_2010']

    # --------------------------- UTILS functions --------------------------

    def _getInputVolume(self):
        if self.inputVolume.get() is None:
            fnVol = self.inputStructure.get().getVolume()
        else:
            fnVol = self.inputVolume.get()
        return fnVol

    def _parseFile(self, fileName):
        with open(fileName) as f:
            line = f.readline()
            while line:
                words = line.strip().split()
                if len(words) > 1:
                    if (words[0] == 'Ramachandran' and words[1] == 'outliers'):
                        self.ramachandranOutliers = Float(words[3])
                    elif (words[0] == 'favored' and words[1] == '='):
                        self.ramachandranFavored = Float(words[2])
                    elif (words[0] == 'Rotamer' and words[1] == 'outliers'):
                        self.rotamerOutliers = Float(words[3])
                    elif (words[0] == 'C-beta' and words[1] == 'deviations'):
                        self.cbetaOutliers = Integer(words[3])
                    elif (words[0] == 'Clashscore' and words[1] == '='):
                        self.clashscore = Float(words[2])
                    elif (words[0] == 'MolProbity' and words[1] == 'score'):
                        self.overallScore = Float(words[3])
                line = f.readline()

