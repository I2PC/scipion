
import os
from pyworkflow.object import String, Float, Integer
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.protocol.params import PointerParam, FloatParam
from pyworkflow.em.headers import adaptFileToCCP4, START
from convert import runPhenixProgram, getProgram, MOLPROBITY
from pyworkflow.utils import magentaStr

class PhenixProtRunRefinementBase(EMProtocol):
    """MolProbity is a Phenix application to validate the geometry of an
atomic structure derived from a cryo-EM density map.
"""
    _label = 'refinementBase'
    _program = ""

    MOLPROBITYOUTFILENAME = 'molprobity.out'
    MOLPROBITYCOOTFILENAME = 'molprobity_coot.py'
    MOLPROBITYPKLFILENAME = 'molprobity.pkl'

    # --------------------------- DEFINE param functions -------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVolume', PointerParam, pointerClass="Volume",
                      label='Input Volume', allowsNull=True,
                      help="Set the starting volume.")
        form.addParam('resolution', FloatParam,
                      label='Resolution (A):',
                      default=3.0,
                      help='Set the resolution of the input volume.')
        form.addParam('inputStructure', PointerParam,
                      pointerClass="PdbFile", allowsNull=False,
                      label='Input atomic structure',
                      help="Set the PDBx/mmCIF to be processed.")

    # --------------------------- STEPS functions --------------------------

    def convertInputStep(self, tmpMapFileName):
        """ convert 3D maps to MRC '.mrc' format
        """
        vol = self._getInputVolume()
        inVolName = vol.getFileName()
        newFn = self._getExtraPath(tmpMapFileName)
        origin = vol.getOrigin(force=True).getShifts()
        sampling = vol.getSamplingRate()
        adaptFileToCCP4(inVolName, newFn, origin, sampling, START)  # ORIGIN


    def _summary(self):
        #  Think on how to update this summary with created PDB
        summary = []
        try:
            summary.append("MolProbity statistics:\n")
            summary.append("Ramachandran outliers: %0.2f %%   (Goal: < 0.2%%)  "\
                           "Ramachandran favored: %0.2f %%   (Goal: > 98%%) " %\
                           (self.ramachandranOutliers.get(),
                            self.ramachandranFavored.get())
                           )
            summary.append("Rotamer outliers:           %0.2f"\
                           " %%   (Goal: < 1%%)    "\
                           "C-beta outliers:              %d"\
                           "            (Goal: 0) " % \
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

    def validateBase(self, program, label):
        errors = []
        # Check that the program exists
        program = getProgram(program)
        if program is None:
            errors.append("Missing variables %s and/or PHENIX_HOME" % label)
        elif not os.path.exists(program):
            errors.append("Binary '%s' does not exists.\n" % program)

        # If there is any error at this point it is related to config variables
        if errors:
            errors.append("Check configuration file: "
                          "~/.config/scipion/scipion.conf")
            errors.append("and set %s and PHENIX_HOME variables "
                          "properly."% label)
            if program is not None:
                errors.append("Current values:")
                errors.append("PHENIX_HOME = %s" % os.environ['PHENIX_HOME'])
                errors.append("%s = %s" % label, program)

        return errors

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
