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
from pyworkflow.em.data import EMObject
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.protocol.params import PointerParam, FloatParam
from pyworkflow.utils.properties import Message
from pyworkflow.em.convert_header.CCP4.convert import adaptFileToCCP4, START
from convert import runPhenixProgram

class PhenixProtRunMolprobity(EMProtocol):
    """MolProbity is a Phenix application to validate the geometry of an
atomic structure derived from a cryo-EM density map.
"""
    _label = 'molprobity: model validation'
    _program = ""
    #_version = VERSION_1_2
    MOLPROBITY = 'molprobity.py'
    MOLPROBITYFILE = 'molprobity.mrc'

    # --------------------------- DEFINE param functions -------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVolume', PointerParam, pointerClass="Volume",
                      label='Input Volume', allowsNull=True,
                      help="Set the starting volume")
        form.addParam('resolution', FloatParam, allowsNull=False,
                      label='Resolution (A):',
                      help='Set the resolution of the input volume')
        form.addParam('inputStructure', PointerParam,
                      pointerClass="PdbFile", allowsNull=False,
                      label='Input atomic structure',
                      help="Set the PDBx/mmCIF to be validated.")

    # --------------------------- INSERT steps functions --------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('runMolprobityStep')

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
        vol = os.path.abspath(self._getTmpPath(self.MOLPROBITYFILE))
        inFileName = vol.getFileName()
        volume = os.path.abspath(inFileName)
        args += "map_file_name=" + volume
        args += "d_min=" + self.resolution.get()
        print "PROGRAM: ", (self.MOLPROBITY + ' ' + args)

        # script with auxiliary files

        self._log.info('Launching: ' + runPhenixProgram(self.MOLPROBITY + ' ' +
                       args))


    # def createOutputStep(self, inVolumes, norVolumesNames, init_counter=1):
    #     """ Copy the PDB structure and register the output object.
    #     """
    #     template = self._getExtraPath(cootPdbTemplateFileName)
    #     counter = init_counter
    #     counter -= 1
    #     while os.path.isfile(template % counter):
    #         pdb = PdbFile()
    #         pdb.setFileName(template % counter)
    #
    #         outputs = {"outputPdb_%04d" % counter: pdb}
    #         self._defineOutputs(**outputs)
    #
    #         # self._defineOutputs(outputPdb=pdb)
    #         self._defineSourceRelation(self.inputPdbFiles, pdb)
    #         # self._defineSourceRelation(self.inputVolumes, self.outputPdb)
    #
    #         for vol in inVolumes:
    #             self._defineSourceRelation(vol, pdb)
    #         counter += 1
    #
    #     if not os.path.isfile(template % 2):  # only the first time get inside
    #         # here
    #         counter = 1
    #         for inVol, norVolName in zip(inVolumes, norVolumesNames):
    #             outVol = Volume()
    #             sampling = inVol.getSamplingRate()
    #             origin = inVol.getOrigin(
    #                 force=True)
    #             outVol.setSamplingRate(sampling)
    #             outVol.setOrigin(origin)
    #
    #             if norVolName.endswith('.mrc'):
    #                 norVolName = norVolName + ":mrc"
    #             outFileName = self._getVolumeFileName(norVolName)
    #             outVol.setFileName(outFileName)
    #             outputs = {"output3DMap_%04d" % counter: outVol}
    #             counter += 1
    #             self._defineOutputs(**outputs)
    #             self._defineSourceRelation(inVol, outVol)
    #     if os.path.isfile(self._getExtraPath('STOPPROTCOL')):
    #         self.setStatus(STATUS_FINISHED)
    #         # NOTE: (ROB) can a derthy way to make an interactive process finish but I do not
    #         # think there is a clean one
    #         self._steps[self.step-1].setInteractive(False)

    # --------------------------- INFO functions ---------------------------
    def _validate(self):
        errors = []
        # Check that the program exists
        #program = getProgram(self.EMRINGER)
        program = runPhenixProgram(self.MOLPROBITY)
        if program is None:
            errors.append("Missing variables MOLPROBITY and/or PHENIX_HOME")
            print "AAAAAAAAAAAAAAA: ", program, "  ", os.path(program)
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
                errors.append("MOLPROBITY= %s" % self.EMRINGER)

        # Check that the input volume exist
        if (not self.imputPdbFile.get().hasVolume()) \
                and self.inputVolume is None:
            errors.append("Error: You should provide a volume.\n")

        return errors

    def _summary(self):
        #  Think on how to update this summary with created PDB
        summary = []
        if self.getOutputsSize() >= 1:
            for key, output in self.iterOutputAttributes(EMObject):
                summary.append("*%s:* \n %s " % (key, output.getObjComment()))
        else:
            summary.append(Message.TEXT_NO_OUTPUT_CO)
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