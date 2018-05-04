# **************************************************************************
# *
# * Authors:     Marta Martinez (mmmtnez@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

from protocol_base import createScriptFile
from pyworkflow.protocol.params import PointerParam, StringParam
import os
from pyworkflow.em.viewers.chimera_utils import \
    runChimeraProgram,\
    chimeraPdbTemplateFileName, chimeraMapTemplateFileName, \
    chimeraScriptFileName, getProgram, sessionFile
from pyworkflow.em.packages.chimera.protocol_base import ChimeraProtBase

class ChimeraProtRestore(ChimeraProtBase):
    """Protocol to perform rigid fit using Chimera.
        Execute command *scipionwrite [model #n] [refmodel #p]
        [saverefmodel 0|1]* from command line in order to transfer fitted
        pdb to scipion. Default values are model=#0,
        refmodel =#1 and saverefmodel 0 (false).
        model refers to the pdb file. refmodel to a 3Dmap"""
    _label = 'chimera restore session'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputProtocol', PointerParam,
                  label="Input protocols", important=True,
                  pointerClass='ChimeraProtOperate, ChimeraProtRigidFit',
                  help="protocol to be reloaded")

        form.addParam('extraCommands', StringParam,
                      default='',
                      condition='False',
                      label='Extra commands for chimera viewer',
                      help="Add extra commands in cmd file. Use for testing")
        form.addSection(label='Help')
        form.addLine('''Execute command *scipionwrite [model #n] [refmodel #p]
        [saverefmodel 0|1]* from command line in order to transfer structures
        and 3D map volumes to SCIPION.
        In the particular case in which you have only a volume and a structure,
        default values are model #2, refmodel #1 and saverefmodel 0 (false).
        Model refers to the PDBx/mmCIF file, refmodel to a 3D map volume.
        If you have several structures and no volumes, you can save
        all of them by executing commands *scipionwrite [model #1]*,
        *scipionwrite [model #2]*, *scipionwrite [model #3]*, and so on.''')

    # --------------------------- INSERT steps functions --------------------

    def prerequisitesStep(self):
        """
        """
        self.parentProt = self.inputProtocol.get()
        self.parentProt.setProject(self.getProject()) # I do not really
        # understand this line

        self.inputVolume = self.parentProt.inputVolume
        self.pdbFileToBeRefined = self.parentProt.pdbFileToBeRefined
        self.inputPdbFiles = self.parentProt.inputPdbFiles

    def runChimeraStep(self):
        #create CMD file
        parentSessionFileName = self.parentProt._getExtraPath(sessionFile)
        sessionFileName = self._getExtraPath(sessionFile)
        f = open(self._getTmpPath(chimeraScriptFileName), "w")
        f.write("from chimera import runCommand\n")
        # create coherent header
        createScriptFile(1,  # model id pdb
                         1,  # model id 3D map
                         self._getExtraPath(chimeraPdbTemplateFileName),
                         self._getExtraPath(chimeraMapTemplateFileName),
                         f,
                         sessionFileName
                         )

        f.write("restoreSession('%s')\n"%parentSessionFileName)

        if len(self.extraCommands.get()) > 2:
            f.write(self.extraCommands.get())
            args = " --nogui --script " + self._getTmpPath(
                chimeraScriptFileName)
        else:
            args = " --script " + self._getTmpPath(chimeraScriptFileName)
        f.close()
        program = getProgram()
        self._log.info('Launching: ' + program + ' ' + args)

        # run in the background
        runChimeraProgram(program, args)

    def createOutput(self):
        super(ChimeraProtRestore, self).createOutput()

    def _validate(self):
        errors = super(ChimeraProtRestore, self)._validate()
        parentProt = self.inputProtocol.get()
        parentProt.setProject(self.getProject()) # I do not really understand
                                                 # this line
        sessionFileName = parentProt._getExtraPath(sessionFile)
        # Check SESSION.py exists
        if (not os.path.exists(sessionFileName)):
            errors.append("Error: No session saved by protocol: %s\n"
                          % parentProt.getObjLabel())

        return errors
