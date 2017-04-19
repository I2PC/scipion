# **************************************************************************
# *
# * Authors:     Grigory Sharov (sharov@igbmc.fr)
# *
# * L'Institut de genetique et de biologie moleculaire et cellulaire (IGBMC)
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

import pyworkflow.utils as pwutils
from pyworkflow import VERSION_1_2
from pyworkflow.utils.properties import Message
from pyworkflow.protocol.params import MultiPointerParam, PointerParam, BooleanParam
from pyworkflow.em.protocol import EMProtocol
from convert import (getProgram, runCCP4Program,
                     chimeraPdbTemplateFileName, chimeraScriptFileName)
from pyworkflow.em.convert import ImageHandler
from pyworkflow.em import PdbFile



class ChimeraProtRigidFit(EMProtocol):
    """Protocol to perform rigid fit using Chimera.
    Execute command scipion_write from command line in order
    to transferm fitted pdb to scipion
"""
    _label = 'Rigid Fit'
    _version = VERSION_1_2

    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVolumes', PointerParam, pointerClass="Volume",
                      label='Input Volume',
                      help="Volume to process")
        form.addParam('pdbFileToBeRefined', PointerParam, pointerClass="PdbFile",
                      label='PDB to be refined',
                      help="PDB file to be refined. This PDB object, after refinement, will be saved")
        form.addParam('inputPdbFiles', MultiPointerParam, pointerClass="PdbFile",
                      label='Other referece PDBs',
                      help="Other PDB files used as reference. This PDB objects will not be saved")
        form.addSection(label='Help')
        form.addLine('''Execute command scipion_write from command line in order
    to transferm fitted pdb to scipion''')
        # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):

        self._insertFunctionStep('prerequisitesStep')
        self._insertFunctionStep('runChimeraStep')
        self._insertFunctionStep('createOutput')

    # --------------------------- STEPS functions ---------------------------------------------------

    def prerequisitesStep(self):
        """
        """
        createScriptFile(0,  #model id pdb
                         1,  #model id 3D map
                         self._getTmpPath(chimeraScriptFileName),
                         self._getExtraPath(chimeraPdbTemplateFileName))

    def runChimeraStep(self):
        #find last created PDB output file
        template = self._getExtraPath(chimeraPdbTemplateFileName)
        counter=1
        while os.path.isfile(template%counter):
            counter += 1

        #if there is not previous output use pdb file form form
        #otherwise use last created pdb file
        if counter == 1:
            pdbFileToBeRefined = self.pdbFileToBeRefined.get().getFileName()
        else:
            pdbFileToBeRefined = template%(counter-1)
            self._log.info("Using last created PDB file named=%s", pdbFileToBeRefined)
        args = ""
        args +=  " --pdb " + pdbFileToBeRefined
        for pdb in self.inputPdbFiles:
            args += " --pdb " + pdb.get().getFileName() # other pdb files
        args +=  " --script " + self._getTmpPath(chimeraScriptFileName) # script wit auxiliary files
        for vol in self.inputVolumes:
            inFileName  = vol.get().getFileName()
            args += " --map " + self._getVolumeFileName(inFileName)
        #_envDict['COOT_PDB_TEMPLATE_FILE_NAME'] = self._getExtraPath(cootPdbTemplateFileName)
        self._log.info('Launching: ' + getProgram(os.environ['COOT']) + ' ' + args)

        #run in the background
        runCCP4Program(getProgram(os.environ['COOT']), args)

        # while creating files
        #if askYesNo(Message.TITLE_SAVE_OUTPUT, Message.LABEL_SAVE_OUTPUT, None):
        self.createOutput(counter)

    def createOutput(self):
        """ Copy the PDB structure and register the output object.
        """
        template = self._getExtraPath(chimeraPdbTemplateFileName)
        counter=1
        while os.path.isfile(template%counter):
            #counter -=1
            pdb = PdbFile()
            pdb.setFileName(template%counter)

            outputs = {"outputPdb_%04d"%counter: pdb}
            self._defineOutputs(**outputs)

            #self._defineOutputs(outputPdb=pdb)
            self._defineSourceRelation(self.inputPdbFiles, pdb)
            #self._defineSourceRelation(self.inputVolumes, self.outputPdb)
            self._defineSourceRelation(self.inputVolumes, pdb)
            counter += 1

    # --------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        # Check that the program exists
        program = getProgram("chimera")
        if program is None:
            errors.append("Missing variable CHIMERA_HOME")
        elif not os.path.exists(program):
            errors.append("Binary '%s' does not exists.\n" % program)

        # If there is any error at this point it is related to config variables
        if errors:
            errors.append("Check configuration file: ~/.config/scipion/scipion.conf")
            errors.append("and set CHIMERA_HOME variables properly.")
            if program is not None:
                errors.append("Current value:")
                errors.append("CHIMERA_HOME = %s" % os.environ['CHIMERA_HOME'])

        return errors

    def _summary(self):
        #Think on how to update this summary with created PDB
        summary = []
        summary.append("Number of input micrographs: %d"
                       % self.getInputMicrographs().getSize())
        if (self.getOutputsSize() > 0):
            summary.append("Number of particles picked: %d" % self.getCoords().getSize())
            summary.append("Particle size: %d px" % self.getCoords().getBoxSize())
            summary.append("Threshold min: %0.2f" % self.threshold)
        else:
            summary.append(Message.TEXT_NO_OUTPUT_CO)
        return summary

    def _methods(self):
        methodsMsgs = []
        methodsMsgs.append("TODO")

        return methodsMsgs

    def _citations(self):
        return ['Emsley_2004']

    # --------------------------- UTILS functions --------------------------------------------------

    def _getVolumeFileName(self, inFileName):
        return os.path.join(self._getExtraPath(''),
                            pwutils.replaceBaseExt(inFileName, 'mrc'))

chimeraScript= '''
def beep(time):
   """I simply do not know how to create a portable beep sound.
      This system call seems to work pretty well if you have sox
      installed"""
   try:
      command = "play --no-show-progress -n synth %f sin 880"%time
      os.system(command)
   except:
      pass

"""load the script in chimera for example chimera --script ScipioChimeraExt/ChimeraExtension.py
from chimea command line type one of the following three options
scipionwrite
scipionwrite model #n
scipionwrite model #n refmodel #p
"""
def cmd_scipionWrite(scipionWrite,args):
  from Midas.midas_text import doExtensionFunc

  def scipionWrite(model="#%d",refmodel="#%d"):
     #get model (pdb) id
     modelId = int(model[1:])# model to write
     refModelId = int(refmodel[1:])# coordenate system refers to this model

     # get actual models
     model    = chimera.openModels.list()[modelId]
     refModel = chimera.openModels.list()[refModelId]

     #write the model
     #wrapper should incorporate this pdb as scipion object
     from Midas import write
     write(model, refModel, "%s")

  doExtensionFunc(scipionWrite,args)

from Midas.midas_text import addCommand
addCommand('scipionwrite', cmd_scipionWrite, help="http://scipion.cnb.csic.es")
'''

def createScriptFile(pdbID, _3DmapId, scriptFile, pdbFile):
    f = open(scriptFile,"w")
    f.write(chimeraScript% (pdbID, _3DmapId, pdbFile))
    f.close()

