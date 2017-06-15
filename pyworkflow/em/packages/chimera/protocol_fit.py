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
from convert import (getProgram, runChimeraProgram,
                     chimeraPdbTemplateFileName,
                     chimeraMapTemplateFileName,
                     chimeraScriptFileName)
from pyworkflow.em import PdbFile
from pyworkflow.em import Volume

#todo add function scipion_change_scale
#TODO:tests coherence betweem metadata and header
#convert to mrc if needed

class ChimeraProtRigidFit(EMProtocol):
    """Protocol to perform rigid fit using Chimera.
        Execute command *scipionwrite [model #n] [refmodel #p] [saverefmodel=0|1]* from command
        line in order to transferm fitted pdb to scipion. Default values are model=#0,
        refmodel =#1 and saverefmodel=0 (false).
        model refers to the pdb file. refmodel to a 3Dmap"""
    _label = 'Rigid Fit'
    _version = VERSION_1_2

    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVolume', PointerParam, pointerClass="Volume",
                      label='Input Volume',
                      help="Volume to process")
        form.addParam('pdbFileToBeRefined', PointerParam, pointerClass="PdbFile",
                      label='PDB to be refined',
                      help="PDB file to be refined. This PDB object, after refinement, will be saved")
        form.addParam('inputPdbFiles', MultiPointerParam, pointerClass="PdbFile",
                      label='Other referece PDBs',
                      help="Other PDB files used as reference. This PDB objects will not be saved")
        form.addSection(label='Help')
        form.addLine('''Execute command *scipionwrite [model #n] [refmodel #p] [saverefmodel=0|1]* from command
        line in order to transferm fitted pdb to scipion. Default values are model=#0,
        refmodel =#1 and saverefmodel=0 (false).
        model refers to the pdb file. refmodel to a 3Dmap''')
        # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):

        self._insertFunctionStep('prerequisitesStep')
        self._insertFunctionStep('runChimeraStep')
        self._insertFunctionStep('createOutput')

    # --------------------------- STEPS functions ---------------------------------------------------

    def prerequisitesStep(self):
        """
        """
        #convert to mrc
        #create coherent header
        createScriptFile(0,  #model id pdb
                         1,  #model id 3D map
                         self._getTmpPath(chimeraScriptFileName),
                         self._getExtraPath(chimeraPdbTemplateFileName),
                         self._getExtraPath(chimeraMapTemplateFileName),
                         )

    def runChimeraStep(self):
        args = ""
        args +=  " --script " + self._getTmpPath(chimeraScriptFileName) + " "
        args += self.pdbFileToBeRefined.get().getFileName() + " "
        args += self.inputVolume.get().getFileName() + " "
        for pdb in self.inputPdbFiles:
            args += " " + pdb.get().getFileName() # other pdb files
        self._log.info('Launching: ' + getProgram() + ' ' + args)

        #run in the background
        runChimeraProgram(getProgram(), args)

    def createOutput(self):
        """ Copy the PDB structure and register the output object.
        """
        pdb = PdbFile()
        pdb.setFileName(self._getExtraPath(chimeraPdbTemplateFileName)%1)
        self._defineOutputs(outputPdb=pdb)
        self._defineSourceRelation(self.inputPdbFiles, pdb)
        self._defineSourceRelation(self.inputVolume, pdb)

        volFileName = self._getExtraPath(chimeraMapTemplateFileName%1)
        if os.path.exists(volFileName):
            vol = Volume()
            vol.setLocation(volFileName)
            #TODO: this sampling rate needs to be modified
            #HORROR: this way of changing sampling is wrong
            #develope a central one
            #the only nice thing is that I know that the file is MRC
            #since it is output and I did create it
            import struct
            f = open(volFileName,'rb')
            s = f.read(13*4)#read header up to angles incluse word 6
            f.close()
            chain = "< 3i i 3i 3i 3f"
            a = struct.unpack(chain, s)
            sampling = a[12]/a[9]
            #end of the HORROR
            #vol.setSamplingRate(self.inputVolume.get().getSamplingRate())
            vol.setSamplingRate(sampling)
            self._defineOutputs(output3Dmap=vol)
            self._defineSourceRelation(self.inputPdbFiles, vol)
            self._defineSourceRelation(self.inputVolume, vol)


    # --------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        # Check that the program exists
        program = getProgram()
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
        summary.append("TO BE DONE")
        if (self.getOutputsSize() > 0):
            summary.append("we have some result")
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

#TODO: add change scale, scipio_write_pdb, scipion_write_map
chimeraScriptHeader= '''
import os
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
scipionwrite model #n refmodel #p saverefmodel 0/1
"""
def cmd_scipionWrite(scipionWrite,args):
  from Midas.midas_text import doExtensionFunc
'''

chimeraScriptMain = '''
  def scipionWrite(model="#%(pdbID)d",refmodel="#%(_3DmapId)d", saverefmodel=0):
     #get model (pdb) id
     modelId = int(model[1:])# model to write
     refModelId = int(refmodel[1:])# coordenate system refers to this model

     # get actual models
     model    = chimera.openModels.list()[modelId]
     refModel = chimera.openModels.list()[refModelId]

     # Save the PDB relative to the volume coordinate system
     # TODO: check if this Will work if the reference is a PDB?
     from Midas import write
     write(model, refModel, "%(pdbFileTemplate)s")
     # alternative way to save  the pdb file using a command
     #run('write relative #1 #0 pdb_path')

     # Save the map if sampling rate has been changed
     if saverefmodel:
         from VolumeData import save_grid_data
         save_grid_data(refModel.data, "%(chimeraMapTemplateFileName)s")
     beep(0.1)

  doExtensionFunc(scipionWrite,args)

from Midas.midas_text import addCommand
addCommand('scipionwrite', cmd_scipionWrite, help="http://scipion.cnb.csic.es")
'''

def createScriptFile(pdbID, _3DmapId, scriptFile, pdbFileTemplate, mapFileTemplate):
    f = open(scriptFile,"w")
    f.write(chimeraScriptHeader)
    d={}
    d['pdbID'] = pdbID
    d['_3DmapId'] = _3DmapId
    d['pdbFileTemplate'] = pdbFileTemplate%1
    d['chimeraMapTemplateFileName'] = mapFileTemplate%1

    f.write(chimeraScriptMain % d)
    f.close()

