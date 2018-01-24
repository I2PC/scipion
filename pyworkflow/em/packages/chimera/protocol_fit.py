# **************************************************************************
# *
# * Authors:     Grigory Sharov (sharov@igbmc.fr)
# *              Marta Martinez (mmmtnez@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
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
from pyworkflow import VERSION_1_2
from pyworkflow.em import PdbFile
from pyworkflow.em import Volume
from pyworkflow.em.convert import ImageHandler
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.em.utils.ccp4_utilities.convert import Ccp4Header
from pyworkflow.em.utils.chimera_utilities.convert import \
    createCoordinateAxisFile, \
    adaptOriginFromCCP4ToChimera, getProgram, runChimeraProgram,\
    chimeraPdbTemplateFileName, chimeraMapTemplateFileName, \
    chimeraScriptFileName
from pyworkflow.protocol.params import MultiPointerParam, PointerParam, \
    StringParam
from pyworkflow.utils.properties import Message


class ChimeraProtRigidFit(EMProtocol):
    """Protocol to perform rigid fit using Chimera.
        Execute command *scipionwrite [model #n] [refmodel #p]
        [saverefmodel 0|1]* from command line in order to transferm fitted
        pdb to scipion. Default values are model=#0,
        refmodel =#1 and saverefmodel 0 (false).
        model refers to the pdb file. refmodel to a 3Dmap"""
    _label = 'chimera rigid fit'
    _version = VERSION_1_2

    # --------------------------- DEFINE param functions --------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVolume', PointerParam, pointerClass="Volume",
                      label='Input Volume', allowsNull=True,
                      help="Volume to process")
        form.addParam('pdbFileToBeRefined', PointerParam,
                      pointerClass="PdbFile",
                      label='PDBx/mmCIF file to be refined',
                      help="PDBx/mmCIF file to be refined. This cif object, "
                           "after refinement, will be saved")
        form.addParam('inputPdbFiles', MultiPointerParam,
                      pointerClass="PdbFile",
                      label='Other referece PDBx/mmCIF files',
                      help="Other PDBx/mmCIF files used as reference. "
                           "This cif objects will not be saved")
        form.addParam('extraCommands', StringParam,
                      default='',
                      condition='False',
                      label='Extra commands for chimera viewer',
                      help="""Add extra commands in cmd file. Use for testing
                      """)
        form.addSection(label='Help')
        form.addLine('''Execute command *scipionwrite [model #n] [refmodel
        #p] [saverefmodel 0|1]* from command
        line in order to transfer fitted cif to scipion. Default values are
        model=#2,
        refmodel =#1 and saverefmodel 0 (false).
        model refers to the cif file. refmodel to a 3Dmap''')

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):

        self._insertFunctionStep('prerequisitesStep')
        self._insertFunctionStep('runChimeraStep')
        self._insertFunctionStep('createOutput')

    # --------------------------- STEPS functions ---------------------------

    def prerequisitesStep(self):
        """
        """
        if self.inputVolume.get() is None:
            fnVol = self.pdbFileToBeRefined.get().getVolume()
            print fnVol, type(fnVol)
            index, fn = fnVol.getLocation()
            print "Volume: Volume associated to atomic structure %s(%d)\n" \
                  % (fn, index)
        else:
            fnVol = self.inputVolume.get()
            print "Volume: Input volume %s\n" % fnVol

    def runChimeraStep(self):
        # building script file including the coordinate axes and the input
        # volume with samplingRate and Origin information
        f = open(self._getTmpPath(chimeraScriptFileName), "w")
        f.write("from chimera import runCommand\n")
        # create coherent header
        createScriptFile(0,  # model id pdb
                         1,  # model id 3D map
                         self._getExtraPath(chimeraPdbTemplateFileName),
                         self._getExtraPath(chimeraMapTemplateFileName),
                         f
                         )
        if self.inputVolume.get() is None:
            _inputVol = self.pdbFileToBeRefined.get().getVolume()
        else:
            _inputVol = self.inputVolume.get()

        # building coordinate axes
        dim = _inputVol.getDim()[0]
        sampling = _inputVol.getSamplingRate()
        tmpFileName = os.path.abspath(self._getTmpPath("axis_input.bild"))
        createCoordinateAxisFile(dim,
                                 bildFileName=tmpFileName,
                                 sampling=sampling)
        f.write("runCommand('open %s')\n" % tmpFileName)

        # input vol with its origin coordinates
        x_input, y_input, z_input = adaptOriginFromCCP4ToChimera(
            _inputVol.getVolOriginAsTuple())
        inputVolFileName = os.path.abspath(ImageHandler.removeFileType(
            _inputVol.getFileName()))
        f.write("runCommand('open %s')\n" % inputVolFileName)
        f.write("runCommand('volume #1 style surface voxelSize %f')\n"
                % _inputVol.getSamplingRate())
        f.write("runCommand('volume #1 origin %0.2f,%0.2f,%0.2f')\n"
                % (x_input, y_input, z_input))
        #f.write("runCommand('volume #1 style surface voxelSize %f origin "
        #        "%0.2f,%0.2f,%0.2f')\n"
        #        % (_inputVol.getSamplingRate(), x_input, y_input, z_input))
        pdbFileToBeRefined = self.pdbFileToBeRefined.get()
        f.write("runCommand('open %s')\n" % os.path.abspath(
            pdbFileToBeRefined.getFileName()))
        if pdbFileToBeRefined.hasOrigin():
            #x, y, z = adaptOriginFromCCP4ToChimera(
            #    pdbFileToBeRefined.getOrigin().getShifts())
            x, y, z = (pdbFileToBeRefined.getOrigin().getShifts())
            print "x, y, z: ", x, y, z
            f.write("runCommand('move %0.2f,%0.2f,%0.2f model #%d coord #0')\n"
                % (x, y, z, 2))
        # other pdb files
        pdbModelCounter=3
        for pdb in self.inputPdbFiles:
            f.write("runCommand('open %s')\n" % os.path.abspath(pdb.get(
            ).getFileName()))
            if pdb.get().hasOrigin():
                # x, y, z = adaptOriginFromCCP4ToChimera(
                #     pdb.get().getOrigin().getShifts())
                x, y, z = pdb.get().getOrigin().getShifts()
                f.write("runCommand('move %0.2f,%0.2f,%0.2f model #%d coord #0')\n"
                        % (x,y,z,pdbModelCounter))

        # run the text:
        if len(self.extraCommands.get()) > 2:
            f.write(self.extraCommands.get())
            args = " --nogui --script " + self._getTmpPath(
                chimeraScriptFileName)
        else:
            args = " --script " + self._getTmpPath(chimeraScriptFileName)

        f.close()

        self._log.info('Launching: ' + getProgram() + ' ' + args)

        # run in the background
        runChimeraProgram(getProgram(), args)

    def createOutput(self):
        """ Copy the PDB structure and register the output object.
        """
        volFileName = self._getExtraPath((chimeraMapTemplateFileName) % 1)
        if os.path.exists(volFileName):
            vol = Volume()
            vol.setLocation(volFileName)

            ccp4header = Ccp4Header(volFileName, readHeader=True)
            sampling = ccp4header.computeSampling()
            vol.setSamplingRate(sampling)
            if self.inputVolume.get() is None:
                origin = self.pdbFileToBeRefined.get().getVolume(). \
                    getOrigin(returnInitIfNone=True)
            else:
                origin = self.inputVolume.get().getOrigin(
                    returnInitIfNone=True)
            vol.setOrigin(origin)
            print "origin.getShifts: ", origin.getShifts()#origin.getShifts:  (37.5, 37.5, 37.5)
            print "ccp4header.getStartAngstrom(sampling): ", ccp4header.getStartAngstrom(
                sampling)#ccp4header.getStartAngstrom(sampling):  (0.0, 0.0, 0.0)

            # ##DELETE THIS
            # ccp4header.setStartAngstrom(origin.getShifts(), sampling)
            # # ccp4header.writeHeader()
            # ###ccp4header.getStartAngstrom(sampling)
            # print "ccp4header.getStartAngstrom(sampling): ", \
            #     ccp4header.getStartAngstrom(sampling)#ccp4header.getStartAngstrom(sampling):  (37.5, 37.5, 37.5)
            # ##

            self._defineOutputs(output3Dmap=vol)
            self._defineSourceRelation(self.inputPdbFiles, vol)
            if self.inputVolume.get() is None:
                self._defineSourceRelation(
                    self.pdbFileToBeRefined.get().getVolume(), vol)
            else:
                self._defineSourceRelation(self.inputVolume.get(), vol)
        else:
            if self.inputVolume.get() is None:
                vol = self.pdbFileToBeRefined.get().getVolume()
            else:
                vol = self.inputVolume.get()

        pdb = PdbFile()
        pdb.setFileName(self._getExtraPath(chimeraPdbTemplateFileName) % 1)
        pdb.setVolume(vol)
        self._defineOutputs(outputPdb=pdb)
        self._defineSourceRelation(self.inputPdbFiles, pdb)
        if self.inputVolume.get() is None:
            self._defineSourceRelation(
                self.pdbFileToBeRefined.get().getVolume(), pdb)
        else:
            self._defineSourceRelation(self.inputVolume.get(), pdb)

    # --------------------------- INFO functions ----------------------------
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
            errors.append("Check configuration file: ~/.config/scipion/"
                          "scipion.conf")
            errors.append("and set CHIMERA_HOME variables properly.")
            if program is not None:
                errors.append("Current value:")
                errors.append("CHIMERA_HOME = %s" % os.environ['CHIMERA_HOME'])

        # Check that the input volume exist
        if (not self.pdbFileToBeRefined.get().hasVolume()) and (
                    self.inputVolume.get() is None):
            errors.append("Error: You should provide a volume.\n")

        return errors

    def _summary(self):
        # Think on how to update this summary with created PDB
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

    # --------------------------- UTILS functions --------------------------

# TODO: add change scale, scipio_write_pdb, scipion_write_map
chimeraScriptHeader = '''
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

"""load the script in chimera for example chimera --script ScipioChimeraExt/
ChimeraExtension.py
from chimera command line type one of the following three options
scipionwrite
scipionwrite model #n
scipionwrite model #n refmodel #p saverefmodel 0/1
"""
def cmd_scipionWrite(scipionWrite,args):
  from Midas.midas_text import doExtensionFunc
'''

chimeraScriptMain = '''
  def scipionWrite(model="#%(pdbID)d",refmodel="#%(_3DmapId)d",
     saverefmodel=0):
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


def createScriptFile(pdbID, _3DmapId,
                     pdbFileTemplate, mapFileTemplate, f):
    f.write(chimeraScriptHeader)
    d = {}
    d['pdbID'] = pdbID
    d['_3DmapId'] = _3DmapId
    d['pdbFileTemplate'] = pdbFileTemplate % 1
    d['chimeraMapTemplateFileName'] = mapFileTemplate % 1

    f.write(chimeraScriptMain % d)
