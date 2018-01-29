# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (jgomez@cnb.csic.es)
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

import os

from pyworkflow.object import String
from pyworkflow.utils.properties import Message
from pyworkflow.utils.path import join, getExt
from pyworkflow.gui.dialog import askYesNo
from pyworkflow.em.protocol import ProtParticlePicking, BooleanParam

import eman2
from pyworkflow.em.packages.eman2.convert import loadJson
from convert import readSetOfCoordinates


class EmanProtBoxing(ProtParticlePicking):
    """ Picks particles in a set of micrographs using eman2 boxer. """
    _label = 'boxer'
        
    def __init__(self, **args):     
        ProtParticlePicking.__init__(self, **args)
        # The following attribute is only for testing
        self.importFolder = String(args.get('importFolder', None))

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self.inputMics = self.inputMicrographs.get()
        micList = [os.path.relpath(mic.getFileName(), self.workingDir.get()) for mic in self.inputMics]

        self._params = {'inputMics': ' '.join(micList)}
        # Launch Boxing GUI
        self._insertFunctionStep('launchBoxingGUIStep', interactive=True)

    def _defineParams(self, form):
        ProtParticlePicking._defineParams(self, form)

        form.addParam('invertY', BooleanParam, default=False,
                      label='Invert Y coordinates',
                      help='In some cases, using dm3 or tiff Y coordinates '
                           'must be flipped. Check output and activate this'
                           ' if needed.')

    #--------------------------- STEPS functions ---------------------------------------------------
    def launchBoxingGUIStep(self):
        # Print the eman version, useful to report bugs
        self.runJob(eman2.getEmanProgram('e2version.py'), '')
        # Program to execute and it arguments
        program = eman2.getEmanProgram("e2boxer.py")
        arguments = "%(inputMics)s"
        # Run the command with formatted parameters
        self._log.info('Launching: ' + program + ' ' + arguments % self._params)
        self.runJob(program, arguments % self._params)

        # Open dialog to request confirmation to create output
        if askYesNo(Message.TITLE_SAVE_OUTPUT, Message.LABEL_SAVE_OUTPUT, None):
            self.check_gauss()
            self._leaveDir()# going back to project dir
            self._createOutput(self.getWorkingDir())

    def check_gauss(self):
        # Function to check if gaussian algorithm was used to pick and is so
        # ask user if she wants to perform an automatic picking for the remaining micrographs
        gaussJsonFile = join("e2boxercache", "gauss_box_DB.json")
        # Check if gauss json file exists and load it
        if os.path.exists(gaussJsonFile):
            jsonGaussDict = loadJson(gaussJsonFile)
            gaussParsDict = None
            micList = [os.path.relpath(mic.getFileName(), self.workingDir.get()) for mic in self.inputMics]
            # Go over the list of input micrographs and see if gaussian was used to pick any of them
            for mic in micList:
                if mic in jsonGaussDict:
                    gaussParsDict = jsonGaussDict[mic]
                    break
            if gaussParsDict is not None:
                # If found ask user if she wats to perform an automatic gaussian picking for the rest of mics
                # if askYesNo(Message.TITLE_PICK_GAUSS, Message.LABEL_PICK_GAUSS, None):
                self._params['boxSize'] = gaussParsDict['boxsize']
                # Run sxprocess.py to store parameters
                program = eman2.getEmanProgram("sxprocess.py")
                argsList = ["'%s'=%s:" %(key, val) for (key, val) in gaussParsDict.iteritems()]
                args = 'demoparms --makedb ' + "".join(argsList)
                # Remove last ":" to avoid error
                args = args[:-1]
                # Run the command with formatted parameters
                self._log.info('Launching: ' + program + ' ' + args)
                self.runJob(program, args)
                # Now run e2boxer.py with stored parameters
                #arguments = "--gauss_autoboxer=demoparms --write_ptcl --boxsize=%(boxSize)s --norm=normalize.ramp.normvar" + arguments
                arguments = "--gauss_autoboxer=demoparms --write_dbbox --boxsize=%(boxSize)s " + "%(inputMics)s"
                program = eman2.getEmanProgram("e2boxer.py")
                self._log.info('Launching: ' + program + ' ' + arguments % self._params)
                self.runJob(program, arguments % self._params)

    #--------------------------- INFO functions ---------------------------------------------------
    def _validate(self):
        errors = []
        eman2.validateVersion(self, errors)
        return errors

    def _warnings(self):
        warnings = []
        firstMic = self.inputMicrographs.get().getFirstItem()
        fnLower = firstMic.getFileName().lower()

        ext = getExt(fnLower)

        if ext in ['.tif', '.dm3'] and not self.invertY.get():
            warnings.append(
                'We have seen a flip in Y when using %s files in EMAN2' % ext)
            warnings.append(
                'The generated coordinates may or may not be valid in Scipion.')
            warnings.append(
                'TIP: Activate "Invert Y coordinates" if you find it wrong.')
        return warnings

    #--------------------------- UTILS functions ---------------------------------------------------
    def _runSteps(self, startIndex):
        # Redefine run to change to workingDir path
        # Change to protocol working directory
        self._enterWorkingDir()
        ProtParticlePicking._runSteps(self, startIndex)
    
    def getFiles(self):
        filePaths = self.inputMicrographs.get().getFiles() | ProtParticlePicking.getFiles(self)
        return filePaths

    def readSetOfCoordinates(self, workingDir, coordSet):
        readSetOfCoordinates(workingDir, self.inputMics, coordSet, self.invertY.get())
