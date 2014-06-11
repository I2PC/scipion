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
# *  e-mail address 'jgomez@cnb.csic.es'
# *
# **************************************************************************

import os

from pyworkflow.em import *  
from pyworkflow.utils import * 
import eman2
from data import *
from convert import readSetOfCoordinates
from pyworkflow.gui.dialog import askYesNo


class EmanProtBoxing(ProtParticlePicking):
    """Protocol to pick particles in a set of micrographs using eman"""
    _label = 'boxer'
        
    def __init__(self, **args):     
        ProtParticlePicking.__init__(self, **args)
        # The following attribute is only for testing
        self.importFolder = String(args.get('importFolder', None))
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputMicrographs', PointerParam, label="Micrographs",
                      pointerClass='SetOfMicrographs',
                      help='Select the SetOfMicrograph ')
        form.addParam('boxSize', IntParam, label='Box size',
                      validators=[Positive])
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self.inputMics = self.inputMicrographs.get()
        micList = [os.path.relpath(mic.getFileName(), self.workingDir.get()) for mic in self.inputMics]

        self._params = {'inputMics': ' '.join(micList), 
                        'boxSize': self.boxSize.get()}      
        # Launch Boxing GUI
        if not self.importFolder.hasValue():
            self._insertFunctionStep('launchBoxingGUIStep', isInteractive=True)
        else: # This is only used for test purposes
            self._insertFunctionStep('_importFromFolderStep')  
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def launchBoxingGUIStep(self):
        # First we go to runs directory (we create if it does not exist)
        #path.makePath("runs")
        # Program to execute and it arguments
        program = eman2.getEmanProgram("e2boxer.py")
        arguments = "%(inputMics)s --boxsize=%(boxSize)i"
        # Run the command with formatted parameters
        self._log.info('Launching: ' + program + ' ' + arguments % self._params)
        self.runJob(program, arguments % self._params)
        
        # Open dialog to request confirmation to create output
        if askYesNo(Message.TITLE_SAVE_OUTPUT, Message.LABEL_SAVE_OUTPUT, None):
            self._createOutput()
        
    def _createOutput(self):
            workDir = self.workingDir.get()
            
            # As we move to workingDir we must leave it. 
            self._leaveWorkingDir()
            micSet = self.inputMics
            
            count = 0;
            for key, output in self.iterOutputAttributes(EMObject):
                count += 1
            suffix = str(count + 1) if count > 0 else ''
            outputName = 'outputCoordinates' + suffix
            
            coordSet = self._createSetOfCoordinates(self.inputMics, suffix)
            readSetOfCoordinates(workDir, self.inputMics, coordSet)
            outputs = {outputName: coordSet}
            self._defineOutputs(**outputs)
            self._defineSourceRelation(micSet, coordSet)
    
    def _importFromFolderStep(self):
        """ This function will copy Xmipp .pos files for
        simulating an particle picking run...this is only
        for testing purposes.
        """
        from pyworkflow.utils.path import copyTree

        print "COPYTREE from %s TO %s" % (self.importFolder.get(), os.getcwd())
        copyTree(self.importFolder.get(), os.getcwd())
        
        self._createOutput()
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _runSteps(self, startIndex):
        # Redefine run to change to workingDir path
        # Change to protocol working directory
        self._enterWorkingDir()
        eman2.loadEnvironment()
        Protocol._runSteps(self, startIndex)
    
    def getFiles(self):
        filePaths = self.inputMicrographs.get().getFiles() | ProtParticlePicking.getFiles(self)
        return filePaths
