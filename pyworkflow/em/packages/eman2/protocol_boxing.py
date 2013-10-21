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



class EmanProtBoxing(ProtParticlePicking):
    """Protocol to pick particles manually of a set of micrographs in the project"""    
    def __init__(self, **args):     
        ProtParticlePicking.__init__(self, **args)
        # The following attribute is only for testing
        self.importFolder = String(args.get('importFolder', None))
        
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputMicrographs', PointerParam, label="Micrographs",
                      pointerClass='SetOfMicrographs',
                      help='Select the SetOfMicrograph ')
        form.addParam('boxSize', IntParam, label='Box size',
                      validators=[Positive])
        
        
    def _runSteps(self, startIndex):
        # Redefine run to change to workingDir path
        # Change to protocol working directory
        self._enterWorkingDir()
        eman2.loadEnvironment()
        Protocol._runSteps(self, startIndex)
        
    def _defineSteps(self):
        self.inputMics = self.inputMicrographs.get()
        micList = [os.path.relpath(mic.getFileName(), self.workingDir.get()) for mic in self.inputMics]
        self._params = {'inputMics': ' '.join(micList), 
                        'boxSize': self.boxSize.get()}      
        # Launch Boxing GUI
        if not self.importFolder.hasValue():
            self._insertFunctionStep('launchBoxingGUI', isInteractive=True)
        else: # This is only used for test purposes
            self._insertFunctionStep('_importFromFolder')  
        # Insert step to create output objects       
        self._insertFunctionStep('createOutput')
    
    def launchBoxingGUI(self):
        # First we go to runs directory (we create if it does not exist)
        #path.makePath("runs")
        # Program to execute and it arguments
        program = eman2.getEmanProgram("e2boxer.py")
        arguments = "%(inputMics)s --boxsize=%(boxSize)i"
        # Run the command with formatted parameters
        self._log.info('Launching: ' + program + ' ' + arguments % self._params)
        self.runJob(None, program, arguments % self._params)
        
    def createOutput(self):
        workDir = self.workingDir.get()
        # As we move to workingDir we must leave it. 
        self._leaveWorkingDir()
        coordSet = self._createSetOfCoordinates()
        coordSet.setMicrographs(self.inputMics)
        readSetOfCoordinates(workDir, self.inputMics, coordSet)
        coordSet.write()
        self._defineOutputs(outputCoordinates=coordSet)
        
        
    def getFiles(self):
        filePaths = self.inputMicrographs.get().getFiles() | ProtParticlePicking.getFiles(self)
        return filePaths
        
    def _importFromFolder(self):
        """ This function will copy Xmipp .pos files for
        simulating an particle picking run...this is only
        for testing purposes.
        """
        from pyworkflow.utils.path import copyTree

        print "COPYTREE from %s TO %s" % (self.importFolder.get(), os.getcwd())
        
        copyTree(self.importFolder.get(), os.getcwd())

