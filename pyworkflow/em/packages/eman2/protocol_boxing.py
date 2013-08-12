'''
Created on Apr 12, 2013
modified by Josue Gomez, July 22, 2013

@author: antonio
'''

import os

from pyworkflow.em import *  
from pyworkflow.utils import * 
import eman2
from data import *
from convert import readSetOfCoordinates

class EmanDefParticlePicking(Form):
    """Create the definition of parameters for
    the Eman Boxing protocol.
    """
    def __init__(self):
        Form.__init__(self)        
        self.addSection(label='Input')
        self.addParam('inputMicrographs', PointerParam, label="Micrographs",
                      pointerClass='SetOfMicrographs',
                      help='Select the SetOfMicrograph ')
        self.addParam('boxSize', IntParam, label='Box size',
                      validators=[Positive])


class EmanProtBoxing(ProtParticlePicking):
    """Protocol to pick particles manually of a set of micrographs in the project"""    
    _definition = EmanDefParticlePicking()
    
    def __init__(self, **args):     
        ProtParticlePicking.__init__(self, **args)
        # The following attribute is only for testing
        self.importFolder = String(args.get('importFolder', None))
        
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
        coodrSet = self._createSetOfCoordinates()
        coordSet.setFileName(fn)
        coordSet.setMicrographs(self.inputMics)
        readSetOfCoordinates(self.inputMics, coordSet)
        coordSet.write()
        self._defineOutputs(outputCoordinates=coordSet)
        # As we move to workingDir we must leave it. 
        self._leaveWorkingDir()    

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

