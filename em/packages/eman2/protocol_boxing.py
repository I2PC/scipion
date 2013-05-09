'''
Created on Apr 12, 2013

@author: antonio
'''

from pyworkflow.em import *  
from pyworkflow.utils import * 
from pyworkflow.em.packages.eman2.data import *
import os

class EmanDefParticlePicking(Form):
    """Create the definition of parameters for
    the Eman Boxing protocol.
    """
    def __init__(self):
        Form.__init__(self)        
        self.addSection(label='Input')
        self.addParam('inputMicrographs', PointerParam, label="Micrographs", pointerClass='SetOfMicrographs')
        self.addParam('boxSize', IntParam, label='Box size')

class EmanProtBoxing(ProtParticlePicking):
    
    _definition = EmanDefParticlePicking()
    
    def __init__(self, **args):        
        ProtParticlePicking.__init__(self, **args)
        
    def _runSteps(self, startIndex):
        # Redefine run to change to workingDir path
        # Change to protocol working directory
        self._enterWorkingDir()
        Protocol._runSteps(self, startIndex)
        
    def _defineSteps(self):
        micList = [os.path.relpath(mic.getFileName(), self.workingDir.get()) for mic in self.inputMicrographs.get()]
        self._params = {'inputMics': ' '.join(micList), 
                        'boxSize': self.boxSize.get()}      
        # Launch Boxing GUI
        self._insertFunctionStep('launchBoxingGUI', isInteractive=True) 
        # Insert step to create output objects       
        self._insertFunctionStep('createOutput')
    
    def launchBoxingGUI(self):
        # First we go to runs directory (we create if it does not exist)
        #path.makePath("runs")
        # Program to execute and it arguments
        program = "e2boxer.py"
        arguments = "%(inputMics)s --gui --boxsize=%(boxSize)i"
        # Run the command with formatted parameters
        self._log.info('Launching... ' + program + ' ' + arguments % self._params)
        runJob(None, program, arguments % self._params)
        
    def createOutput(self):
        # Get the box size store in Eman db
        self._params['boxSize'] = self.__getBoxingBoxSize()
        program = "pwd; e2boxer.py"
        arguments = "%(inputMics)s --boxsize=%(boxSize)i --write_dbbox"
        self._log.info('Creating output... ' + program + ' ' + arguments % self._params)
        runJob(None, program, arguments % self._params) 
        # As we move to workingDir we must leave it. 
        self._leaveWorkingDir()      
        # Create the SetOfCoordinates object on the database 
        self.outputCoordinates = EmanSetOfCoordinates(filename=self.workingDir.get())
        self.outputCoordinates.setBoxSize(self._params['boxSize'])
        self.outputCoordinates.setMicrographs(self.inputMicrographs.get())
        self._defineOutputs(outputCoordinates=self.outputCoordinates) 
    
    def __getBoxingBoxSize(self):
        """ Recover the box size from EMAN Berkeley data base. """        
        command = "e2bdb.py -D bdb:emboxerbase"
        pipe = os.popen(command)
        stOutput = pipe.readlines()
        pipe.close()
        auxBoxSize = None
        for line in stOutput:
            if ("box_size" in line):
                auxBoxSize = int(line.split(" : ")[1])
        if auxBoxSize is None:
            raise Exception("Error getting the stored boxsize with command: " + command) 
        return auxBoxSize
    
    def getFiles(self):
        files = self.inputMicrographs.get().getFiles() | self.outputCoordinates.getFiles()
        return files

    