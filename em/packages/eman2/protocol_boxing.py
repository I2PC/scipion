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
        # The following attribute is only for testing
        self.importFolder = args.get('importFolder', None)
        
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
        if self.importFolder is None:
            self._insertFunctionStep('launchBoxingGUI', isInteractive=True)
        else: # This is only used for test purposes
            self._insertFunctionStep('_importFromFolder')  
        # Insert step to create output objects       
        self._insertFunctionStep('createOutput')
    
    def launchBoxingGUI(self):
        # First we go to runs directory (we create if it does not exist)
        #path.makePath("runs")
        # Program to execute and it arguments
        program = "e2boxer.py"
        arguments = "%(inputMics)s --gui --boxsize=%(boxSize)i"
        # Run the command with formatted parameters
        self._log.info('Launching: ' + program + ' ' + arguments % self._params)
        self.runJob(None, program, arguments % self._params)
        
    def createOutput(self):
        # Get the box size store in Eman db
        self._params['boxSize'] = int(self.__getEmanParamValue('box_size'))
        program = "pwd; e2boxer.py"
        arguments = "%(inputMics)s --boxsize=%(boxSize)i --write_dbbox"
        self._log.info('Creating output: ' + program + ' ' + arguments % self._params)
        self.runJob(None, program, arguments % self._params) 
        # As we move to workingDir we must leave it. 
        self._leaveWorkingDir()      
        # Create the SetOfCoordinates object on the database 
        self.outputCoordinates = EmanSetOfCoordinates(filename=self.workingDir.get())
        self.outputCoordinates.setBoxSize(self._params['boxSize'])
        self.outputCoordinates.setMicrographs(self.inputMicrographs.get())
        particlesWritten = bool(self.__getEmanParamValue('write_particles'))
        if particlesWritten:
            print 'siiii tenemos particulas'
            self.outputImages = EmanSetOfImages(filename=self.workingDir.get())
            
        self._defineOutputs(outputCoordinates=self.outputCoordinates) 
    
    def __getEmanParamValue(self, paramName):
        """ Recover a parameter value from EMAN Berkeley data base. """        
        command = "e2bdb.py -D bdb:emboxerbase"
        pipe = os.popen(command)
        stOutput = pipe.readlines()
        pipe.close()
        auxValue = None
        for line in stOutput:
            if (paramName in line):
                auxValue = line.split(" : ")[1]
        if auxValue is None:
            raise Exception("Error getting the stored paramter with command: " + command) 
        return auxValue
    
    def getFiles(self):
        filePaths = self.inputMicrographs.get().getFiles() | ProtParticlePicking.getFiles(self)
        return filePaths
      
    def _importFromFolder(self):
        """ This function will copy Eman .box files for
        simulating an particle picking run...this is only
        for testing purposes.
        """
        from pyworkflow.utils.path import copyTree

        copyTree(self.importFolder, os.getcwd())


