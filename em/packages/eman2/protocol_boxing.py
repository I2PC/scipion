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
    the Eman Boxing protocol"""
    def __init__(self):
        Form.__init__(self)
        
        self.addSection(label='Input')
        self.addParam('inputMicrographs', PointerParam, label="Micrographs", pointerClass='SetOfMicrographs')
        self.addParam('boxSize', IntParam, label='Box size')

class EmanProtBoxing(ProtParticlePicking):
    
    _definition = EmanDefParticlePicking()
    
    def __init__(self, **args):
        
        ProtParticlePicking.__init__(self, **args)
        
    def _defineSteps(self):
        self._params = {'inputMics': ' '.join([os.path.relpath(mic.getFileName(), self.workingDir.get()) for mic in self.inputMicrographs.get()]),
                       'boxSize': self.boxSize.get()}      
        # Launch Boxing GUI
        self._insertFunctionStep('launchBoxingGUI', isInteractive=True) 
        # Insert step to create output objects       
        self._insertFunctionStep('createOutput')
    
    def launchBoxingGUI(self):
        # First we go to runs directory (we create if it does not exist)
        #path.makePath("runs")
        os.chdir(self.workingDir.get())
        # Program to execute and it arguments
        program = "e2boxer.py"
        arguments = "%(inputMics)s --gui --boxsize=%(boxSize)i"
        # Run the command with formatted parameters
        runJob(None, program, arguments % self._params)
        
    def createOutput(self):
        # Generate .box files with picked coordinates.
        newBoxSize =  self.__getBoxingBoxSize()
        self._params['boxSize'] = newBoxSize
        runJob(None, "e2boxer.py", 
                     "%(inputMics)s --boxsize=%(boxSize)d --write_dbbox" % self._params)   
        # Create the SetOfCoordinates object on the database        
        self.outputCoordinates = EmanSetOfCoordinates()
        self.outputCoordinates.setBoxSize(self._params['boxSize'])
        self.outputCoordinates.setMicrographs(self.inputMicrographs.get())
        self._defineOutputs(outputCoordinates=self.outputCoordinates) 
    
    def __getBoxingBoxSize(self):
        """ Recover the box size from EMAN Berkeley data base. """
        command = "e2bdb.py -D bdb:emboxerbase"
        pipe = os.popen(command)
        stOutput = pipe.readlines()
        pipe.close()
        for line in stOutput:
            if ("box_size" in line):
                auxBoxSize = line.split(" : ")[1]
        return int(auxBoxSize)

    