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
        self.addParam('boxSize', IntParam, label='Box size', validators=[Positive()])

class EmanProtBoxing(ProtParticlePicking):
    
    _definition = EmanDefParticlePicking()
    
    def __init__(self, **args):     
        ProtParticlePicking.__init__(self, **args)
        # The following attribute is only for testing
        self.importFolder = String(args.get('importFolder', None))
        
        
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
        program = "e2boxer.py"
        arguments = "%(inputMics)s --gui --boxsize=%(boxSize)i"
        # Run the command with formatted parameters
        self._log.info('Launching: ' + program + ' ' + arguments % self._params)
        self.runJob(None, program, arguments % self._params)
        
    def createOutput(self):
        # Get the box size store in Eman db
        self._params['boxSize'] = int(EmanDbd.getEmanParamValue('box_size'))
        if not self.importFolder.hasValue():
            program = "pwd; e2boxer.py"
            arguments = "%(inputMics)s --boxsize=%(boxSize)i --write_dbbox"
            self._log.info('Creating output: ' + program + ' ' + arguments % self._params)
            self.runJob(None, program, arguments % self._params) 
  
        # Create the SetOfCoordinates object on the database 
        outputCoordinates = EmanSetOfCoordinates(filename=self.workingDir.get())
        outputCoordinates.setBoxSize(self._params['boxSize'])
        outputCoordinates.setMicrographs(self.inputMicrographs.get())
        particlesWritten = bool(EmanDbd.getEmanParamValue('write_particles'))
        particlesFormat = str(EmanDbd.getEmanParamValue('format'))
        # As we move to workingDir we must leave it. 
        self._leaveWorkingDir()    
        
        self._defineOutputs(outputCoordinates=outputCoordinates) 
        
        if particlesWritten:
            print 'siiii tenemos particulas con format %s ' % particlesFormat.strip()
            outputParticles = EmanSetOfParticles(filename=self.workingDir.get(), format=particlesFormat.strip())
            self._defineOutputs(outputParticles=outputParticles) 

    
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

