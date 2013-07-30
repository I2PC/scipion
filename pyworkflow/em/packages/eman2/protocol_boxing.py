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
        loadEnvironment()
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
        arguments = "%(inputMics)s --boxsize=%(boxSize)i"
        # Run the command with formatted parameters
        self._log.info('Launching: ' + program + ' ' + arguments % self._params)
        self.runJob(None, program, arguments % self._params)
    
    def _createSetOfCoordinates(self, size):
        mics = self.inputMicrographs.get()
        # Create a json file with the coordinates file for each micrograph
        boxId = 0L
        jsonDict = {}
        for mic in mics:
            micId = mic.getId()
            micFnroot = removeBaseExt(mic.getFileName()) + '_info.json'
            micPosFn = self._getRelPath("info", micFnroot)
            if exists(micPosFn):
                jsonDict[micId] = micPosFn
                jsonPosDict = loadJson(micPosFn)
                boxes = jsonPosDict["boxes"]
                listbox = []
                for box in boxes:
                    boxId += 1
                    listbox.append(boxId)
                jsonPosDict["coordId"] = listbox
                writeJson(jsonPosDict, micPosFn)
        jsoncoordsFn = self._getRelPath('scipion_micrographs_coordinates.json')
        writeJson(jsonDict, jsoncoordsFn)
        coords = EmanSetOfCoordinates(filename=jsoncoordsFn)
        coords.setMicrographs(mics)
        coords.boxSize.set(size)
        return coords
        
    def createOutput(self):
        # Get the box size store in Eman json
        jsonBoxDict = loadJson(self._getRelPath("e2boxercache", 'base.json'))
        size = int(jsonBoxDict["box_size"])
        
#         if not self.importFolder.hasValue():
#             program = "e2boxer.py"
#             arguments = "%(inputMics)s --boxsize=%(boxSize)i"
#             self._log.info('Creating output: ' + program + ' ' + arguments % self._params)
#             self.runJob(None, program, arguments % self._params) 
  
        # Create the SetOfCoordinates object on the database 
        #TODO: Create a json file with pairs micrographId, jsonPosFile and point EmanSetOFCoordinates to it
        
        outputCoordinates = self._createSetOfCoordinates(size)
        
        # TODO: Now particlesWritten will come from form paramter and if yes e2boxer will need to be executed again
#        particlesWritten = bool(EmanDbd.getEmanParamValue('write_particles'))
        # TODO: No need to retreive format anymore
#        particlesFormat = str(EmanDbd.getEmanParamValue('format'))
        # As we move to workingDir we must leave it. 
        self._leaveWorkingDir()    
        self._defineOutputs(outputCoordinates=outputCoordinates) 
        
#         if particlesWritten:
#             # TODO: GEnerate lst with e2buildsets or maybe a unique hdf
#             print 'siiii tenemos particulas con format %s ' % particlesFormat.strip()
#             outputParticles = EmanSetOfParticles(filename=self.workingDir.get(), format=particlesFormat.strip())
#             self._defineOutputs(outputParticles=outputParticles) 
#     
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

