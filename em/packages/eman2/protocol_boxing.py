'''
Created on Apr 12, 2013

@author: antonio
'''

from pyworkflow.em import *  
from pyworkflow.utils import * 

def defineBoxing():
    """Create the definition of parameters for
    the Eman Boxing protocol"""
    f = Form()
    
    f.addSection(label='Input')
    f.addParam('inputMicrographs', PointerParam, label="Micrographs", pointerClass='SetOfMicrographs')
    f.addParam('boxSize', IntParam, label='Box size')
    return f

class EmanProtBoxing():
    
    _definition = defineBoxing()
    
    def __init__(self, **args):
        
        Protocol.__init__(self, **args)
        
    def defineSteps(self):
        self.inputMics = self.inputMicrographs.get() 
        
        micrographListParam = ""
        
        for mic in self.inputMics:
            fn = mic.getFileName()
            micrographListParam += fn + " "
        
        self.params = {'inputMics': micrographListParam,
                       'boxSize': self.boxSize.get()}
        
        self.insertStepsForBoxing()
        
    def insertStepsForBoxing(self):
        # Boxing
        self._insertOneStep("e2boxing.py",
                            "%(inputMic)s --gui --boxsize=%(boxSize)i ")
        
        # Generate coordinates
        # User can change box size
        newBoxSize =  self.__getBoxingBoxSize()
        self._insertOneStep("e2boxing.py",
                            "%(inputMic)s --boxsize=" + newBoxSize +" --write_dbbox") 
                    
        
    def _insertOneStep(self,  program, arguments):        
        # Insert the command with the formatted parameters
        self.insertRunJobStep(program, arguments % self.params)
        
    def createOutput(self):
        self.__createBoxingCoordinatesMetadata()   
    
    def __getBoxingBoxSize(self):
        #os.system("cd " + directory)
        command = "e2bdb.py -D bdb:emboxerbase"
        pipe = os.popen(command)
        stOutput = pipe.readlines()
        pipe.close()
        for line in stOutput:
            if ("box_size" in line):
                auxBoxSize = line.split(" : ")[1]
        return int(auxBoxSize)
    
    def __createBoxingCoordinatesMetadata(self):
        command = "ls *.box"
        pipe = os.popen(command)
        stOutput = pipe.readlines()
        pipe.close()
        boxMetadataFile=open("boxMetaData.mdbox","w")
        for line in stOutput:
            boxMetadataFile.write(line + "\n")
        boxMetadataFile.close()