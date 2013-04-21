'''
Created on Apr 12, 2013

@author: antonio
'''

from pyworkflow.em import *  
from pyworkflow.utils import * 
from pyworkflow.em.packages.eman2.data import *

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
        self.params = {'inputMics': ' '.join([mic.getFileName() for mic in self.inputMicrographs.get()]),
                       'boxSize': self.boxSize.get()}      
        # Insert the step with the command and the formatted parameters
        self._insertRunJobStep("e2boxer.py",  "%(inputMics)s --gui --boxsize=%(boxSize)i" % self.params)
        # Insert step to create output objects       
        self._insertFunctionStep('createOutput')
        
    def createOutput(self):
        # Generate .box files with picked coordinates.
        newBoxSize =  self.__getBoxingBoxSize()
        self.params['boxSize'] = newBoxSize
        runJob(None, "e2boxer.py", 
                     "%(inputMics)s --boxsize=%(boxSize)d --write_dbbox" % self.params) 
        # Create metadata file with all .box files referenced
        metaDataFilePath = self.__createBoxingCoordinatesMetadata()  
        # Create the SetOfCoordinates object on the database        
        self.outputCoordinates = EmanSetOfCoordinates()
        self.outputCoordinates.setBoxSize(self.params['boxSize'])
        self.outputCoordinates.setMicrographs(self.inputMicrographs.get())
        self.outputCoordinates.setFileName(metaDataFilePath)
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
    
    def __createBoxingCoordinatesMetadata(self):
        """ Create a metadata registering all coordinates files. """
        from glob import glob
        files = glob("*.box")
        if len(files) == 0:
            raise Exception('Eman Particle Picking:There is not .box files')
        boxMetadataFile=open("boxMetaData.mdbox","w")
        for line in files:
            boxMetadataFile.write(line + "\n")            
        boxMetadataFile.close()
        return self._getPath("boxMetaData.mdbox")

    