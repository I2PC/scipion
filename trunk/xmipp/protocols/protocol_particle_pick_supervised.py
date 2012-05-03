#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# General script for Xmipp-based particle picking
#
# Author: Carlos Oscar Sorzano, August, 2011
#
#------------------------------------------------------------------------------------------------

from config_protocols import protDict
from protlib_base import *
from protlib_utils import runJob
from protlib_filesystem import createLink, copyFile
import xmipp
from glob import glob
from os.path import exists, join

from protocol_particle_pick import getPosFiles, validateMicrographs,\
    launchParticlePickingGUI, PM_READONLY, PM_SUPERVISED, countParticles

# Create a GUI automatically from a selfile of micrographs
class ProtParticlePickingSupervised(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.particle_pick_supervised.name, scriptname, project)
        self.Import = "from protocol_particle_pick_supervised import *"
        pickingProt = self.getProtocolFromRunName(self.PickingRun)
        self.pickingDir = pickingProt.WorkingDir
        self.inputMicrographs = pickingProt.inputMicrographs
        self.micrographs = self.getEquivalentFilename(pickingProt, self.inputMicrographs)
        self.pickingProt = pickingProt

    def defineSteps(self):
        self.insertStep('copyFile', verifyfiles=[self.micrographs], source=self.inputMicrographs, dest=self.micrographs)
        filesToCopy = getPosFiles(self.pickingProt)
        filesToCopy += [self.pickingProt.getFilename('families'), 
                        self.pickingProt.getFilename('macros'),
                        self.inputMicrographs]
        for pos in filesToCopy:
            posDest = self.getEquivalentFilename(self.pickingProt, pos)
            self.insertStep('copyFile', verifyfiles=[posDest], source=pos, dest=posDest)
        self.insertStep('createLink2', filename="acquisition_info.xmd",dirSrc=self.pickingDir, dirDest=self.WorkingDir)
        modeWithArgs = PM_SUPERVISED + " %(NumberOfThreads)d %(Fast)s %(InCore)s" % self.ParamsDict
        self.insertStep('launchParticlePickingGUI',execution_mode=SqliteDb.EXEC_ALWAYS,
                           InputMicrographs=self.micrographs, WorkingDir=self.WorkingDir,
                           PickingMode=modeWithArgs, Memory=self.Memory)       
        
    def summary(self):
        md = xmipp.MetaData(self.micrographs)
        micrographs, particles, familiesDict = countParticles(self, 'auto')
        suffix = "micrographs"        
        summary = ["Manual picking RUN: <%s> " % self.PickingRun,
                   "Input picking: [%s] with <%u> %s" % (self.pickingDir, md.size(), suffix)]        

        return summary
    
    def validate(self):
        return validateMicrographs(self.inputMicrographs)
    
    def visualize(self):
        launchParticlePickingGUI(None, self.micrographs, self.WorkingDir, PM_READONLY)

# Main
#     
if __name__ == '__main__':
    protocolMain(ProtParticlePickingSupervised)
