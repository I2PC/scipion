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
from protlib_xmipp import redStr
from protlib_particles import countParticles
import xmipp
from glob import glob
from os.path import exists, join

# Picking modes
PM_MANUAL, PM_READONLY, PM_REVIEW = ('manual', 'readonly', 'review')

# Create a GUI automatically from a selfile of micrographs
class ProtParticlePicking(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.particle_pick.name, scriptname, project)
        self.Import = "from protocol_particle_pick import *"
        self.setPreviousRun(self.ImportRun)
        self.importDir = self.PrevRun.WorkingDir
        self.inputProperty('TiltPairs')
        self.MicrographsMd = self.PrevRun.getFilename('micrographs')
        self.Input['micrographs'] = self.MicrographsMd
        self.inputFilename('microscope', 'acquisition')
        self.MicrographsMd = self.getEquivalentFilename(self.PrevRun, self.MicrographsMd)

    def defineSteps(self):
        self.insertStep("createDir",verifyfiles=[self.ExtraDir],path=self.ExtraDir)
        self.insertImportOfFiles([self.Input['micrographs'], self.Input['acquisition']])
        if getattr(self, 'LaunchGUI', True):
            self.insertStep('launchParticlePickingGUI',execution_mode=SqliteDb.EXEC_ALWAYS,
                           InputMicrographs=self.MicrographsMd, ExtraDir=self.ExtraDir,
                           TiltPairs=self.TiltPairs, Memory=self.Memory)       
    
    def createFilenameTemplates(self):
        return {
                 'pos': join('%(ExtraDir)s', '%(micrograph)s.pos'),
                 'templates': join('%(ExtraDir)s', 'templates.stk'),
#      
                'config': join('%(ExtraDir)s','config.xmd')
                }  
          
    def summary(self):
        micrographs, particles = countParticles(self.ExtraDir)
        if self.TiltPairs: 
            suffix = "tilt pairs"
            items = "pairs"
            micrographs /= 2
            particles /= 2
        else: 
            suffix = "micrographs"
            items = "particles"        
        from protlib_xmipp import getMdSize
        size = getMdSize(self.MicrographsMd)
        summary = ["Input: [%s] with <%u> %s" % (self.importDir, size, suffix),         
                   "Number of %(items)s manually picked: <%(particles)d> (from <%(micrographs)d> micrographs)" % locals()]
        md=xmipp.MetaData(self.MicrographsMd)
        if not md.containsLabel(xmipp.MDL_CTF_MODEL):
            summary.append(redStr("There is no CTF information in the input micrographs: "))
            summary.append("[%s]"%self.MicrographsMd)
        
        return summary
    
    def validate(self):
        errors = []
        if not exists(self.Input['micrographs']):
            errors.append("Cannot find input micrographs: \n   <%s>" % inputMicrographs)
        return errors
    
    def visualize(self):
        launchParticlePickingGUI(None, self.MicrographsMd, self.ExtraDir, PM_READONLY, self.TiltPairs)

def launchParticlePickingGUI(log, InputMicrographs, ExtraDir, PickingMode=PM_MANUAL,
                             TiltPairs=False, Memory=2):
    ''' Utility function to launch the Particle Picking application '''
    args = "-i %(InputMicrographs)s -o %(ExtraDir)s --mode %(PickingMode)s --memory %(Memory)dg"
   
    if TiltPairs:
        program = "xmipp_micrograph_tiltpair_picking"
    else:
        program = "xmipp_micrograph_particle_picking"
    runJob(log, program, args % locals(), RunInBackground=False)
#		
# Main
#     
if __name__ == '__main__':
    protocolMain(ProtParticlePicking)
