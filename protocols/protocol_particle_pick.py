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
        if self.TiltPairs:
            self.MicrographsMd = self.PrevRun.getFilename('tilted_pairs')
        else:
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
        micrographs, particles,  = countParticles(self)
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
        
        return summary
    
    def validate(self):
        return validateMicrographs(self.Input['micrographs'], self.TiltPairs)
    
    def visualize(self):
        launchParticlePickingGUI(None, self.MicrographsMd, self.ExtraDir, PM_READONLY, self.TiltPairs)


def getPosFiles(prot, pattern=''):
    '''Return the .pos files of this picking protocol'''
    return glob(os.path.join(prot.ExtraDir,'*%s.pos' % pattern))

def getTemplateFiles(prot, pattern=''):
    '''Return the .pos files of this picking protocol'''
    return glob(os.path.join(prot.ExtraDir,'*%s_templates.stk' % pattern))


def validateMicrographs(inputMicrographs, tiltPairs=False):
    ''' Validate the existence of input micrographs metadata file 
    and also each of the micrographs, return an error list if some 
    file is missing '''
    # Check that there is a valid list of micrographs
    if not exists(inputMicrographs):
        return ["Cannot find input micrographs: \n   <%s>" % inputMicrographs]
    # Check that all micrographs exist
    errors = []
    md = xmipp.MetaData(inputMicrographs)
    missingMicrographs = []        
    
    def checkMicrograph(label, objId): # Check if micrograph exists
        micrograph = md.getValue(label, objId)
        if not xmippExists(micrograph):
            missingMicrographs.append(micrograph)
    for objId in md:
        checkMicrograph(xmipp.MDL_MICROGRAPH, objId)
        if tiltPairs:
            checkMicrograph(xmipp.MDL_MICROGRAPH_TILTED, objId)
    
    if len(missingMicrographs):
        print missingMicrographs
        errors.append("Cannot find the following micrographs: " + "\n".join(missingMicrographs))
    return errors    

def countParticles(prot, pattern=''):
    '''Return the number of picked micrographs and particles '''
    particles = 0
    micrographs = 0
    
    for posfile in getPosFiles(prot, pattern):
        pos_particles = 0
        
        blocks = xmipp.getBlocksInMetaDataFile(posfile) 
        if 'particles' in blocks: #FIXME: also read supervised
            md = xmipp.MetaData("particles@%(posfile)s" % locals());
            md.removeDisabled();
            pos_particles += md.size()
               
            if pos_particles > 0:
                particles += pos_particles
                micrographs += 1
    return micrographs, particles

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
