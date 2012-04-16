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

# Create a GUI automatically from a selfile of micrographs
class ProtParticlePicking(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.particle_pick.name, scriptname, project)
        self.Import = "from protocol_particle_pick import *"
        importProt = self.getProtocolFromRunName(self.ImportRun)
        self.TiltPairs = os.path.exists(os.path.join(importProt.WorkingDir,"tilted_pairs.xmd"))
        if self.TiltPairs:
            self.inputMicrographs = importProt.getFilename('tiltedPairs')
        else:
            self.inputMicrographs = importProt.getFilename('micrographs')
        self.micrographs = self.getEquivalentFilename(importProt, self.inputMicrographs)

    def defineSteps(self):
        self.insertStep('copyFile', source=self.inputMicrographs, dest=self.micrographs)
        fnAcquisition=self.workingDirPath("acquisition_info.xmd")
        self.insertStep('copyAcquisitionInfo',verifyfiles=[fnAcquisition],
                        source=self.inputMicrographs,dest=fnAcquisition)
        self.insertStep('launchParticlePickingGUI',execution_mode=SqliteDb.EXEC_ALWAYS,
                           MicrographSelfile=self.micrographs, WorkingDir=self.WorkingDir,
                           TiltPairs=self.TiltPairs,
                           AutomaticPicking=self.AutomaticPicking, NumberOfThreads=self.NumberOfThreads,
                           Fast=self.Fast, InCore=self.InCore)       

    def createFilenameTemplates(self):
        return {
                'families': join('%(WorkingDir)s', 'families.xmd')
                }
        
    def summary(self):
        md = xmipp.MetaData(self.micrographs)
        if self.TiltPairs: 
            suffix = "tilt pairs"
        else: 
            suffix = "micrographs"        
        
        summary=["Input: <%s> with <%u> %s" % (self.inputMicrographs, md.size(), suffix)]        
        total_manual = 0
        N_manual = 0
        total_auto = 0
        N_auto = 0
        Nblock = {}
        for posfile in glob(self.workingDirPath('*.pos')):
            blockList = xmipp.getBlocksInMetaDataFile(posfile)
            if 'families' in blockList:
                blockList.remove('families')
            particles = 0
            for block in blockList:
                md = xmipp.MetaData("%(block)s@%(posfile)s" % locals());
                md.removeDisabled();
                Nparticles = md.size()
                particles += Nparticles
                if block in Nblock.keys():
                    Nblock[block] += Nparticles
                else:
                    Nblock[block] = Nparticles
            if particles > 0:
                if 'auto' in posfile:
                    total_auto += particles
                    N_auto += 1
                else:
                    total_manual += particles
                    N_manual += 1
        if self.TiltPairs:
            summary.append("Number of pairs manually picked: <%d> (from <%d> micrographs)" % (total_manual/2, N_manual/2))
        else:
            summary.append("Number of particles manually picked: <%d> (from <%d> micrographs)" % (total_manual, N_manual))
        if N_auto > 0:
            summary.append("Number of particles automatically picked: <%d> (from <%d> micrographs)" % (total_auto, N_auto))
        fnFamilies = self.getFilename('families')
        if exists(fnFamilies):
            md = xmipp.MetaData(fnFamilies)
            Nfamilies = md.size()
            if Nfamilies > 1:
                summary.append("Number of families: <%u>" % Nfamilies)
        for block in Nblock.keys():
            if self.TiltPairs:
                summary.append("Family <%s>: <%u> pairs" % (block, Nblock[block]/2))
            else:
                summary.append("Family <%s>: <%u> particles" % (block, Nblock[block]))
        return summary
    
    def validate(self):
        # Check that there is a valid list of micrographs
        if not exists(self.inputMicrographs):
            return ["Cannot find input micrographs: \n   <%s>" % self.inputMicrographs]
        # Check that all micrographs exist
        errors = []
        md = xmipp.MetaData(self.inputMicrographs)
        missingMicrographs = []        
        def checkMicrograph(label, objId): # Check if micrograph exists
            micrograph = md.getValue(label, objId)
            if not exists(micrograph):
                missingMicrographs.append(micrograph)
        for objId in md:
            checkMicrograph(xmipp.MDL_MICROGRAPH, objId)
            if self.TiltPairs:
                checkMicrograph(xmipp.MDL_MICROGRAPH_TILTED, objId)
        
        if len(missingMicrographs):
            errors.append("Cannot find the following micrographs: " + "\n".join(missingMicrographs))
                
        # Check that automatic particle picking is not for tilted
        if self.TiltPairs and self.AutomaticPicking:
            errors.append("Automatic particle picking cannot be done on tilt pairs")
        
        return errors
    
    def visualize(self):
        launchParticlePickingGUI(None, self.micrographs, self.WorkingDir, self.TiltPairs, ReadOnly=True)

def copyAcquisitionInfo(log,source,dest):
    MD=xmipp.MetaData()
    MD.read("acquisition_info@"+source)
    MD.write("acquisition_info@"+dest)

# Execute protocol in the working directory
def launchParticlePickingGUI(log,MicrographSelfile,WorkingDir,
                             TiltPairs=False,
                             AutomaticPicking=False, NumberOfThreads=1, Fast=True, InCore=False, ReadOnly=False):
    if TiltPairs:
        if ReadOnly:
            mode = "readonly"
        else:
            mode = "manual"
        args="-i %(MicrographSelfile)s -o %(WorkingDir)s --mode %(mode)s" % locals()
        runJob(log,"xmipp_micrograph_tiltpair_picking", args, RunInBackground=True)
    else:
        if ReadOnly:
            mode = "readonly"
        elif AutomaticPicking:
            mode = "supervised %(NumberOfThreads)d %(Fast)s %(InCore)s"
        else:
            mode = "manual"
        params = "-i %(MicrographSelfile)s -o %(WorkingDir)s --mode " + mode 
        runJob(log,"xmipp_micrograph_particle_picking", params % locals(), RunInBackground=True)

#		
# Main
#     
if __name__ == '__main__':
    protocolMain(ProtParticlePicking)
