#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# General script for Xmipp-based particle picking
#
# Author: Carlos Oscar Sorzano, August, 2011
#
#------------------------------------------------------------------------------------------------

from config_protocols import protDict
from protlib_base import *
from protlib_utils import runJob
from protlib_filesystem import createLink
import xmipp
import glob

# Create a GUI automatically from a selfile of micrographs
class ProtParticlePicking(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.particle_pick.name, scriptname, project)
        self.Import="from protocol_particle_pick import *"
        self.importRunname=self.ImportRun.replace(protDict.import_micrographs.name,"")
        if self.importRunname[0]=="_":
            self.importRunname=self.importRunname[1:]
        self.micrographSelfile = os.path.join(protDict.import_micrographs.dir,self.importRunname, "micrographs.sel")

    def defineSteps(self):
        self.Db.insertStep('createLink',execution_mode=SqliteDb.EXEC_MAINLOOP,
                           source=self.micrographSelfile,dest=os.path.join(self.WorkingDir,"micrographs.sel"))
        self.Db.insertStep('launchParticlePickingGUI',execution_mode=SqliteDb.EXEC_ALWAYS,
                           MicrographSelfile=self.micrographSelfile, WorkingDir=self.WorkingDir,
                           AutomaticPicking=self.AutomaticPicking, NumberOfThreads=self.NumberOfThreads,
                           Fast=self.Fast,InCore=self.InCore)       

    def summary(self):
        summary = []
        
        mD=xmipp.MetaData(self.micrographSelfile)
        isPairList = mD.containsLabel(xmipp.MDL_ASSOCIATED_IMAGE1) and not xmipp.FileName(self.micrographSelfile).isStar1()

        if isPairList:
            summary=["Input: "+self.micrographSelfile+" with "+str(mD.size())+" tilt pairs"]
        else:
            summary=["Input: "+self.micrographSelfile+" with "+str(mD.size())+" micrographs"]
        
        total_manual = 0
        N_manual = 0
        for posfile in glob.glob(self.WorkingDir+"/*.pos"):
            blockList=xmipp.getBlocksInMetaDataFile(posfile)
            manual=0
            for block in blockList:
                mD=xmipp.MetaData(posfile);
                manual+=mD.size()
            if manual>0:
                total_manual+=manual
                N_manual+=1
        msg="Number of particles picked: %d (from %d micrographs" % (total_manual, N_manual)
        fnFamilies=os.path.join(self.WorkingDir,"families.xmd")
        if os.path.exists(fnFamilies):
            mD=xmipp.MetaData(fnFamilies)
            Nfamilies=mD.size()
            if Nfamilies>1:
                msg+=" and %d families"
        msg+=")"
        summary.append(msg)
        return summary
    
    def validate(self):
        errors = []
        
        # Check that there is a valid list of micrographs
        if not os.path.exists(self.micrographSelfile)>0:
            errors.append("Cannot find "+self.micrographSelfile)
            return errors
        
        # Check that all micrographs exist
        NnotFound=0
        message=""
        mD=xmipp.MetaData(self.micrographSelfile)
        isPairList = mD.containsLabel(xmipp.MDL_ASSOCIATED_IMAGE1) and not xmipp.FileName(self.micrographSelfile).isStar1()
        for id in mD:
             micrograph = mD.getValue(xmipp.MDL_IMAGE,id)
             if not os.path.exists(micrograph):
                message+="  "+micrograph
                NnotFound=NnotFound+1
             if isPairList:
                 micrograph = mD.getValue(xmipp.MDL_ASSOCIATED_IMAGE1,id)
                 if not os.path.exists(micrograph):
                     message+="  "+micrograph
                     NnotFound=NnotFound+1
        
        if NnotFound>0:
            errors.append("Cannot find the following micrographs: "+message)
                
        # Check that automatic particle picking is not for tilted
        if isPairList and AutomaticPicking:
            errors.append("Automatic particle picking cannot be done on tilt pairs")
        
        return errors
    
    def visualize(self):
        launchParticlePickingGUI(None, self.micrographSelfile, self.WorkingDir)

# Execute protocol in the working directory
def launchParticlePickingGUI(log,MicrographSelfile,WorkingDir,
                             AutomaticPicking=False,NumberOfThreads=1,Fast=True,InCore=False):
    params="-i %s -o %s"%(MicrographSelfile,WorkingDir)
    if AutomaticPicking:
        params+=" --auto %d %s %s"%(NumberOfThreads,Fast,InCore)
    runJob(log,"xmipp_micrograph_particle_picking",params,RunInBackground=True)

#		
# Main
#     
if __name__ == '__main__':
    protocolMain(ProtParticlePicking)
