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
import xmipp

# Create a GUI automatically from a selfile of micrographs
class ProtParticlePicking(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.particle_pick.name, scriptname, project)
        self.Import="from protocol_particle_pick import *"
        self.preprocessingRunname=self.PreprocessingRun.replace(protDict.preprocess_micrographs.name,"")
        if self.preprocessingRunname[0]=="_":
            self.preprocessingRunname=self.preprocessingRunname[1:]
        self.micrographSelfile = os.path.join(protDict.preprocess_micrographs.dir,self.preprocessingRunname, "micrographs.sel")

    def defineSteps(self):
        print "Defining steps"
        self.Db.insertStep('launchParticlePickingGUI',execution_mode=SqliteDb.EXEC_ALWAYS,
                           MicrographSelfile=self.micrographSelfile, WorkingDir=self.WorkingDir,
                           AutomaticPicking=self.AutomaticPicking, NumberOfThreads=self.NumberOfThreads)       

    def summary(self):
        summary = []

        mD=xmipp.MetaData(self.micrographSelfile)
        isPairList = mD.containsLabel(xmipp.MDL_ASSOCIATED_IMAGE1) and not xmipp.FileName(self.micrographSelfile).isStar1()

        if self.isPairList:
            summary=["Input: "+self.micrographSelfile+" (Tilt pairs)"]
        else:
            summary=["Input: "+self.micrographSelfile]
        
        total_manual = 0
        total_auto = 0
        N_manual = 0
        N_auto = 0
        for id in mD:
             micrograph = mD.getValue(xmipp.MDL_IMAGE,id)
             manual=CountPicked(self.WorkingDir,micrograph,"Common")
             if manual>0:
                 total_manual+=manual
                 N_manual+=1
             if AutomaticPicking:
                 auto=CountPicked(self.WorkingDir,micrograph,"Common.auto")
                 if auto>0:
                     total_auto+=auto
                     N_auto+=1
        summary.append("# Manually picked: %d (from %d micrographs)" % (total_manual, N_manual))
        if AutomaticPicking:
            summary.append("# Automatically picked: %d (from %d micrographs) " % (total_auto, N_auto))
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

def CountPicked(WorkingDir,micrograph,label):
    posfile=WorkingDir+"/"+str(micrograph)+'.'+label+'.pos'
    if os.path.exists(posfile):
        mD=xmipp.MetaData(posfile);
        return mD.size()
    return 0
    
# Execute protocol in the working directory
def launchParticlePickingGUI(log,MicrographSelfile,WorkingDir,AutomaticPicking,NumberOfThreads):
    print "in launch"
    params="-i %s -o %s"%(MicrographSelfile,WorkingDir)
    if AutomaticPicking:
        params+=" -auto"
        if NumberOfThreads>1:
            params+=" -thr %d"%NumberOfThreads
    runJob(log,"xmipp_micrograph_particle_picking",params,RunInBackground=True)

#		
# Main
#     
if __name__ == '__main__':
    protocolMain(ProtParticlePicking)
