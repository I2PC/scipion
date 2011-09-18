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
from protlib_filesystem import createLink, deleteFiles
from xmipp import MetaData, MD_APPEND,MDL_IMAGE, MDL_PICKING_FAMILY, \
                  MDL_PICKING_MICROGRAPH_FAMILY_STATE, MDL_PICKING_PARTICLE_SIZE,\
                  MDL_ENABLED, MDL_COST
import glob

# Create a GUI automatically from a selfile of micrographs
class ProtParticlePickingAuto(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.particle_pick_auto.name, scriptname, project)
        self.Import="from protocol_particle_pick_auto import *"
        self.pickingDir=getWorkingDirFromRunName(self.PickingRun)
        self.micrographSelfile = os.path.join(self.pickingDir, "micrographs.sel")
        self.familiesFile=os.path.join(self.pickingDir, "families.xmd")

        # Families to process
        self.familiesForAuto=[]
        self.particleSizeForAuto=[]
        mD=MetaData(self.familiesFile)
        for id in mD:
            family=mD.getValue(MDL_PICKING_FAMILY,id)
            if os.path.exists(os.path.join(self.pickingDir,family+"_training.txt")):
                self.familiesForAuto.append(family)
                particleSize=mD.getValue(MDL_PICKING_PARTICLE_SIZE,id)
                self.particleSizeForAuto.append(particleSize)

    def defineSteps(self):
        self.Db.insertStep('createLink',execution_mode=SqliteDb.EXEC_MAINLOOP,
                           source=self.micrographSelfile,dest=self.workingDirPath("micrographs.sel"))
        self.Db.insertStep('createLink',execution_mode=SqliteDb.EXEC_MAINLOOP,
                           source=self.familiesFile,dest=self.workingDirPath("families.xmd"))
        for familyIdx in range(len(self.familiesForAuto)):
            family=self.familiesForAuto[familyIdx]
            modelRoot=self.workingDirPath(family)
            self.Db.insertStep('createLink',execution_mode=SqliteDb.EXEC_MAINLOOP,
                               source=os.path.join(self.pickingDir,family+"_training.txt"),dest=modelRoot+"_training.txt")
            self.Db.insertStep('createLink',execution_mode=SqliteDb.EXEC_MAINLOOP,
                               source=os.path.join(self.pickingDir,family+"_mask.xmp"),dest=modelRoot+"_mask.xmp")

        idMPI=self.Db.insertStep('runStepGapsMpi',passDb=True, script=self.scriptName, NumberOfMpi=self.NumberOfMpi)
        mDmicrographs=MetaData(self.micrographSelfile)
        for familyIdx in range(len(self.familiesForAuto)):
            family=self.familiesForAuto[familyIdx]
            particleSize=self.particleSizeForAuto[familyIdx]
            modelRoot=self.workingDirPath(family)
            for id in mDmicrographs:
                # Get Micrograph name
                micrographFullPath=mDmicrographs.getValue(MDL_IMAGE,id)
                micrographName=os.path.split(os.path.split(micrographFullPath)[0])[1]
                proceed=False
                fnPos=os.path.join(self.pickingDir,micrographName+".pos")
                if not os.path.exists(fnPos):
                    proceed=True
                else:
                    mDaux=MetaData("families@"+fnPos)
                    proceed=True
                    for idaux in mDaux:
                        familyAux=mDaux.getValue(MDL_PICKING_FAMILY,idaux)
                        if family==familyAux:
                            state=mDaux.getValue(MDL_PICKING_MICROGRAPH_FAMILY_STATE,idaux)
                            if state!="Available":
                                proceed=False
                                break
                if proceed:
                    script=getScriptFromRunName(self.PickingRun)
                    protPicking=getProtocolFromModule(script, self.project)
                    id=self.Db.insertStep('autoPick',parent_step_id=XmippProjectDb.FIRST_STEP,
                                          execution_mode=SqliteDb.EXEC_GAP,WorkingDir=self.WorkingDir,
                                          ModelRoot=modelRoot,MicrographFullPath=micrographFullPath,MicrographName=micrographName,
                                          ParticleSize=particleSize,Fast=protPicking.Fast,InCore=protPicking.InCore)
        parent_id=idMPI
        for familyIdx in range(len(self.familiesForAuto)):
            family=self.familiesForAuto[familyIdx]
            parent_id=self.Db.insertStep('gatherResults',parent_step_id=parent_id,Family=family,WorkingDir=self.WorkingDir,PickingDir=self.pickingDir)
        self.Db.insertStep('deleteTempFiles',parent_step_id=parent_id,WorkingDir=self.WorkingDir)

    def summary(self):
        summary = []
        summary.append("Manual picking run: "+self.PickingRun)
        summary.append("Automatic picking of the following models: "+",".join(self.familiesForAuto))
        autoFiles = glob.glob(self.workingDirPath("*auto.pos"))
        Nparticles=0
        Nmicrographs=len(autoFiles)
        first=True
        maxTime=0
        minTime=0
        for file in autoFiles:
            mD=MetaData(file)
            Nparticles+=mD.size()
            if first:
                first=False
                minTime=os.stat(file).st_mtime
                maxTime=minTime
            else:
                t=os.stat(file).st_mtime
                minTime=min(minTime,t)
                maxTime=max(maxTime,t)
        summary.append("%d particles automatically picked from %d micrographs in %d minutes"%(Nparticles,Nmicrographs,int((maxTime-minTime)/60.0)))
        return summary
    
    def validate(self):
        errors = []
        return errors
    
    def visualize(self):
        mD=MetaData(self.familiesFile)
        for file in glob.glob(self.workingDirPath("*extract_list.xmd")):
            params="-i %s -o %s --mode review %s"%(self.micrographSelfile,self.WorkingDir,file)
            runJob(None,"xmipp_micrograph_particle_picking",params,RunInBackground=True)

def autoPick(log,WorkingDir,ModelRoot,MicrographFullPath,MicrographName,ParticleSize,Fast,InCore):
    args="-i "+MicrographFullPath+" --particleSize "+str(ParticleSize)+" --model "+ModelRoot+\
         " --outputRoot "+os.path.join(WorkingDir,MicrographName)+" --mode autoselect"
    runJob(log,"xmipp_micrograph_automatic_picking",args)

def gatherResults(log,Family,WorkingDir,PickingDir):
    mD=MetaData(os.path.join(WorkingDir,"micrographs.sel"))
    mDpos=MetaData()
    mDposAux=MetaData()
    fnExtractList=os.path.join(WorkingDir,Family+"_extract_list.xmd")
    for id in mD:
        micrographFullName=mD.getValue(MDL_IMAGE,id)
        micrographName=os.path.split(os.path.split(micrographFullName)[0])[1]
        
        fnManual=os.path.join(PickingDir,micrographName+".pos")
        fnAuto1=os.path.join(PickingDir,micrographName+"_auto.pos")
        fnAuto2=os.path.join(WorkingDir,micrographName+"_auto.pos")
        
        mDpos.clear()
        if os.path.exists(fnManual):
            try:            
                mDpos.read(Family+"@"+fnManual)
                mDpos.setValueCol(MDL_COST,2.0)
            except:
                pass
        if os.path.exists(fnAuto1):
            try:
                mDposAux.read(Family+"@"+fnAuto1)
                mDposAux.removeDisabled();
                mDposAux.removeLabel(MDL_ENABLED)
                mDpos.unionAll(mDposAux) 
            except:
                pass
        if os.path.exists(fnAuto2):
            try:
                mDposAux.read(Family+"@"+fnAuto2)
                mDposAux.removeDisabled();
                mDposAux.removeLabel(MDL_ENABLED)
                mDpos.unionAll(mDposAux) 
            except:
                pass
        # Append alphanumeric prefix to help identifying the block 
        mDpos.write("mic_"+micrographName+"@"+fnExtractList,MD_APPEND)

def deleteTempFiles(log,WorkingDir):
    deleteFiles(log,glob.glob(os.path.join(WorkingDir,"*_auto.pos")),True)

#		
# Main
#     
if __name__ == '__main__':
    protocolMain(ProtParticlePicking)
