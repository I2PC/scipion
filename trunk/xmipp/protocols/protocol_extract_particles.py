#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based pre-processing of single-particles: 
#  - phase flipping
#  - extraction of particles
#  - normalization
#  - sort_junk

# Author: Carlos Oscar, August 2011
#
from protlib_base import *
import xmipp
import glob
import os
from protlib_utils import runJob
from protlib_filesystem import deleteFile

class ProtExtractParticles(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.extract_particles.name, scriptname, project)
        self.Import = 'from protocol_extract_particles import *'
        # COSS: Falta incluir una ejecucion de automaticos
        # COSS: Falta tiltpairs
        self.pickingRunname=self.PickingRun.replace(protDict.particle_pick.name,"")
        if self.pickingRunname[0]=="_":
            self.pickingRunname=self.pickingRunname[1:]
        self.pickingDir= os.path.join(protDict.particle_pick.dir,self.pickingRunname)
        self.familyFile = os.path.join(self.pickingDir, "families.xmd")

    def defineSteps(self):
        families=xmipp.MetaData(self.familyFile)
        self.familyList=[]
        for id in families:
            familyName=families.getValue(xmipp.MDL_PICKING_FAMILY,id)
            particleSize=families.getValue(xmipp.MDL_PICKING_PARTICLE_SIZE,id)
            self.familyList.append((familyName,particleSize))
            familyDir=os.path.join(self.WorkingDir,familyName)
            self.Db.insertStep('createDir',verifyfiles=[familyDir],path=familyDir)
        self.Db.insertStep('copyFile',source=os.path.join(self.pickingDir,"micrographs.sel"),dest=os.path.join(self.WorkingDir,"micrographs.sel"))
        
        idMPI=self.Db.insertStep('runStepGapsMpi',passDb=True, script=self.scriptName, NumberOfMpi=self.NumberOfMpi)
        verifyFiles=[]        
        mD=xmipp.MetaData(os.path.join(self.pickingDir,"micrographs.sel"))
        containsCTF=mD.containsLabel(xmipp.MDL_CTFMODEL)
        for id in mD:
            micrograph=mD.getValue(xmipp.MDL_IMAGE,id)
            micrographName=os.path.split(os.path.split(micrograph)[0])[1]
            posFile=os.path.join(self.pickingDir,micrographName+".pos")
            if os.path.exists(posFile):
                parent_id=XmippProjectDb.FIRST_STEP
                micrographToExtract=micrograph
                if containsCTF:
                    ctf=mD.getValue(xmipp.MDL_CTFMODEL,id)
                else:
                    ctf=None
                if self.DoFlip:
                    micrographToExtract=os.path.join(self.TmpDir,micrographName+"_flipped.xmp")
                    parent_id=self.Db.insertStep('phaseFlip',parent_step_id=XmippProjectDb.FIRST_STEP,
                                                 verifyfiles=[micrographToExtract],execution_mode=SqliteDb.EXEC_GAP,
                                                 micrograph=micrograph,ctf=ctf,fnOut=micrographToExtract)
                (tasks,outputFiles)=self.whichTasks(posFile,micrographName)
                verifyFiles+=outputFiles
                parent_id=self.Db.insertStep('extractTasks',verifyfiles=outputFiles,execution_mode=SqliteDb.EXEC_GAP,parent_step_id=parent_id,
                                      micrograph=micrograph,ctf=ctf,micrographToExtract=micrographToExtract,tasks=tasks,
                                      doNorm=self.DoNorm,doLog=self.DoLog,doInvert=self.DoInvert,
                                      bgRadius=self.BackGroundRadius, doRemoveDust=self.DoRemoveDust,
                                      dustRemovalThreshold=self.DustRemovalThreshold)
                if self.DoFlip:
                    self.Db.insertStep('deleteFile',execution_mode=SqliteDb.EXEC_GAP,parent_step_id=parent_id,filename=micrographToExtract,verbose=True)
        self.Db.updateVerifyFiles(idMPI,verifyFiles)
        id=self.Db.insertStep('gatherSelfiles',parent_step_id=idMPI,
                           WorkingDir=self.WorkingDir,familyList=self.familyList)
        idMPI=self.Db.insertStep('runStepGapsMpi',passDb=True, script=self.scriptName, NumberOfMpi=self.NumberOfMpi)
        verifyFiles=[]        
        for family in self.familyList:
            selfileRoot=os.path.join(self.WorkingDir,family[0])
            fnOut=selfileRoot+"_sorted.sel"
            self.Db.insertStep('sortImageInFamily',parent_step_id=XmippProjectDb.FIRST_STEP,execution_mode=SqliteDb.EXEC_GAP,
                               verifyfiles=[fnOut],selfileRoot=selfileRoot)
            verifyFiles.append(fnOut)
        self.Db.updateVerifyFiles(idMPI,verifyFiles)
        self.Db.insertStep('avgZscore',parent_step_id=idMPI,
                           WorkingDir=self.WorkingDir,familyList=self.familyList,micrographSelfile=os.path.join(self.WorkingDir,"micrographs.sel"))
                
    def validate(self):
        errors = []
        if not os.path.exists(self.familyFile):
            errors.append("Cannot find "+self.familyFile)
        fnMicrographs=os.path.join(self.pickingDir,"micrographs.sel")
        if not os.path.exists(fnMicrographs):
            errors.append("Cannot find "+fnMicrographs)
        else:
            mD=xmipp.MetaData(fnMicrographs)
            if self.DoFlip and not mD.containsLabel(xmipp.MDL_CTFMODEL):
                errors.append(fnMicrographs+" does not contain CTF information for phase flipping")
        return errors

    def summary(self):
        message=[]
        families=xmipp.MetaData(self.familyFile)
        message.append(str(families.size())+" different families: ")
        for id in families:
            familyName=families.getValue(xmipp.MDL_PICKING_FAMILY,id)
            particleSize=families.getValue(xmipp.MDL_PICKING_PARTICLE_SIZE,id)
            selfile=os.path.join(self.WorkingDir,familyName+".sel")
            if os.path.exists(selfile):
                mD=xmipp.MetaData(selfile)
                message.append(familyName+" ("+str(mD.size())+" particles, size="+str(particleSize)+")")
            
        return message

    def whichTasks(self,posFile,micrographName):
        fileList=[]
        tasks=[]
        for family in self.familyList:
            familyName=family[0]
            particleSize=family[1]
            blockName="%s@%s"%(familyName,posFile)
            mD=xmipp.MetaData(blockName)
            if mD.size()>0:
                fnOut=os.path.join(self.WorkingDir,familyName,micrographName+".stk")
                tasks.append((blockName,fnOut,particleSize))
                fileList.append(fnOut)
        return (tasks,fileList)

    def visualize(self):
        summaryFile=os.path.join(self.WorkingDir,"micrographs.sel")
        if not os.path.exists(summaryFile):
            import tkMessageBox
            tkMessageBox.showerror("Error", "There is no result yet")                    
        if self.DoShowSummary:
            os.system("xmipp_visualize_preprocessing_micrographj -i "+summaryFile+" --memory 2048m &")
        if self.DoShowFamilies:
            fileList=""
            families=xmipp.MetaData(self.familyFile)
            for id in families:
                familyName=families.getValue(xmipp.MDL_PICKING_FAMILY,id)
                fileList+=os.path.join(self.WorkingDir,familyName+".sel")+" "
            if fileList!="":
                os.system("xmipp_showj -i "+fileList+" --memory 1024m &")

def phaseFlip(log,micrograph,ctf,fnOut):
    runJob(log,"xmipp_ctf_phase_flip"," -i "+micrograph+" --ctf "+ctf+" -o "+fnOut)

def extractTasks(log,micrograph,ctf,micrographToExtract,tasks,
                    doNorm,doLog,doInvert,bgRadius,doRemoveDust,dustRemovalThreshold):
    for task in tasks:
        (blockName,fnOut,particleSize)=task
        extract(log,micrographToExtract,blockName,particleSize,fnOut,
                doNorm,doLog,doInvert,bgRadius,doRemoveDust,dustRemovalThreshold)

        # Substitute the micrograph name if it comes from the flipped version
        # Add information about the ctf if available
        if micrograph!=micrographToExtract or ctf!=None:
            selfile=fnOut.replace(".stk",".sel")
            mD=xmipp.MetaData(selfile)
            if micrograph!=micrographToExtract:
                mD.setValueCol(xmipp.MDL_MICROGRAPH,micrograph)
            if ctf!=None:
                mD.setValueCol(xmipp.MDL_CTFMODEL,ctf)
            mD.write(selfile)

def extract(log,micrographToExtract,blockName,particleSize,fnOut,
            doNorm,doLog,doInvert,bgRadius,doRemoveDust,dustRemovalThreshold):
    # Extract 
    rootname=os.path.splitext(fnOut)[0]
    arguments="-i "+micrographToExtract+" --pos "+blockName+" --oroot "+rootname+" --Xdim "+str(particleSize)
    if doInvert:
        arguments+=" --invert"
    if doLog:
        arguments+=" --log"
    runJob(log,"xmipp_micrograph_scissor",arguments)
    
    # Normalize 
    if doNorm:
        if bgRadius==0:
            bgRadius=int(particleSize/2)
        arguments="-i "+fnOut+' --method Ramp --background circle '+str(bgRadius)
        if doRemoveDust:
            arguments+=' --thr_black_dust -' + str(dustRemovalThreshold)+' --thr_white_dust ' + str(dustRemovalThreshold)
        runJob(log,"xmipp_transform_normalize",arguments)

def gatherSelfiles(log,WorkingDir,familyList):
    for family in familyList:
        familyName=family[0]
        selfiles=glob.glob(os.path.join(WorkingDir,familyName)+"/*.sel")
        selfiles.sort()
        familySelfile=xmipp.MetaData()
        for selfile in selfiles:
            mD=xmipp.MetaData(selfile)
            familySelfile.unionAll(mD)
        familySelfile.write(os.path.join(WorkingDir,familyName+".sel"))

def sortImageInFamily(log,selfileRoot):
    runJob(log,"xmipp_image_sort_by_statistics","-i "+selfileRoot+".sel --multivariate --addToInput -o "+selfileRoot+"_sorted")

def avgZscore(log,WorkingDir,familyList,micrographSelfile):
    allParticles=xmipp.MetaData()
    for family in familyList:
        allParticles.unionAll(xmipp.MetaData(os.path.join(WorkingDir,family[0]+".sel")))
    mDavgZscore=xmipp.MetaData()
    mDavgZscore.aggregate(allParticles, xmipp.AGGR_AVG, xmipp.MDL_MICROGRAPH, xmipp.MDL_ZSCORE, xmipp.MDL_ZSCORE)
    oldMicrographsSel=xmipp.MetaData(micrographSelfile)
    oldMicrographsSel.removeLabel(xmipp.MDL_ZSCORE)
    newMicrographsSel=xmipp.MetaData()
    # Make copy of metadata because removeLabel leaves the values in the table
    newMicrographsSel.join(xmipp.MetaData(oldMicrographsSel),mDavgZscore,xmipp.MDL_IMAGE,xmipp.MDL_MICROGRAPH,xmipp.LEFT)
    newMicrographsSel.write(micrographSelfile)

