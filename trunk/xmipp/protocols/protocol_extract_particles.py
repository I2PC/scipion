#!/usr/bin/env xmipp_python
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
from xmipp import MetaData, AGGR_AVG, LEFT, MD_APPEND, MDL_IMAGE, MDL_CTFMODEL, MDL_MICROGRAPH, MDL_ZSCORE, getBlocksInMetaDataFile
import glob
import os
from protlib_utils import runJob
from protlib_filesystem import deleteFile, deleteDir, createLink, copyFile

class ProtExtractParticles(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.extract_particles.name, scriptname, project)
        self.Import = 'from protocol_extract_particles import *'
        # COSS: Falta tiltpairs
        self.pickingDir= getWorkingDirFromRunName(self.PickingRun)

    def defineSteps(self):
        fnMicrographsSel=os.path.join(self.pickingDir,"micrographs.sel")
        self.Db.insertStep('copyFile',source=fnMicrographsSel,dest=self.workingDirPath("micrographs.sel"))
        mD=MetaData(fnMicrographsSel)
        self.containsCTF=mD.containsLabel(MDL_CTFMODEL)

        fnExtractList=os.path.join(self.pickingDir,self.Family+"_extract_list.xmd")
        if os.path.exists(fnExtractList):
            self.Db.insertStep("createLink",source=fnExtractList,dest=self.workingDirPath(self.Family+"_extract_list.xmd"))
            micrographs=getBlocksInMetaDataFile(fnExtractList)
        else:
            self.Db.insertStep("createExtractList",Family=self.Family,fnMicrographsSel=fnMicrographsSel,pickingDir=self.pickingDir,
                               fnExtractList=self.workingDirPath(self.Family+"_extract_list.xmd"))
            fnExtractList=self.workingDirPath(self.Family+"_extract_list.xmd")
            micrographs=self.createBlocksInExtractFile(fnMicrographsSel)
        
        idMPI=self.Db.insertStep('runStepGapsMpi',passDb=True, script=self.scriptName, NumberOfMpi=self.NumberOfMpi)

        for micrograph in micrographs:
            parent_id=XmippProjectDb.FIRST_STEP
            micrographName=micrograph[4:] # Remove "mic_" from the name
            originalMicrograph,ctf=self.getMicrographInfo(micrographName, mD)
            micrographToExtract=originalMicrograph
            if self.DoFlip:
                micrographFlipped=os.path.join(self.TmpDir,micrographName+"_flipped.xmp")
                parent_id=self.Db.insertStep('phaseFlip',parent_step_id=XmippProjectDb.FIRST_STEP,
                                             verifyfiles=[micrographFlipped],execution_mode=SqliteDb.EXEC_GAP,
                                             micrograph=micrographToExtract,ctf=ctf,fnOut=micrographFlipped)
                micrographToExtract=micrographFlipped
            fnOut=self.workingDirPath(micrographName+".stk")
            parent_id=self.Db.insertStep('extractParticles',execution_mode=SqliteDb.EXEC_GAP,parent_step_id=parent_id,
                                  WorkingDir=self.WorkingDir,
                                  micrographName=micrographName,ctf=ctf,
                                  originalMicrograph=originalMicrograph,micrographToExtract=micrographToExtract,
                                  fnExtractList=fnExtractList,particleSize=self.ParticleSize,
                                  doFlip=self.DoFlip,doNorm=self.DoNorm,doLog=self.DoLog,doInvert=self.DoInvert,
                                  bgRadius=self.BackGroundRadius, doRemoveDust=self.DoRemoveDust,
                                  dustRemovalThreshold=self.DustRemovalThreshold)
            if self.DoFlip:
                self.Db.insertStep('deleteFile',execution_mode=SqliteDb.EXEC_GAP,parent_step_id=parent_id,filename=micrographToExtract,verbose=True)
        self.Db.insertStep('gatherSelfiles',parent_step_id=idMPI,WorkingDir=self.WorkingDir,family=self.Family)

        selfileRoot=self.workingDirPath(self.Family)
        fnOut=selfileRoot+"_sorted.xmd"
        self.Db.insertStep('sortImagesInFamily',verifyfiles=[fnOut],selfileRoot=selfileRoot)
        self.Db.insertStep('avgZscore',WorkingDir=self.WorkingDir,family=self.Family,
                           micrographSelfile=self.workingDirPath("micrographs.sel"))
                
    def validate(self):
        errors = []
        if self.pickingDir:
            fnMicrographs=os.path.join(self.pickingDir,"micrographs.sel")
            if not os.path.exists(fnMicrographs):
                errors.append("Cannot find "+fnMicrographs)
            else:
                mD=MetaData(fnMicrographs)
                if self.DoFlip and not mD.containsLabel(MDL_CTFMODEL):
                    errors.append(fnMicrographs+" does not contain CTF information for phase flipping")
        else:
            errors.append("Picking run is not valid")
        return errors

    def summary(self):
        message=[]
        message.append("Picking "+self.Family+" with size "+str(self.ParticleSize))
        selfile=self.workingDirPath(self.Family+".sel")
        if os.path.exists(selfile):
            mD=MetaData(selfile)
            message.append(str(mD.size())+" particles extracted")            
        return message

    def getMicrographInfo(self,micrograph,mD):
        for id in mD:
            micrographFullName=mD.getValue(MDL_IMAGE,id)
            if micrograph in micrographFullName:
                if self.containsCTF:
                    ctfFile=mD.getValue(MDL_CTFMODEL,id)
                else:
                    ctfFile=None
                return (micrographFullName,ctfFile)

    def visualize(self):
        summaryFile=self.workingDirPath("micrographs.sel")
        if not os.path.exists(summaryFile):
            import tkMessageBox
            tkMessageBox.showerror("Error", "There is no result yet")                    
        os.system("xmipp_visualize_preprocessing_micrographj -i "+summaryFile+" --memory 2048m &")
        fnSelFile=self.workingDirPath(self.Family+".sel")
        if os.path.exists(fnSelFile):
            os.system("xmipp_showj -i "+fnSelFile+" --memory 1024m &")
    
    def createBlocksInExtractFile(self,fnMicrographsSel):
        mD=MetaData(fnMicrographsSel)
        blocks=[]
        for id in mD:
            micrographFullName=mD.getValue(MDL_IMAGE,id)
            micrographName=os.path.split(os.path.split(micrographFullName)[0])[1]
            blocks.append("mic_"+micrographName)
        return blocks

def createExtractList(log,Family,fnMicrographsSel,pickingDir,fnExtractList):
    mD=MetaData(fnMicrographsSel)
    mDpos=MetaData()
    mDposAux=MetaData()
    for id in mD:
        micrographFullName=mD.getValue(MDL_IMAGE,id)
        micrographName=os.path.split(os.path.split(micrographFullName)[0])[1]
        
        fnManual=os.path.join(pickingDir,micrographName+".pos")
        fnAuto1=os.path.join(pickingDir,micrographName+"_auto.pos")
        
        mDpos.clear()
        if os.path.exists(fnManual):
            try:            
                mDpos.read(Family+"@"+fnManual)
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
        # Append alphanumeric prefix to help identifying the block 
        mDpos.write("mic_"+micrographName+"@"+fnExtractList,MD_APPEND)

def phaseFlip(log,micrograph,ctf,fnOut):
    runJob(log,"xmipp_ctf_phase_flip"," -i "+micrograph+" --ctf "+ctf+" -o "+fnOut)

def extractParticles(log,WorkingDir,micrographName,ctf,originalMicrograph,micrographToExtract,fnExtractList,
                     particleSize,doFlip,doNorm,doLog,doInvert,bgRadius,doRemoveDust,dustRemovalThreshold):
    fnBlock="mic_"+micrographName+"@"+fnExtractList
    mD=MetaData(fnBlock)
    if mD.size()==0:
        return
    
    # Extract 
    rootname=os.path.join(WorkingDir,micrographName)
    arguments="-i "+micrographToExtract+" --pos "+fnBlock+" --oroot "+rootname+" --Xdim "+str(particleSize)
    if doInvert:
        arguments+=" --invert"
    if doLog:
        arguments+=" --log"
    runJob(log,"xmipp_micrograph_scissor",arguments)
    
    # Normalize 
    if doNorm:
        if bgRadius==0:
            bgRadius=int(particleSize/2)
        arguments="-i "+rootname+'.stk --method Ramp --background circle '+str(bgRadius)
        if doRemoveDust:
            arguments+=' --thr_black_dust -' + str(dustRemovalThreshold)+' --thr_white_dust ' + str(dustRemovalThreshold)
        runJob(log,"xmipp_transform_normalize",arguments)

    # Substitute the micrograph name if it comes from the flipped version
    # Add information about the ctf if available
    if originalMicrograph!=micrographToExtract or ctf!=None:
        mD=MetaData(rootname+".sel")
        if originalMicrograph!=micrographToExtract:
            mD.setValueCol(MDL_MICROGRAPH,originalMicrograph)
        if ctf!=None:
            mD.setValueCol(MDL_CTFMODEL,ctf)
        mD.write(selfile)

def gatherSelfiles(log,WorkingDir,family):
    stackFiles=glob.glob(os.path.join(WorkingDir,"*.stk"))
    stackFiles.sort()
    familySelfile=MetaData()
    for stackFile in stackFiles:
        selfile=stackFile.replace(".stk",".sel")
        mD=MetaData(selfile)
        familySelfile.unionAll(mD)
    familySelfile.write(os.path.join(WorkingDir,family+".sel"))

def sortImagesInFamily(log,selfileRoot):
    mD=MetaData(selfileRoot+".sel")
    if mD.size()>0:
        runJob(log,"xmipp_image_sort_by_statistics","-i "+selfileRoot+".sel --multivariate --addToInput -o "+selfileRoot+"_sorted")
    else:
        createLink(log,selfileRoot+".sel",selfileRoot+"_sorted.xmd")

def avgZscore(log,WorkingDir,family,micrographSelfile):
    allParticles=MetaData(os.path.join(WorkingDir,family+".sel"))
    mDavgZscore=MetaData()
    mDavgZscore.aggregate(allParticles, AGGR_AVG, MDL_MICROGRAPH, MDL_ZSCORE, MDL_ZSCORE)
    oldMicrographsSel=MetaData(micrographSelfile)
    oldMicrographsSel.removeLabel(MDL_ZSCORE)
    newMicrographsSel=MetaData()
    # Make copy of metadata because removeLabel leaves the values in the table
    newMicrographsSel.join(MetaData(oldMicrographsSel),mDavgZscore,MDL_IMAGE,MDL_MICROGRAPH,LEFT)
    newMicrographsSel.write(micrographSelfile)

