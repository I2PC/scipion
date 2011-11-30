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
from xmipp import MetaData, AGGR_AVG, LEFT, MD_APPEND, MDL_IMAGE, MDL_CTFMODEL, MDL_MICROGRAPH, MDL_MICROGRAPH_TILTED, \
                  MDL_ZSCORE, MDL_IMAGE_TILTED, getBlocksInMetaDataFile
import glob
import os
from protlib_utils import runJob
from protlib_filesystem import deleteFile, deleteDir, createLink, copyFile

class ProtExtractParticles(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.extract_particles.name, scriptname, project)
        self.Import = 'from protocol_extract_particles import *'
        self.pickingDir= getWorkingDirFromRunName(self.PickingRun)

        protPicking = getProtocolFromModule(getScriptFromRunName(self.PickingRun),self.project)
        self.fnMicrographsSel=protPicking.getFilename("micrographs")
        self.tiltPairs =self.fnMicrographsSel.endswith("tilted_pairs.xmd")
        if self.tiltPairs:
            self.FilenamesDict["micrographs"]=self.workingDirPath('tilted_pairs.xmd')

    def createFilenameTemplates(self):
        return {
                'micrographs': self.workingDirPath('micrographs.xmd'),
                'extractList': "%s_extract_list.xmd"%self.Family,
                'family': "%s.xmd"%self.Family
                }

    def defineSteps(self):
        self.Db.insertStep('copyFile',source=self.fnMicrographsSel,dest=self.getFilename("micrographs"))
        mD=MetaData(self.fnMicrographsSel)
        self.containsCTF=mD.containsLabel(MDL_CTFMODEL)

        destFnExtractList=self.workingDirPath(self.getFilename("extractList"))
        if self.tiltPairs:
            self.Db.insertStep("createExtractListTiltPairs",family=self.Family,
                               fnMicrographsSel=self.fnMicrographsSel,pickingDir=self.pickingDir,
                               fnExtractList=fnExtractList)
            micrographs=self.createBlocksInExtractFile(self.fnMicrographsSel)
        else:
            srcFnExtractList=os.path.join(self.pickingDir,self.getFilename("extractList"))
            if os.path.exists(srcFnExtractList):
                self.Db.insertStep("createLink",source=srcFnExtractList,dest=destFnExtractList)
                micrographs=getBlocksInMetaDataFile(srcFnExtractList)
            else:
                self.Db.insertStep("createExtractList",Family=self.Family,fnMicrographsSel=self.fnMicrographsSel,pickingDir=self.pickingDir,
                                   fnExtractList=destFnExtractList)
                micrographs=self.createBlocksInExtractFile(self.fnMicrographsSel)
        
        for micrograph in micrographs:
            parent_id=XmippProjectDb.FIRST_STEP
            micrographName=micrograph[4:] # Remove "mic_" from the name
            originalMicrograph,ctf=self.getMicrographInfo(micrographName, mD)
            micrographToExtract=originalMicrograph
            if self.DoFlip:
                micrographFlipped=os.path.join(self.TmpDir,micrographName+"_flipped.xmp")
                parent_id=self.insertParallelStep('phaseFlip',parent_step_id=parent_id,
                                                        verifyfiles=[micrographFlipped],micrograph=micrographToExtract,
                                                        ctf=ctf,fnOut=micrographFlipped)
                micrographToExtract=micrographFlipped
            fnOut=self.workingDirPath(micrographName+".stk")
            parent_id=self.insertParallelStep('extractParticles',parent_step_id=parent_id,
                                  WorkingDir=self.WorkingDir,
                                  micrographName=micrographName,ctf=ctf,
                                  originalMicrograph=originalMicrograph,micrographToExtract=micrographToExtract,
                                  fnExtractList=destFnExtractList,particleSize=self.ParticleSize,
                                  doFlip=self.DoFlip,doNorm=self.DoNorm,doLog=self.DoLog,doInvert=self.DoInvert,
                                  bgRadius=self.BackGroundRadius, doRemoveDust=self.DoRemoveDust,
                                  dustRemovalThreshold=self.DustRemovalThreshold)
            if self.DoFlip:
                self.Db.insertStep('deleteFile',execution_mode=SqliteDb.EXEC_PARALLEL,parent_step_id=parent_id,filename=micrographToExtract,verbose=True)

        if self.tiltPairs:
            self.Db.insertStep('gatherTiltPairSelfiles',family=self.Family,
                               WorkingDir=self.WorkingDir,fnMicrographsSel=self.getFilename("micrographs"))
        else:
            self.Db.insertStep('gatherSelfiles',WorkingDir=self.WorkingDir,family=self.Family)
            selfileRoot=self.workingDirPath(self.Family)
            fnOut=selfileRoot+"_sorted.xmd"
            self.Db.insertStep('sortImagesInFamily',verifyfiles=[fnOut],selfileRoot=selfileRoot)
            self.Db.insertStep('avgZscore',WorkingDir=self.WorkingDir,family=self.Family,
                               micrographSelfile=self.getFilename("micrographs"))
                
    def validate(self):
        errors = []
        if self.pickingDir:
            if not os.path.exists(self.fnMicrographsSel):
                errors.append("Cannot find "+self.fnMicrographsSel)
            else:
                mD=MetaData(self.fnMicrographsSel)
                if self.DoFlip and not mD.containsLabel(MDL_CTFMODEL):
                    errors.append(self.fnMicrographsSel+" does not contain CTF information for phase flipping")
        else:
            errors.append("Picking run is not valid")
        return errors

    def summary(self):
        message=[]
        selfile=self.getFilename("family")
        if self.tiltPairs:
            message.append("Picking tilt pairs with size "+str(self.ParticleSize))
        else:
            message.append("Picking "+self.Family+" with size "+str(self.ParticleSize))
        if os.path.exists(selfile):
            mD=MetaData(selfile)
            if self.tiltPairs:
                message.append(str(mD.size())+" particle pairs extracted")            
            else:
                message.append(str(mD.size())+" particles extracted")            
        return message

    def getMicrographInfo(self,micrograph,mD):
        for id in mD:
            micrographFullName=mD.getValue(MDL_MICROGRAPH,id)
            micrographName=os.path.splitext(os.path.split(micrographFullName)[1])[0]
            if micrograph==micrographName:
                if self.containsCTF:
                    ctfFile=mD.getValue(MDL_CTFMODEL,id)
                else:
                    ctfFile=None
                return (micrographFullName,ctfFile)
            if self.tiltPairs:
                micrographFullName=mD.getValue(MDL_MICROGRAPH_TILTED,id)
                micrographName=os.path.splitext(os.path.split(micrographFullName)[1])[0]
                if micrograph==micrographName:
                    return (micrographFullName,None)

    def visualize(self):
        selfile=self.getFilename("family")
        if not os.path.exists(selfile):
            import tkMessageBox
            tkMessageBox.showerror("Error", "There is no result yet")
        if not self.tiltPairs:
            summaryFile=self.getFilename("micrographs")
            if os.path.exists(summaryFile):                    
                os.system("xmipp_visualize_preprocessing_micrographj -i "+summaryFile+" --memory 2048m &")
        if os.path.exists(selfile):
            from protlib_utils import runShowJ
            runShowJ(selfile, memory="1024m")
    
    def createBlocksInExtractFile(self,fnMicrographsSel):
        mD=MetaData(fnMicrographsSel)
        blocks=[]
        tiltPairs=mD.containsLabel(MDL_MICROGRAPH_TILTED)
        for id in mD:
            micrographFullName=mD.getValue(MDL_MICROGRAPH,id)
            micrographName=os.path.splitext(os.path.split(micrographFullName)[1])[0]
            blocks.append("mic_"+micrographName)
            if tiltPairs:
                micrographFullName=mD.getValue(MDL_MICROGRAPH_TILTED,id)
                micrographName=os.path.splitext(os.path.split(micrographFullName)[1])[0]
                blocks.append("mic_"+micrographName)
        return blocks

def createExtractList(log,Family,fnMicrographsSel,pickingDir,fnExtractList):
    mD=MetaData(fnMicrographsSel)
    mDpos=MetaData()
    mDposAux=MetaData()
    for id in mD:
        micrographFullName=mD.getValue(MDL_MICROGRAPH,id)
        micrographName=os.path.splitext(os.path.split(micrographFullName)[1])[0]
        
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

def createExtractListTiltPairs(log,family,fnMicrographsSel,pickingDir,fnExtractList):
    mD=MetaData(fnMicrographsSel)
    mDUntiltedPos=MetaData()
    mDTiltedPos=MetaData()
    for id in mD:
        micrographFullName=mD.getValue(MDL_MICROGRAPH,id)
        untiltedMicrographName=os.path.splitext(os.path.split(micrographFullName)[1])[0]
        micrographFullName=mD.getValue(MDL_MICROGRAPH_TILTED,id)
        tiltedMicrographName=os.path.splitext(os.path.split(micrographFullName)[1])[0]
        
        fnUntilted=os.path.join(pickingDir,untiltedMicrographName+".pos")
        fnTilted=os.path.join(pickingDir,tiltedMicrographName+".pos")

        mDUntiltedPos.clear()
        mDTiltedPos.clear()
        if os.path.exists(fnUntilted) and os.path.exists(fnTilted):
            try:            
                mDUntiltedPos.read("%s@%s"%(family,fnUntilted))
                mDTiltedPos.read("%s@%s"%(family,fnTilted))
            except:
                pass
        # Append alphanumeric prefix to help identifying the block 
        mDUntiltedPos.write("mic_"+untiltedMicrographName+"@"+fnExtractList,MD_APPEND)
        mDTiltedPos.write("mic_"+tiltedMicrographName+"@"+fnExtractList,MD_APPEND)

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
        selfile=rootname+".xmd"
        mD=MetaData(selfile)
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
        selfile=stackFile.replace(".stk",".xmd")
        mD=MetaData(selfile)
        familySelfile.unionAll(mD)
    familySelfile.write(os.path.join(WorkingDir,family+".xmd"))

def gatherTiltPairSelfiles(log,family,WorkingDir,fnMicrographsSel):
    MDpairs=MetaData(fnMicrographsSel)
    untiltedSelfile=MetaData()
    tiltedSelfile=MetaData()
    for id in MDpairs:
        micrographFullName=MDpairs.getValue(MDL_MICROGRAPH,id)
        untiltedMicrographName=os.path.splitext(os.path.split(micrographFullName)[1])[0]
        fnUntiltedSel=os.path.join(WorkingDir,untiltedMicrographName+".xmd")
        if not os.path.exists(fnUntiltedSel):
            continue
        mDuntilted=MetaData(fnUntiltedSel)
        untiltedSelfile.unionAll(mDuntilted)
        
        micrographFullName=MDpairs.getValue(MDL_MICROGRAPH_TILTED,id)
        tiltedMicrographName=os.path.splitext(os.path.split(micrographFullName)[1])[0]
        mDtilted=MetaData(os.path.join(WorkingDir,tiltedMicrographName+".xmd"))
        tiltedSelfile.unionAll(mDtilted)
        
    untiltedSelfile.write(os.path.join(WorkingDir,"%s_untilted.xmd"%family))
    tiltedSelfile.write(os.path.join(WorkingDir,"%s_tilted.xmd"%family))

    pairsSelfile=MetaData()
    pairsSelfile.setColumnValues(MDL_IMAGE,untiltedSelfile.getColumnValues(MDL_IMAGE))
    pairsSelfile.setColumnValues(MDL_IMAGE_TILTED,tiltedSelfile.getColumnValues(MDL_IMAGE))
    pairsSelfile.write(os.path.join(WorkingDir,"%s.xmd"%family))

def sortImagesInFamily(log,selfileRoot):
    mD=MetaData(selfileRoot+".xmd")
    if mD.size()>0:
        runJob(log,"xmipp_image_sort_by_statistics","-i "+selfileRoot+".xmd --multivariate --addToInput -o "+selfileRoot+"_sorted")
    else:
        createLink(log,selfileRoot+".xmd",selfileRoot+"_sorted.xmd")

def avgZscore(log,WorkingDir,family,micrographSelfile):
    allParticles=MetaData(os.path.join(WorkingDir,family+".xmd"))
    mDavgZscore=MetaData()
    mDavgZscore.aggregate(allParticles, AGGR_AVG, MDL_MICROGRAPH, MDL_ZSCORE, MDL_ZSCORE)
    oldMicrographsSel=MetaData(micrographSelfile)
    oldMicrographsSel.removeLabel(MDL_ZSCORE)
    newMicrographsSel=MetaData()
    # Make copy of metadata because removeLabel leaves the values in the table
    newMicrographsSel.join(MetaData(oldMicrographsSel),mDavgZscore,MDL_IMAGE,MDL_MICROGRAPH,LEFT)
    newMicrographsSel.write(micrographSelfile)

