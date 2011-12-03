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
from os.path import exists
from protlib_utils import runJob
from protlib_filesystem import createLink

class ProtExtractParticles(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.extract_particles.name, scriptname, project)
        self.Import = 'from protocol_extract_particles import *'
        # Take some parameter from previous picking protocol run
        pickingProt = self.getProtocolFromRunName(self.PickingRun)
        self.TiltPairs = pickingProt.TiltPairs
        self.pickingDir = pickingProt.WorkingDir
        if self.TiltPairs:
            self.pickingMicrographs = pickingProt.getFilename('tiltedPairs')
        else:
            self.pickingMicrographs = pickingProt.getFilename("micrographs")
        self.micrographs = self.getEquivalentFilename(pickingProt, self.pickingMicrographs)

    def createFilenameTemplates(self):
        return {
                'extractList': self.workingDirPath("%(Family)s_extract_list.xmd"),
                'family': self.workingDirPath("%(Family)s.xmd")
                }

    def defineSteps(self):
        self.insertStep('copyFile', source=self.pickingMicrographs, dest=self.micrographs)
        md = MetaData(self.pickingMicrographs)
        self.containsCTF = md.containsLabel(MDL_CTFMODEL)

        destFnExtractList = self.getFilename("extractList")
        if self.TiltPairs:
            self.insertStep("createExtractListTiltPairs", family=self.Family,
                               fnMicrographsSel=self.pickingMicrographs, pickingDir=self.pickingDir,
                               fnExtractList=destFnExtractList)
            micrographs = self.createBlocksInExtractFile(self.pickingMicrographs)
        else:
            srcFnExtractList = os.path.join(self.pickingDir,self.getFilename("extractList"))
            if exists(srcFnExtractList):
                self.insertStep("createLink",source=srcFnExtractList,dest=destFnExtractList)
                micrographs = getBlocksInMetaDataFile(srcFnExtractList)
            else:
                self.insertStep("createExtractList",Family=self.Family,fnMicrographsSel=self.pickingMicrographs,pickingDir=self.pickingDir,
                                   fnExtractList=destFnExtractList)
                micrographs = self.createBlocksInExtractFile(self.pickingMicrographs)
        
        for micrograph in micrographs:
            parent_id = XmippProjectDb.FIRST_STEP
            micrographName = micrograph[4:] # Remove "mic_" from the name
            originalMicrograph, ctf = self.getMicrographInfo(micrographName, md)
            micrographToExtract = originalMicrograph
            if self.DoFlip:
                micrographFlipped=os.path.join(self.TmpDir,micrographName+"_flipped.xmp")
                parent_id=self.insertParallelStep('phaseFlip',parent_step_id=parent_id,
                                                        verifyfiles=[micrographFlipped],micrograph=micrographToExtract,
                                                        ctf=ctf,fnOut=micrographFlipped)
                micrographToExtract=micrographFlipped
            fnOut = self.workingDirPath(micrographName + ".stk")
            parent_id=self.insertParallelStep('extractParticles', parent_step_id=parent_id,
                                  WorkingDir=self.WorkingDir,
                                  micrographName=micrographName,ctf=ctf,
                                  originalMicrograph=originalMicrograph,micrographToExtract=micrographToExtract,
                                  fnExtractList=destFnExtractList,particleSize=self.ParticleSize,
                                  doFlip=self.DoFlip,doNorm=self.DoNorm,doLog=self.DoLog,doInvert=self.DoInvert,
                                  bgRadius=self.BackGroundRadius, doRemoveDust=self.DoRemoveDust,
                                  dustRemovalThreshold=self.DustRemovalThreshold)
            if self.DoFlip:
                self.Db.insertStep('deleteFile',execution_mode=SqliteDb.EXEC_PARALLEL,parent_step_id=parent_id,filename=micrographToExtract,verbose=True)

        if self.TiltPairs:
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
            if not exists(self.pickingMicrographs):
                errors.append("Cannot find "+self.pickingMicrographs)
            else:
                md=MetaData(self.pickingMicrographs)
                if self.DoFlip and not md.containsLabel(MDL_CTFMODEL):
                    errors.append(self.pickingMicrographs+" does not contain CTF information for phase flipping")
        else:
            errors.append("Picking run is not valid")
        return errors

    def summary(self):
        message=[]
        familyFn = self.getFilename("family")
        if self.TiltPairs:
            part1 = "tilt pairs"
            part2 = "particle pairs"
        else:
            part1 = "family <%s>" % self.Family
            part2 = "particles"
        message.append("Picking %s with size  <%d>" % (part1, self.ParticleSize))
    
        if exists(familyFn):
            md = MetaData(familyFn)
            message.append("%d %s extracted" % (md.size(), part2))
        return message

    def getMicrographInfo(self,micrograph,md):
        for id in md:
            micrographFullName=md.getValue(MDL_MICROGRAPH,id)
            micrographName=os.path.splitext(os.path.split(micrographFullName)[1])[0]
            if micrograph==micrographName:
                if self.containsCTF:
                    ctfFile=md.getValue(MDL_CTFMODEL,id)
                else:
                    ctfFile=None
                return (micrographFullName,ctfFile)
            if self.TiltPairs:
                micrographFullName=md.getValue(MDL_MICROGRAPH_TILTED,id)
                micrographName=os.path.splitext(os.path.split(micrographFullName)[1])[0]
                if micrograph==micrographName:
                    return (micrographFullName,None)

    def visualize(self):
        selfile=self.getFilename("family")
        if not exists(selfile):
            import tkMessageBox
            tkMessageBox.showerror("Error", "There is no result yet")
        if not self.TiltPairs:
            summaryFile=self.getFilename("micrographs")
            if exists(summaryFile):                    
                os.system("xmipp_visualize_preprocessing_micrographj -i "+summaryFile+" --memory 2048m &")
        if exists(selfile):
            from protlib_utils import runShowJ
            runShowJ(selfile, memory="1024m")
    
    def createBlocksInExtractFile(self,fnMicrographsSel):
        md=MetaData(fnMicrographsSel)
        blocks=[]
        tiltPairs=md.containsLabel(MDL_MICROGRAPH_TILTED)
        for id in md:
            micrographFullName=md.getValue(MDL_MICROGRAPH,id)
            micrographName=os.path.splitext(os.path.split(micrographFullName)[1])[0]
            blocks.append("mic_"+micrographName)
            if tiltPairs:
                micrographFullName=md.getValue(MDL_MICROGRAPH_TILTED,id)
                micrographName=os.path.splitext(os.path.split(micrographFullName)[1])[0]
                blocks.append("mic_"+micrographName)
        return blocks

def createExtractList(log,Family,fnMicrographsSel,pickingDir,fnExtractList):
    md=MetaData(fnMicrographsSel)
    mdpos=MetaData()
    mdposAux=MetaData()
    for id in md:
        micrographFullName=md.getValue(MDL_MICROGRAPH,id)
        micrographName=os.path.splitext(os.path.split(micrographFullName)[1])[0]
        
        fnManual=os.path.join(pickingDir,micrographName+".pos")
        fnAuto1=os.path.join(pickingDir,micrographName+"_auto.pos")
        
        mdpos.clear()
        if exists(fnManual):
            try:            
                mdpos.read(Family+"@"+fnManual)
            except:
                pass
        if exists(fnAuto1):
            try:
                mdposAux.read(Family+"@"+fnAuto1)
                mdposAux.removeDisabled();
                mdposAux.removeLabel(MDL_ENABLED)
                mdpos.unionAll(mdposAux) 
            except:
                pass
        # Append alphanumeric prefix to help identifying the block 
        mdpos.write("mic_"+micrographName+"@"+fnExtractList,MD_APPEND)

def createExtractListTiltPairs(log,family,fnMicrographsSel,pickingDir,fnExtractList):
    md=MetaData(fnMicrographsSel)
    mdUntiltedPos=MetaData()
    mdTiltedPos=MetaData()
    for id in md:
        micrographFullName=md.getValue(MDL_MICROGRAPH,id)
        untiltedMicrographName=os.path.splitext(os.path.split(micrographFullName)[1])[0]
        micrographFullName=md.getValue(MDL_MICROGRAPH_TILTED,id)
        tiltedMicrographName=os.path.splitext(os.path.split(micrographFullName)[1])[0]
        
        fnUntilted=os.path.join(pickingDir,untiltedMicrographName+".pos")
        fnTilted=os.path.join(pickingDir,tiltedMicrographName+".pos")

        mdUntiltedPos.clear()
        mdTiltedPos.clear()
        if exists(fnUntilted) and exists(fnTilted):
            try:            
                mdUntiltedPos.read("%s@%s"%(family,fnUntilted))
                mdTiltedPos.read("%s@%s"%(family,fnTilted))
            except:
                pass
        # Append alphanumeric prefix to help identifying the block 
        mdUntiltedPos.write("mic_"+untiltedMicrographName+"@"+fnExtractList,MD_APPEND)
        mdTiltedPos.write("mic_"+tiltedMicrographName+"@"+fnExtractList,MD_APPEND)

def phaseFlip(log,micrograph,ctf,fnOut):
    runJob(log,"xmipp_ctf_phase_flip"," -i "+micrograph+" --ctf "+ctf+" -o "+fnOut)

def extractParticles(log,WorkingDir,micrographName,ctf,originalMicrograph,micrographToExtract,fnExtractList,
                     particleSize,doFlip,doNorm,doLog,doInvert,bgRadius,doRemoveDust,dustRemovalThreshold):
    fnBlock="mic_"+micrographName+"@"+fnExtractList
    md=MetaData(fnBlock)
    if md.size()==0:
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
        md=MetaData(selfile)
        if originalMicrograph!=micrographToExtract:
            md.setValueCol(MDL_MICROGRAPH,originalMicrograph)
        if ctf!=None:
            md.setValueCol(MDL_CTFMODEL,ctf)
        md.write(selfile)

def gatherSelfiles(log,WorkingDir,family):
    stackFiles=glob.glob(os.path.join(WorkingDir,"*.stk"))
    stackFiles.sort()
    familySelfile=MetaData()
    for stackFile in stackFiles:
        selfile=stackFile.replace(".stk",".xmd")
        md=MetaData(selfile)
        familySelfile.unionAll(md)
    familySelfile.write(os.path.join(WorkingDir,family+".xmd"))

def gatherTiltPairSelfiles(log,family,WorkingDir,fnMicrographsSel):
    MDpairs=MetaData(fnMicrographsSel)
    untiltedSelfile=MetaData()
    tiltedSelfile=MetaData()
    for id in MDpairs:
        micrographFullName=MDpairs.getValue(MDL_MICROGRAPH,id)
        untiltedMicrographName=os.path.splitext(os.path.split(micrographFullName)[1])[0]
        fnUntiltedSel=os.path.join(WorkingDir,untiltedMicrographName+".xmd")
        if not exists(fnUntiltedSel):
            continue
        mduntilted=MetaData(fnUntiltedSel)
        untiltedSelfile.unionAll(mduntilted)
        
        micrographFullName=MDpairs.getValue(MDL_MICROGRAPH_TILTED,id)
        tiltedMicrographName=os.path.splitext(os.path.split(micrographFullName)[1])[0]
        mdtilted=MetaData(os.path.join(WorkingDir,tiltedMicrographName+".xmd"))
        tiltedSelfile.unionAll(mdtilted)
        
    untiltedSelfile.write(os.path.join(WorkingDir,"%s_untilted.xmd"%family))
    tiltedSelfile.write(os.path.join(WorkingDir,"%s_tilted.xmd"%family))

    pairsSelfile=MetaData()
    pairsSelfile.setColumnValues(MDL_IMAGE,untiltedSelfile.getColumnValues(MDL_IMAGE))
    pairsSelfile.setColumnValues(MDL_IMAGE_TILTED,tiltedSelfile.getColumnValues(MDL_IMAGE))
    pairsSelfile.write(os.path.join(WorkingDir,"%s.xmd"%family))

def sortImagesInFamily(log,selfileRoot):
    md=MetaData(selfileRoot+".xmd")
    if md.size()>0:
        runJob(log,"xmipp_image_sort_by_statistics","-i "+selfileRoot+".xmd --multivariate --addToInput -o "+selfileRoot+"_sorted")
    else:
        createLink(log,selfileRoot+".xmd",selfileRoot+"_sorted.xmd")

def avgZscore(log,WorkingDir,family,micrographSelfile):
    allParticles=MetaData(os.path.join(WorkingDir,family+".xmd"))
    mdavgZscore=MetaData()
    mdavgZscore.aggregate(allParticles, AGGR_AVG, MDL_MICROGRAPH, MDL_ZSCORE, MDL_ZSCORE)
    oldMicrographsSel=MetaData(micrographSelfile)
    oldMicrographsSel.removeLabel(MDL_ZSCORE)
    newMicrographsSel=MetaData()
    # Make copy of metadata because removeLabel leaves the values in the table
    newMicrographsSel.join(MetaData(oldMicrographsSel),mdavgZscore,MDL_IMAGE,MDL_MICROGRAPH,LEFT)
    newMicrographsSel.write(micrographSelfile)

