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
                  MDL_ZSCORE, MDL_IMAGE_TILTED, MDL_ENABLED, getBlocksInMetaDataFile
import glob
from os.path import exists
from protlib_utils import runJob, printLog
from protlib_filesystem import createLink, removeFilenameExt
from protlib_gui_ext import showError

# The dictionary with specific filename templates 
# is defined here to allow use of it outside the protocol
_mic_block = 'mic_%(micName)s'
_templateDict = {
        # This templates are relative to a micrographDir
        'mic_block': _mic_block,
        'mic_block_fn': _mic_block+'@%(fn)s',
        'sorted': '%(root)s_sorted.xmd'
        }

def _getFilename(key, **args):
    return _templateDict[key] % args

class ProtExtractParticles(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.extract_particles.name, scriptname, project)
        self.Import = 'from protocol_extract_particles import *'
        # Take some parameter from previous picking protocol run
        pickingProt = self.getProtocolFromRunName(self.PickingRun)
        if hasattr(pickingProt,'TiltPairs'):
            self.TiltPairs = pickingProt.TiltPairs
        else:
            self.TiltPairs = False
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
                               fnMicrographs=self.pickingMicrographs, pickingDir=self.pickingDir,
                               fnExtractList=destFnExtractList)
            micrographs = self.createBlocksInExtractFile(self.pickingMicrographs)
        else:
            srcFnExtractList = join(self.pickingDir, self.getFilename("extractList"))
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
                micrographFlipped=join(self.TmpDir,micrographName+"_flipped.xmp")
                parent_id=self.insertParallelStep('phaseFlip',parent_step_id=parent_id,
                                                        verifyfiles=[micrographFlipped],micrograph=micrographToExtract,
                                                        ctf=ctf,fnOut=micrographFlipped)
                micrographToExtract=micrographFlipped
            fnOut = self.workingDirPath(micrographName + ".stk")
            parent_id = self.insertParallelStep('extractParticles', parent_step_id=parent_id,
                                  WorkingDir=self.WorkingDir,
                                  micrographName=micrographName,ctf=ctf,
                                  originalMicrograph=originalMicrograph,micrographToExtract=micrographToExtract,
                                  fnExtractList=destFnExtractList,particleSize=self.ParticleSize,
                                  doFlip=self.DoFlip,doNorm=self.DoNorm,doLog=self.DoLog,doInvert=self.DoInvert,
                                  bgRadius=self.BackGroundRadius, doRemoveDust=self.DoRemoveDust,
                                  dustRemovalThreshold=self.DustRemovalThreshold)
            if self.DoFlip:
                self.insertParallelStep('deleteFile', parent_step_id=parent_id,filename=micrographToExtract,verbose=True)

        if self.TiltPairs:
            self.insertStep('gatherTiltPairSelfiles',family=self.Family,
                               WorkingDir=self.WorkingDir, fnMicrographs=self.micrographs)
        else:
            self.insertStep('gatherSelfiles',WorkingDir=self.WorkingDir,family=self.Family)
            selfileRoot=self.workingDirPath(self.Family)
            fnOut = _getFilename('sorted', root=selfileRoot)
            self.insertStep('sortImagesInFamily',verifyfiles=[fnOut],selfileRoot=selfileRoot)
            self.insertStep('avgZscore',WorkingDir=self.WorkingDir,family=self.Family,
                               micrographSelfile=self.micrographs)
                
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

    def getMicrographInfo(self, micrograph, md):
        for objId in md:
            micFullName = md.getValue(MDL_MICROGRAPH,objId)
            micName = baseWithoutExt(micFullName)
            if micrograph==micName:
                if self.containsCTF:
                    ctfFile = md.getValue(MDL_CTFMODEL,objId)
                else:
                    ctfFile = None
                return (micFullName, ctfFile)
            if self.TiltPairs:
                micFullName = md.getValue(MDL_MICROGRAPH_TILTED, objId)
                micName = baseWithoutExt(micFullName)
                if micrograph == micName:
                    return (micFullName, None)

    def visualize(self):
        selfile = self.getFilename("family")
        if not exists(selfile):
            showError("Error", "There is no result yet")
        else:
            from protlib_utils import runShowJ
            if self.TiltPairs:
                runShowJ(selfile,extraParams="--mode metadata --render")
            else:
                runShowJ(selfile)
    
    def createBlocksInExtractFile(self,fnMicrographsSel):
        md = MetaData(fnMicrographsSel)
        blocks = []
        tiltPairs = md.containsLabel(MDL_MICROGRAPH_TILTED)
        def addBlock(label, objId): 
            micName = baseWithoutExt(md.getValue(label, objId)) 
            blocks.append(_getFilename('mic_block', micName=micName))
        for objId in md:
            addBlock(MDL_MICROGRAPH, objId)            
            if tiltPairs:
                addBlock(MDL_MICROGRAPH_TILTED, objId)
        return blocks

def baseWithoutExt(filename):
    return removeFilenameExt(os.path.basename(filename))

def createExtractList(log,Family,fnMicrographsSel,pickingDir,fnExtractList):
    md=MetaData(fnMicrographsSel)
    mdpos=MetaData()
    mdposAux=MetaData()
    for objId in md:
        micName = baseWithoutExt(md.getValue(MDL_MICROGRAPH, objId))
        
        fnManual = join(pickingDir, micName + ".pos")
        fnAuto1 = join(pickingDir, micName + "_auto.pos")
        
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
        fn = _getFilename('mic_block_fn', micName=micName, fn=fnExtractList)
        mdpos.write(fn, MD_APPEND)

def createExtractListTiltPairs(log, family, fnMicrographs, pickingDir, fnExtractList):
    md = MetaData(fnMicrographs)
    mdUntiltedPos = MetaData()
    mdTiltedPos = MetaData()
    for objId in md:
        umicName = baseWithoutExt(md.getValue(MDL_MICROGRAPH, objId))
        tmicName =baseWithoutExt(md.getValue(MDL_MICROGRAPH_TILTED, objId))
        
        fnUntilted = join(pickingDir,umicName + ".pos")
        fnTilted = join(pickingDir,tmicName + ".pos")
        print "Looking for "+fnUntilted+" "+fnTilted

        mdUntiltedPos.clear()
        mdTiltedPos.clear()
        if exists(fnUntilted) and exists(fnTilted):
            try:            
                mdUntiltedPos.read("%s@%s"%(family,fnUntilted))
                mdTiltedPos.read("%s@%s"%(family,fnTilted))
            except:
                pass
        # Append alphanumeric prefix to help identifying the block 
        fn = _getFilename('mic_block_fn', micName=umicName, fn=fnExtractList)
        print "Writing untilted to "+fn
        mdUntiltedPos.write(fn, MD_APPEND)
        fn =_getFilename('mic_block_fn', micName=tmicName, fn=fnExtractList)
        print "Writing tilted to "+fn
        mdTiltedPos.write(fn, MD_APPEND)

def phaseFlip(log, micrograph, ctf, fnOut):
    runJob(log,"xmipp_ctf_phase_flip"," -i %(micrograph)s --ctf %(ctf)s -o %(fnOut)" % locals())

def extractParticles(log,WorkingDir,micrographName, ctf,originalMicrograph, micrographToExtract, fnExtractList,
                     particleSize, doFlip, doNorm, doLog, doInvert, bgRadius, doRemoveDust, dustRemovalThreshold):
    print "---->Extracting",micrographName
    fnBlock = _getFilename('mic_block_fn', micName=micrographName, fn=fnExtractList)
    md = MetaData(fnBlock)
    printLog( 'Metadata block: %s' % fnBlock, log)
    if md.isEmpty():
        printLog( "EMPTY", log)
        printLog('Metadata: %s is empty, returning without extracting' % fnBlock, log)
        return
    
    # Extract 
    rootname = join(WorkingDir, micrographName)
    arguments="-i "+micrographToExtract+" --pos "+fnBlock+" --oroot "+rootname+" --Xdim "+str(particleSize)
    if doInvert:
        arguments += " --invert"
    if doLog:
        arguments += " --log"
    runJob(log,"xmipp_micrograph_scissor", arguments)
    
    # Normalize 
    if doNorm:
        if bgRadius == 0:
            bgRadius = int(particleSize/2)
        arguments = "-i "+rootname+'.stk --method Ramp --background circle '+str(bgRadius)
        if doRemoveDust:
            arguments+=' --thr_black_dust -' + str(dustRemovalThreshold)+' --thr_white_dust ' + str(dustRemovalThreshold)
        runJob(log,"xmipp_transform_normalize",arguments)

    # Substitute the micrograph name if it comes from the flipped version
    # Add information about the ctf if available
    if originalMicrograph != micrographToExtract or ctf != None:
        selfile = rootname + ".xmd"
        md = MetaData(selfile)
        if originalMicrograph != micrographToExtract:
            md.setValueCol(MDL_MICROGRAPH, originalMicrograph)
        if ctf != None:
            md.setValueCol(MDL_CTFMODEL,ctf)
        md.write(selfile)

def gatherSelfiles(log, WorkingDir, family):
    stackFiles = glob.glob(join(WorkingDir,"*.stk"))
    stackFiles.sort()
    familySelfile=MetaData()
    for stackFile in stackFiles:
        selfile=stackFile.replace(".stk",".xmd")
        md=MetaData(selfile)
        familySelfile.unionAll(md)
    familySelfile.write(join(WorkingDir,family + ".xmd"))

def gatherTiltPairSelfiles(log, family, WorkingDir, fnMicrographs):
    mdPairs = MetaData(fnMicrographs)
    mdUntilted = MetaData()
    mdTilted = MetaData()
    for objId in mdPairs:
        umicName = baseWithoutExt(mdPairs.getValue(MDL_MICROGRAPH, objId))        
        fnUntilted = join(WorkingDir, umicName + ".xmd")
        # Check if there are picked particles in this micrographs
        if exists(fnUntilted):
            mdUntilted.unionAll(MetaData(fnUntilted))            
            tmicName = baseWithoutExt(mdPairs.getValue(MDL_MICROGRAPH_TILTED,objId))
            mdTilted.unionAll(MetaData(join(WorkingDir, tmicName + ".xmd")))
        
    mdUntilted.write(join(WorkingDir, "%s_untilted.xmd" % family))
    mdTilted.write(join(WorkingDir, "%s_tilted.xmd" % family))

    mdTiltedPairs = MetaData()
    mdTiltedPairs.setColumnValues(MDL_IMAGE, mdUntilted.getColumnValues(MDL_IMAGE))
    mdTiltedPairs.setColumnValues(MDL_IMAGE_TILTED, mdTilted.getColumnValues(MDL_IMAGE))
    mdTiltedPairs.write(join(WorkingDir,"%s.xmd" % family))

def sortImagesInFamily(log, selfileRoot):
    fn = selfileRoot + '.xmd'
    fnSorted = _getFilename('sorted', root=selfileRoot)
    md = MetaData(fn)
    if not md.isEmpty():
        runJob(log, "xmipp_image_sort_by_statistics","-i %(fn)s --multivariate --addToInput -o %(fnSorted)s" % locals())
    else:
        createLink(log, fn, fnSorted)

def avgZscore(log,WorkingDir,family,micrographSelfile):
    allParticles = MetaData(join(WorkingDir, family + ".xmd"))
    mdavgZscore = MetaData()
    mdavgZscore.aggregate(allParticles, AGGR_AVG, MDL_MICROGRAPH, MDL_ZSCORE, MDL_ZSCORE)
    oldMicrographsSel = MetaData(micrographSelfile)
    oldMicrographsSel.removeLabel(MDL_ZSCORE)
    newMicrographsSel = MetaData()
    # Make copy of metadata because removeLabel leaves the values in the table
    newMicrographsSel.join(oldMicrographsSel,mdavgZscore,MDL_MICROGRAPH,MDL_MICROGRAPH,LEFT)
    newMicrographsSel.write(micrographSelfile)

