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
from xmipp import MetaData, AGGR_AVG, LEFT, MD_APPEND, MDL_IMAGE, MDL_CTF_MODEL, MDL_MICROGRAPH, MDL_MICROGRAPH_TILTED, \
                  MDL_MICROGRAPH_ORIGINAL, MDL_MICROGRAPH_TILTED_ORIGINAL, MDL_ZSCORE, MDL_IMAGE_TILTED, MDL_ENABLED, \
                  MDL_SAMPLINGRATE, MDL_SAMPLINGRATE_ORIGINAL, getBlocksInMetaDataFile
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
        'extract_list':  join('%(ExtraDir)s', "extract_list.xmd")      
        }

def _getFilename(key, **args):
    return _templateDict[key] % args

class DownsamplingMode:
    SameAsPicking, SameAsOriginal, NewDownsample = range(3)
    
class ProtExtractParticles(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.extract_particles.name, scriptname, project)
        self.Import = 'from protocol_extract_particles import *'
        # Take some parameter from previous picking protocol run
        self.setPreviousRun(self.PickingRun)
        self.pickingDir = self.PrevRun.WorkingDir
        self.inputFilename('acquisition')
        self.inputProperty('TiltPairs', 'MicrographsMd')
        self.micrographs = self.getEquivalentFilename(self.PrevRun, self.MicrographsMd)
        
    def setSamplingMode(self):
        md = MetaData(self.PrevRun.getFilename('acquisition'))
        objId = md.firstObject()
        self.TsOriginal= self.TsInput = md.getValue(MDL_SAMPLINGRATE, objId)        
        if md.containsLabel(MDL_SAMPLINGRATE_ORIGINAL):
            self.TsOriginal = md.getValue(MDL_SAMPLINGRATE_ORIGINAL, objId)

        if self.DownsampleType == "same as picking":
            self.TsFinal = md.getValue(MDL_SAMPLINGRATE, objId)
        elif self.DownsampleType == "original":
            self.TsFinal = self.TsOriginal
        else:
            self.TsFinal = self.TsOriginal*self.DownsampleFactor
        
        if abs(self.TsFinal-self.TsInput)<0.001:
            self.downsamplingMode=DownsamplingMode.SameAsPicking
        elif abs(self.TsFinal-self.TsOriginal)<0.001:
            self.downsamplingMode=DownsamplingMode.SameAsOriginal
        else:
            self.downsamplingMode=DownsamplingMode.NewDownsample

    def defineSteps(self):
        self.setSamplingMode()
        filesToCopy = [self.MicrographsMd, self.Input['acquisition']]
        self.insertImportOfFiles(filesToCopy, copy=True)
        self.insertStep('createDir',path=self.ExtraDir)

        # Update sampling rate in 'acquisition_info.xmd' if necessary
        if self.downsamplingMode != DownsamplingMode.SameAsPicking:
            self.insertStep('updateSampling', SamplingMd=self.getFilename('acquisition'), 
                            TsOriginal=self.TsOriginal, Ts=self.TsFinal, TsMode=self.downsamplingMode)
        
        md = MetaData(self.MicrographsMd)
        self.containsCTF = md.containsLabel(MDL_CTF_MODEL)
        # Create or look for the extract list
        destFnExtractList = _getFilename('extract_list', ExtraDir=self.ExtraDir)
        if self.TiltPairs:
            self.insertStep("createExtractListTiltPairs", family=self.Family,
                               fnMicrographs=self.MicrographsMd, pickingDir=self.pickingDir,
                               fnExtractList=destFnExtractList)
            micrographs = self.createBlocksInExtractFile(self.MicrographsMd)
        else:
            srcFnExtractList = self.PrevRun.getFilename('extract_list')
            if exists(srcFnExtractList):
                self.insertCopyFile(srcFnExtractList, destFnExtractList)
                micrographs = getBlocksInMetaDataFile(srcFnExtractList)
            else:
                self.insertStep("createExtractList",Family=self.Family,fnMicrographsSel=self.MicrographsMd,pickingDir=self.pickingDir,
                                   fnExtractList=destFnExtractList)
                micrographs = self.createBlocksInExtractFile(self.MicrographsMd)

        # Process each micrograph                        
        for micrograph in micrographs:
            parent_id = XmippProjectDb.FIRST_STEP
            micrographName = micrograph[4:] # Remove "mic_" from the name
            fullMicrographName, originalMicrograph, ctf = self.getMicrographInfo(micrographName, md)
            
            # Downsampling?
            if self.downsamplingMode==DownsamplingMode.SameAsPicking:
                micrographToExtract = fullMicrographName
            elif self.downsamplingMode==DownsamplingMode.SameAsOriginal:
                micrographToExtract = originalMicrograph
            else:
                micrographDownsampled = self.tmpPath(micrographName+"_downsampled.xmp")
                parent_id=self.insertParallelRunJobStep("xmipp_transform_downsample", "-i %s -o %s --step %f --method fourier" % \
                                                        (originalMicrograph,micrographDownsampled,self.DownsampleFactor),
                                                        parent_step_id=parent_id)
                micrographToExtract=micrographDownsampled

            # Removing dust?
            if self.DoRemoveDust:
                micrographNoDust = self.tmpPath(micrographName+"_noDust.xmp")
                threshold=self.DustRemovalThreshold
                args=" -i %(micrographToExtract)s -o %(micrographNoDust)s --bad_pixels outliers %(threshold)f" % locals()
                parent_id=self.insertParallelRunJobStep("xmipp_transform_filter", args, parent_step_id=parent_id)
                micrographToExtract=micrographNoDust
            
            # Flipping?
            if self.DoFlip:
                micrographFlipped = self.tmpPath(micrographName+"_flipped.xmp")
                args=" -i %(micrographToExtract)s --ctf %(ctf)s -o %(micrographFlipped)s" % locals()
                if self.downsamplingMode!=DownsamplingMode.SameAsOriginal:
                    args+=" --downsampling "+str(self.TsFinal/self.TsOriginal)
                parent_id=self.insertParallelRunJobStep("xmipp_ctf_phase_flip", args, parent_step_id=parent_id)
                micrographToExtract=micrographFlipped
            
            # Actually extract
            fnOut = self.workingDirPath(micrographName + ".stk")
            parent_id = self.insertParallelStep('extractParticles', parent_step_id=parent_id,
                                  ExtraDir=self.ExtraDir,
                                  micrographName=micrographName,ctf=ctf,
                                  fullMicrographName=fullMicrographName,originalMicrograph=originalMicrograph,
                                  micrographToExtract=micrographToExtract,
                                  TsFinal=self.TsFinal, TsInput=self.TsInput, downsamplingMode=self.downsamplingMode,
                                  fnExtractList=destFnExtractList,particleSize=self.ParticleSize,
                                  doFlip=self.DoFlip,doNorm=self.DoNorm,doInvert=self.DoInvert,
                                  bgRadius=self.BackGroundRadius)
            if self.DoRemoveDust:
                self.insertParallelStep('deleteFile', parent_step_id=parent_id,filename=micrographNoDust,verbose=True)
            if self.downsamplingMode==DownsamplingMode.NewDownsample:
                self.insertParallelStep('deleteFile', parent_step_id=parent_id,filename=micrographDownsampled,verbose=True)
            if self.DoFlip:
                self.insertParallelStep('deleteFile', parent_step_id=parent_id,filename=micrographFlipped,verbose=True)

        # Gather results
        if self.TiltPairs:
            self.insertStep('gatherTiltPairSelfiles',
                               WorkingDir=self.WorkingDir, ExtraDir=self.ExtraDir, fnMicrographs=self.micrographs)
        else:
            self.insertStep('gatherSelfiles',WorkingDir=self.WorkingDir,ExtraDir=self.ExtraDir)
            self.insertStep('sortImagesInFamily',WorkingDir=self.WorkingDir)
            self.insertStep('avgZscore',WorkingDir=self.WorkingDir,
                               micrographSelfile=self.micrographs)

    def validate(self):
        errors = []
        if self.pickingDir:
            if not exists(self.MicrographsMd):
                errors.append("Cannot find: \n<%s>" % self.MicrographsMd)
            else:
                md = MetaData(self.MicrographsMd)
                if self.DoFlip and not md.containsLabel(MDL_CTF_MODEL):
                    errors.append("Micrographs metadata: <%s>\n does not contain CTF information for phase flipping" % self.MicrographsMd)
        else:
            errors.append("Picking run is not valid")
        if self.DownsampleFactor<1:
            errors.append("Downsampling factor must be >=1")
        if self.TiltPairs and self.DoFlip:
            errors.append("Phase cannot be corrected on tilt pairs")
        return errors

    def summary(self):
        self.setSamplingMode()
        message=[]
        familyFn = self.getFilename("family", family=self.Family)
        if self.TiltPairs:
            part1 = "tilt pairs"
            part2 = "particle pairs"
        else:
            part1 = "family <%s>" % self.Family
            part2 = "particles"
        message.append("Picking %s with size  <%d> from [%s]" % (part1, self.ParticleSize,self.pickingDir))

        msg="Extraction sampling rate <%3.2f> (" % (self.TsFinal)
        if self.downsamplingMode==DownsamplingMode.SameAsPicking:
            msg+="same as the one used in picking)"
        elif self.downsamplingMode==DownsamplingMode.SameAsOriginal:
            msg+="same as the original micrographs)"
        else:
            msg+="new sampling)"
        message.append(msg)
    
        if exists(familyFn):
            md = MetaData(familyFn)
            message.append("%d %s extracted" % (md.size(), part2))
        return message

    def getMicrographInfo(self, micrograph, md):
        for objId in md:
            micFullName = md.getValue(MDL_MICROGRAPH,objId)
            micName = removeBasenameExt(micFullName)
            if micrograph==micName:
                if md.containsLabel(MDL_MICROGRAPH_ORIGINAL):
                    micOriginalName = md.getValue(MDL_MICROGRAPH_ORIGINAL,objId)
                elif md.containsLabel(MDL_MICROGRAPH):
                    micOriginalName = md.getValue(MDL_MICROGRAPH,objId)
                else:
                    micOriginalName = None
                if self.containsCTF:
                    ctfFile = md.getValue(MDL_CTF_MODEL, objId)
                else:
                    ctfFile = None
                return (micFullName, micOriginalName, ctfFile)
            if self.TiltPairs:
                micFullName = md.getValue(MDL_MICROGRAPH_TILTED, objId)
                micName = removeBasenameExt(micFullName)
                if micrograph == micName:
                    if md.containsLabel(MDL_MICROGRAPH_TILTED_ORIGINAL):
                        micOriginalName=md.getValue(MDL_MICROGRAPH_TILTED_ORIGINAL,objId)
                    elif md.containsLabel(MDL_MICROGRAPH_TILTED):
                        micOriginalName=md.getValue(MDL_MICROGRAPH_TILTED, objId)
                    else:
                        micOriginalName=None
                    return (micFullName, micOriginalName, None)

    def visualize(self):
        selfile = self.getFilename("family", family=self.Family)
        if not exists(selfile):
            showError("Error", "There is no result yet")
        else:
            from protlib_utils import runShowJ
            if self.TiltPairs:
                runShowJ(selfile,extraParams="--mode metadata --render first")
            else:
                runShowJ(selfile)
                
            selfileRoot = self.workingDirPath(self.Family)
            fnMD = selfileRoot+".xmd"
            if exists(fnMD):
                from protlib_gui_figure import XmippPlotter
                from xmipp import MDL_ZSCORE
                MD=MetaData(fnMD)
                if MD.containsLabel(MDL_ZSCORE):
                    #MD.sort(MDL_ZSCORE)
                    xplotter = XmippPlotter(windowTitle="Zscore particles sorting")
                    xplotter.createSubPlot("Particle sorting", "Particle number", "Zscore")
                    xplotter.plotMd(MD, False, mdLabelY=MDL_ZSCORE)
                    xplotter.show()
    
    def createBlocksInExtractFile(self,fnMicrographsSel):
        md = MetaData(fnMicrographsSel)
        blocks = []
        tiltPairs = md.containsLabel(MDL_MICROGRAPH_TILTED)
        def addBlock(label, objId): 
            micName = removeBasenameExt(md.getValue(label, objId)) 
            blocks.append(_getFilename('mic_block', micName=micName))
        for objId in md:
            addBlock(MDL_MICROGRAPH, objId)            
            if tiltPairs:
                addBlock(MDL_MICROGRAPH_TILTED, objId)
        return blocks

def updateSampling(log, SamplingMd, TsOriginal, Ts, TsMode):
    md = MetaData(SamplingMd)
    objId = md.firstObject()
    
    if TsMode == DownsamplingMode.SameAsOriginal:
        md.removeLabel(MDL_SAMPLINGRATE_ORIGINAL)
    elif TsMode == DownsamplingMode.NewDownsample:
        md.setValue(MDL_SAMPLINGRATE_ORIGINAL, TsOriginal, objId)
        
    md.setValue(MDL_SAMPLINGRATE, Ts, objId)
    md.write(SamplingMd)

def createExtractList(log,Family,fnMicrographsSel,pickingDir,fnExtractList):
    md=MetaData(fnMicrographsSel)
    mdpos=MetaData()
    mdposAux=MetaData()
    for objId in md:
        micName = removeBasenameExt(md.getValue(MDL_MICROGRAPH, objId))
        
        fnManual = join(pickingDir, "extra", micName + ".pos")
        fnAuto1 = join(pickingDir, "extra", micName + "_auto.pos")
        
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
        umicName = removeBasenameExt(md.getValue(MDL_MICROGRAPH, objId))
        tmicName =removeBasenameExt(md.getValue(MDL_MICROGRAPH_TILTED, objId))
        
        fnUntilted = join(pickingDir,"extra",umicName + ".pos")
        fnTilted = join(pickingDir,"extra",tmicName + ".pos")

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
        mdUntiltedPos.write(fn, MD_APPEND)
        fn =_getFilename('mic_block_fn', micName=tmicName, fn=fnExtractList)
        mdTiltedPos.write(fn, MD_APPEND)

def extractParticles(log,ExtraDir,micrographName, ctf, fullMicrographName, originalMicrograph, micrographToExtract,
                     TsFinal, TsInput, downsamplingMode,
                     fnExtractList, particleSize, doFlip, doNorm, doInvert, bgRadius):
    
    
    fnBlock = _getFilename('mic_block_fn', micName=micrographName, fn=fnExtractList)
    md = MetaData(fnBlock)
    printLog( 'Metadata block: %s' % fnBlock, log)
    if md.isEmpty():
        printLog( "EMPTY", log)
        printLog('Metadata: %s is empty, returning without extracting' % fnBlock, log)
        return
    
    # Extract 
    rootname = join(ExtraDir, micrographName)
    arguments="-i "+micrographToExtract+" --pos "+fnBlock+" -o "+rootname+" --Xdim "+str(particleSize)
    if abs(TsFinal-TsInput)>0.001:
        arguments+=" --downsampling "+str(TsFinal/TsInput)
    if doInvert:
        arguments += " --invert"
    try:
        runJob(log,"xmipp_micrograph_scissor", arguments)
    except OSError, e:
       print ("Execution failed %s, command: %s" % (e, 'xmipp_micrograph_scissor'))
    # Normalize 
    if doNorm:
        if bgRadius == 0:
            bgRadius = int(particleSize/2)
        arguments = "-i "+rootname+'.stk --method Ramp --background circle '+str(bgRadius)
        runJob(log,"xmipp_transform_normalize",arguments)

    # Substitute the micrograph name if it comes from the flipped version
    # Add information about the ctf if available

    if fullMicrographName != micrographToExtract or ctf != None:
        selfile = rootname + ".xmd"
        md = MetaData(selfile)
        if downsamplingMode==DownsamplingMode.SameAsOriginal:
            md.setValueCol(MDL_MICROGRAPH, originalMicrograph)
        else:
            md.setValueCol(MDL_MICROGRAPH, fullMicrographName)
        if downsamplingMode==DownsamplingMode.NewDownsample:
            downsamplingFactor=TsFinal/TsInput
            md.operate("Xcoor=Xcoor*%f"%downsamplingFactor)
            md.operate("Ycoor=Ycoor*%f"%downsamplingFactor)
        if ctf != None:
            md.setValueCol(MDL_CTF_MODEL,ctf)
        md.write(selfile)

def gatherSelfiles(log, WorkingDir, ExtraDir):
    stackFiles = glob.glob(join(ExtraDir,"*.stk"))
    stackFiles.sort()
    familySelfile=MetaData()
    for stackFile in stackFiles:
        selfile=stackFile.replace(".stk",".xmd")
        md=MetaData(selfile)
        familySelfile.unionAll(md)
    familySelfile.write(join(WorkingDir,"images.xmd"))

def gatherTiltPairSelfiles(log, WorkingDir, ExtraDir, fnMicrographs):
    mdPairs = MetaData(fnMicrographs)
    mdUntilted = MetaData()
    mdTilted = MetaData()
    for objId in mdPairs:
        umicName = removeBasenameExt(mdPairs.getValue(MDL_MICROGRAPH, objId))        
        fnUntilted = join(ExtraDir, umicName + ".xmd")
        # Check if there are picked particles in this micrographs
        if exists(fnUntilted):
            mdUntilted.unionAll(MetaData(fnUntilted))            
            tmicName = removeBasenameExt(mdPairs.getValue(MDL_MICROGRAPH_TILTED,objId))
            mdTilted.unionAll(MetaData(join(ExtraDir, tmicName + ".xmd")))
        
    mdUntilted.write(join(WorkingDir, "images_untilted.xmd"))
    mdTilted.write(join(WorkingDir, "images_tilted.xmd"))

    mdTiltedPairs = MetaData()
    mdTiltedPairs.setColumnValues(MDL_IMAGE, mdUntilted.getColumnValues(MDL_IMAGE))
    mdTiltedPairs.setColumnValues(MDL_IMAGE_TILTED, mdTilted.getColumnValues(MDL_IMAGE))
    mdTiltedPairs.write(join(WorkingDir,"images.xmd"))

def sortImagesInFamily(log, WorkingDir):
    fn = os.path.join(WorkingDir, 'images.xmd')
    md = MetaData(fn)
    if not md.isEmpty():
        runJob(log, "xmipp_image_sort_by_statistics","-i %(fn)s --multivariate --addToInput" % locals())
        md.read(fn) # Should have ZScore label after runJob
        md.sort(MDL_ZSCORE)
        md.write(fn)
 
def avgZscore(log,WorkingDir,micrographSelfile):
    allParticles = MetaData(join(WorkingDir, "images.xmd"))
    mdavgZscore = MetaData()
    mdavgZscore.aggregate(allParticles, AGGR_AVG, MDL_MICROGRAPH, MDL_ZSCORE, MDL_ZSCORE)
    oldMicrographsSel = MetaData(micrographSelfile)
    oldMicrographsSel.removeLabel(MDL_ZSCORE)
    newMicrographsSel = MetaData()
    # Make copy of metadata because removeLabel leaves the values in the table
    newMicrographsSel.join(oldMicrographsSel,mdavgZscore,MDL_MICROGRAPH,MDL_MICROGRAPH,LEFT)
    newMicrographsSel.write(micrographSelfile)

