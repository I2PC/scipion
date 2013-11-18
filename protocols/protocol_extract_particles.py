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
from protlib_particles import *
from xmipp import MetaData, MetaDataInfo, AGGR_AVG, LEFT, MD_APPEND, MDL_IMAGE, MDL_CTF_MODEL, MDL_MICROGRAPH, MDL_MICROGRAPH_TILTED, \
                  MDL_MICROGRAPH_ORIGINAL, MDL_MICROGRAPH_TILTED_ORIGINAL, MDL_ZSCORE, MDL_IMAGE_TILTED, MDL_ENABLED, \
                  MDL_SAMPLINGRATE, MDL_SAMPLINGRATE_ORIGINAL, getBlocksInMetaDataFile
import glob
from os.path import exists
from protlib_utils import runJob, printLog
from protlib_filesystem import createLink, removeFilenameExt, replaceBasenameExt
from protlib_gui_ext import showError

class DownsamplingMode:
    SameAsPicking, SameAsOriginal, NewDownsample = range(3)
    
class ProtExtractParticles(ProtParticlesBase):
    def __init__(self, scriptname, project):
        ProtParticlesBase.__init__(self, protDict.extract_particles.name, scriptname, project)
        self.Import += 'from protocol_extract_particles import *'
        # Take some parameter from previous picking protocol run
        self.setPreviousRun(self.PickingRun)
        self.pickingDir = self.PrevRun.WorkingDir
        self.inputFilename('acquisition','micrographs')
        self.inputProperty('TiltPairs')
        self.MicrographsMd=self.Input['micrographs']
        self.micrographs = self.getEquivalentFilename(self.PrevRun, self.Input['micrographs'])
        self.RejectionMethod = getattr(self, 'RejectionMethod', 'none')
        
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
        
        if abs(self.TsFinal-self.TsInput) < 0.001:
            self.downsamplingMode=DownsamplingMode.SameAsPicking
        elif abs(self.TsFinal-self.TsOriginal) < 0.001:
            self.downsamplingMode=DownsamplingMode.SameAsOriginal
        else:
            self.downsamplingMode=DownsamplingMode.NewDownsample

    def defineSteps(self):
        self.setSamplingMode()
        filesToCopy = [ self.Input['micrographs'], self.Input['acquisition']]
        self.insertImportOfFiles(filesToCopy, copy=True)
        self.insertStep('createDir',path=self.ExtraDir)

        # Update sampling rate in 'acquisition_info.xmd' if necessary
        if self.downsamplingMode != DownsamplingMode.SameAsPicking:
            self.insertStep('updateSampling', SamplingMd=self.getFilename('acquisition'), 
                            TsOriginal=self.TsOriginal, Ts=self.TsFinal, TsMode=self.downsamplingMode)
        
        md = MetaData( self.Input['micrographs'])
        md.removeDisabled()
        self.containsCTF = md.containsLabel(MDL_CTF_MODEL)

        # Process each micrograph                        
        for id in md:
            micrograph,micrographOriginal,ctf=self.getMicrographInfo(md, id)
            micrographName=removeBasenameExt(micrograph)
            micrographPos=self.PrevRun.extraPath(micrographName+".pos")
            if not os.path.exists(micrographPos):
                continue
            
            parent_id = XmippProjectDb.FIRST_STEP
            
            # Downsampling?
            if self.downsamplingMode==DownsamplingMode.SameAsPicking:
                micrographToExtract = micrograph
            elif self.downsamplingMode==DownsamplingMode.SameAsOriginal:
                micrographToExtract = micrographOriginal
            else:
                micrographDownsampled = self.tmpPath(micrographName+"_downsampled.xmp")
                parent_id=self.insertParallelRunJobStep("xmipp_transform_downsample", "-i %s -o %s --step %f --method fourier" % \
                                                        (micrographOriginal,micrographDownsampled,self.DownsampleFactor),
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
            parent_id = self.insertParallelStep('extractParticles', parent_step_id=parent_id,
                                  ExtraDir=self.ExtraDir,micrographName=micrographName,micrographPos=micrographPos,ctf=ctf,
                                  fullMicrographName=micrograph,micrographOriginal=micrographOriginal,
                                  micrographToExtract=micrographToExtract,
                                  TsFinal=self.TsFinal, TsInput=self.TsInput, downsamplingMode=self.downsamplingMode,
                                  particleSize=self.ParticleSize,doFlip=self.DoFlip,doInvert=self.DoInvert,
                                  doNorm=self.DoNorm, normType=self.NormType, bgRadius=self.BackGroundRadius)
            if self.DoRemoveDust:
                self.insertParallelStep('deleteFile', parent_step_id=parent_id,filename=micrographNoDust,verbose=True)
            if self.downsamplingMode==DownsamplingMode.NewDownsample:
                self.insertParallelStep('deleteFile', parent_step_id=parent_id,filename=micrographDownsampled,verbose=True)
            if self.DoFlip:
                self.insertParallelStep('deleteFile', parent_step_id=parent_id,filename=micrographFlipped,verbose=True)

        # Create results metadata
        ImagesFn = self.getFilename('images')
        if self.TiltPairs:
            self.insertStep('createTiltPairsImagesMd', verifyfiles=[ImagesFn], WorkingDir=self.WorkingDir, ExtraDir=self.ExtraDir, 
                            fnMicrographs=self.micrographs)
        else:
            self.insertStep('createImagesMd', verifyfiles=[ImagesFn], ImagesFn=ImagesFn, ExtraDir=self.ExtraDir)
            self.insertStep('sortImages',ImagesFn=ImagesFn,rejectionMethod=self.RejectionMethod, maxZscore=self.MaxZscore, percentage=self.Percentage)
            self.insertStep('avgZscore',WorkingDir=self.WorkingDir, micrographSelfile=self.micrographs)

    def validate(self):
        errors = []
        if self.pickingDir:
            if exists(self.MicrographsMd):
                md = MetaData(self.MicrographsMd)
                if self.DoFlip and not md.containsLabel(MDL_CTF_MODEL):
                    errors.append("Micrographs metadata: <%s>\n does not contain CTF information for phase flipping" % self.MicrographsMd)
        else:
            errors.append("Picking run is not valid")
        if self.DownsampleFactor<1:
            errors.append("Downsampling factor must be >=1")
        if self.TiltPairs and self.DoFlip:
            errors.append("Phase cannot be corrected on tilt pairs")
        if self.RejectionMethod=='MaxZscore' and self.MaxZscore<0:
            errors.append("MaxZscore must be positive")
        elif self.RejectionMethod=='Percentage' and (self.Percentage<0 or self.Percentage>100):
            errors.append("Percentage must be between 0 and 100")
        if self.ParticleSize<=0:
            errors.append("Size has to be a positive amount")
        return errors

    def summary(self):
        self.setSamplingMode()
        message = ["Picking run: [%s]" % self.pickingDir]
        if self.TiltPairs:
            part1 = "Tilt pairs"
            part2 = "Particle pairs"
        else:
            part1 = "Particles"
            part2 = "Particles"
        message.append("%s with size  <%d>" % (part1, self.ParticleSize))

        msg = "Sampling rate: <%3.2f> (" % (self.TsFinal)
        if self.downsamplingMode==DownsamplingMode.SameAsPicking:
            msg += "same as the one used in picking)"
        elif self.downsamplingMode==DownsamplingMode.SameAsOriginal:
            msg += "same as the original micrographs)"
        else:
            msg += "different sampling)"
        message.append(msg)
    
        if self.RejectionMethod=='MaxZscore':
            message.append('Rejecting: Zscore>'+str(self.MaxZscore))
        elif self.RejectionMethod=='Percentage':
            message.append('Rejecting: '+str(self.Percentage)+"%")

        fn = self.getFilename("images")
        if exists(fn):
            _, _, _, _, size = MetaDataInfo(fn)
            message.append("%s extracted: <%d>" % (part2, size))
        return message

    def getMicrographInfo(self, md, id):
        micrograph=md.getValue(MDL_MICROGRAPH,id)
        micrographOriginal=md.getValue(MDL_MICROGRAPH_ORIGINAL,id)
        if micrographOriginal is None:
            micrographOriginal=micrograph
        ctf=md.getValue(MDL_CTF_MODEL,id)
        return [micrograph,micrographOriginal,ctf]

def updateSampling(log, SamplingMd, TsOriginal, Ts, TsMode):
    md = MetaData(SamplingMd)
    objId = md.firstObject()
    
    if TsMode == DownsamplingMode.SameAsOriginal:
        md.removeLabel(MDL_SAMPLINGRATE_ORIGINAL)
    elif TsMode == DownsamplingMode.NewDownsample:
        md.setValue(MDL_SAMPLINGRATE_ORIGINAL, TsOriginal, objId)
        
    md.setValue(MDL_SAMPLINGRATE, Ts, objId)
    md.write(SamplingMd)

def extractParticles(log,ExtraDir,micrographName, micrographPos, ctf, fullMicrographName, micrographOriginal, micrographToExtract,
                     TsFinal, TsInput, downsamplingMode, particleSize, doFlip, doInvert, doNorm, normType, bgRadius):
    md=readPosCoordinates(micrographPos)
    fnTmpPos=os.path.join(ExtraDir,"tmp_%s.pos"%micrographName)
    md.write(fnTmpPos)

    # Extract 
    rootname = join(ExtraDir, micrographName)
    arguments="-i %(micrographToExtract)s --pos %(fnTmpPos)s -o %(rootname)s --Xdim %(particleSize)d" % locals()
    if abs(TsFinal-TsInput)>0.001:
        arguments+=" --downsampling "+str(TsFinal/TsInput)
    if doInvert:
        arguments += " --invert"
    runJob(log,"xmipp_micrograph_scissor", arguments)
    # Normalize 
    if doNorm:
        runNormalize(log, rootname+'.stk',normType, bgRadius, 1)
    deleteFile(log,fnTmpPos)        

    # Substitute the micrograph name if it comes from the flipped version
    # Add information about the ctf if available
    if fullMicrographName != micrographToExtract or ctf != None:
        selfile = rootname + ".xmd"
        md = MetaData(selfile)
        if downsamplingMode==DownsamplingMode.SameAsOriginal:
            md.setValueCol(MDL_MICROGRAPH, micrographOriginal)
        else:
            md.setValueCol(MDL_MICROGRAPH, fullMicrographName)
        if downsamplingMode==DownsamplingMode.NewDownsample:
            downsamplingFactor=TsFinal/TsInput
            md.operate("Xcoor=Xcoor*%f"%downsamplingFactor)
            md.operate("Ycoor=Ycoor*%f"%downsamplingFactor)
        if ctf != None:
            md.setValueCol(MDL_CTF_MODEL,ctf)
        md.write(selfile)

def createImagesMd(log, ImagesFn, ExtraDir):
    stackFiles = glob.glob(join(ExtraDir,"*.stk"))
    stackFiles.sort()
    imagesMd = MetaData()
    for stack in stackFiles:
        fn = stack.replace(".stk",".xmd")
        md = MetaData(fn)
        imagesMd.unionAll(md)
    imagesMd.addItemId()
    imagesMd.write(ImagesFn)

def createTiltPairsImagesMd(log, WorkingDir, ExtraDir, fnMicrographs):
    mdPairs = MetaData("micrographPairs@"+fnMicrographs)
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
        
    fn = getImagesFilename(WorkingDir)    
    mdUntilted.addItemId()
    mdTilted.addItemId()
    mdUntilted.write(getProtocolFilename('images_untilted', WorkingDir=WorkingDir))
    mdTilted.write(getProtocolFilename('images_tilted', WorkingDir=WorkingDir))
    mdTiltedPairs = MetaData()
    mdTiltedPairs.setColumnValues(MDL_IMAGE, mdUntilted.getColumnValues(MDL_IMAGE))
    mdTiltedPairs.setColumnValues(MDL_IMAGE_TILTED, mdTilted.getColumnValues(MDL_IMAGE))
    mdTiltedPairs.addItemId()
    mdTiltedPairs.write(fn)
 
def avgZscore(log,WorkingDir,micrographSelfile):
    allParticles = MetaData(getImagesFilename(WorkingDir))
    mdavgZscore = MetaData()
    mdavgZscore.aggregate(allParticles, AGGR_AVG, MDL_MICROGRAPH, MDL_ZSCORE, MDL_ZSCORE)
    oldMicrographsSel = MetaData(micrographSelfile)
    oldMicrographsSel.removeLabel(MDL_ZSCORE)
    newMicrographsSel = MetaData()
    # Make copy of metadata because removeLabel leaves the values in the table
    newMicrographsSel.join(oldMicrographsSel,mdavgZscore,MDL_MICROGRAPH,MDL_MICROGRAPH,LEFT)
    newMicrographsSel.write(micrographSelfile)
