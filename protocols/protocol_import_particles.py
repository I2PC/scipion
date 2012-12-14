#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based importing single-particles: 
#  - normalization
#  - sort_junk

# Author: Carlos Oscar, August 2011
#
from protlib_base import *
from protlib_particles import *
import xmipp
from xmipp import MetaData, FileName, ImgSize, MDL_IMAGE, MDL_SAMPLINGRATE
import glob
import os
from os.path import relpath, dirname
from protlib_utils import runJob
from protlib_filesystem import deleteFile, createDir, copyFile, fixPath, \
findProjectInPathTree, xmippRelpath, splitFilename


class ProtImportParticles(ProtParticlesBase):
    def __init__(self, scriptname, project):
        ProtParticlesBase.__init__(self, protDict.import_particles.name, scriptname, project)
        self.Import += 'from protocol_import_particles import *'

    def defineSteps(self):
        fnOut = self.getFilename('acquisition')
        self.insertStep('createAcquisitionMd',verifyfiles=[fnOut],samplingRate=self.SamplingRate, fnOut=fnOut)
        
        fnOut = self.getFilename('micrographs')
        self.insertStep('createEmptyMicrographSel',verifyfiles=[fnOut],fnOut=fnOut)
        
        fnOut = self.getFilename('images')
        self.insertStep('importImages', verifyfiles=[fnOut],
                           InputFile=self.InputFile, WorkingDir=self.WorkingDir, DoCopy=self.DoCopy,
                           ImportAll=self.ImportAll, SubsetMode=self.SubsetMode, Nsubset=self.Nsubset)
        
        if self.DoInvert and self.ImportAll:
            self.insertStep('invert', ImagesMd=fnOut,Nproc=self.NumberOfMpi)
        
        if self.DoRemoveDust and self.ImportAll:
            self.insertStep('removeDust', ImagesMd=fnOut,threshold=self.DustRemovalThreshold,Nproc=self.NumberOfMpi)
        
        if self.DoNorm and self.ImportAll:
            self.insertStep('normalize', ImagesMd=fnOut,bgRadius=self.BackGroundRadius,Nproc=self.NumberOfMpi)
        
        if self.DoSort:
            self.insertStep('sortImages', verifyfiles=[fnOut], fn=fnOut)
            
    def validate(self):
        errors = []
        inputExt=os.path.splitext(self.InputFile)[1]
        if not inputExt in ['.mrc','.stk','.sel','.xmd','.hed','.img', '.ctfdat']:
            errors.append("Input file must be stack or metadata (valid extensions are .mrc, .stk, .sel, .xmd, .hed, .img, .ctfdat")
        else:
            if inputExt in ['.sel', '.xmd', '.ctfdat']:
                md = MetaData(self.InputFile)
                if not md.containsLabel(MDL_IMAGE):
                    errors.append("Cannot find label for images in the input file")
        return errors

    def summary(self):
        message=[]
        message.append("Stack imported from: [%s]"%self.InputFile)
        if self.DoCopy:
            message.append("Copied into [%s]"%self.WorkingDir)
        steps=[]
        if self.DoInvert:
            steps.append("Constrast inversion")
        if self.DoRemoveDust:
            steps.append("Dust removal")
        if self.DoNorm:
            steps.append("Ramp normalization")
        if len(steps)>0:
            message.append("Steps applied: "+",".join(steps))
            
        return message
       

def createAcquisitionMd(log, samplingRate, fnOut):
    md = MetaData()
    md.setValue(MDL_SAMPLINGRATE, float(samplingRate), md.addObject())
    md.write(fnOut)

def createEmptyMicrographSel(log, fnOut):
    md = MetaData()
    md.setValue(MDL_IMAGE,"ImportedImages", md.addObject())
    md.write(fnOut)

def writeImagesMd(log, md, ImportAll, SubsetMode, Nsubset, imagesFn, imagesStk, DoCopy):
    if not ImportAll:
        if SubsetMode=="Random particles":
            mdaux=MetaData()
            mdaux.randomize(md)
        else:
            mdaux=MetaData(md)
        md.selectPart(mdaux, 0, Nsubset)
    md.write(imagesFn)
    if DoCopy:
        runJob(log,"xmipp_image_convert","-i %(imagesFn)s -o %(imagesStk)s" % locals())
        md = MetaData(imagesStk)
        md.write(imagesFn)      
   
def importImages(log, InputFile, WorkingDir, DoCopy, ImportAll, SubsetMode, Nsubset):
    imagesFn = getImagesFilename(WorkingDir)
    fnInput = FileName(InputFile)    
    md = MetaData(InputFile)
    
    if fnInput.isMetaData():        
        inputRelativePath = dirname(relpath(InputFile, '.'))
        projectPath = findProjectInPathTree(InputFile)
        for id in md:
            imgFn = md.getValue(MDL_IMAGE, id)
            imgNo, imgFn = splitFilename(imgFn)
            imgFn = xmippRelpath(fixPath(imgFn, projectPath, inputRelativePath, '.'))
            md.setValue(MDL_IMAGE, "%s@%s" % (imgNo, imgFn), id)
        outExt = '.stk'
    else:
        outExt = '.%s' % fnInput.getExtension()
    imagesStk = imagesFn.replace('.xmd', outExt)
    writeImagesMd(log, md, ImportAll, SubsetMode, Nsubset, imagesFn, imagesStk, DoCopy)
       

def invert(log,ImagesMd,Nproc):
    runJob(log,'xmipp_image_operate','-i %(ImagesMd)s --mult -1' % locals(),Nproc)

def removeDust(log,ImagesMd,threshold,Nproc):
    runJob(log,'xmipp_transform_filter','-i %(ImagesMd)s --bad_pixels outliers %(threshold)f' % locals(),Nproc)

def normalize(log,ImagesMd,bgRadius,Nproc):
    if bgRadius <= 0:
        particleSize = ImgSize(ImagesMd)[0]
        bgRadius = int(particleSize/2)
    runJob(log,"xmipp_transform_normalize", '-i %(ImagesMd)s --method Ramp --background circle %(bgRadius)d' % locals(),Nproc)
