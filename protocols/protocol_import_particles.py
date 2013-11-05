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
from xmipp import MetaData, FileName, MDL_IMAGE, MDL_SAMPLINGRATE
import glob
import os
from os.path import relpath, dirname
from protlib_utils import runJob
from protlib_filesystem import deleteFile, createDir, copyFile, fixPath, \
findProjectPath, xmippRelpath, splitFilename

class ProtImportParticles(ProtParticlesBase):
    def __init__(self, scriptname, project):
        ProtParticlesBase.__init__(self, protDict.import_particles.name, scriptname, project)
        self.Import += 'from protocol_import_particles import *'

    def defineSteps(self):
        if self.DoInvert or self.DoRemoveDust or self.DoNorm:
            self.DoCopy = True 
        fnOut = self.getFilename('acquisition')
        self.insertStep('createAcquisitionMd',verifyfiles=[fnOut],samplingRate=self.SamplingRate, fnOut=fnOut)
        
        fnOut = self.getFilename('micrographs')
        self.insertStep('createEmptyMicrographSel',verifyfiles=[fnOut],fnOut=fnOut)
        
        fnOut = self.getFilename('images')
        self.insertStep('importImages', verifyfiles=[fnOut],
                           InputFile=self.InputFile, WorkingDir=self.WorkingDir, DoCopy=self.DoCopy,
                           ImportAll=self.ImportAll, SubsetMode=self.SubsetMode, Nsubset=self.Nsubset)
        
        if self.DoInvert:
            self.insertStep('invert', ImagesMd=fnOut,Nproc=self.NumberOfMpi)
        
        if self.DoRemoveDust:
            self.insertStep('removeDust', ImagesMd=fnOut,threshold=self.DustRemovalThreshold,Nproc=self.NumberOfMpi)
        
        if self.DoNorm:
            self.insertStep('runNormalize', stack=fnOut, normType=self.NormType, 
                            bgRadius=self.BackGroundRadius, Nproc=self.NumberOfMpi)
        
        if self.DoSort:
            self.insertStep('sortImages', verifyfiles=[fnOut], ImagesFn=fnOut)
            
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
        messages = []
        messages.append("Import from: [%s]" % self.InputFile)
        self.addBasicPreprocessSteps()
        self.addStepsSummary(messages)
        
        if self.DoCopy:
            messages.append("Copied into [%s]"%self.WorkingDir)
       
        return messages

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
        md = MetaData(imagesFn)
        mdStk = MetaData(imagesStk)
        md.merge(mdStk) # Change MDL_IMAGE to point to new stack
        md.write(imagesFn)      
   
def importImages(log, InputFile, WorkingDir, DoCopy, ImportAll, SubsetMode, Nsubset):
    imagesFn = getImagesFilename(WorkingDir)
    fnInput = FileName(InputFile)    
    md = MetaData(InputFile)
    
    if fnInput.isMetaData():        
        inputRelativePath = dirname(relpath(InputFile, '.'))
        projectPath = findProjectPath(InputFile)
        for id in md:
            imgFn = md.getValue(MDL_IMAGE, id)
            imgNo, imgFn = splitFilename(imgFn)
            imgFn = xmippRelpath(fixPath(imgFn, projectPath, inputRelativePath, '.'))
            if imgNo != None:
                imgFn = "%s@%s" % (imgNo, imgFn)
            md.setValue(MDL_IMAGE, imgFn, id)

            imgFn = md.getValue(MDL_MICROGRAPH, id)
            imgNo, imgFn = splitFilename(imgFn)
            imgFn = xmippRelpath(fixPath(imgFn, projectPath, inputRelativePath, '.'))
            if imgNo != None:
                imgFn = "%s@%s" % (imgNo, imgFn)
            md.setValue(MDL_MICROGRAPH, imgFn, id)
            
            imgFn = md.getValue(MDL_CTF_PARAMS, id)
            imgFn = xmippRelpath(fixPath(imgFn, projectPath, inputRelativePath, '.'))
            md.setValue(MDL_CTF_PARAMS, imgFn, id)
            
        outExt = '.stk'
    else:
        outExt = '.%s' % fnInput.getExtension()
    imagesStk = imagesFn.replace('.xmd', outExt)
    writeImagesMd(log, md, ImportAll, SubsetMode, Nsubset, imagesFn, imagesStk, DoCopy)
       

