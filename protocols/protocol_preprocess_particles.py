#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based pre-processing of single-particles: 

# Author: Carlos Oscar, August 2011
#
from protlib_base import *
from protlib_particles import *
from xmipp import MetaData, MetaDataInfo, FileName, MDL_SAMPLINGRATE, MDL_IMAGE, MDL_ZSCORE, MDL_MICROGRAPH, MDL_XCOOR, MDL_YCOOR, MDL_IMAGE_ORIGINAL, \
    MDL_IMAGE_TILTED, INNER_JOIN
import os
from os.path import exists, split, splitext
from protlib_utils import runJob, runShowJ
from protlib_filesystem import deleteFile,findAcquisitionInfo, createLink, removeBasenameExt
import glob
from protlib_gui_ext import showError

class ProtPreprocessParticles(ProtParticlesBase):
    def __init__(self, scriptname, project):
        ProtParticlesBase.__init__(self, protDict.preprocess_particles.name, scriptname, project)
        self.Import = 'from protocol_preprocess_particles import *'
        (self.InputDir,file)=split(self.InSelFile)
        baseFile = removeBasenameExt(file)
        self.OutStack = self.getFilename('images_stk')
        self.OutMetadata = self.getFilename('images')
        self.TiltPair = exists(getProtocolFilename('tilted_pairs', WorkingDir=self.InputDir)) and self.InSelFile.find('images.xmd')!=-1
        if self.TiltPair:
            self.TiltedPairs=os.path.join(self.InputDir,'tilted_pairs.xmd')
            self.InputUntilted=self.InSelFile.replace(".xmd","_untilted.xmd")
            self.InputTilted=self.InSelFile.replace(".xmd","_tilted.xmd")
            self.OutUntiltedStack=self.getFilename('images_untilted_stk')
            self.OutUntiltedMetadata=self.getFilename('images_untilted')
            self.OutTiltedStack=self.getFilename('images_tilted_stk')
            self.OutTiltedMetadata=self.getFilename('images_tilted')
            self.OutStack=self.OutUntiltedStack

    def defineSteps(self):
        if self.TiltPair:
            fnTiltedPairs=self.workingDirPath('tilted_pairs.xmd')
            self.insertStep('copyFile',verifyfiles=[fnTiltedPairs],source=self.TiltedPairs,dest=fnTiltedPairs)
            self.insertStep('copyImages',verifyfiles=[self.OutUntiltedMetadata, self.OutUntiltedStack],
                               InputFile=self.InputUntilted, OutputStack=self.OutUntiltedStack)
            self.insertStep('copyImages',verifyfiles=[self.OutTiltedMetadata, self.OutTiltedStack],
                               InputFile=self.InputTilted, OutputStack=self.OutTiltedStack)
            fnImages=self.getFilename('images')
            self.insertStep('copyFile',verifyfiles=[fnImages],source=self.InSelFile,dest=fnImages)
            self.insertStep('joinMetaDatas',InputFile=fnImages,UntiltedMetadata=self.OutUntiltedMetadata,TiltedMetadata=self.OutTiltedMetadata)
        else:
            self.insertStep('copyImages',verifyfiles=[self.OutMetadata, self.OutStack],
                               InputFile=self.InSelFile, OutputStack=self.OutStack)
        self.insertStep('createAcquisition',InputFile=self.InSelFile, WorkingDir=self.WorkingDir,
                        DoResize=self.DoResize, NewSize=self.NewSize)
        
        # Apply threshold if selected
        if self.DoThreshold:
            self.insertStep('runThreshold',stack=self.OutStack,selectionMode=self.SelectionMode,threshold=self.Threshold,
                            substituteBy=self.SubstituteBy, substituteValue=self.SubstituteValue, Nproc=self.NumberOfMpi)
            if self.TiltPair:
                self.insertStep('runFourierFilter',stack=self.OutTiltedStack,freq_low=self.Freq_low,freq_high=self.Freq_high,freq_decay=self.Freq_decay,Nproc=self.NumberOfMpi)

        # Apply filters if selected
        if self.DoFourier:
            self.insertStep('runFourierFilter',stack=self.OutStack,freq_low=self.Freq_low,freq_high=self.Freq_high,freq_decay=self.Freq_decay,Nproc=self.NumberOfMpi)
            if self.TiltPair:
                self.insertStep('runFourierFilter',stack=self.OutTiltedStack,freq_low=self.Freq_low,freq_high=self.Freq_high,freq_decay=self.Freq_decay,Nproc=self.NumberOfMpi)
        if self.DoGaussian:
            self.insertStep('runGaussianFilter',stack=self.OutStack,freq_sigma=self.Freq_sigma, Nproc=self.NumberOfMpi)
            if self.TiltPair:
                self.insertStep('runGaussianFilter',stack=self.OutTiltedStack,freq_sigma=self.Freq_sigma, Nproc=self.NumberOfMpi)
        if self.DoRealGaussian:
            self.insertStep('runRealGaussianFilter',stack=self.OutStack,real_sigma=self.Real_sigma, Nproc=self.NumberOfMpi)
            if self.TiltPair:
                self.insertStep('runRealGaussianFilter',stack=self.OutTiltedStack,real_sigma=self.Real_sigma, Nproc=self.NumberOfMpi)
            
        # Apply mask
        if self.DoMask:
            if self.Substitute == "value":
                self.Substitute = str(self.SubstituteValue)
            params = "--substitute %s --mask %s " % (self.Substitute, self.MaskType)
            if self.MaskType == 'raised_cosine':
                params += "-%d -%d" % (self.MaskRadius, self.MaskRadius + self.MaskRadiusOuter)
            elif self.MaskType == 'circular':
                params += '-%d' % self.MaskRadius
            else: # from file:
                params += self.MaskFile
            self.insertRunJobStep("xmipp_transform_mask", "-i %s"%self.OutStack+" "+params)
            if self.TiltPair:
                self.insertRunJobStep("xmipp_transform_mask", "-i %s"%self.OutTiltedStack+" "+params)
            
        # Resize images
        if self.DoCrop:
            self.insertStep('runCrop',stack=self.OutStack, cropSize=self.CropSize, tmpStack=self.tmpPath('tmpCrop.stk'))
            if self.TiltPair:
                self.insertStep('runCrop',stack=self.OutTiltedStack, cropSize=self.CropSize, tmpStack=self.tmpPath('tmpCrop.stk'))
        if self.DoResize:
            self.insertStep('runResize',stack=self.OutStack,new_size=self.NewSize,Nproc=self.NumberOfMpi)
            if self.TiltPair:
                self.insertStep('runResize',stack=self.OutTiltedStack,new_size=self.NewSize,Nproc=self.NumberOfMpi)
        
    def validate(self):
        errors = []
        if self.DoResize:
            if self.NewSize < 1:
                errors.append("New size for scale have not correctly set")
            if self.TiltPair:
                errors.append("Cannot scale particles extracted from tilt pairs. Re-extract the particles at a different sampling rate, instead.")
            if self.DoCrop and self.OutputSize > self.NewSize:
                errors.append("Crop output size cannot be greater than resize output size")
        return errors
        
    def summary(self):
        messages = []        
        messages.append("Input images: [%s]" % self.InSelFile)

        self.addSummaryStep(self.DoFourier, "Fourier filter: freq_low = %(Freq_low)f freq_high = %(Freq_high)f freq_decay = %(Freq_decay)f")            
        self.addSummaryStep(self.DoGaussian, "Gaussian filter: freq_sigma = %(Freq_sigma)f")
        self.addSummaryStep(self.DoMask, "Mask: mask file = %(MaskFile)s substituted value = %(Substitute)s")
        self.addSummaryStep(self.DoResize, "Resize: NewSize = %(NewSize)d")            
        self.addSummaryStep(self.DoCrop, "Crop: CropSize = %(CropSize)d")
        self.addStepsSummary(messages)
        
        messages.append("Output: [%s]" % self.OutMetadata)
        if self.TiltPair:
            messages.append("Operations applied to both, untilted and tilted, images")
        
        return messages


def createAcquisition(log,InputFile,WorkingDir,DoResize,NewSize):
    fnAcqIn = findAcquisitionInfo(InputFile)
    if not fnAcqIn is None:
        fnAcqOut = getProtocolFilename('acquisition', WorkingDir=WorkingDir)
              
        if not DoResize:
            createLink(log, fnAcqIn, fnAcqOut)
        else:    
            md = MetaData(fnAcqIn)
            id = md.firstObject()
            Ts = md.getValue(MDL_SAMPLINGRATE, id)
            (Xdim, Ydim, Zdim, Ndim, _) = MetaDataInfo(InputFile)
            downsampling = float(Xdim)/NewSize;
            md.setValue(MDL_SAMPLINGRATE,Ts*downsampling,id)
            md.write(getProtocolFilename('acquisition', WorkingDir=WorkingDir))

def copyImages(log,InputFile,OutputStack):
    fnOutputMetadata=OutputStack.replace(".stk",".xmd")
    runJob(log,"xmipp_image_convert","-i %(InputFile)s -o %(OutputStack)s --track_origin --save_metadata_stack %(fnOutputMetadata)s --keep_input_columns" % locals())
    mDstack = MetaData(fnOutputMetadata)
    mDstack.removeLabel(MDL_ZSCORE)
    mDstack.write(fnOutputMetadata)

def joinMetaDatas(log,InputFile,UntiltedMetadata,TiltedMetadata):    
    mDInput=MetaData(InputFile)
    mDstack = MetaData(UntiltedMetadata)
    mDstack.removeLabel(MDL_MICROGRAPH)
    mDstack.removeLabel(MDL_XCOOR)
    mDstack.removeLabel(MDL_YCOOR)
    mdPart1=MetaData()
    mdPart1.join(mDstack,mDInput,MDL_IMAGE_ORIGINAL,MDL_IMAGE,INNER_JOIN)
    mdPart1.removeLabel(MDL_IMAGE_ORIGINAL)

    mDstack = MetaData(TiltedMetadata)
    mDstack.removeLabel(MDL_MICROGRAPH)
    mDstack.removeLabel(MDL_XCOOR)
    mDstack.removeLabel(MDL_YCOOR)
    mDstack.renameColumn(MDL_IMAGE,MDL_IMAGE_TILTED)
    mdJoined=MetaData()
    mdJoined.join(mDstack,mdPart1,MDL_IMAGE_ORIGINAL,MDL_IMAGE_TILTED,INNER_JOIN)
    mdJoined.removeLabel(MDL_IMAGE_ORIGINAL)
    mdJoined.write(InputFile)
