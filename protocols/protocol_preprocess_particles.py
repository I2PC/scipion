#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based pre-processing of single-particles: 

# Author: Carlos Oscar, August 2011
#
from protlib_base import *
import xmipp
import os
from protlib_utils import runJob, runShowJ
from protlib_filesystem import deleteFile,linkAcquisitionInfo, moveFile
import glob
from protlib_gui_ext import showError

class ProtPreprocessParticles(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.preprocess_particles.name, scriptname, project)
        self.Import = 'from protocol_preprocess_particles import *'
        (self.InputDir,file)=os.path.split(self.InSelFile)
        baseFile=os.path.splitext(file)[0]
        self.OutStack=self.workingDirPath(baseFile+".stk")
        self.OutMetadata=self.workingDirPath(baseFile+".xmd")
        self.TiltPair=os.path.exists(os.path.join(self.InputDir,'tilted_pairs.xmd'))

    def defineSteps(self):
        self.Db.insertStep('copyFiles',verifyfiles=[self.OutStack],InputFile=self.InSelFile,OutputStack=self.OutStack,OutputMetadata=self.OutMetadata)
        self.Db.insertStep('createAcquisition',InputFile=self.InSelFile,WorkingDir=self.WorkingDir,DoScale=self.DoScale,NewSize=self.NewSize)
        if self.DoScale:
            self.Db.insertStep('doScale',stack=self.OutStack,new_size=self.NewSize,Nproc=self.NumberOfMpi)
        if self.DoFourier:
            self.Db.insertStep('doFourier',stack=self.OutStack,freq_low=self.Freq_low,freq_high=self.Freq_high,freq_decay=self.Freq_decay,Nproc=self.NumberOfMpi)
        if self.DoGaussian:
            self.Db.insertStep('doGaussian',stack=self.OutStack,freq_sigma=self.Freq_sigma, Nproc=self.NumberOfMpi)
        if self.DoCrop:
            self.Db.insertStep('doCrop',stack=self.OutStack, cropSize=self.CropSize, tmpStack=self.tmpPath('tmpCrop.stk'))
        if self.DoRemoveDust:
            self.Db.insertStep('doRemoveDust',stack=self.OutStack,threshold=self.DustRemovalThreshold,Nproc=self.NumberOfMpi)
        if self.DoNorm:
            self.Db.insertStep('doNorm',stack=self.OutStack,normType=self.NormType,bgRadius=self.BackGroundRadius,Nproc=self.NumberOfMpi)
        if self.DoMask:
            if self.Substitute == "value":
                self.Substitute = str(self.SubstituteValue)
            params = "-i %s --substitute %s --mask %s " % (self.OutStack, self.Substitute, self.MaskType)
            if self.MaskType == 'raised_cosine':
                params += "-%d -%d" % (self.MaskRadius, self.MaskRadius + self.MaskRadiusOuter)
            elif self.MaskType == 'circular':
                params += '-%d' % self.MaskRadius
            else: # from file:
                params += self.MaskFile
            self.insertRunJobStep("xmipp_transform_mask", params)
        if self.TiltPair:
            if self.InSelFile.find("_untilted")!=-1:
                fnBase=os.path.split(self.InSelFile)[1]
                fnFamily=fnBase.replace("_untilted","")
                self.Db.insertStep('translateTiltPair',verifyfiles=[self.workingDirPath(fnFamily)],
                                   WorkingDir=self.WorkingDir,InputDir=self.InputDir,FnFamily=fnFamily,OutStack=self.OutStack)
                fnTilted=fnBase.replace("_untilted","_tilted")
                self.Db.insertStep('copyFile',verifyfiles=[self.workingDirPath(fnTilted)],source=self.InSelFile.replace('_untilted','_tilted'),
                                   dest=self.workingDirPath(fnTilted))
                self.Db.insertStep('copyFile',verifyfiles=[self.workingDirPath('tilted_pairs.xmd')],source=os.path.join(self.InputDir,'tilted_pairs.xmd'),
                                   dest=self.workingDirPath('tilted_pairs.xmd'))
        
    def validate(self):
        errors = []
        if self.DoScale and self.NewSize<=0:
            errors.append("New size for scale have not correctly set")
        if self.DoScale and self.TiltPair:
            errors.append("Cannot scale particles extracted from tilt pairs. Re-extract the particles at a different sampling rate, instead.")
        return errors

    def setStepMessage(self, stepMsg):
        step = len(self.messages) # assumed just one message before calling this function
        msg = "Step %d -> " % step
        msg += stepMsg % self.ParamsDict
        self.messages.append(msg)            
        
    def summary(self):
        self.messages = []
        self.messages.append("Steps applied to [%s]" % self.InSelFile)
        
        if self.DoScale:
            self.setStepMessage("Scale applied: new_size = %(NewSize)d")            
        if self.DoFourier:
            self.setStepMessage("Fourier filter applied: freq_low = %(Freq_low)f freq_high = %(Freq_high)f freq_decay = %(Freq_decay)f")            
        if self.DoGaussian:
            self.setStepMessage("Gaussian filter applied: freq_sigma = %(Freq_sigma)f")
        if self.DoCrop:
            self.setStepMessage("Crop applied: cropSize = %(CropSize)d")
        if self.DoRemoveDust:
            self.setStepMessage("Dust removal filter applied: threshold = %(DustRemovalThreshold)f")
        if self.DoNorm:
            if self.NormType == "OldXmipp":
                self.setStepMessage("Normalization applied: type = %(NormType)s")
            else:
                self.setStepMessage("Normalization applied: type = %(NormType)s backgroundRadius = %(BackGroundRadius)d")
        if self.DoMask:
            self.setStepMessage("Mask applied: mask file = %(MaskFile)s substituted value = %(Substitute)s")
        self.messages.append("Output: [%s]" % self.OutMetadata)
        return self.messages

    def visualize(self):
        if not os.path.exists(self.OutStack):
            showError("Error", "There is no result yet")
        else:
            runShowJ(self.OutMetadata)

def createAcquisition(log,InputFile,WorkingDir,DoScale,NewSize):
    if not DoScale:
        linkAcquisitionInfo("acquisition_info.xmd", InputFile, WorkingDir)
    else:
        dirSrc=os.path.dirname(InputFile)
        fnAcquisitionIn=os.path.join(dirSrc,"acquisition_info.xmd")
        if os.path.exists(fnAcquisitionIn):
            MD=xmipp.MetaData(fnAcquisitionIn)
            id=MD.firstObject()
            Ts=MD.getValue(xmipp.MDL_SAMPLINGRATE,id)
            (Xdim, Ydim, Zdim, Ndim) = xmipp.ImgSize(InputFile)
            downsampling=float(Xdim)/NewSize;
            MD.setValue(xmipp.MDL_SAMPLINGRATE,Ts*downsampling,id)
            MD.write(os.path.join(WorkingDir,"acquisition_info.xmd"))

def copyFiles(log,InputFile,OutputStack,OutputMetadata):
    runJob(log,"xmipp_image_convert","-i "+InputFile+" -o "+OutputStack)
    mDstack=xmipp.MetaData(OutputStack)
    if xmipp.FileName.isMetaData(xmipp.FileName(InputFile)):
        mDin=xmipp.MetaData(InputFile)
        mDin.removeDisabled()
        mDin.removeLabel(xmipp.MDL_IMAGE)
        mDin.removeLabel(xmipp.MDL_ZSCORE)
        mDaux=xmipp.MetaData(mDin)
        mDstack.merge(mDaux)
    mDstack.write(OutputMetadata)

def doScale(log,stack,new_size,Nproc):
    runJob(log,"xmipp_image_resize","-i %(stack)s --fourier %(new_size)d"%locals(),Nproc)

def doFourier(log,stack,freq_low,freq_high,freq_decay,Nproc):
    if freq_low==0:
        runJob(log,"xmipp_transform_filter","-i %(stack)s --fourier low_pass %(freq_high)f %(freq_decay)f"%locals(),Nproc)
    elif freq_high==0.5:
        runJob(log,"xmipp_transform_filter","-i %(stack)s --fourier high_pass %(freq_low)f %(freq_decay)f"%locals(),Nproc)
    else:
        runJob(log,"xmipp_transform_filter","-i %(stack)s --fourier band_pass %(freq_low)f %(freq_high)f %(freq_decay)f"%locals(),Nproc)

def doGaussian(log, stack, freq_sigma, Nproc):
    runJob(log,"xmipp_transform_filter","-i %(stack)s --fourier gaussian %(freq_sigma)f"%locals(),Nproc)

def doCrop(log, stack, cropSize, tmpStack):
    runJob(log,"xmipp_transform_window","-i %(stack)s --size %(cropSize)d -o %(tmpStack)s" % locals())
    moveFile(log, tmpStack, stack)
    
def doRemoveDust(log,stack,threshold,Nproc):
    runJob(log,"xmipp_transform_filter","-i %(stack)s --bad_pixels outliers %(threshold)f"%locals(),Nproc)

def doNorm(log,stack,normType,bgRadius,Nproc):
    if bgRadius <= 0:
        particleSize = xmipp.ImgSize(stack)[0]
        bgRadius=int(particleSize/2)
    if normType=="OldXmipp":
        runJob(log,"xmipp_transform_normalize","-i %(stack)s --method OldXmipp"%locals(),Nproc)
    elif normType=="NewXmipp":
        runJob(log,"xmipp_transform_normalize","-i %(stack)s --method NewXmipp --background circle %(bgRadius)d"%locals(),Nproc)
    else:
        runJob(log,"xmipp_transform_normalize","-i %(stack)s --method Ramp --background circle %(bgRadius)d"%locals(),Nproc)

def doMask(log,stack,maskFile,substitute,Nproc):
    runJob(log,"xmipp_transform_mask","-i %(stack)s --mask binary_file %(maskFile)s %(substitute)s"%locals(),Nproc)

def translateTiltPair(log,WorkingDir,InputDir,FnFamily,OutStack):
    MDoutStack=xmipp.MetaData(OutStack)
    MDfamily=xmipp.MetaData(os.path.join(InputDir,FnFamily))
    MDfamily.merge(MDoutStack)
    MDfamily.write(os.path.join(WorkingDir,FnFamily))

