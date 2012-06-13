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
from protlib_filesystem import deleteFile,linkAcquisitionInfoIfPresent
import glob
from protlib_gui_ext import showError

class ProtPreprocessParticles(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.preprocess_particles.name, scriptname, project)
        self.Import = 'from protocol_preprocess_particles import *'
        file=os.path.split(self.InSelFile)[1]
        self.OutSelFile=self.workingDirPath(os.path.splitext(file)[0]+".stk")

    def defineSteps(self):
        self.Db.insertStep('copyFiles',verifyfiles=[self.OutSelFile],InputFile=self.InSelFile,OutputFile=self.OutSelFile)
        self.Db.insertStep('createAcquisition',InputFile=self.InSelFile,WorkingDir=self.WorkingDir,DoScale=self.DoScale,NewSize=self.NewSize)
        if self.DoScale:
            self.Db.insertStep('doScale',stack=self.OutSelFile,new_size=self.NewSize,Nproc=self.NumberOfMpi)
        if self.DoFourier:
            self.Db.insertStep('doFourier',stack=self.OutSelFile,freq_low=self.Freq_low,freq_high=self.Freq_high,freq_decay=self.Freq_decay,Nproc=self.NumberOfMpi)
        if self.DoGaussian:
            self.Db.insertStep('doGaussian',stack=self.OutSelFile,freq_sigma=self.Freq_sigma,Nproc=self.NumberOfMpi)
        if self.DoRemoveDust:
            self.Db.insertStep('doRemoveDust',stack=self.OutSelFile,threshold=self.DustRemovalThreshold,Nproc=self.NumberOfMpi)
        if self.DoNorm:
            self.Db.insertStep('doNorm',stack=self.OutSelFile,normType=self.NormType,bgRadius=self.BackGroundRadius,Nproc=self.NumberOfMpi)
        if self.DoMask:
            self.Db.insertStep('doMask',stack=self.OutSelFile,maskFile=self.MaskFile,substitute=self.Substitute,Nproc=self.NumberOfMpi)
        
    def validate(self):
        errors = []
        if self.DoScale and self.NewSize<=0:
            errors.append("New size for scale have not correctly set")
        return errors

    def summary(self):
        message=[]
        step=1
        message.append("Steps applied to [%s]"%self.InSelFile)
        if self.DoScale:
            message.append("Step %d -> scale applied: new_size=%d"%(step,self.NewSize))
            step+=1
        if self.DoFourier:
            message.append("Step %d -> Fourier filter applied: freq_low=%f freq_high=%f freq_decay=%f"%(step,self.Freq_low,self.Freq_high,self.Freq_decay))
            step+=1
        if self.DoGaussian:
            message.append("Step %d -> Gaussian filter applied: freq_sigma=%f"%(step,self.Freq_sigma))
            step+=1
        if self.DoRemoveDust:
            message.append("Step %d -> Dust removal filter applied: threshold=%f"%(step,self.DustRemovalThreshold))
            step+=1
        if self.DoNorm:
            if self.NormType=="OldXmipp":
                message.append("Step %d -> Normalization applied: type=%s"%(step,self.NormType))
            else:
                message.append("Step %d -> Normalization applied: type=%s backgroundRadius=%d"%(step,self.NormType,self.BackGroundRadius))
            step+=1
        if self.DoMask:
            message.append("Step %d -> Mask applied: mask file=%f substituted value=%f"%(step,self.MaskFile,self.Substitute))
            step+=1
        return message

    def visualize(self):
        if not os.path.exists(self.OutSelFile):
            showError("Error", "There is no result yet")
        else:
            runShowJ(self.OutSelFile)

def createAcquisition(log,InputFile,WorkingDir,DoScale,NewSize):
    if not DoScale:
        linkAcquisitionInfoIfPresent("acquisition_info.xmd", InputFile, WorkingDir)
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

def copyFiles(log,InputFile,OutputFile):
    runJob(log,"xmipp_image_convert","-i "+InputFile+" -o "+OutputFile)
    if xmipp.FileName.isMetaData(xmipp.FileName(InputFile)):
        mDin=xmipp.MetaData(InputFile)
        mDin.removeLabel(xmipp.MDL_IMAGE)
        mDin.removeLabel(xmipp.MDL_ZSCORE)
        mDaux=xmipp.MetaData(mDin)
        mDstack=xmipp.MetaData(OutputFile)
        mDaux.merge(mDstack)
        mDaux.write(os.path.splitext(OutputFile)[0]+".xmd")

def doScale(log,stack,new_size,Nproc):
    runJob(log,"xmipp_transform_geometry","-i %(stack)s --scale dim %(new_size)d"%locals(),Nproc)

def doFourier(log,stack,freq_low,freq_high,freq_decay,Nproc):
    if freq_low==0:
        runJob(log,"xmipp_transform_filter","-i %(stack)s --fourier low_pass %(freq_high)f %(freq_decay)f"%locals(),Nproc)
    elif freq_high==0.5:
        runJob(log,"xmipp_transform_filter","-i %(stack)s --fourier high_pass %(freq_low)f %(freq_decay)f"%locals(),Nproc)
    else:
        runJob(log,"xmipp_transform_filter","-i %(stack)s --fourier band_pass %(freq_low)f %(freq_high)f %(freq_decay)f"%locals(),Nproc)

def doGaussian(log,stack,freq_sigma,Nproc):
    runJob(log,"xmipp_transform_filter","-i %(stack)s --fourier gaussian %(freq_sigma)f"%locals(),Nproc)

def doRemoveDust(log,stack,threshold,Nproc):
    runJob(log,"xmipp_transform_filter","-i %(stack)s --bad_pixels outliers %(threshold)f"%locals(),Nproc)

def doNorm(log,stack,normType,bgRadius,Nproc):
    if bgRadius==0:
        particleSize=xmipp.ImgSize(stack)[0]
        bgRadius=int(particleSize/2)
    if normType=="OldXmipp":
        runJob(log,"xmipp_transform_normalize","-i %(stack)s --method OldXmipp"%locals(),Nproc)
    elif normType=="NewXmipp":
        runJob(log,"xmipp_transform_normalize","-i %(stack)s --method NewXmipp --background circle %(bgRadius)d"%locals(),Nproc)
    else:
        runJob(log,"xmipp_transform_normalize","-i %(stack)s --method Ramp --background circle %(bgRadius)d"%locals(),Nproc)

def doMask(log,stack,maskFile,substitute,Nproc):
    runJob(log,"xmipp_transform_mask","-i %(stack)s --mask binary_file %(maskFile)s %(substitute)s"%locals(),Nproc)

