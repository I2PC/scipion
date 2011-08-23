#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based pre-processing of single-particles: 

# Author: Carlos Oscar, August 2011
#
from protlib_base import *
import xmipp
import os
from protlib_utils import runJob
from protlib_filesystem import deleteFile
import glob

class ProtPreprocessParticles(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.preprocess_particles.name, scriptname, project)
        self.Import = 'from protocol_preprocess_particles import *'
        self.extractionDir=getWorkingDirFromRunName(self.ExtractionRun)

    def defineSteps(self):
        fnOut=os.path.join(self.WorkingDir,"micrographs.sel")
        self.Db.insertStep('createLink',verifyfiles=[fnOut],source=os.path.join(self.extractionDir,"micrographs.sel"),dest=fnOut)

        fnIn=os.path.join(self.extractionDir,self.Family)
        fnOut=os.path.join(self.WorkingDir,self.Family)
        self.Db.insertStep('copyDir',verifyfiles=[fnOut],source=fnIn,dest=fnOut)
        fnIn+=".sel"
        fnOut+=".sel"
        self.Db.insertStep('transposeMetadata',verifyfiles=[fnOut],source=fnIn,sourceWorkingDir=self.extractionDir,
                           destWorkingDir=self.WorkingDir,dest=fnOut)

        if self.DoFourier:
            self.Db.insertStep('doFourier',stack=fnOut,freq_low=self.Freq_low,freq_high=self.Freq_high,freq_decay=self.Freq_decay,
                               Nproc=self.NumberOfMpi)
        if self.DoGaussian:
            self.Db.insertStep('doGaussian',stack=fnOut,freq_sigma=self.Freq_sigma,Nproc=self.NumberOfMpi)
        if self.DoRemoveDust:
            self.Db.insertStep('doRemoveDust',stack=fnOut,threshold=self.DustRemovalThreshold,Nproc=self.NumberOfMpi)
        if self.DoNorm:
            self.Db.insertStep('doNorm',stack=fnOut,normType=self.NormType,bgRadius=self.BackGroundRadius,Nproc=self.NumberOfMpi)
        
    def validate(self):
        errors = []
        fnSel=os.path.join(self.extractionDir,self.Family+".sel")
        if not os.path.exists(fnSel):
            errors.append("Cannot find "+fnSel)
        return errors

    def summary(self):
        message=[]
        message.append("Still to do")
        return message


    def visualize(self):
        selfile=os.path.join(self.WorkingDir,self.Family+".sel")
        if not os.path.exists(selfile):
            import tkMessageBox
            tkMessageBox.showerror("Error", "There is no result yet")                    
        os.system("xmipp_showj -i "+selfile+" --memory 1024m &")

def transposeMetadata(log,source,sourceWorkingDir,destWorkingDir,dest):
    mD=xmipp.MetaData(source)
    mD.removeLabel(xmipp.MDL_ZSCORE)
    mD.operate("image=replace(image,'%s','%s')"%(sourceWorkingDir,destWorkingDir))
    mD.write(dest)

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
    if normType=="OldXmipp":
        runJob(log,"xmipp_transform_normalize","-i %(stack)s --method OldXmipp"%locals(),Nproc)
    elif normType=="NewXmipp":
        runJob(log,"xmipp_transform_normalize","-i %(stack)s --method NewXmipp --background circle %(bgRadius)d"%locals(),Nproc)
    else:
        runJob(log,"xmipp_transform_normalize","-i %(stack)s --method Ramp --background circle %(bgRadius)d"%locals(),Nproc)
