#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based pre-processing of volumes 

# Author: Carlos Oscar, August 2013
#
from protlib_base import *
import os
from xmipp import createEmptyFile
from os.path import exists, split, splitext
from protlib_utils import runJob, runShowJ
from protlib_filesystem import copyFile
import glob
from protlib_gui_ext import showError

class ProtCreateVolumeMask(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.create_volume_mask.name, scriptname, project)
        self.Import = 'from protocol_create_volume_mask import *'
        self.OutMask=self.workingDirPath("mask.vol")

    def defineSteps(self):
        # Create mask
        if self.MaskSource=="Volume":
            if self.VolumeOperation=="Threshold":
                self.insertStep("threshold",WorkingDir=self.WorkingDir,InModel=self.InModel,Threshold=self.Threshold)
            elif self.VolumeOperation=="Segment":
                self.insertStep("segment",WorkingDir=self.WorkingDir,InModel=self.InModel,SegmentationType=self.SegmentationType,
                                SegmentationMass=self.SegmentationMass,Ts=self.Ts)
        elif self.MaskSource=="Geometry":
            self.insertStep("createEmptyVolume",Mask=self.OutMask,MaskSize=self.MaskSize)
            if self.SimpleOperation=="Sphere":
                self.insertStep("createSphere",Mask=self.OutMask,Radius=self.Radius)
            elif self.SimpleOperation=="Box":
                self.insertStep("createBox",Mask=self.OutMask,BoxSize=self.BoxSize)
            elif self.SimpleOperation=="Crown":
                self.insertStep("createCrown",Mask=self.OutMask,InnerRadius=self.InnerRadius,OuterRadius=self.OuterRadius)
            elif self.SimpleOperation=="Cylinder":
                self.insertStep("createCylinder",Mask=self.OutMask,Radius=self.Radius,Height=self.Height)
            elif self.SimpleOperation=="Gaussian":
                self.insertStep("createGaussian",Mask=self.OutMask,Sigma=self.Sigma)
            elif self.SimpleOperation=="Raised cosine":
                self.insertStep("createRaisedCosine",Mask=self.OutMask,InnerRadius=self.InnerRadius,OuterRadius=self.OuterRadius)
            elif self.SimpleOperation=="Raised crown":
                self.insertStep("createRaisedCrown",Mask=self.OutMask,InnerRadius=self.InnerRadius,OuterRadius=self.OuterRadius,
                                PixelWidth=self.PixelWidth)
        elif self.MaskSource=="Binary mask file":
            self.insertStep("copyFile",source=self.BinaryMask,dest=self.OutMask)
        
        # Postprocess mask
        if self.DoSmall:
            self.insertStep("removeSmall",Mask=self.OutMask,SmallSize=self.SmallSize)
        if self.DoBig:
            self.insertStep("keepLargest",Mask=self.OutMask)
        if self.DoSymmetrize:
            self.insertStep("symmetrize",Mask=self.OutMask,Symmetry=self.Symmetry)
        if self.DoMorphological:
            self.insertStep("morphology",Mask=self.OutMask,MorphologicalOperation=self.MorphologicalOperation,
                            ElementSize=self.ElementSize)
        if self.DoInvert:
            self.insertStep("invert",Mask=self.OutMask)        
        if self.DoSmooth:
            self.insertStep("convolve",Mask=self.OutMask,SigmaConvolution=self.SigmaConvolution)
        
        # Apply Mask
        if self.MaskSource=="Volume" and self.ApplyMask:
            self.insertStep("applyMask",Mask=self.OutMask,WorkingDir=self.WorkingDir,InModel=self.InModel)

    def summary(self):
        messages = []      
        messages.append("Output: [%s]" % self.OutMask)
        messages.append("<Mask creation>")
        if self.MaskSource=="Binary mask file":
            messages.append("   Input: [%s]"%self.BinaryMask)
        elif self.MaskSource=="Volume":
            messages.append("   Input: [%s]"%self.InModel)
            if self.VolumeOperation=="Threshold":
                messages.append("   Thresholding %f"%self.Threshold)
            elif self.VolumeOperation=="Segment":
                if self.SegmentationType=="Automatic":
                    messages.append("   Automatically segmented")
                else:
                    m="   Segmented to a mass of "
                    if self.SegmentationType=="Voxel mass":
                        m+="%d voxels"%(int(SegmentationMass))
                    elif self.SegmentationType=="Aminoacid mass":
                        m+="%d aminoacids"%(int(SegmentationMass))
                    elif self.SegmentationType=="Dalton mass":
                        m+="%d daltons"%(int(SegmentationMass))
                    messages.append(m)
        elif self.MaskSource=="Geometry":
            if self.SimpleOperation=="Sphere":
                messages.append("   Sphere of radius %d"%int(self.Radius))
            elif self.SimpleOperation=="Box":
                messages.append("   Box of size %d"%int(self.BoxSize))
            elif self.SimpleOperation=="Crown":
                messages.append("   Crown between %d and %d"%(int(self.InnerRadius),int(self.OuterRadius)))
            elif self.SimpleOperation=="Cylinder":
                messages.append("   Cylinder of radius %f and height %f"%(int(self.Radius),int(self.Height)))
            elif self.SimpleOperation=="Gaussian":
                messages.append("   Gaussian of sigma %f"%(self.Sigma))
            elif self.SimpleOperation=="Raised cosine":
                messages.append("   Raised cosine between %f and %f"%(self.InnerRadius,self.OuterRadius))
            elif self.SimpleOperation=="Raised crown":
                messages.append("   Raised crown between %f and %f (width=%f)"%(self.InnerRadius,self.OuterRadius,self.PixelWidth))
        messages.append("<Mask processing>")
        if self.DoSmall:
            messages.append("   Removing components smaller than %d"%(int(self.SmallSize)))
        if self.DoBig:
            messages.append("   Keeping largest component")
        if self.DoSymmetrize:
            messages.append("   Symmetrized %s"%self.Symmetry)
        if self.DoMorphological:
            messages.append("   Morphological operation: %s"%self.MorphologicalOperation)
        if self.DoInvert:
            messages.append("   Inverted")
        if self.DoSmooth:
            messages.append("   Smoothed (sigma=%f)"%float(self.SigmaConvolution))
        return messages

    def visualize(self):
        from protlib_utils import runShowJ
        if os.path.exists(self.OutMask):
            runShowJ(self.OutMask)
        fnMasked=self.workingDirPath("volume_masked.vol")
        if os.path.exists(fnMasked):
            runShowJ(fnMasked)

def threshold(log,WorkingDir,InModel,Threshold):
    runJob(log,"xmipp_transform_threshold","-i %s -o %s/mask.vol --select below %f --substitute binarize"%(InModel,WorkingDir,float(Threshold)))

def segment(log,WorkingDir,InModel,SegmentationType,SegmentationMass,Ts):
    fnMask="%s/mask.vol"%WorkingDir
    args="-i %s -o %s --method "%(InModel,fnMask)
    if SegmentationType=="Voxel mass":
        args+="voxel_mass %d"%(int(SegmentationMass))
    elif SegmentationType=="Aminoacid mass":
        args+="aa_mass %d %f"%(int(SegmentationMass),float(Ts))
    elif SegmentationType=="Dalton mass":
        args+="dalton_mass %d %f"%(int(SegmentationMass),float(Ts))
    else:
        args+="otsu"
    runJob(log,"xmipp_volume_segment",args)

def createEmptyVolume(log,Mask,MaskSize):
    createEmptyFile(Mask,MaskSize,MaskSize,MaskSize)

def createSphere(log,Mask,Radius):
    runJob(log,"xmipp_transform_mask","-i %s --mask circular %d --create_mask %s"%(Mask,-int(Radius),Mask))

def createBox(log,Mask,BoxSize):
    b=int(BoxSize)
    runJob(log,"xmipp_transform_mask","-i %s --mask rectangular %d %d %d --create_mask %s"%(Mask,-b,-b,-b,Mask))

def createCrown(log,Mask,InnerRadius,OuterRadius):
    runJob(log,"xmipp_transform_mask","-i %s --mask crown %d %d --create_mask %s"%(Mask,-int(InnerRadius),-int(OuterRadius),Mask))

def createCylinder(log,Mask,Radius,Height):
    runJob(log,"xmipp_transform_mask","-i %s --mask cylinder %d %d --create_mask %s"%(Mask,-int(Radius),-int(Height),Mask))    

def createGaussian(log,Mask,Sigma):
    runJob(log,"xmipp_transform_mask","-i %s --mask gaussian %f --create_mask %s"%(Mask,-float(Sigma),Mask))    

def createRaisedCosine(log,Mask,InnerRadius,OuterRadius):
    runJob(log,"xmipp_transform_mask","-i %s --mask raised_cosine %d %d --create_mask %s"%(Mask,-int(InnerRadius),-int(OuterRadius),Mask))    

def createRaisedCrown(log,Mask,InnerRadius,OuterRadius,PixelWidth):
    runJob(log,"xmipp_transform_mask","-i %s --mask raised_crown %d %d %d --create_mask %s"%(Mask,-int(InnerRadius),-int(OuterRadius),
                                                                                             int(PixelWidth),Mask))    

def keepLargest(log,Mask):
    runJob(log,"xmipp_transform_morphology","-i %s --binaryOperation keepBiggest"%Mask)

def removeSmall(log,Mask,SmallSize):
    runJob(log,"xmipp_transform_morphology","-i %s --binaryOperation removeSmall %d"%(Mask,int(SmallSize)))

def symmetrize(log,Mask,Symmetry):
    if Symmetry!='c1':
        args="-i %s --sym %s --dont_wrap"%(Mask,Symmetry)
        runJob(log,"xmipp_transform_symmetrize",args)

def morphology(log,Mask,MorphologicalOperation,ElementSize):
    runJob(log,"xmipp_transform_morphology","-i %s --binaryOperation %s --size %d"%(Mask,MorphologicalOperation,int(ElementSize)))

def invert(log,Mask):
    runJob(log,"xmipp_image_operate","-i %s --mult -1"%Mask)
    runJob(log,"xmipp_image_operate","-i %s --plus  1"%Mask)

def convolve(log,Mask,SigmaConvolution):
    runJob(log,"xmipp_transform_filter","-i %s --fourier real_gaussian %f"%(Mask,float(SigmaConvolution)))

def applyMask(log,Mask,WorkingDir,InModel):
    runJob(log,"xmipp_transform_mask","-i %s --mask binary_file %s -o %s"%(InModel,Mask,os.path.join(WorkingDir,"volume_masked.vol")))
