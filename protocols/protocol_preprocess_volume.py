#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based pre-processing of volumes 

# Author: Carlos Oscar, August 2013
#
from protlib_base import *
import os
from xmipp import MetaData, MDL_SAMPLINGRATE
from os.path import exists, split, splitext
from protlib_utils import runJob, runShowJ
from protlib_filesystem import copyFile
import glob
from protlib_gui_ext import showError

class ProtPreprocessVolumes(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.preprocess_volume.name, scriptname, project)
        self.Import = 'from protocol_preprocess_volume import *'
        self.OutModel=self.workingDirPath("volume.vol")

    def defineSteps(self):
        self.insertStep('createAcquisition',WorkingDir=self.WorkingDir,Ts=self.FinalTs)
        if self.Action=="convert from PDB":
            self.insertStep("convertFromPDB",InModel=self.InModel,WorkingDir=self.WorkingDir,Ts=self.FinalTs,Size=self.FinalSize)
        else:
            if self.InitialTs!=self.FinalTs or self.FinalSize!=-1:
                self.insertStep("changeSamplingRateAndOrBox",InModel=self.InModel,WorkingDir=self.WorkingDir,InitialTs=self.InitialTs,
                                FinalTs=self.FinalTs,Size=self.FinalSize)
            else:
                self.insertStep("copyFile",source=self.InModel,dest=self.OutModel)

        if self.DoChangehand:
            self.insertStep("changeHand",WorkingDir=self.WorkingDir)
        
        if self.DoRandomize:
            self.insertStep("randomize",WorkingDir=self.WorkingDir,Ts=self.FinalTs,MaxResolution=self.MaxResolutionRandomize)
        
        if self.DoFilter:
            self.insertStep("filter",WorkingDir=self.WorkingDir,Ts=self.FinalTs,MaxResolution=self.MaxResolution)
        
        if self.DoSymmetrize:
            self.insertStep("symmetrize",WorkingDir=self.WorkingDir,Symmetry=self.Symmetry, SymmetryAggregation=self.SymmetryAggregation)
        
        if self.DoMask:
            self.insertStep("mask",WorkingDir=self.WorkingDir,MaskRadius=self.MaskRadius)
        
        if self.DoAdjust:
            self.insertStep("adjust",WorkingDir=self.WorkingDir,SetOfImages=self.SetOfImages)

        if self.DoNormalize:
            self.insertStep("normalize",WorkingDir=self.WorkingDir,MaskRadius=self.MaskRadiusNormalize)
        
        if self.DoThreshold:
            self.insertStep("threshold",WorkingDir=self.WorkingDir,Threshold=self.Threshold)
        
        if self.DoSegment:
            self.insertStep("segment",WorkingDir=self.WorkingDir,SegmentationType=self.SegmentationType,SegmentationMass=self.SegmentationMass,
                            Ts=self.FinalTs)
        
    def validate(self):
        errors = []
        if self.Action=="preprocess density volume":
            if self.InitialTs<0:
                errors.append("Initial sampling rate must be provided")
            if self.FinalTs<0:
                errors.append("Initial sampling rate must be provided")
        maxFreq=2.0*self.FinalTs
        if self.DoRandomize and self.MaxResolutionRandomize<maxFreq:
            errors.append("Phase randomization cannot be performed beyond %f A (Nyquist)"%maxFreq)
        if self.DoFilter and self.MaxResolutionRandomize<maxFreq:
            errors.append("Low pass filtering cannot be performed beyond %f A (Nyquist)"%maxFreq)
        return errors
        
    def summary(self):
        messages = []      
        messages.append("Input model: [%s]" % self.InModel)
        messages.append("Output: [%s]" % self.OutModel)
        messages.append("Operations: ")
        if self.Action=="preprocess density volume" and self.InitialTs!=self.FinalTs:
            messages.append("   Sampling rate changed from %f to %f"%(float(self.InitialTs),float(self.FinalTs)))
        if self.FinalSize>0:
            messages.append("   Volume boxed to %dx%dx%d voxels"%(int(self.FinalSize),int(self.FinalSize),int(self.FinalSize)))
        if self.DoChangehand:
            messages.append("   Hand changed")
        if self.DoRandomize:
            messages.append("   Phases randomized beyond %f A"%float(self.MaxResolutionRandomize))
        if self.DoFilter:
            messages.append("   Filtered to %f A"%float(self.MaxResolution))
        if self.DoSymmetrize:
            messages.append("   Symmetrized %s"%self.Symmetry)
        if self.DoMask:
            if self.MaskRadius>0:
                messages.append("   Masked within a sphere of radius %d"%self.MaskRadius)
            else:
                messages.append("   Masked within the maximal sphere fitting in its box")
        if self.DoAdjust:
            messages.append("   Gray values adjusted to fit [%s]"%self.SetOfImages)
        if self.DoNormalize:
            if self.MaskRadiusNormalize>0:
                messages.append("   Normalized with background beyond radius %d"%self.MaskRadiusNormalize)
            else:
                messages.append("   Normalized with background beyond the maximal fitting sphere")
        if self.DoThreshold:
            messages.append("   Thresholded below %f"%(float(self.Threshold)))
        if self.DoSegment:
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
        return messages

    def visualize(self):
        from protlib_utils import runShowJ
        fnVolume=self.workingDirPath('volume.vol')
        if os.path.exists(fnVolume):
            runShowJ(fnVolume)

def createAcquisition(log,WorkingDir,Ts):
    md=MetaData()
    id=md.addObject()
    md.setValue(MDL_SAMPLINGRATE,float(Ts),id)
    md.write(os.path.join(WorkingDir,"acquisition_info.xmd"))

def convertFromPDB(log,InModel,WorkingDir,Ts,Size):
    args="-i %s -o %s/volume --centerPDB"%(InModel,WorkingDir)
    if Size>0:
        args+=" --size %d"%(int(Size))
    if Ts>4:
        args+=" --poor_Gaussian"
    runJob(log,"xmipp_volume_from_pdb",args)

def changeSamplingRateAndOrBox(log,InModel,WorkingDir,InitialTs,FinalTs,Size):
    input=InModel
    fnVol="%s/volume.vol"%WorkingDir
    if InitialTs!=FinalTs:
        runJob(log,"xmipp_transform_geometry","-i %s -o %s --scale %f --dont_wrap"%(input,fnVol,InitialTs/FinalTs))
        input=fnVol
    if Size>0:
        runJob(log,"xmipp_transform_window","-i %s -o %s --size %d"%(input,fnVol,int(Size)))

def randomize(log,WorkingDir,Ts,MaxResolution):
    runJob(log,"xmipp_transform_randomize_phases","-i %s/volume.vol --freq continuous %f %f"%(WorkingDir,float(MaxResolution),float(Ts)))

def filter(log,WorkingDir,Ts,MaxResolution):
    runJob(log,"xmipp_transform_filter","-i %s/volume.vol --fourier low_pass %f --sampling %f"%(WorkingDir,float(MaxResolution),float(Ts)))

def mask(log,WorkingDir,MaskRadius):
    fnVol='%s/volume.vol'%WorkingDir
    if MaskRadius==-1:
        (Xdim, Ydim, Zdim, Ndim) = getImageSize(fnVol)
        MaskRadius=Xdim/2
    runJob(log,"xmipp_transform_mask","-i %s --mask circular %d"%(fnVol,-int(MaskRadius)))

def symmetrize(log,WorkingDir,Symmetry,SymmetryAggregation):
    if Symmetry!='c1':
        args="-i %s/volume.vol --sym %s"%(WorkingDir,Symmetry)
        if SymmetryAggregation=="sum":
            args+=" --sum"
        runJob(log,"xmipp_transform_symmetrize",args)

def adjust(log,WorkingDir,SetOfImages):
    runJob(log,"xmipp_adjust_volume_grey_levels","-i %s/volume.vol -m %s"%(WorkingDir,SetOfImages))

def normalize(log,WorkingDir,MaskRadius):
    fnVol='%s/volume.vol'%WorkingDir
    if MaskRadius==-1:
        (Xdim, Ydim, Zdim, Ndim) = getImageSize(fnVol)
        MaskRadius=Xdim/2
    runJob(log,"xmipp_transform_normalize","-i %s --method NewXmipp --background circle %d"%(fnVol,int(MaskRadius)))

def threshold(log,WorkingDir,Threshold):
    runJob(log,"xmipp_transform_threshold","-i %s/volume.vol --select below %f --substitute value 0"%(WorkingDir,float(Threshold)))

def changeHand(log,WorkingDir):
    runJob(log,"xmipp_transform_mirror","-i %s/volume.vol --flipX"%(WorkingDir))

def segment(log,WorkingDir,SegmentationType,SegmentationMass,Ts):
    fnIn="%s/volume.vol"%WorkingDir
    fnMask="%s/mask.vol"%WorkingDir
    args="-i %s -o %s --method "%(fnIn,fnMask)
    if SegmentationType=="Voxel mass":
        args+="voxel_mass %d"%(int(SegmentationMass))
    elif SegmentationType=="Aminoacid mass":
        args+="aa_mass %d %f"%(int(SegmentationMass),float(Ts))
    elif SegmentationType=="Dalton mass":
        args+="dalton_mass %d %f"%(int(SegmentationMass),float(Ts))
    else:
        args+="otsu"
    runJob(log,"xmipp_volume_segment",args)
    if os.path.exists(fnMask):
        runJob(log,"xmipp_transform_mask","-i %s --mask binary_file %s"%(fnIn,fnMask))
