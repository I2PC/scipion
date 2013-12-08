#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based pre-processing of volumes 

# Author: Carlos Oscar, August 2013
#
from protlib_base import *
import os
from xmipp import MetaData, MDL_SAMPLINGRATE, MetaDataInfo
from os.path import exists, split, splitext
from protlib_utils import runJob, runShowJ
from protlib_filesystem import copyFile
import glob
from math import floor
from protlib_gui_ext import showError

class ProtPreprocessVolumes(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.preprocess_volume.name, scriptname, project)
        self.Import = 'from protocol_preprocess_volume import *'
        self.mdIn=MetaData(self.InModel)
        if self.mdIn.size()==1:
            self.OutModel=self.workingDirPath("volume.vol")
            self.singleVolume=True
        else:
            self.OutModel=self.workingDirPath("volumes.stk")
            self.singleVolume=False

    def defineSteps(self):
        self.insertStep('createAcquisition',WorkingDir=self.WorkingDir,Ts=self.FinalTs)
        if self.InitialTs!=self.FinalTs or self.FinalSize!=-1:
            if self.InitialTs!=self.FinalTs and self.FinalSize==-1:
                x, _, _, _, _ = MetaDataInfo(self.mdIn)    
                self.FinalSize=floor(x/(self.FinalTs/self.InitialTs))
            self.insertStep("changeSamplingRateAndOrBox",InModel=self.InModel,OutModel=self.OutModel,
                            SingleVolume=self.singleVolume, InitialTs=self.InitialTs,FinalTs=self.FinalTs,Size=self.FinalSize)
        else:
            if self.singleVolume:
                self.insertStep('copyFile',source=self.InModel,dest=self.OutModel)
            else:
                self.insertRunJobStep('xmipp_image_convert', "-i %s -o %s --save_metadata_stack --track_origin"%(self.InModel,self.OutModel),
                                      verifyFiles=[self.OutModel])

        if self.DoChangehand:
            self.insertStep("changeHand",OutModel=self.OutModel)
        
        if self.DoRandomize:
            self.insertStep("randomize",OutModel=self.OutModel,Ts=self.FinalTs,MaxResolution=self.MaxResolutionRandomize)
        
        if self.DoFilter:
            self.insertStep("filter",OutModel=self.OutModel,Ts=self.FinalTs,MaxResolution=self.MaxResolution)
        
        if self.DoSymmetrize:
            self.insertStep("symmetrize",OutModel=self.OutModel,Symmetry=self.Symmetry, SymmetryAggregation=self.SymmetryAggregation)
        
        if self.DoMask:
            self.insertStep("mask",OutModel=self.OutModel,MaskRadius=self.MaskRadius)
        
        if self.DoAdjust:
            self.insertStep("adjust",OutModel=self.OutModel,SetOfImages=self.SetOfImages)

        if self.DoNormalize:
            self.insertStep("normalize",OutModel=self.OutModel,MaskRadius=self.MaskRadiusNormalize)
        
        if self.DoThreshold:
            self.insertStep("threshold",OutModel=self.OutModel,Threshold=self.Threshold)
        
        if self.DoSegment:
            self.insertStep("segment",OutModel=self.OutModel,SegmentationType=self.SegmentationType,SegmentationMass=self.SegmentationMass,
                            Ts=self.FinalTs)
        
    def validate(self):
        errors = []
        if self.InitialTs<0:
            errors.append("Initial sampling rate must be provided")
        if self.FinalTs<0:
            errors.append("Initial sampling rate must be provided")
        maxFreq=2.0*self.FinalTs
        if self.DoRandomize and self.MaxResolutionRandomize<maxFreq:
            errors.append("Phase randomization cannot be performed beyond %f A (Nyquist)"%maxFreq)
        if self.DoFilter and self.MaxResolutionRandomize<maxFreq:
            errors.append("Low pass filtering cannot be performed beyond %f A (Nyquist)"%maxFreq)
        if self.DoAdjust and not self.singleVolume:
            errors.append("Gray adjusting is meant only for single volumes")
        if self.DoSegment and not self.singleVolume:
            errors.append("Segmentation is meant only for single volumes")
        return errors
        
    def summary(self):
        messages = []      
        messages.append("Input model: [%s]" % self.InModel)
        messages.append("Output: [%s]" % self.OutModel)
        messages.append("Operations: ")
        if self.InitialTs!=self.FinalTs:
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

    def papers(self):
        papers=[]
        if self.DoNormalize:
            papers.append('Sorzano, Ultramic (2004) [http://www.ncbi.nlm.nih.gov/pubmed/15450658]')
        if self.InitialTs!=self.FinalTs:
            papers.append('Sorzano, IEEE WISP (2009) [http://ieeexplore.ieee.org/xpl/login.jsp?arnumber=5286563]')
        if self.DoRandomize:
            papers.append('Chen, Ultramic (2013) [http://www.ncbi.nlm.nih.gov/pubmed/23872039]')
        return papers

    def visualize(self):
        from protlib_utils import runShowJ
        if os.path.exists(self.OutModel):
            runShowJ(self.OutModel)

def createAcquisition(log,WorkingDir,Ts):
    md=MetaData()
    id=md.addObject()
    md.setValue(MDL_SAMPLINGRATE,float(Ts),id)
    md.write(os.path.join(WorkingDir,"acquisition_info.xmd"))

def window(log,input,OutModel,SingleVolume,Size):
    if Size>0:
        args="-i %s --size %d"%(input,int(Size))
        if not SingleVolume:
            args+=" --save_metadata_stack"
        if input!=OutModel:
            args+=" -o %s"%OutModel
        runJob(log,"xmipp_transform_window",args)

def scale(log,input,OutModel,SingleVolume,scale):
    args="-i %s -o %s --scale %f --dont_wrap"%(input,OutModel,scale)
    if not SingleVolume:
        args+=" --save_metadata_stack"
    runJob(log,"xmipp_transform_geometry",args)

def changeSamplingRateAndOrBox(log,InModel,OutModel,SingleVolume,InitialTs,FinalTs,Size):
    if InitialTs==FinalTs:
        window(log,InModel,OutModel,SingleVolume,Size)
    elif InitialTs<FinalTs:
        scale(log,InModel,OutModel,SingleVolume,InitialTs/FinalTs)
        window(log,OutModel,OutModel,SingleVolume,Size)
    else:
        window(log,InModel,OutModel,SingleVolume,Size)
        scale(log,OutModel,OutModel,SingleVolume,InitialTs/FinalTs)

def randomize(log,OutModel,Ts,MaxResolution):
    runJob(log,"xmipp_transform_randomize_phases","-i %s --freq continuous %f %f"%(OutModel,float(MaxResolution),float(Ts)))

def filter(log,OutModel,Ts,MaxResolution):
    runJob(log,"xmipp_transform_filter","-i %s --fourier low_pass %f --sampling %f"%(OutModel,float(MaxResolution),float(Ts)))

def mask(log,OutModel,MaskRadius):
    if MaskRadius==-1:
        md=MetaData(OutModel)
        Xdim, _, _, _, _ = MetaDataInfo(md)    
        MaskRadius=Xdim/2
    runJob(log,"xmipp_transform_mask","-i %s --mask circular %d"%(OutModel,-int(MaskRadius)))

def symmetrize(log,OutModel,Symmetry,SymmetryAggregation):
    if Symmetry!='c1':
        args="-i %s --sym %s"%(OutModel,Symmetry)
        if SymmetryAggregation=="sum":
            args+=" --sum"
        runJob(log,"xmipp_transform_symmetrize",args)

def adjust(log,OutModel,SetOfImages):
    runJob(log,"xmipp_adjust_volume_grey_levels","-i %s -m %s"%(OutModel,SetOfImages))

def normalize(log,OutModel,MaskRadius):
    if MaskRadius==-1:
        md=MetaData(OutModel)
        Xdim, _, _, _, _ = MetaDataInfo(md)    
        MaskRadius=Xdim/2
    runJob(log,"xmipp_transform_normalize","-i %s --method NewXmipp --background circle %d"%(OutModel,int(MaskRadius)))

def threshold(log,OutModel,Threshold):
    runJob(log,"xmipp_transform_threshold","-i %s --select below %f --substitute value 0"%(OutModel,float(Threshold)))

def changeHand(log,OutModel):
    runJob(log,"xmipp_transform_mirror","-i %s --flipX"%OutModel)

def segment(log,OutModel,SegmentationType,SegmentationMass,Ts):
    WorkingDir=os.path.split(OutModel)[0]
    fnMask=os.path.join(WorkingDir,"mask.vol")
    args="-i %s -o %s --method "%(OutModel,fnMask)
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
        runJob(log,"xmipp_transform_mask","-i %s --mask binary_file %s"%(OutModel,fnMask))
