#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using hierarchical clustering principles
# Author: Carlos Oscar Sanchez Sorzano, August 2011
#

import glob,os,re,sys,shutil,time
from protlib_base import *
from config_protocols import protDict
from protlib_xmipp import getMdSize
from protlib_utils import runJob, getListFromRangeString
from protlib_filesystem import createLink, deleteFile, linkAcquisitionInfo, getXmippPath
from xmipp import MetaData, MD_APPEND, MDL_IMAGE

class ProtAlignVolume(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.align_volume.name, scriptname, project)
        self.Import = 'from protocol_align_volume import *'    

    def defineSteps(self):
        self.Db.insertStep("linkAcquisitionInfo",InputFile=self.InputVolume,dirDest=self.WorkingDir)
        fnOut=self.workingDirPath("volume_aligned.vol")
        self.Db.insertStep("copyFile",source=self.InputVolume,dest=fnOut)
        args="--i1 %s --i2 %s --apply"%(self.ReferenceVolume,fnOut)
        if self.DoMask:
            if self.MaskType=="circular":
                args+=" --mask circular -%d"%self.MaskRadius
            else:
                args+=" --mask binary_file %s"%self.MaskFile
        if self.AlignmentMethod=="Fast Fourier":
            args+=" --frm"
        elif self.AlignmentMethod=="Local":
            args+=" --local --rot %f %f 1 --tilt %f %f 1 --psi %f %f 1 -x %f %f 1 -y %f %f 1 -z %f %f 1 --scale %f %f 0.005"%\
               (self.RotCurrent,self.RotCurrent,self.TiltCurrent,self.TiltCurrent,self.PsiCurrent,self.PsiCurrent,\
                self.ShiftXCurrent,self.ShiftXCurrent,self.ShiftYCurrent,self.ShiftYCurrent,\
                self.ShiftZCurrent,self.ShiftZCurrent, self.ScaleCurrent, self.ScaleCurrent)
        else:
            args+=" --rot %f %f %f --tilt %f %f %f --psi %f %f %f -x %f %f %f -y %f %f %f -z %f %f %f --scale %f %f %f"%\
               (self.Rot0,self.RotF,self.RotStep,self.Tilt0,self.TiltF,self.TiltStep,self.Psi0,self.PsiF,self.PsiStep,\
                self.X0,self.XF,self.XStep,self.Y0,self.YF,self.YStep,self.Z0,self.ZF,self.ZStep,
                self.Scale0,self.ScaleF,self.ScaleStep)
        self.insertRunJobStep("xmipp_volume_align", args)
        
        if self.AlignmentMethod=="Exhaustive+Local":
            args="--i1 %s --i2 %s --apply --local"%(self.ReferenceVolume,fnOut)
            self.insertRunJobStep("xmipp_volume_align", args)
    
    def summary(self):
        message=[]
        message.append("Reference volume: [%s] "%self.ReferenceVolume)
        message.append("Input volume: [%s] "%self.InputVolume)
        message.append("Alignment method: %s"%self.AlignmentMethod)
        return message
    
    def papers(self):
        papers=[]
        if self.AlignmentMethod=="Fast Fourier":
            papers.append('Chen, JSB (2013) [http://www.ncbi.nlm.nih.gov/pubmed/23523719]')
        return papers

    def validate(self):
        errors = []
        fnShAlignment=getXmippPath('lib/python2.7/site-packages/sh_alignment')
        if not os.path.exists(fnShAlignment) and self.AlignmentMethod=="Fast Fourier":
            errors.append("Fast Fourier can only be used if Xmipp has been installed with FRM support")
        return errors
    
    def visualize(self):
        from protlib_utils import runShowJ
        fnOut=self.workingDirPath("volume_aligned.vol")
        if os.path.exists(fnOut):
            runShowJ(fnOut)
