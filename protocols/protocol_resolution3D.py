#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based pre-processing of volumes 

# Author: Carlos Oscar, August 2013
#
from protlib_base import *
import os
from protlib_utils import runJob, runShowJ
from protlib_filesystem import linkAcquisitionInfo
import glob
from xmipp import MetaData, MDL_SAMPLINGRATE, MDL_RESOLUTION_FREQREAL, MDL_RESOLUTION_FRC
from protlib_gui_ext import showError

class ProtResolution3D(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.resolution3D.name, scriptname, project)
        self.Import = 'from protocol_resolution3D import *'

    def defineSteps(self):
        self.insertStep('linkAcquisitionInfo',InputFile=self.InputVol, dirDest=self.WorkingDir)
        self.insertStep('fsc',ReferenceVol=self.ReferenceVol,InputVol=self.InputVol,WorkingDir=self.WorkingDir)
        
    def validate(self):
        errors = []
        return errors
        
    def summary(self):
        messages = []      
        messages.append("Reference volume: [%s]" % self.ReferenceVol)
        messages.append("Input     volume: [%s]" % self.InputVol)
        fnFSC=self.workingDirPath("fsc.xmd")
        if os.path.exists(fnFSC):
            md = MetaData(fnFSC)
            resolution = [md.getValue(MDL_RESOLUTION_FREQREAL, id) for id in md]
            frc = [md.getValue(MDL_RESOLUTION_FRC, id) for id in md]
            for i in range(len(frc)):
                if frc[i]<0.5:
                    messages.append("Resolution FSC(0.5) = %f"%resolution[i])
                    break
            for i in range(len(frc)):
                if frc[i]<0.143:
                    messages.append("Resolution FSC(0.143)= %f"%resolution[i])
                    break
        return messages

    def visualize(self):
        from protlib_gui_figure import XmippPlotter
        fnFSC=self.workingDirPath("fsc.xmd")
        if os.path.exists(fnFSC):
            os.system('xmipp_metadata_plot -i %s -x resolutionFreqFourier -y resolutionFRC --title "Fourier Shell Correlation" --xtitle "1/Angstroms"'%fnFSC) 

def fsc(log,ReferenceVol,InputVol,WorkingDir):
    fnAcquisition=os.path.join(WorkingDir,'acquisition_info.xmd')
    if os.path.exists(fnAcquisition):
        md=MetaData(fnAcquisition)
        Ts=md.getValue(MDL_SAMPLINGRATE,md.firstObject())
    else:
        Ts=1.0

    fnOut=os.path.join(WorkingDir,'fsc.xmd')
    args="--ref %s -i %s -o %s --sampling_rate %f"%(ReferenceVol,InputVol,fnOut,float(Ts))
    runJob(log,"xmipp_resolution_fsc",args)
