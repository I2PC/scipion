#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based pre-processing of volumes 

# Author: Carlos Oscar, August 2013
#
from protlib_base import *
import os
from protlib_utils import runJob, runShowJ, grepFirst
from protlib_filesystem import findAcquisitionInfo
import glob
from xmipp import MetaData, MDL_SAMPLINGRATE, MDL_RESOLUTION_FREQREAL, MDL_RESOLUTION_FRC
from protlib_gui_ext import showError
from os.path import dirname, basename

class ProtResolution3D(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.resolution3D.name, scriptname, project)
        self.Import = 'from protocol_resolution3D import *'

    def defineSteps(self):
        if self.DoFSC:
            self.insertStep('fsc',ReferenceVol=self.ReferenceVol,InputVol=self.InputVol,WorkingDir=self.WorkingDir)
        if self.DoStructureFactor:
            self.insertStep('structureFactor',Structure=self.workingDirPath('structureFactor.xmd'),InputVol=self.InputVol)
        if self.DoSSNR or self.DoVSSNR:
            self.insertStep('createDir',path=self.ExtraDir)
            projmatchProtocol=self.project.getProtocolFromFile(self.InputVol)
            reconstructionCmd=grepFirst("Logs/%s.log"%projmatchProtocol.getExtendedRunName(),"-o "+self.InputVol)
            i0=reconstructionCmd.find("Running command: ")
            iF=reconstructionCmd.find(" --xmipp_protocol_script")
            trueCmd=reconstructionCmd[(i0+24):iF]
            self.insertStep('createNoisyVolume',WorkingDir=self.WorkingDir,cmd=trueCmd)
        if self.DoSSNR:
            self.insertStep('calculateSSNR',WorkingDir=self.WorkingDir,cmd=trueCmd)
        if self.DoVSSNR:
            self.insertStep('calculateVSSNR',WorkingDir=self.WorkingDir,cmd=trueCmd)

    def validate(self):
        errors = []
        if self.DoSSNR or self.DoVSSNR:
            if not self.InputVol.startswith(protDict.projmatch.dir):
                errors.append("SSNR and VSSNR are meant only for volumes coming out of a projection matching protocol")
            fnVol=basename(self.InputVol)
            import re
            plainReconstruction=re.compile("reconstruction_Ref3D_...\.vol")
            if re.match(plainReconstruction,fnVol) is None:
                errors.append("SSNR and VSSNR are meant only for volumes after reconstruction (no filtering). Their filename is of the kind reconstruction_Rec3D_???.vol")
        return errors
        
    def summary(self):
        messages = []      
        messages.append("Input     volume: [%s]" % self.InputVol)
        if self.DoFSC:
            messages.append("Reference volume: [%s]" % self.ReferenceVol)
            fnFSC=self.workingDirPath("fsc.xmd")
            if os.path.exists(fnFSC):
                md = MetaData(fnFSC)
                resolution = [md.getValue(MDL_RESOLUTION_FREQREAL, id) for id in md]
                frc = [md.getValue(MDL_RESOLUTION_FRC, id) for id in md]
                for i in range(len(frc)):
                    if frc[i]<0.5:
                        messages.append("Resolution FSC(=0.5) = %f"%resolution[i])
                        break
                for i in range(len(frc)):
                    if frc[i]<0.143:
                        messages.append("Resolution FSC(=0.143)= %f"%resolution[i])
                        break
        if self.DoSSNR:
            fnSSNR=self.workingDirPath("ssnr.xmd")
            if os.path.exists(fnSSNR):
                md = MetaData(fnSSNR)
                resolution = [md.getValue(MDL_RESOLUTION_FREQREAL, id) for id in md]
                ssnr = [md.getValue(MDL_RESOLUTION_SSNR, id) for id in md]
                for i in range(len(ssnr)):
                    if ssnr[i]<1:
                        messages.append("Resolution SSNR(=1) = %f"%resolution[i])
                        break
        return messages

    def visualize(self):
        fnFSC=self.workingDirPath("fsc.xmd")
        if os.path.exists(fnFSC):
            if self.DisplayFSC and self.DoFSC:
                os.system('xmipp_metadata_plot -i %s -x resolutionFreqFourier -y resolutionFRC --title "Fourier Shell Correlation" --xtitle "1/Angstroms" &'%fnFSC) 
            if self.DisplayDPR and self.DoFSC:
                os.system('xmipp_metadata_plot -i %s -x resolutionFreqFourier -y resolutionDPR --title "Differential Phase Residual" --xtitle "1/Angstroms" &'%fnFSC) 
        fnStructure=self.workingDirPath("structureFactor.xmd")
        if os.path.exists(fnStructure):
            if self.DisplayStructureFactor and self.DoStructureFactor:
                os.system('xmipp_metadata_plot -i %s -x resolutionFreqFourier -y resolutionLogStructure --title "Structure factor" --xtitle "Frequency (1/A)" --ytitle "Log(StructureFactor)" &'%fnStructure)
            if self.DisplayGuinier:
                if self.UseMatlab:
                    os.system("matlab -r \"xmipp_show_structure_factor(\'"+self.WorkingDir+"\')\"")
                else:
                    os.system('xmipp_metadata_plot -i %s -x resolutionFreqFourier2 -y resolutionLogStructure --title "Guinier plot" --xtitle "Frequency (1/A^2)" --ytitle "Log(StructureFactor)" &'%fnStructure)
        fnSSNR=self.workingDirPath("ssnr.xmd")
        if os.path.exists(fnSSNR):
            if self.DisplaySSNR and self.DoSSNR:
                os.system('xmipp_metadata_plot -i %s -x resolutionFreqFourier -y resolutionSSNR --title "Spectral SNR" --xtitle "Frequency (1/A)" --ytitle "SNR" &'%fnSSNR)
        fnVSSNR=self.workingDirPath("vssnr.vol")
        if os.path.exists(fnVSSNR):
            if self.DisplayVSSNR and self.DoVSSNR:
                runShowJ(fnVSSNR)

def getTs(InputVol):
    fnAcquisition=findAcquisitionInfo(InputVol)
    if os.path.exists(fnAcquisition):
        md=MetaData(fnAcquisition)
        Ts=md.getValue(MDL_SAMPLINGRATE,md.firstObject())
    else:
        print("Cannot find acquisition_info.xmd. Using sampling rate of 1 Angstrom/pixel")
        Ts=1.0
    return Ts

def fsc(log,ReferenceVol,InputVol,WorkingDir):
    Ts=getTs(InputVol)
    fnOut=os.path.join(WorkingDir,'fsc.xmd')
    args="--ref %s -i %s -o %s --sampling_rate %f --do_dpr"%(ReferenceVol,InputVol,fnOut,float(Ts))
    runJob(log,"xmipp_resolution_fsc",args)

def structureFactor(log,Structure,InputVol):
    Ts=getTs(InputVol)
    runJob(log,"xmipp_volume_structure_factor","-i %s -o %s --sampling %f"%(InputVol,Structure,float(Ts)))

def createNoisyVolume(log,WorkingDir,cmd):
    tokens=cmd.split(' ')
    fnMetadata=tokens[tokens.index('-i')+1]
    fnOut=os.path.join(WorkingDir,'extra/noisyImages')
    runJob(log,'xmipp_image_operate','-i %s -o %s.stk --reset --save_metadata_stack'%(fnMetadata,fnOut))
    runJob(log,'xmipp_transform_add_noise','-i %s.stk --type gaussian 1'%fnOut)
    tokens[tokens.index('-i')+1]=fnOut+".xmd"
    tokens[tokens.index('-o')+1]=os.path.join(WorkingDir,'extra/noisyVolume.vol')
    newArgs=' '.join(tokens[1:])
    runJob(log,tokens[0],newArgs)

def getSSNRParams(cmd,WorkingDir):
    tokens=cmd.split(' ')
    fnSignalVolume=tokens[tokens.index('-o')+1]
    fnSignalSel=tokens[tokens.index('-i')+1]
    fnNoiseVolume=os.path.join(WorkingDir,'extra/noisyVolume.vol')
    fnNoiseSel=os.path.join(WorkingDir,'extra/noisyImages.xmd')
    Ts=getTs(fnSignalVolume)
    return [fnSignalVolume,fnSignalSel,fnNoiseVolume,fnNoiseSel,Ts]

def calculateSSNR(log,WorkingDir,cmd):
    fnSignalVolume,fnSignalSel,fnNoiseVolume,fnNoiseSel,Ts=getSSNRParams(cmd,WorkingDir)
    args="--signal "+fnSignalVolume+" --sel_signal "+fnSignalSel+" --noise "+fnNoiseVolume+" --sel_noise "+fnNoiseSel+\
         " -o "+os.path.join(WorkingDir,"ssnr.xmd")+" --sampling_rate "+str(Ts)
    runJob(log,'xmipp_resolution_ssnr',args)

def calculateVSSNR(log,WorkingDir,cmd):
    fnSignalVolume,fnSignalSel,fnNoiseVolume,fnNoiseSel,Ts=getSSNRParams(cmd,WorkingDir)
    args="--signal "+fnSignalVolume+" --sel_signal "+fnSignalSel+" --noise "+fnNoiseVolume+" --sel_noise "+fnNoiseSel+\
         " -o "+os.path.join(WorkingDir,"ssnr.xmd")+" --sampling_rate "+str(Ts)+" --gen_VSSNR --VSSNR "+os.path.join(WorkingDir,"vssnr.vol")
    runJob(log,'xmipp_resolution_ssnr',args)
