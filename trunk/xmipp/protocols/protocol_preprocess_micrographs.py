#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# General script for Xmipp-based pre-processing of micrographs: 
#  - downsampling
#  - power spectral density (PSD) and CTF estimation on the micrograph
#
# For each micrograph given, this script will perform 
# the requested operations below.
# For each micrograph a subdirectory will be created
#
# Author: Carlos Oscar Sorzano, July 2011


import glob, os, shutil, sys
from config_protocols import protDict
from protlib_base import *
from protlib_utils import which, runJob
from protlib_filesystem import *
import xmipp

#FIXME: IMPLEMENTATION SHOULD BE REVISED
class ProtPreprocessMicrographs(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.preprocess_micrographs.key, scriptname, project)
        self.Import="from protocol_preprocess_micrographs import launchPreprocessMicrographsBatch, gatherResults"

        # Check if ctffind executable is available
        self.CtffindExec =  which('ctffind3.exe')
        self.DoCtffind = self.CtffindExec != ''

    def defineActions(self):
        for filename in glob.glob(self.DirMicrographs + '/' + self.ExtMicrographs):
            # Get the shortname and extension
            (filepath, micrographName)=os.path.split(filename)
            (shortname, extension)=os.path.splitext(micrographName)
            micrographDir=os.path.join(self.WorkingDir,shortname)                    
            finalname='micrograph'
            if self.DoPreprocess and (not self.Stddev == -1 or not self.Crop == -1 or not self.Down == 1):
                finalname += ".mrc"
            else:
                finalname += extension
                
            if self.SamplingRateMode=="From image":
                AngPix=self.SamplingRate
            else:
                AngPix=(10000. * self.ScannedPixelSize * self.Down) / self.Magnification
            
            # Insert actions in the database
            self.Db.insertAction('createDir',path=micrographDir)
            id=self.Db.insertAction('preprocessMicrograph',[os.path.join(micrographDir,"micrograph."+extension)],None,
                                    False,False,micrograph=filename,micrographDir=micrographDir,DoPreprocess=self.DoPreprocess,
                                    Crop=self.Crop,Stddev=self.Stddev,Down=self.Down)
            if self.DoCtfEstimate:
                self.Db.insertAction('estimateCtfXmipp',[os.path.join(micrographDir,"xmipp_ctf.ctfparam")],id,False,False,
                                     micrograph=finalname,micrographDir=micrographDir,Voltage=self.Voltage,
                                     SphericalAberration=self.SphericalAberration,AngPix=AngPix,
                                     AmplitudeContrast=self.AmplitudeContrast,LowResolCutoff=self.LowResolCutoff,
                                     HighResolCutoff=self.HighResolCutoff,
                                     MinFocus=self.MinFocus,MaxFocus=self.MaxFocus,WinSize=self.WinSize)
                if self.DoCtffind:
                    self.Db.insertAction('estimateCtfCtffind',[os.path.join(micrographDir,"ctffind.ctfparam")],id,False,False,
                                         CtffindExec=self.CtffindExec,micrograph=finalname,micrographDir=micrographDir,
                                         Voltage=self.Voltage,SphericalAberration=self.SphericalAberration,
                                         AngPix=AngPix,Magnification=self.Magnification,AmplitudeContrast=self.AmplitudeContrast,
                                         LowResolCutoff=self.LowResolCutoff,
                                         HighResolCutoff=self.HighResolCutoff,MinFocus=self.MinFocus,MaxFocus=self.MaxFocus,
                                         StepFocus=self.StepFocus,WinSize=self.WinSize)

        # Launch all the external actions
        self.Db.insertAction('launchPreprocessMicrographsBatch',WorkingDir=self.WorkingDir,protocolHeader=self.scriptName)
        
        # Gather results after external actions
        self.Db.insertAction('gatherResults',WorkingDir=self.WorkingDir,DirMicrographs=self.DirMicrographs,
                             ExtMicrographs=self.ExtMicrographs, DoCtfEstimate=self.DoCtfEstimate,DoCtffind=self.DoCtffind)
               
    def validate(self):
        errors = []

        # Check that there are any micrograph to process
        listOfMicrographs=glob.glob(self.DirMicrographs + '/' + self.ExtMicrographs)
        if len(listOfMicrographs) == 0:
            errors.append("There are no micrographs to process in ") + DirMicrographs + '/' + ExtMicrographs
        
        # Check that Q0 is negative
        if self.AmplitudeContrast>0:
            errors.append("Q0 should be negative ")
    
        return errors

    def summary(self):
        message=[]
        listOfMicrographs=glob.glob(self.DirMicrographs + '/' + self.ExtMicrographs)
        message.append("Preprocessing of %d micrographs from %s" % (len(listOfMicrographs), self.DirMicrographs))
        if self.DoCtfEstimate:
            msg="CTF estimated with Xmipp"
            if self.DoCtffind:
                msg+=" and validated with CTFFIND3"
            message.append(msg)
        return message
    
    def visualize(self):
        summaryFile=os.path.join(self.WorkingDir,"micrographs.sel")
        if os.path.exists(summaryFile):
            os.system("xmipp_visualize_preprocessing_micrographj -i "+summaryFile+" --mem 2048m &")
        else:
            import tkMessageBox
            tkMessageBox.showerror("Error", "There is no result yet")        
    
def launchPreprocessMicrographsBatch(log,WorkingDir,protocolHeader):
    pass

def gatherResults(log,WorkingDir,DirMicrographs,ExtMicrographs,DoCtfEstimate,DoCtffind):
    MD=xmipp.MetaData()
    for filename in glob.glob(DirMicrographs + '/' + ExtMicrographs):
        (filepath, micrographName)=os.path.split(filename)
        (shortname, extension)=os.path.splitext(micrographName)
        micrographDir=WorkingDir+"/"+shortname

        objId=MD.addObject()
        MD.setValue(xmipp.MDL_IMAGE, filename,objId)
        if DoCtfEstimate:
            MD.setValue(xmipp.MDL_PSD,               os.path.join(micrographDir,"xmipp_ctf.psd"),objId)
            MD.setValue(xmipp.MDL_PSD_ENHANCED,      os.path.join(micrographDir,"xmipp_ctf_psd_enhanced.xmp"),objId)
            MD.setValue(xmipp.MDL_CTFMODEL,          os.path.join(micrographDir,"xmipp_ctf.ctfparam"),objId)
            MD.setValue(xmipp.MDL_ASSOCIATED_IMAGE1, os.path.join(micrographDir,"xmipp_ctf_ctfmodel_quadrant.xmp"),objId)
            MD.setValue(xmipp.MDL_ASSOCIATED_IMAGE2, os.path.join(micrographDir,"xmipp_ctf_ctfmodel_halfplane.xmp"),objId)
            if DoCtffind:
                MD.setValue(xmipp.MDL_CTFMODEL2, os.path.join(micrographDir,"ctffind.ctfparam"),objId)
                MD.setValue(xmipp.MDL_ASSOCIATED_IMAGE3, os.path.join(micrographDir,"ctffind_spectrum.mrc"),objId)
    MD.sort(xmipp.MDL_IMAGE);
    MD.write(os.path.join(WorkingDir,"micrographs.sel"))

    # CTF Quality control
    if DoCtfEstimate:
        runJob(log,"xmipp_ctf_sort_psds","-i "+WorkingDir + "/micrographs.sel")
