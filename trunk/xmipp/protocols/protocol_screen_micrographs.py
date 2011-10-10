#!/usr/bin/env xmipp_python
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

import glob
from config_protocols import protDict
from protlib_base import *
from protlib_utils import which, runJob
from protlib_filesystem import *
import xmipp

class ProtScreenMicrographs(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.screen_micrographs.name, scriptname, project)
        self.Import="from protocol_screen_micrographs import *"
        self.CtffindExec =  which('ctffind3.exe')
        self.ContinueAtStep=1

    def defineSteps(self):
        CtfFindActions=[]
        idMPI=self.Db.insertStep('runStepGapsMpi',passDb=True, script=self.scriptName, NumberOfMpi=self.NumberOfMpi)
        verifyFiles=[]        
        for filename in glob.glob(os.path.join(self.DirMicrographs, self.ExtMicrographs)):
            # Get the shortname and extension
            micrographName = os.path.split(filename)[1]
            (shortname, extension)=os.path.splitext(micrographName)
            micrographDir=self.workingDirPath(shortname)                    
            finalname='micrograph'
            if self.DoPreprocess and (not self.Stddev == -1 and not self.Crop == -1 and not self.Down == 1):
                finalname += ".mrc"
            else:
                finalname += extension
                
            if self.SamplingRateMode=="From image":
                AngPix=self.SamplingRate
            else:
                AngPix=(10000. * self.ScannedPixelSize * self.Down) / self.Magnification
            
            # Insert actions in the database
            id=self.Db.insertStep('createDir',path=micrographDir,parent_step_id=XmippProjectDb.FIRST_STEP,execution_mode=SqliteDb.EXEC_GAP)
            fnOut=os.path.join(micrographDir,"micrograph"+extension)
            verifyFiles.append(fnOut)
            id=self.Db.insertStep('preprocessMicrograph',#COSS verifyfiles=[fnOut],
                                    parent_step_id=id, execution_mode=SqliteDb.EXEC_GAP,
                                    micrograph=filename,micrographDir=micrographDir,DoPreprocess=self.DoPreprocess,
                                    Crop=self.Crop,Stddev=self.Stddev,Down=self.Down)
            if self.DoCtfEstimate:
                fnOut=os.path.join(micrographDir,"xmipp_ctf.ctfparam")
                verifyFiles.append(fnOut)
                self.Db.insertStep('estimateCtfXmipp',verifyfiles=[fnOut],
                                     parent_step_id=id,execution_mode=SqliteDb.EXEC_GAP,
                                     micrograph=finalname,micrographDir=micrographDir,Voltage=self.Voltage,
                                     SphericalAberration=self.SphericalAberration,AngPix=AngPix,
                                     AmplitudeContrast=self.AmplitudeContrast,LowResolCutoff=self.LowResolCutoff,
                                     HighResolCutoff=self.HighResolCutoff,
                                     MinFocus=self.MinFocus,MaxFocus=self.MaxFocus,WinSize=self.WinSize)
                if self.DoCtffind:
                    fnOut=os.path.join(micrographDir,"ctffind.ctfparam")
                    verifyFiles.append(fnOut)
                    CtfFindActions.append([dict(verifyfiles=[fnOut],
                                         execution_mode=SqliteDb.EXEC_GAP,
                                         CtffindExec=self.CtffindExec,micrograph=finalname,micrographDir=micrographDir,
                                         tmpDir=self.TmpDir,
                                         Voltage=self.Voltage,SphericalAberration=self.SphericalAberration,
                                         AngPix=AngPix,Magnification=self.Magnification,AmplitudeContrast=self.AmplitudeContrast,
                                         LowResolCutoff=self.LowResolCutoff,
                                         HighResolCutoff=self.HighResolCutoff,MinFocus=self.MinFocus,MaxFocus=self.MaxFocus,
                                         StepFocus=self.StepFocus,WinSize=self.WinSize)])
        for action in CtfFindActions:
            action["parent_step_id"]=id # This makes all ctffinds to go after the last preprocessing
            self.Db.insertStep('estimateCtfCtffind',action)
        self.Db.updateVerifyFiles(idMPI,verifyFiles)
        
        # Gather results after external actions
        self.Db.insertStep('gatherResults',verifyfiles=[self.workingDirPath("micrographs.sel")],
                           parent_step_id=idMPI,
                           WorkingDir=self.WorkingDir,DirMicrographs=self.DirMicrographs,
                           ExtMicrographs=self.ExtMicrographs, DoCtfEstimate=self.DoCtfEstimate,DoCtffind=self.DoCtffind)
               
    def validate(self):
        errors = []

        # Check that there are any micrograph to process
        listStr = self.DirMicrographs + '/' + self.ExtMicrographs
        listOfMicrographs = glob.glob(listStr)
        if len(listOfMicrographs) == 0:
            errors.append("There are no micrographs to process in " + listStr)
        
        # Check that Q0 is negative
        if self.AmplitudeContrast>0:
            errors.append("Q0 should be negative ")
        
        # Check CTFFIND is available
        if self.DoCtffind and self.CtffindExec=="":
            errors.append("Cannot locate ctffind3.exe")
    
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
        summaryFile=self.workingDirPath("micrographs.sel")
        if os.path.exists(summaryFile):
            os.system("xmipp_visualize_preprocessing_micrographj -i "+summaryFile+" --memory 2048m &")
        else:
            summaryFile=os.path.join(self.TmpDir,"micrographs.sel")
            buildSummaryMetadata(self.WorkingDir,self.DoCtfEstimate,self.DoCtffind,summaryFile)
            if os.path.exists(summaryFile):
                os.system("xmipp_visualize_preprocessing_micrographj -i "+summaryFile+" --memory 2048m &")
            else:
                import tkMessageBox
                tkMessageBox.showerror("Error", "There is no result yet")        
    
def preprocessMicrograph(log,micrograph,micrographDir,DoPreprocess,Crop,Stddev,Down):
    # Decide name after preprocessing
    finalname=os.path.join(micrographDir,'micrograph')
    (filepath, micrographName)=os.path.split(os.path.relpath(micrograph))
    (shortname2, extension)=os.path.splitext(micrographName)        
    if DoPreprocess and (not Stddev == -1 or not Crop == -1 or not Down == 1):
        finalname += ".mrc"
    else:
        finalname += extension
        if not os.path.exists(finalname):
            relpath = os.path.relpath(micrograph, micrographDir)
            runJob(log,"ln",'-s %(relpath)s %(finalname)s' % locals())
            if micrograph.endswith(".raw"):
                runJob(log,"ln",'-s %(relpath)s.inf %(finalname)s.inf' % locals())
            return
    if not DoPreprocess:
        return
    
    # Crop
    iname=micrograph
    if not Crop == -1:
        runJob(log,"xmipp_transform_window"," -i %(iname)s -o %(finalname)s --crop %(Crop)d -v 0" % locals())
        iname=finalname
    
    # Remove bad pixels
    if not Stddev == -1:
        params = " -i %(iname)s --bad_pixels outliers %(Stddev)f -v 0" % locals()
        if not iname == finalname:
            params += " -o " + finalname
            iname=finalname
        runJob(log,"xmipp_transform_filter",params)
    
    # Downsample
    if not Down == 1:
        tmpFile=os.path.join(micrographDir,"tmp.mrc")
        runJob(log,"xmipp_transform_downsample","-i %(iname)s -o %(tmpFile)s --step %(Down)f --method fourier" % locals())
        runJob(log,"mv", "-f %(tmpFile)s %(finalname)s" % locals())

def estimateCtfXmipp(log,micrograph,micrographDir,Voltage,SphericalAberration,AngPix,
                     AmplitudeContrast,LowResolCutoff,HighResolCutoff,MinFocus,MaxFocus,WinSize):
    params="--micrograph "+os.path.join(micrographDir,micrograph)+" --oroot "+os.path.join(micrographDir,"xmipp_ctf")+\
           " --kV "+str(Voltage)+\
           " --Cs "+str(SphericalAberration)+\
           " --sampling_rate "+str(AngPix)+\
           " --ctfmodelSize 256"+\
           " --Q0 "+str(AmplitudeContrast)+\
           " --min_freq "+str(LowResolCutoff)+\
           " --max_freq "+str(HighResolCutoff)+\
           " --pieceDim "+str(WinSize)+\
           " --defocus_range "+str((MaxFocus-MinFocus)*10000/2)+\
           " --defocusU "+str(-(MaxFocus+MinFocus)*10000/2)
    runJob(log,"xmipp_ctf_estimate_from_micrograph",params)

def estimateCtfCtffind(log,CtffindExec,micrograph,micrographDir,tmpDir,Voltage,SphericalAberration,AngPix,Magnification,
                       AmplitudeContrast,LowResolCutoff,HighResolCutoff,MinFocus,MaxFocus,StepFocus,WinSize):
    # Convert image to MRC
    if not micrograph.endswith('.mrc'):
        import uuid
        deleteTempMicrograph=True
        mrcMicrograph= os.path.join(tmpDir,os.path.splitext(micrograph)[0]+"_"+str(uuid.uuid4())+'.mrc')
        runJob(log,'xmipp_image_convert','-i ' + micrograph + ' -o ' + mrcMicrograph + ' -v 0')
    else:
        deleteTempMicrograph=False
        mrcMicrograph = micrograph;

    # Prepare parameters for CTFTILT
    params = '  << eof > ' + micrographDir + '/ctffind.log\n'
    params += os.path.join(micrographDir,mrcMicrograph) + "\n"
    params += micrographDir + '/ctffind_spectrum.mrc\n'
    params += str(SphericalAberration) + ',' + \
              str(Voltage) + ',' + \
              str(AmplitudeContrast) + ',' + \
              str(Magnification) + ',' + \
              str(Magnification*AngPix*1e-4) + "\n"
    params += str(WinSize) + ',' + \
              str(AngPix / LowResolCutoff) + ',' + \
              str(AngPix / HighResolCutoff) + ',' + \
              str(MinFocus*10000) + ',' + \
              str(MaxFocus*10000) + ',' + \
              str(StepFocus*10000) + "\n"
    runJob(log, "export NATIVEMTZ=kk ; "+CtffindExec,params)

    fnOut=micrographDir + '/ctffind.ctfparam'

    # Remove temporary files
    if deleteTempMicrograph:
        deleteFile(log, mrcMicrograph)

    # Pick values from ctffind
    ctffindLog=os.path.join(micrographDir,'ctffind.log')
    if not os.path.exists(ctffindLog):
        raise xmipp.XmippError("Cannot find "+ctffindLog)
    
    # Effectively pickup results
    fh=open(ctffindLog, 'r')
    lines=fh.readlines()
    fh.close()
    DF1=0.
    DF2=0.
    Angle=0.
    found=False
    for i in range(len(lines)):
        if not (lines[i].find('Final Values') == -1):
            words=lines[i].split()
            DF1=float(words[0])
            DF2=float(words[1])
            Angle=float(words[2])
            found=True
            break
    if not found:
        raise xmipp.XmippError("Cannot find defocus values in "+ctffindLog)
    
    # Generate Xmipp .ctfparam file:
    MD=xmipp.MetaData()
    MD.setColumnFormat(False)
    objId=MD.addObject()
    MD.setValue(xmipp.MDL_CTF_SAMPLING_RATE, float(AngPix), objId)
    MD.setValue(xmipp.MDL_CTF_VOLTAGE,       float(Voltage), objId)
    MD.setValue(xmipp.MDL_CTF_DEFOCUSU,      float(-DF2), objId)
    MD.setValue(xmipp.MDL_CTF_DEFOCUSV,      float(-DF1), objId)
    MD.setValue(xmipp.MDL_CTF_DEFOCUS_ANGLE, float(Angle), objId)
    MD.setValue(xmipp.MDL_CTF_CS,            float(SphericalAberration), objId)
    MD.setValue(xmipp.MDL_CTF_Q0,            float(-AmplitudeContrast), objId)
    MD.setValue(xmipp.MDL_CTF_K,             1.0, objId)
    MD.write(fnOut)

def buildSummaryMetadata(WorkingDir,DoCtfEstimate,DoCtffind,summaryFile):
    MD=xmipp.MetaData()
    for filename in glob.glob(WorkingDir + '/*/micrograph.*'):
        micrographDir=os.path.dirname(filename)
        objId=MD.addObject()
        MD.setValue(xmipp.MDL_IMAGE,filename,objId)
        if DoCtfEstimate:
            ctfparam=os.path.join(micrographDir,"xmipp_ctf.ctfparam")
            if os.path.exists(ctfparam):
                MD.setValue(xmipp.MDL_PSD,               os.path.join(micrographDir,"xmipp_ctf.psd"),objId)
                MD.setValue(xmipp.MDL_PSD_ENHANCED,      os.path.join(micrographDir,"xmipp_ctf_enhanced_psd.xmp"),objId)
                MD.setValue(xmipp.MDL_CTFMODEL,          ctfparam,objId)
                MD.setValue(xmipp.MDL_ASSOCIATED_IMAGE1, os.path.join(micrographDir,"xmipp_ctf_ctfmodel_quadrant.xmp"),objId)
                MD.setValue(xmipp.MDL_ASSOCIATED_IMAGE2, os.path.join(micrographDir,"xmipp_ctf_ctfmodel_halfplane.xmp"),objId)
            else:
                MD.setValue(xmipp.MDL_PSD,               "NA",objId)
                MD.setValue(xmipp.MDL_PSD_ENHANCED,      "NA",objId)
                MD.setValue(xmipp.MDL_CTFMODEL,          "NA",objId)
                MD.setValue(xmipp.MDL_ASSOCIATED_IMAGE1, "NA",objId)
                MD.setValue(xmipp.MDL_ASSOCIATED_IMAGE2, "NA",objId)
            if DoCtffind:
                ctffindCTF=os.path.join(micrographDir,"ctffind.ctfparam")
                if os.path.exists(ctffindCTF):
                    MD.setValue(xmipp.MDL_CTFMODEL2,         ctffindCTF,objId)
                    MD.setValue(xmipp.MDL_ASSOCIATED_IMAGE3, os.path.join(micrographDir,"ctffind_spectrum.mrc"),objId)
                else:
                    MD.setValue(xmipp.MDL_CTFMODEL2,         "NA",objId)
                    MD.setValue(xmipp.MDL_ASSOCIATED_IMAGE3, "NA",objId)
    if MD.size()!=0:
        MD.sort(xmipp.MDL_IMAGE);
        MD.write(summaryFile)

def gatherResults(log,WorkingDir,DirMicrographs,ExtMicrographs,DoCtfEstimate,DoCtffind):
    summaryFile=os.path.join(WorkingDir,"micrographs.sel")
    buildSummaryMetadata(WorkingDir,DoCtfEstimate,DoCtffind,summaryFile)
    if DoCtfEstimate:
        runJob(log,"xmipp_ctf_sort_psds","-i "+summaryFile)
