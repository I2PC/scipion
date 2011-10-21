#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Author: Carlos Oscar Sorzano, October 2011

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
        self.importDir=getWorkingDirFromRunName(self.ImportRun)
        self.CtffindExec =  which('ctffind3.exe')

    def defineSteps(self):
        # Read Microscope parameters
        MD=xmipp.MetaData(os.path.join(self.importDir,"microscope.xmd"))
        id=MD.firstObject()
        Voltage=MD.getValue(xmipp.MDL_CTF_VOLTAGE,id)
        SphericalAberration=MD.getValue(xmipp.MDL_CTF_CS,id)
        AngPix=MD.getValue(xmipp.MDL_CTF_SAMPLING_RATE,id)
        Magnification=MD.getValue(xmipp.MDL_MAGNIFICATION,id)

        # Create verifyFiles for the MPI and output directories
        selfile=os.path.join(self.importDir,"micrographs.sel")
        MD=xmipp.MetaData(selfile)
        mDict={}        
        verifyFiles=[]
        for id in MD:
            inputFile=MD.getValue(xmipp.MDL_IMAGE,id)
            micrographName = os.path.split(inputFile)[1]
            (shortname, extension)=os.path.splitext(micrographName)
            micrographDir=self.workingDirPath(shortname)                    
            mDict[inputFile]=(micrographName,shortname,micrographDir)
            
            id=self.Db.insertStep('createDir',verifyfiles=[micrographDir],path=micrographDir)

            verifyFiles.append(os.path.join(micrographDir,"xmipp_ctf.ctfparam"))
            if self.DoCtffind:
                verifyFiles.append(os.path.join(micrographDir,"ctffind.ctfparam"))
        idMPI=self.insertRunMpiGapsStep()
        
        # Now the estimation actions
        CtfFindActions=[]
        for id in MD:
            inputFile=MD.getValue(xmipp.MDL_IMAGE,id)
            (micrographName,shortname,micrographDir)=mDict[inputFile]
            
            # Downsample if necessary
            parent_id=XmippProjectDb.FIRST_STEP
            if self.Down!=1:
                finalname=os.path.join(self.TmpDir,shortname+"_tmp.mrc")
                parent_id=self.insertRunJobGapStep("xmipp_transform_downsample",
                                                   "-i %s -o %s --step %f --method fourier" % (inputFile,finalname,self.Down),
                                                   [finalname])
            else:
                finalname=inputFile
            
            # CTF estimation with Xmipp
            self.insertRunJobGapStep('xmipp_ctf_estimate_from_micrograph',
                                     "--micrograph "+finalname+\
                                     " --oroot "+os.path.join(micrographDir,"xmipp_ctf")+\
                                     " --kV "+str(Voltage)+\
                                     " --Cs "+str(SphericalAberration)+\
                                     " --sampling_rate "+str(AngPix)+\
                                     " --ctfmodelSize 256"+\
                                     " --Q0 "+str(self.AmplitudeContrast)+\
                                     " --min_freq "+str(self.LowResolCutoff)+\
                                     " --max_freq "+str(self.HighResolCutoff)+\
                                     " --pieceDim "+str(self.WinSize)+\
                                     " --defocus_range "+str((self.MaxFocus-self.MinFocus)*10000/2)+\
                                     " --defocusU "+str(-(self.MaxFocus+self.MinFocus)*10000/2),
                                     verifyfiles=[os.path.join(micrographDir,"xmipp_ctf.ctfparam")],parent_step_id=parent_id)

            # CTF estimation with Ctffind
            if self.DoCtffind:
                CtfFindActions.append([dict(verifyfiles=[os.path.join(micrographDir,"ctffind.ctfparam")],
                                     execution_mode=SqliteDb.EXEC_GAP,
                                     parent_step_id=parent_id,
                                     CtffindExec=self.CtffindExec,micrograph=finalname,micrographDir=micrographDir,
                                     tmpDir=self.TmpDir,
                                     Voltage=Voltage,SphericalAberration=SphericalAberration,
                                     AngPix=AngPix,Magnification=Magnification,AmplitudeContrast=self.AmplitudeContrast,
                                     LowResolCutoff=self.LowResolCutoff,HighResolCutoff=self.HighResolCutoff,
                                     MinFocus=self.MinFocus,MaxFocus=self.MaxFocus,
                                     StepFocus=self.StepFocus,WinSize=self.WinSize)])
        for action in CtfFindActions:
            self.insertStep('estimateCtfCtffind',action)
        
        # Gather results after external actions
        self.insertStep('gatherResults',verifyfiles=[self.workingDirPath("micrographs.sel")],
                           parent_step_id=idMPI,
                           TmpDir=self.TmpDir,
                           WorkingDir=self.WorkingDir,
                           Selfile=selfile,
                           DoCtffind=self.DoCtffind)
               
    def validate(self):
        errors = []

        # Check that there are any micrograph to process
        fnSel=os.path.join(self.importDir,"micrographs.sel")
        if not os.path.exists(fnSel):
            errors.append("Cannot find micrographs.sel in "+self.importDir)
        else:
            MD=xmipp.MetaData(os.path.join(self.importDir,"micrographs.sel"))
            if MD.size()==0:
                errors.append("No micrographs to process")
        fnMic=os.path.join(self.importDir,"microscope.xmd")
        if not os.path.exists(fnMic):
            errors.append("Cannot find microscope.xmd in "+self.importDir)
        
        # Check that Q0 is negative
        if self.AmplitudeContrast>0:
            errors.append("Q0 should be negative ")
        
        # Check CTFFIND is available
        if self.DoCtffind and self.CtffindExec=="":
            errors.append("Cannot locate ctffind3.exe")
    
        return errors

    def summary(self):
        message=[]
        fnSel=os.path.join(self.importDir,"micrographs.sel")
        MD=xmipp.MetaData(fnSel)
        message.append("CTF screening of %d micrographs from %s" % (MD.size(), self.importDir))
        if self.DoCtffind:
            message.append("CTF validated with CTFFIND3")
        return message
    
    def visualize(self):
        summaryFile=self.workingDirPath("micrographs.sel")
        if os.path.exists(summaryFile):
            os.system("xmipp_visualize_preprocessing_micrographj -i "+summaryFile+" --memory 2048m &")
        else:
            summaryFile=os.path.join(self.TmpDir,"micrographs.sel")
            selfile=os.path.join(self.importDir,"micrographs.sel")
            buildSummaryMetadata(self.WorkingDir,self.DoCtffind,selfile,summaryFile)
            if os.path.exists(summaryFile):
                os.system("xmipp_visualize_preprocessing_micrographj -i "+summaryFile+" --memory 2048m &")
            else:
                import tkMessageBox
                tkMessageBox.showerror("Error", "There is no result yet")        
    
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

def buildSummaryMetadata(WorkingDir,DoCtffind,Selfile,summaryFile):
    MD=xmipp.MetaData()
    MDmicrographs=xmipp.MetaData(Selfile)
    for id in MDmicrographs:
        inputFile=MDmicrographs.getValue(xmipp.MDL_IMAGE,id)
        micrographName = os.path.split(inputFile)[1]
        (shortname, extension)=os.path.splitext(micrographName)
        micrographDir=os.path.join(WorkingDir,shortname)                    
        
        objId=MD.addObject()
        MD.setValue(xmipp.MDL_IMAGE,inputFile,objId)
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

def gatherResults(log,TmpDir,WorkingDir,Selfile,DoCtffind):
    summaryFile=os.path.join(WorkingDir,"micrographs.sel")
    buildSummaryMetadata(WorkingDir,DoCtffind,Selfile,summaryFile)
    runJob(log,"xmipp_ctf_sort_psds","-i "+summaryFile)
