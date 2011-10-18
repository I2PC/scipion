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

class ProtImportMicrographs(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.import_micrographs.name, scriptname, project)
        self.Import="from protocol_import_micrographs import *"

    def defineSteps(self):
        # Create microscope
        self.insertCreateMicroscope()

        # Decide name after preprocessing
        fileDict={}
        self.actualDoPreprocess=self.DoPreprocess and (self.Stddev != -1 or self.Crop != -1 or self.Down != 1)
        for filename in glob.glob(os.path.join(self.DirMicrographs, self.ExtMicrographs)):
            (filepath, micrographName)=os.path.split(filename)
            if self.actualDoPreprocess:
                (finalname, extension)=os.path.splitext(micrographName)
                finalname += ".mrc"
            else:
                finalname = micrographName
            fileDict[filename]=self.workingDirPath(finalname)

        # Preprocess
        idMPI=self.insertRunMpiGapsStep(fileDict.values())
        for filename in glob.glob(os.path.join(self.DirMicrographs, self.ExtMicrographs)):
            self.insertPreprocessStep(filename,fileDict[filename])
        
        # Gather results
        self.insertStep('gatherResults',verifyfiles=[self.workingDirPath("micrographs.sel")],
                        parent_step_id=idMPI,WorkingDir=self.WorkingDir)

    def validate(self):
        errors = []

        # Check that there are any micrograph to process
        listStr = self.DirMicrographs + '/' + self.ExtMicrographs
        listOfMicrographs = glob.glob(listStr)
        if len(listOfMicrographs) == 0:
            errors.append("There are no micrographs to process in " + listStr)

        return errors

    def summary(self):
        message=[]
        listOfMicrographs=glob.glob(self.DirMicrographs + '/' + self.ExtMicrographs)
        message.append("Import of %d micrographs from %s" % (len(listOfMicrographs), self.DirMicrographs))
        return message
    
    def visualize(self):
        summaryFile=self.workingDirPath("micrographs.sel")
        if os.path.exists(summaryFile):
            os.system("xmipp_visualize_preprocessing_micrographj -i "+summaryFile+" --memory 2048m &")

    def insertCreateMicroscope(self):    
        if self.SamplingRateMode=="From image":
            AngPix=self.SamplingRate
        else:
            AngPix=(10000. * self.ScannedPixelSize * self.Down) / self.Magnification
        fnOut=self.workingDirPath("microscope.xmd")
        self.insertStep("createMicroscope", verifyfiles = [fnOut], fnOut=fnOut, Voltage=self.Voltage,
                        SphericalAberration=self.SphericalAberration,SamplingRate=AngPix)
            
    def insertPreprocessStep(self,micrograph,finalname):
        previousId=XmippProjectDb.FIRST_STEP
        if not self.actualDoPreprocess:
            if not os.path.exists(finalname):
                previousId=self.insertStep("createLink",verifyfiles=[finalname],
                                           parent_step_id=previousId, execution_mode=SqliteDb.EXEC_GAP,
                                           source=micrograph,dest=finalname)
                if finalname.endswith(".raw"):
                    previousId=self.insertStep("createLink",verifyfiles=[finalname+".inf"],
                                               parent_step_id=previousId, execution_mode=SqliteDb.EXEC_GAP,
                                               source=micrograph+".inf",dest=finalname+".inf")
        else:    
            # Crop
            iname=micrograph
            if self.Crop != -1:
                previousId=self.insertStep("runJob",programname="xmipp_transform_window",
                                           params=" -i %s -o %s --crop %d -v 0" %(iname,finalname,self.Crop),
                                           verifyfiles = [finalname], parent_step_id=previousId, execution_mode=SqliteDb.EXEC_GAP)
                iname=finalname
            
            # Remove bad pixels
            if self.Stddev != -1:
                params = " -i %s --bad_pixels outliers %f -v 0" % (iname,self.Stddev)
                if not iname == finalname:
                    params += " -o " + finalname
                    iname=finalname
                previousId=self.insertStep("runJob",programname="xmipp_transform_filter",
                                           params=params, verifyfiles = [finalname], parent_step_id=previousId, execution_mode=SqliteDb.EXEC_GAP)
            
            # Downsample
            if self.Down != 1:
                tmpFile=finalname+"_tmp.mrc"
                previousId=self.insertStep("runJob",programname="xmipp_transform_downsample",
                                           params="-i %s -o %s --step %f --method fourier" % (iname,tmpFile,self.Down),
                                           verifyfiles = [tmpFile], parent_step_id=previousId, execution_mode=SqliteDb.EXEC_GAP)
                self.insertStep("renameFile",verifyfiles=[finalname],
                                parent_step_id=previousId, execution_mode=SqliteDb.EXEC_GAP,
                                source=tmpFile,dest=finalname),

def createMicroscope(log,fnOut,Voltage,SphericalAberration,SamplingRate):
    MD=xmipp.MetaData()
    MD.setColumnFormat(False)
    id=MD.addObject()
    MD.setValue(xmipp.MDL_CTF_VOLTAGE,float(Voltage),id)    
    MD.setValue(xmipp.MDL_CTF_CS,float(SphericalAberration),id)    
    MD.setValue(xmipp.MDL_CTF_SAMPLING_RATE,float(SamplingRate),id)
    MD.write(fnOut)    

def gatherResults(log,WorkingDir):
    summaryFile=os.path.join(WorkingDir,"micrographs.sel")
    MD=xmipp.MetaData()
    for filename in glob.glob(WorkingDir+"/*"):
        if filename.endswith(".xmd") or filename.endswith("tmp") or filename.endswith(".sel"):
            continue
        objId=MD.addObject()
        MD.setValue(xmipp.MDL_IMAGE,filename,objId)
    if MD.size()!=0:
        MD.sort(xmipp.MDL_IMAGE);
        MD.write(summaryFile)
