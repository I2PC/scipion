#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# General script for Xmipp-based pre-processing of micrographs
# Author: Carlos Oscar Sorzano, July 2011

from glob import glob
import os
from os.path import join
from config_protocols import protDict
from protlib_base import XmippProtocol
#from protlib_filesystem import *
import xmipp

class ProtImportMicrographs(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.import_micrographs.name, scriptname, project)
        self.Import="from protocol_import_micrographs import *"

    def createFilenameTemplates(self):
        return {
                'micrographs': self.workingDirPath('micrographs.sel'),
                'micrographsPattern': join(self.DirMicrographs, self.ExtMicrographs)
                }
        
    def defineSteps(self):
        # Create microscope
        self.insertCreateMicroscope()

        # Decide name after preprocessing
        fileDict={}
        self.actualDoPreprocess = self.DoPreprocess and (self.Stddev != -1 or self.Crop != -1 or self.Down != 1)
        for filename in glob(join(self.DirMicrographs, self.ExtMicrographs)):
            (filepath, micrographName) = os.path.split(filename)
            if self.actualDoPreprocess:
                (finalname, extension) = os.path.splitext(micrographName)
                finalname += ".mrc"
            else:
                finalname = micrographName
            fileDict[filename]=self.workingDirPath(finalname)

        # Preprocess
        #idMPI=self.insertRunMpiGapsStep(fileDict.values())
        for filename in glob(join(self.DirMicrographs, self.ExtMicrographs)):
            self.insertPreprocessStep(filename,fileDict[filename])
        
        # Gather results
        summaryFile = self.getFilename('micrographs')
        self.insertStep('gatherResults',verifyfiles=[summaryFile],
                        WorkingDir=self.WorkingDir, summaryFile=summaryFile)

    def validate(self):
        errors = []
        # Check that there are any micrograph to process
        if len(self.getMicrographs()) == 0:
            errors.append("There are no micrographs to process in " + self.getFilename('micrographsPattern'))
        return errors

    def summary(self):
        message = []
        micrographs = self.getMicrographs()
        message.append("Import of <%d> micrographs from <%s>" % (len(micrographs), self.DirMicrographs))
        return message
    
    def getMicrographs(self):
        ''' Return a list with micrographs in WorkingDir'''
        return glob(self.getFilename('micrographsPattern'))
    
    def visualize(self):
        summaryFile = self.getFilename('micrographs')
        if os.path.exists(summaryFile):
            os.system("xmipp_visualize_preprocessing_micrographj -i %s --memory 2048m &" % summaryFile)

    def insertCreateMicroscope(self):    
        if self.SamplingRateMode == "From image":
            AngPix = self.SamplingRate
        else:
            AngPix = (10000. * self.ScannedPixelSize * self.Down) / self.Magnification
        fnOut = self.workingDirPath("microscope.xmd")
        self.insertStep("createMicroscope", verifyfiles = [fnOut], fnOut=fnOut, Voltage=self.Voltage,
                        SphericalAberration=self.SphericalAberration,SamplingRate=AngPix,
                        Magnification=self.Magnification)
            
    def insertPreprocessStep(self,micrograph,finalname):
        from protlib_sql.XmippProjectDb import FIRST_STEP
        previousId = FIRST_STEP
        if not self.actualDoPreprocess:
            if not os.path.exists(finalname):
                previousId = self.insertParallelStep('createLink', verifyfiles=[finalname], source=micrograph, dest=finalname)
                if finalname.endswith(".raw"):
                    previousId = self.insertParallelStep("createLink",verifyfiles=[finalname+".inf"], parent_step_id=previousId,
                                               source=micrograph+".inf", dest=finalname+".inf")
        else:    
            # Crop
            iname=micrograph
            if self.Crop != -1:
                previousId = self.insertParallelRunJobStep("xmipp_transform_window", " -i %s -o %s --crop %d -v 0" %(iname,finalname,self.Crop),
                                                           verifyfiles=[finalname])
                iname=finalname
            # Remove bad pixels
            if self.Stddev != -1:
                params = " -i %s --bad_pixels outliers %f -v 0" % (iname, self.Stddev)
                if not iname == finalname:
                    params += " -o " + finalname
                    iname = finalname
                previousId = self.insertParallelRunJobStep("xmipp_transform_filter", params, verifyfiles=[finalname], parent_step_id=previousId)
            
            # Downsample
            if self.Down != 1:
                tmpFile = finalname + "_tmp.mrc"
                previousId = self.insertParallelRunJobStep("xmipp_transform_downsample", "-i %s -o %s --step %f --method fourier" % (iname,tmpFile,self.Down),
                                           verifyfiles = [tmpFile], parent_step_id=previousId)
                self.insertParallelStep("renameFile",verifyfiles=[finalname], parent_step_id=previousId, 
                                source=tmpFile, dest=finalname)

def createMicroscope(log,fnOut,Voltage,SphericalAberration,SamplingRate,Magnification):
    MD = xmipp.MetaData()
    MD.setColumnFormat(False)
    id = MD.addObject()
    MD.setValue(xmipp.MDL_CTF_VOLTAGE,float(Voltage),id)    
    MD.setValue(xmipp.MDL_CTF_CS,float(SphericalAberration),id)    
    MD.setValue(xmipp.MDL_CTF_SAMPLING_RATE,float(SamplingRate),id)
    MD.setValue(xmipp.MDL_MAGNIFICATION,float(Magnification),id)
    MD.write(fnOut)    

def gatherResults(log, WorkingDir, summaryFile):
    MD = xmipp.MetaData()
    for filename in glob(join(WorkingDir, "*")):
        if filename.endswith(".xmd") or filename.endswith("tmp") or filename.endswith(".sel"):
            continue
        objId = MD.addObject()
        MD.setValue(xmipp.MDL_IMAGE, filename, objId)
    if MD.size() != 0:
        MD.sort(xmipp.MDL_IMAGE);
        MD.write(summaryFile)
