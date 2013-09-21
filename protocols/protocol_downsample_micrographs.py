#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# General script for Xmipp-based pre-processing of micrographs
# Author: Carlos Oscar Sorzano, July 2011

from glob import glob
from protlib_base import *
import xmipp
from protlib_filesystem import replaceBasenameExt, createLink, createDir
from protlib_utils import runJob


class ProtDownsampleMicrographs(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.downsample_micrographs.name, scriptname, project)
        self.Import = "from protocol_downsample_micrographs import *"
        self.setPreviousRun(self.ImportRun) 
        self.importDir = self.PrevRun.WorkingDir
        self.inputFilename('microscope', 'micrographs', 'acquisition')
        self.inputProperty('TiltPairs', 'MicrographsMd')
        
    def defineSteps(self):
        self.insertImportOfFiles([self.Input['microscope']])
        #self.insertStep('createLink2', filename="microscope.xmd",dirSrc=self.importDir,dirDest=self.WorkingDir)
        self.insertStep("createDir",verifyfiles=[self.ExtraDir],path=self.ExtraDir)
        self.insertStep("changeSamplingRate",fnIn=self.Input['acquisition'], fnOut=self.getFilename('acquisition'),
                        downsampleFactor=self.DownsampleFactor)

        MD = xmipp.MetaData(self.Input['micrographs'])
        previousId = XmippProjectDb.FIRST_STEP
        IOTable = {}
        MD.removeDisabled();
        for i in MD:
            fnMicrograph = MD.getValue(xmipp.MDL_MICROGRAPH,i)
            fnOut = os.path.join(self.ExtraDir,replaceBasenameExt(fnMicrograph, '.mrc'))
            IOTable[fnMicrograph] = fnOut
            self.insertParallelStep("doDownsample", verifyfiles=[fnOut], parent_step_id=XmippProjectDb.FIRST_STEP, 
                            fnMicrograph=fnMicrograph, fnOut=fnOut, downsampleFactor=self.DownsampleFactor)
        self.insertStep("gatherResults",WorkingDir=self.WorkingDir,ImportDir=self.importDir,IOTable=IOTable,
                        downsampleFactor=self.DownsampleFactor)

    def validate(self):
        errors = []
        if self.DownsampleFactor<=1:
            errors.append("Downsampling must be larger than 1");
        return errors

    def summary(self):
        message = []
        message.append("Downsampling of micrographs from [%s] by a factor <%3.2f>" % (self.importDir,self.DownsampleFactor))
        return message

    def visualize(self):
        summaryFile = self.getFilename('micrographs')
        if os.path.exists(summaryFile):
            from protlib_utils import runShowJ
            runShowJ(summaryFile)

def doDownsample(log,fnMicrograph,fnOut,downsampleFactor):
    runJob(log,"xmipp_transform_downsample", "-i %s -o %s --step %f --method fourier" % (fnMicrograph,fnOut,downsampleFactor))

def changeSamplingRate(log,fnIn,fnOut,downsampleFactor):
    MD=xmipp.MetaData()
    MD.read(fnIn)
    i=MD.firstObject()
    Ts=MD.getValue(xmipp.MDL_SAMPLINGRATE,i)
    MD.setValue(xmipp.MDL_SAMPLINGRATE,Ts*downsampleFactor,i)
    MD.setValue(xmipp.MDL_SAMPLINGRATE_ORIGINAL,Ts,i)
    MD.write(fnOut)

def convertMetaData(fnIn,fnOut,blockname,IOTable,downsampleFactor,tiltPairs):
    MD=xmipp.MetaData(fnIn)
    MD.removeDisabled()
    for i in MD:
        fnMicrograph=MD.getValue(xmipp.MDL_MICROGRAPH,i);
        MD.setValue(xmipp.MDL_MICROGRAPH,IOTable[fnMicrograph],i)
        MD.setValue(xmipp.MDL_MICROGRAPH_ORIGINAL,fnMicrograph,i)
        if tiltPairs:
            fnMicrographTilted=MD.getValue(xmipp.MDL_MICROGRAPH_TILTED,i);
            MD.setValue(xmipp.MDL_MICROGRAPH_TILTED,IOTable[fnMicrographTilted],i)
            MD.setValue(xmipp.MDL_MICROGRAPH_TILTED_ORIGINAL,fnMicrographTilted,i)
    MD.write(blockname+"@"+fnOut)


def gatherResults(log, WorkingDir, ImportDir,IOTable,downsampleFactor):
    convertMetaData(os.path.join(ImportDir,"micrographs.xmd"),os.path.join(WorkingDir,"micrographs.xmd"),
                    "micrographs",IOTable,downsampleFactor,False)
    
    fnTilted=os.path.join(ImportDir,"tilted_pairs.xmd")
    if os.path.exists(fnTilted):
        convertMetaData(fnTilted,os.path.join(WorkingDir,"tilted_pairs.xmd"),
                        "micrographPairs",IOTable,downsampleFactor,True)
    