#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using hierarchical clustering principles
# Author: Carlos Oscar Sanchez Sorzano, August 2011
#

import glob,os,sys,shutil,time
from protlib_base import *
from config_protocols import protDict
from protlib_utils import runJob
from protlib_filesystem import renameFile, deleteFile
from xmipp import MetaData, Image

class ProtCL2DAlignment(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.cl2d_alignment.name, scriptname, project)
        self.Import = 'from protocol_cl2d_align import *'    

    def defineSteps(self):
        self.Db.insertStep('cl2d',verifyfiles=[self.workingDirPath("results_images.xmd")],
                           Selfile=self.InSelFile,WorkingDir=self.WorkingDir,ReferenceImage=self.ReferenceImage,
                           MaxShift=self.MaxShift,NumberOfIterations=self.NumberOfIterations,Nproc=self.NumberOfMpi)
        self.Db.insertStep('gatherResults',
                           verifyfiles=[self.workingDirPath("average.xmp"),self.workingDirPath("alignment.xmd")],
                           WorkingDir=self.WorkingDir)
    
    def summary(self):
        message=["Alignment of "+self.InSelFile]
        fnAlignment=self.workingDirPath("results_level_00_classes.xmd")
        if os.path.exists(fnAlignment):
            mD=MetaData("info@"+fnAlignment)
            date=time.ctime(os.path.getmtime(fnAlignment))
            message.append("Iteration "+str(mD.size())+" at "+date)
        return message
    
    def validate(self):
        errors=[]
        if self.ReferenceImage!="":
            if not os.path.exists(self.ReferenceImage):
                errors.append("Cannot find the file "+self.ReferenceImage)
        return errors
    
    def visualize(self):
        if self.getRunState()==SqliteDb.RUN_FINISHED:
            os.system("xmipp_showj -i %s %s &"
                      %(self.workingDirPath("average.xmp"),self.workingDirPath("alignment.xmd")))
        else:
            if os.path.exists(self.workingDirPath("results_level_00_classes.stk")):
                os.system("xmipp_showj -i %s %s &"
                          %(self.workingDirPath("results_level_00_classes.stk"),self.workingDirPath("results_level_00_classes.xmd")))
    
def cl2d(log,Selfile,WorkingDir,ReferenceImage,MaxShift,NumberOfIterations,Nproc):
    params= '-i '+str(Selfile)+' --oroot '+WorkingDir+'/results --nref 1 --iter '+str(NumberOfIterations)+" --maxShift "+str(MaxShift)
    if ReferenceImage!="":
        params+=" --ref0 "+ReferenceImage
    else:
        params+=" --nref0 1"
    runJob(log,"xmipp_classify_CL2D",params,Nproc)

def gatherResults(log,WorkingDir):
    renameFile(log, os.path.join(WorkingDir,"results_images.xmd"), os.path.join(WorkingDir,"alignment.xmd"))
    fnStack=os.path.join(WorkingDir,"results_level_00_classes.stk")
    I=Image("1@"+fnStack)
    I.write(os.path.join(WorkingDir,"average.xmp"))
    deleteFile(log,fnStack)
    deleteFile(log,os.path.join(WorkingDir,"results_level_00_classes.xmd"))
