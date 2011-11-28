#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using hierarchical clustering principles
# Author: Carlos Oscar Sanchez Sorzano, August 2011
#

import glob, os, sys, shutil,time
from os.path import exists
from protlib_base import *
from config_protocols import protDict
from protlib_filesystem import renameFile, deleteFile
from xmipp import MetaData, Image

class ProtCL2DAlignment(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.cl2d_align.name, scriptname, project)
        self.Import = 'from protocol_cl2d_align import *'    

    def defineSteps(self):
        self.insertCl2dStep()
        self.Db.insertStep('gatherResults',
                           verifyfiles=[self.workingDirPath("average.xmp"),self.workingDirPath("alignment.xmd")],
                           WorkingDir=self.WorkingDir)
    
    def summary(self):
        message=["Alignment of "+self.InSelFile]
        fnAlignment=self.workingDirPath("results_level_00_classes.xmd")
        if exists(fnAlignment):
            mD=MetaData("info@"+fnAlignment)
            date=time.ctime(os.path.getmtime(fnAlignment))
            message.append("Iteration "+str(mD.size())+" at "+date)
        return message
    
    def validate(self):
        errors = []
        if self.ReferenceImage != "":
            if not exists(self.ReferenceImage):
                errors.append("Cannot find the file "+self.ReferenceImage)
        return errors
    
    def visualize(self):
        from protlib_utils import runShowJ
        if self.getRunState() == SqliteDb.RUN_FINISHED:
            runShowJ("%s %s" % (self.workingDirPath("average.xmp"),self.workingDirPath("alignment.xmd")))
        else:
            if exists(self.workingDirPath("results_level_00_classes.stk")):
                runShowJ("%s %s" % (self.workingDirPath("results_level_00_classes.stk"),self.workingDirPath("results_level_00_classes.xmd")))
    
    def insertCl2dStep(self):
    #log,Selfile,WorkingDir,ReferenceImage,MaxShift,NumberOfIterations,Nproc):
        params= '-i %s --oroot %s/results --nref 1 --iter %d --maxShift %d' % 
                (self.InSelFile, self.WorkingDir, self.NumberOfIterations, self.MaxShift)
                
        if self.ReferenceImage!="":
            params += " --ref0 " + self.ReferenceImage
        else:
            params += " --nref0 1"
        self.insertRunJobStep("xmipp_classify_CL2D", params, [self.workingDirPath("results_images.xmd"))

def gatherResults(log, WorkingDir):
    wdPath = lambda path: os.path.join(WorkingDir, path)
    renameFile(log, wdPath("results_images.xmd"), wdPath("alignment.xmd"))
    fnStack=wdPath("results_level_00_classes.stk")
    I=Image("1@"+fnStack)
    I.write(wdPath("average.xmp"))
    deleteFile(log,fnStack)
    deleteFile(log,wdPath("results_level_00_classes.xmd"))
