#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based pre-processing of single-particles: 
#  - phase flipping
#  - extraction of particles
#  - normalization
#  - sort_junk

# Author: Carlos Oscar, August 2011
#
from protlib_base import *
import glob
import os
from protlib_utils import runJob, runShowJ
from protlib_filesystem import createLink, linkAcquisitionInfo

class ProtScreenParticles(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.screen_particles.name, scriptname, project)
        self.Import = 'from protocol_screen_particles import *'
        ext=os.path.splitext(self.InputFile)[1]
        self.addToSelf=True;
        if ext==".xmd" or ext==".sel": 
            self.outputFile=self.workingDirPath(os.path.basename(self.InputFile))
        else:
            baseFile=os.path.splitext(os.path.basename(self.InputFile))[0]
            self.outputFile=self.workingDirPath(baseFile+".xmd")
            self.addToSelf=False;

    def defineSteps(self):
        self.Db.insertStep("linkAcquisitionInfo",InputFile=self.InputFile,dirDest=self.WorkingDir)
        if self.addToSelf:
            self.Db.insertStep("copyFile",verifyfiles=[self.outputFile],source=self.InputFile,dest=self.outputFile)
        self.Db.insertStep('sortImages',inputFile=self.InputFile,outputFile=self.outputFile,addToSelf=self.addToSelf)
                
    def summary(self):
        message=[]
        message.append("Screening of [%s]" % self.InputFile)
        return message
    
    def validate(self):
        errors=[]
        ext=os.path.splitext(self.InputFile)[1]
        if (ext!=".xmd" and ext!=".sel" and ext!=".stk" and ext!=".mrc" and ext!=".hed" and ext!=".img"):
            errors.append("This protocol is designed for stacks")
        return errors

    def visualize(self):
        if not os.path.exists(self.outputFile):
            from protlib_gui_ext import showWarning
            showWarning("Error", "There is no result yet",self.master)
        else:   
            runShowJ(self.outputFile)                                     

def sortImages(log, inputFile, outputFile, addToSelf):
    from xmipp import MetaDataInfo
    (Xdim, Ydim, Zdim, Ndim, _) = MetaDataInfo(inputFile)
    if Ndim > 0:
        if addToSelf: 
            runJob(log,"xmipp_image_sort_by_statistics","-i "+outputFile+" --addToInput")
        else:
            runJob(log,"xmipp_image_sort_by_statistics","-i "+inputFile+" -o "+outputFile)

