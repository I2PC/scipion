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
        self.Db.insertStep('sortImages',inputFile=self.InputFile,outputFile=self.outputFile,addToSelf=self.addToSelf,
                           rejectionMethod=self.RejectionMethod, maxZscore=self.MaxZscore, percentage=self.Percentage)
                
    def summary(self):
        message=[]
        message.append("Screening of [%s]" % self.InputFile)
        if self.RejectionMethod=='MaxZscore':
            message.append('Rejecting: Zscore>'+str(self.MaxZscore))
        elif self.RejectionMethod=='Percentage':
            message.append('Rejecting: '+str(self.Percentage)+"%")
        return message
    
    def validate(self):
        errors=[]
        ext=os.path.splitext(self.InputFile)[1]
        if (ext!=".xmd" and ext!=".sel" and ext!=".stk" and ext!=".mrc" and ext!=".hed" and ext!=".img"):
            errors.append("This protocol is designed for stacks")
        if self.RejectionMethod=='MaxZscore' and self.MaxZscore<0:
            errors.append("MaxZscore must be positive")
        elif self.RejectionMethod=='Percentage' and (self.Percentage<0 or self.Percentage>100):
            errors.append("Percentage must be between 0 and 100")
        return errors

    def visualize(self):
        if not os.path.exists(self.outputFile):
            from protlib_gui_ext import showWarning
            showWarning("Error", "There is no result yet",self.master)
        else:   
            runShowJ(self.outputFile)                                     

def sortImages(log, inputFile, outputFile, addToSelf, rejectionMethod, maxZscore, percentage):
    from xmipp import MetaDataInfo
    (Xdim, Ydim, Zdim, Ndim, _) = MetaDataInfo(inputFile)
    args=""
    if rejectionMethod=='MaxZscore':
        args+=" --zcut "+str(maxZscore)
    elif rejectionMethod=='Percentage':
        args+=" --percent "+str(percentage)
    if Ndim > 0:
        if addToSelf: 
            runJob(log,"xmipp_image_sort_by_statistics","-i "+outputFile+" --addToInput"+args)
        else:
            runJob(log,"xmipp_image_sort_by_statistics","-i "+inputFile+" -o "+outputFile+args)

