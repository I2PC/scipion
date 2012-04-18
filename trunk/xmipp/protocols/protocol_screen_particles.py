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
from protlib_filesystem import deleteFile, deleteDir, createLink, copyFile

class ProtScreenParticles(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.screen_particles.name, scriptname, project)
        self.Import = 'from protocol_screen_particles import *'

    def defineSteps(self):
        outputFile=self.workingDirPath("sorted.xmd")
        self.Db.insertStep("copyAcquisition",inputFile=self.InputFile,WorkingDir=self.WorkingDir)
        self.Db.insertStep('sortImages',verifyfiles=[outputFile],inputFile=self.InputFile, outputFile=outputFile)
                
    def validate(self):
        return []

    def summary(self):
        message=[]
        message.append("Screening of [%s]" % self.InputFile)
        return message

    def visualize(self):
        summaryFile = self.workingDirPath("sorted.xmd")
        if not os.path.exists(summaryFile):
            from protlib_gui_ext import showWarning
            showWarning("Error", "There is no result yet")
        else:   
            runShowJ(summaryFile)                                     

def copyAcquisition(log,inputFile,WorkingDir):
    inputFileDir=os.path.dirname(inputFile)
    fnAcquistionIn=os.path.join(inputFileDir,"acquisition_info.xmd")
    if os.path.exists(fnAcquistionIn):
        fnAcquistionOut=os.path.join(WorkingDir,"acquisition_info.xmd")
        createLink(log, fnAcquistionIn, fnAcquistionOut)

def sortImages(log, inputFile, outputFile):
    from xmipp import MetaData
    mD = MetaData(inputFile)
    if mD.size() > 0:
        runJob(log,"xmipp_image_sort_by_statistics","-i "+inputFile+" --multivariate -o "+outputFile)
    else:
        createLink(log,inputFile,outputFile)

