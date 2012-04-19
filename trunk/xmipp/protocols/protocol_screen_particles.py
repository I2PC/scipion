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
from protlib_filesystem import createLink, linkAcquisitionInfoIfPresent

class ProtScreenParticles(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.screen_particles.name, scriptname, project)
        self.Import = 'from protocol_screen_particles import *'
        self.outputFile=self.workingDirPath(os.path.basename(self.InputFile))

    def defineSteps(self):
        self.Db.insertStep("copyFile",verifyfiles=[self.outputFile],source=self.InputFile,dest=self.outputFile)
        self.Db.insertStep("linkAcquisitionInfoIfPresent",InputFile=self.InputFile,dirDest=self.WorkingDir)
        self.Db.insertStep('sortImages',inputFile=self.outputFile)
                
    def summary(self):
        message=[]
        message.append("Screening of [%s]" % self.InputFile)
        return message

    def visualize(self):
        if not os.path.exists(self.outputFile):
            from protlib_gui_ext import showWarning
            showWarning("Error", "There is no result yet")
        else:   
            runShowJ(self.outputFile)                                     

def sortImages(log, inputFile):
    from xmipp import ImgSize
    (Xdim, Ydim, Zdim, Ndim) = ImgSize(inputFile)
    if Ndim > 0:
        runJob(log,"xmipp_image_sort_by_statistics","-i "+inputFile+" --multivariate --addToInput")

