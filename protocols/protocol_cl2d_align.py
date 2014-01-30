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
from protlib_filesystem import renameFile, xmippExists
from protlib_utils import runJob
from xmipp import MetaData, Image

class ProtCL2DAlignment(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.cl2d_align.name, scriptname, project)
        self.Import = 'from protocol_cl2d_align import *'    

    def defineSteps(self):
        self.Db.insertStep("linkAcquisitionInfo",InputFile=self.InSelFile,dirDest=self.WorkingDir)
        self.insertCl2dStep()
        self.Db.insertStep('gatherResults',
                           verifyfiles=[self.workingDirPath("average.xmp"),self.workingDirPath("images.xmd")],
                           WorkingDir=self.WorkingDir)
    
    def summary(self):
        message=["Alignment of "+self.InSelFile]
        fnAlignment=self.workingDirPath("level_00/level_classes.xmd")
        if exists(fnAlignment):
            from protlib_xmipp import getMdSize
            size = getMdSize("info@"+fnAlignment)
            date = time.ctime(os.path.getmtime(fnAlignment))
            message.append("Iteration %d at %s" % (size, date))
        return message
    
    def papers(self):
        papers=[]
        papers.append('Sorzano, JSB (2010) [http://www.ncbi.nlm.nih.gov/pubmed/20362059]')
        return papers

    def validate(self):
        errors = []
        if self.ReferenceImage != "":
            if not xmippExists(self.ReferenceImage):
                errors.append("Cannot find the file "+self.ReferenceImage)
        return errors
    
    def visualize(self):
        from protlib_utils import runShowJ
        if self.getRunState() == SqliteDb.RUN_FINISHED:
            runShowJ("%s %s" % (self.workingDirPath("average.xmp"),self.workingDirPath("images.xmd")))
        else:
            if exists(self.workingDirPath("level_00/level_classes.stk")):
                runShowJ("%s %s" % (self.workingDirPath("level_00/level_classes.stk"),self.workingDirPath("level_00/level_classes.xmd")))
    
    def insertCl2dStep(self):
        params= '-i %s --oroot level --odir %s --nref 1 --iter %d --maxShift %d' % \
                (self.InSelFile, self.WorkingDir, self.NumberOfIterations, self.MaxShift)
                
        if self.ReferenceImage!="":
            params += " --ref0 " + self.ReferenceImage
        else:
            params += " --nref0 1"
        self.insertRunJobStep("xmipp_classify_CL2D", params, [self.workingDirPath("images.xmd")])

def gatherResults(log, WorkingDir):
    wdPath = lambda path: os.path.join(WorkingDir, path)
    fnStack=wdPath("level_00/level_classes.stk")
    I=Image("1@"+fnStack)
    I.write(wdPath("average.xmp"))
    runJob(log,"rm","-rf "+wdPath("level_00"))
