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
        self.Db.insertStep('cl2d',verifyfiles=[os.path.join(self.WorkingDir,"results_images.xmd")],
                           Selfile=self.InSelFile,WorkingDir=self.WorkingDir,
                           NumberOfIterations=self.NumberOfIterations,Nproc=self.NumberOfMpi)
        self.Db.insertStep('gatherResults',
                           verifyfiles=[os.path.join(self.WorkingDir,"average.xmp"),os.path.join(self.WorkingDir,"alignment.xmd")],
                           WorkingDir=self.WorkingDir)
    
    def summary(self):
        return ["Alignment of "+self.InSelFile]
        # Agnadir la iteracion por la que va
    
    def validate(self):
        return []
    
    def visualize(self):
        levelFiles=glob.glob(os.path.join(self.WorkingDir,"results_level_*.xmd"))
        if levelFiles:
            levelFiles.sort()
            if self.DoShowLast:
                lastLevelFile=levelFiles[-1]
                os.system("xmipp_metadata_viewerj -i "+lastLevelFile+"&")
            else:
                listOfLevels=getRangeValuesFromString(LevelsToShow)
                files=""
                for level in listOfLevels:
                    files+=os.path.join(self.WorkingDir,"results_level_%02d_classes.xmd"%level)
                if files!="":
                    os.system("xmipp_metadata_viewerj -i "+files+" &")
    
def cl2d(log,Selfile,WorkingDir,NumberOfIterations,Nproc):
    params= '-i '+str(Selfile)+' --oroot '+WorkingDir+'/results --nref 1 --nref0 1'+\
            ' --iter '+str(NumberOfIterations)
    runJob(log,"xmipp_classify_CL2D",params,Nproc)

def gatherResults(log,WorkingDir):
    renameFile(log, os.path.join(WorkingDir,"results_images.xmd"), os.path.join(WorkingDir,"alignment.xmd"))
    fnStack=os.path.join(WorkingDir,"results_level_00_classes.stk")
    I=Image("1@"+fnStack)
    I.write(os.path.join(WorkingDir,"average.xmp"))
    deleteFile(log,fnStack)
    deleteFile(log,os.path.join(WorkingDir,"results_level_00_classes.xmd"))
