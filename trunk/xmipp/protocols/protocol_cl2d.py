#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using hierarchical clustering principles
# Author: Carlos Oscar Sanchez Sorzano, August 2011
#

import glob,os,sys,shutil,time
from protlib_base import *
from config_protocols import protDict
from protlib_utils import runJob, getRangeValuesFromString
from protlib_filesystem import createLink, deleteFile
from xmipp import MetaData, MD_APPEND

class ProtCL2D(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.cl2d.name, scriptname, project)
        self.Import = 'from protocol_cl2d import *'    

    def defineSteps(self):
        fnImages=self.workingDirPath("results_images.xmd")
        self.Db.insertStep('cl2d',verifyfiles=[fnImages],Selfile=self.InSelFile,WorkingDir=self.WorkingDir,
                           NumberOfReferences=self.NumberOfReferences,NumberOfInitialReferences=self.NumberOfInitialReferences,
                           NumberOfIterations=self.NumberOfIterations,ComparisonMethod=self.ComparisonMethod,
                           ClusteringMethod=self.ClusteringMethod,AdditionalParameters=self.AdditionalParameters,
                           Nproc=self.NumberOfMpi)
        self.Db.insertStep('evaluateClasses',WorkingDir=self.WorkingDir,pattern="results_level_??_classes.xmd")
        if self.NumberOfReferences>self.NumberOfInitialReferences:
            self.Db.insertStep('core_analysis',verifyfiles=[self.workingDirPath("results_level_00_classes_core.xmd")],
                               WorkingDir=self.WorkingDir,thZscore=self.thZscore,thPCAZscore=self.thPCAZscore,Nproc=self.NumberOfMpi)
            self.Db.insertStep('evaluateClasses',WorkingDir=self.WorkingDir,pattern="results_level_??_classes_core.xmd")
            self.Db.insertStep('stable_core_analysis',
                               verifyfiles=[self.workingDirPath("results_level_%02d_classes_stable_core.xmd"%(self.Tolerance+1))],
                               WorkingDir=self.WorkingDir,tolerance=self.Tolerance,Nproc=self.NumberOfMpi)
            self.Db.insertStep('evaluateClasses',WorkingDir=self.WorkingDir,pattern="results_level_??_classes_stable_core.xmd")
        self.Db.insertStep('sortClasses',WorkingDir=self.WorkingDir,Nproc=self.NumberOfMpi)
    
    def summary(self):
        message=[]
        message.append("Classification of "+self.InSelFile+" into "+str(self.NumberOfReferences)+" classes")
        levelFiles=glob.glob(self.WorkingDir+"/results_level_??_classes.xmd")
        if not levelFiles:
            message.append("No class file has been generated")
        else:
            import re, time
            levelFiles.sort()
            lastLevelFile=levelFiles[-1]
            date=time.ctime(os.path.getmtime(lastLevelFile))
            lastLevel=int(re.search('level_(\d\d)',lastLevelFile).group(1))
            mD=MetaData("info@"+lastLevelFile)
            iteration=mD.size()
            try:
                mD=MetaData("classes@"+lastLevelFile)
                message.append("Last iteration is %s from level %d with %d classes (at %s)"%(iteration,lastLevel,mD.size(),date))
            except:
                pass
        return message
    
    def validate(self):
        errors = []
        if self.NumberOfInitialReferences>self.NumberOfReferences:
            errors.append("The number of initial classes cannot be larger than the number of final classes")
        if self.Tolerance<0:
            errors.append("Tolerance must be larger than 0")
        return errors
    
    def visualize(self):
        if self.WhatToShow=="Classes":
            levelFiles=glob.glob(self.workingDirPath("results_level_??_classes.xmd"))
        elif self.WhatToShow=="Class Cores":
            levelFiles=glob.glob(self.workingDirPath("results_level_??_classes_core.xmd"))
        elif self.WhatToShow=="Class Stable Cores":
            levelFiles=glob.glob(self.workingDirPath("results_level_??_classes_stable_core.xmd"))
        if levelFiles:
            levelFiles.sort()
            if self.DoShowLast:
                lastLevelFile=levelFiles[-1]
                os.system("xmipp_metadata_viewerj -i "+lastLevelFile+"&")
            else:
                listOfLevels=getRangeValuesFromString(self.LevelsToShow)
                files=""
                for level in listOfLevels:
                    files+=levelFiles[level]+" "
                if files!="":
                    print "xmipp_metadata_viewerj -i "+files+" &"
                    os.system("xmipp_metadata_viewerj -i "+files+" &")
    
def cl2d(log,Selfile,WorkingDir,NumberOfReferences,NumberOfInitialReferences,NumberOfIterations,
         ComparisonMethod,ClusteringMethod,AdditionalParameters,Nproc):
    params= '-i '+str(Selfile)+' --oroot '+WorkingDir+'/results '+\
            ' --nref '+str(NumberOfReferences)+\
            ' --nref0 '+str(NumberOfInitialReferences)+\
            ' --iter '+str(NumberOfIterations)+\
            ' '+AdditionalParameters
    if ComparisonMethod=='correlation':
        params+= ' --distance correlation '
    if ClusteringMethod=='classical':
        params+= ' --classicalMultiref '

    runJob(log,"xmipp_classify_CL2D",params,Nproc)
    levelFiles=glob.glob(WorkingDir+"/results_level_*.xmd")
    if not levelFiles:
        levelFiles.sort()
        lastLevelFile=levelFiles[-1]
        mD = MetaData("classes@"+lastLevelFile)
        if mD.size()==NumberOfReferences:
            createLink(log, lastLevelFile, WorkingDir+"results_classes.sel")

def sortClasses(log,WorkingDir,Nproc):
    import re
    for file in glob.glob(os.path.join(WorkingDir,"results_level_??_classes.xmd")):
        level=int(re.search('level_(\d\d)',file).group(1))
        fnRoot=os.path.join(WorkingDir,"results_level_%02d_classes_sorted"%level)
        params= "-i classes@"+file+" --oroot "+fnRoot
        runJob(log,"xmipp_image_sort",params,Nproc)
        mD=MetaData(fnRoot+".xmd")
        mD.write("classes_sorted@"+file,MD_APPEND)
        deleteFile(log,fnRoot+".xmd")

def core_analysis(log,WorkingDir,thZscore,thPCAZscore,Nproc):
    params= "-i "+WorkingDir+'/results --computeCore '+str(thZscore)+' '+str(thPCAZscore)
    runJob(log,"xmipp_classify_CL2D_core_analysis",params,Nproc)

def stable_core_analysis(log,WorkingDir,tolerance,Nproc):
    params= "-i "+WorkingDir+'/results --computeStableCore '+str(tolerance)
    runJob(log,"xmipp_classify_CL2D_core_analysis",params,Nproc)

def evaluateClasses(log,WorkingDir,pattern):
    levelFiles=glob.glob(os.path.join(WorkingDir,pattern))
    for file in levelFiles:
        runJob(log,"xmipp_classify_evaluate_classes","-i "+file)
