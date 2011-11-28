#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using hierarchical clustering principles
# Author: Carlos Oscar Sanchez Sorzano, August 2011
#

import glob,os,re,sys,shutil,time
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
        self.insertCl2dStep()
        self.Db.insertStep('evaluateClasses',WorkingDir=self.WorkingDir,subset="classes")
        if self.NumberOfReferences > self.NumberOfInitialReferences:
            # core analysis
            params= "-i %(WorkingDir)s/results --computeCore %(thZscore)f %(thPCAZscore)f" % self.ParamsDict
            self.insertRunJobStep("xmipp_classify_CL2D_core_analysis",params,
                      [self.workingDirPath("results_level_00_classes_core.xmd")])
            # evaluate classes 
            self.Db.insertStep('evaluateClasses',WorkingDir=self.WorkingDir,subset="classes_core")
            # stable core analysis
            params= "-i %(WorkingDir)s/results --computeStableCore %(Tolerance)f" % self.ParamsDict
            self.insertRunJobStep("xmipp_classify_CL2D_core_analysis", params)
            # evaluate classes again
            self.Db.insertStep('evaluateClasses',WorkingDir=self.WorkingDir,subset="classes_stable_core")
        self.Db.insertStep('sortClasses',WorkingDir=self.WorkingDir,Nproc=self.NumberOfMpi)
    
    def summary(self):
        message=[]
        message.append("Classification of "+self.InSelFile+" into "+str(self.NumberOfReferences)+" classes")
        levelFiles=glob.glob(self.WorkingDir+"/results_level_??_classes.xmd")
        if not levelFiles:
            message.append("No class file has been generated")
        else:
            import time
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
        from protlib_utils import runShowJ
        fnSubset=""
        if self.WhatToShow=="Classes":
            fnSubset="classes"
        elif self.WhatToShow=="Class Cores":
            fnSubset="classes_core"
        elif self.WhatToShow=="Class Stable Cores":
            fnSubset="classes_stable_core"
        levelFiles=glob.glob(self.workingDirPath("results_level_??_%s.xmd"%fnSubset))
        if levelFiles:
            levelFiles.sort()
            if self.DoShowLast:
                lastLevelFile=levelFiles[-1]
                runShowJ(lastLevelFile)
            else:
                listOfLevels = getRangeValuesFromString(self.LevelsToShow)
                files = " ".join([levelFiles[level] for level in listOfLevels])
                if files!="":
                    runShowJ(files)
        if self.DoShowHierarchy:
            fnHierarchy=self.workingDirPath(fnSubset+"_hierarchy.txt")
            if os.path.exists(fnHierarchy):
                from protlib_gui_ext import showTextfileViewer
                showTextfileViewer(fnHierarchy,[fnHierarchy])
                
    def insertCl2dStep(self):
        params= '-i %(InSelFile)s --oroot %(WorkingDir)s/results '+\
                ' --nref %(NumberOfReferences)d'+\
                ' --nref0 %(NumberOfInitialReferences)d'+\
                ' --iter %(NumberOfIterations)d'+\
                ' %(AdditionalParameters)s'
        if self.ComparisonMethod=='correlation':
            params+= ' --distance correlation '
        if self.ClusteringMethod=='classical':
            params+= ' --classicalMultiref '
    
        self.insertRunJobStep("xmipp_classify_CL2D", params % self.ParamsDict, 
                              [self.workingDirPath("results_images.xmd")])
        self.insertStep('postCl2d', WorkingDir=self.WorkingDir, 
                        NumberOfReferences=self.NumberOfReferences)
    
def postCl2d(log, WorkingDir, NumberOfReferences):
    levelFiles=glob.glob(WorkingDir+"/results_level_*.xmd")
    if not levelFiles:
        levelFiles.sort()
        lastLevelFile=levelFiles[-1]
        mD = MetaData("classes@"+lastLevelFile)
        if mD.size()==NumberOfReferences:
            createLink(log, lastLevelFile, WorkingDir+"results_classes.sel")

def sortClasses(log,WorkingDir,Nproc):
    for filename in glob.glob(os.path.join(WorkingDir,"results_level_??_classes.xmd")):
        level=int(re.search('level_(\d\d)',filename).group(1))
        fnRoot=os.path.join(WorkingDir,"results_level_%02d_classes_sorted"%level)
        params= "-i classes@"+filename+" --oroot "+fnRoot
        runJob(log,"xmipp_image_sort",params,Nproc)
        mD=MetaData(fnRoot+".xmd")
        mD.write("classes_sorted@"+filename,MD_APPEND)
        deleteFile(log,fnRoot+".xmd")

def evaluateClasses(log,WorkingDir,subset):
    levelFiles=glob.glob(os.path.join(WorkingDir,"results_level_??_%s.xmd"%subset))
    levelFiles.sort()
    for filename in levelFiles:
        runJob(log,"xmipp_classify_evaluate_classes","-i "+filename)
        level=int(re.search('level_(\d\d)',filename).group(1))
        if level>0:
            previousFile=os.path.join(WorkingDir,"results_level_%02d_%s.xmd"%(level-1,subset))
            if os.path.exists(previousFile):
                fnOut=os.path.join(WorkingDir,subset+"_hierarchy.txt")
                args="--i1 %s --i2 %s -o %s"%(previousFile,filename,fnOut)
                if os.path.exists(fnOut):
                    args+=" --append"
                runJob(log,"xmipp_classify_compare_classes",args)
