#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using hierarchical clustering principles
# Author: Carlos Oscar Sanchez Sorzano, August 2011
#

import glob,os,re,sys,shutil,time
from protlib_base import *
from config_protocols import protDict
from protlib_xmipp import getMdSize
from protlib_utils import runJob, getListFromRangeString
from protlib_filesystem import createLink, deleteFile, linkAcquisitionInfo
from xmipp import MetaData, MD_APPEND, MDL_IMAGE

class ProtCL2D(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.cl2d.name, scriptname, project)
        self.Import = 'from protocol_cl2d import *'    

    def createFilenameTemplates(self):
        return {
            'hierarchy': "%(ExtraDir)s/classes%(subset)s_hierarchy.txt"
            }

    def defineSteps(self):
        self.insertStep('createDir',path=self.ExtraDir)
        self.Db.insertStep("linkAcquisitionInfo",InputFile=self.InSelFile,dirDest=self.WorkingDir)
        self.insertCl2dStep()
        self.Db.insertStep('evaluateClasses',WorkingDir=self.WorkingDir,ExtraDir=self.ExtraDir,subset="")
        self.Db.insertStep('sortClasses',ExtraDir=self.ExtraDir,Nproc=self.NumberOfMpi,suffix="")

        if self.NumberOfReferences > self.NumberOfInitialReferences:
            # core analysis
            params= "--dir %(ExtraDir)s --root level --computeCore %(thZscore)f %(thPCAZscore)f" % self.ParamsDict
            self.insertRunJobStep("xmipp_classify_CL2D_core_analysis",params,
                      [self.workingDirPath("extra/level_00/level_classes_core.xmd")])
            # evaluate classes 
            self.Db.insertStep('evaluateClasses',WorkingDir=self.WorkingDir,ExtraDir=self.ExtraDir,subset="_core")
            self.Db.insertStep('sortClasses',ExtraDir=self.ExtraDir,Nproc=self.NumberOfMpi,suffix="_core")
            # stable core analysis
            params= "--dir %(ExtraDir)s --root level --computeStableCore %(Tolerance)d" % self.ParamsDict
            self.insertRunJobStep("xmipp_classify_CL2D_core_analysis", params)
            # evaluate classes again
            self.Db.insertStep('evaluateClasses',WorkingDir=self.WorkingDir,ExtraDir=self.ExtraDir,subset="_stable_core")
            self.Db.insertStep('sortClasses',ExtraDir=self.ExtraDir,Nproc=self.NumberOfMpi,suffix="_stable_core")

        self.Db.insertStep("postEvaluation",WorkingDir=self.WorkingDir)
    
    def summary(self):
        message=[]
        message.append(("Classification of [%s]"%self.InSelFile)+" into "+str(self.NumberOfReferences)+" classes")
        levelFiles=glob.glob(self.WorkingDir+"/extra/level_??/level_classes.xmd")
        if not levelFiles:
            message.append("No class file has been generated")
        else:
            import time
            levelFiles.sort()
            lastLevelFile=levelFiles[-1]
            date=time.ctime(os.path.getmtime(lastLevelFile))
            lastLevel=int(re.search('level_(\d\d)',lastLevelFile).group(1))
            iteration = getMdSize("info@"+lastLevelFile)
            try:
                size = getMdSize("classes@"+lastLevelFile)
                message.append("Last iteration is %s from level %d with %d classes (at %s)"%(iteration,lastLevel, size,date))
            except:
                pass
        return message
    
    def papers(self):
        papers=[]
        papers.append('Sorzano, JSB (2010) [http://www.ncbi.nlm.nih.gov/pubmed/20362059]')
        return papers

    def validate(self):
        errors = []
        if self.NumberOfInitialReferences>self.NumberOfReferences:
            errors.append("The number of initial classes cannot be larger than the number of final classes")
        mD=MetaData(self.InSelFile)
        if not mD.containsLabel(MDL_IMAGE):
            errors.append("%s does not contain a column of images"%self.InSelFile)
        if self.Tolerance<0:
            errors.append("Tolerance must be larger than 0")
        return errors
    
    def visualize(self):
        from protlib_utils import runShowJ
        fnSubset=""
        if self.WhatToShow=="Classes":
            fnSubset=""
        elif self.WhatToShow=="Class Cores":
            fnSubset="_core"
        elif self.WhatToShow=="Class Stable Cores":
            fnSubset="_stable_core"
        levelFiles=glob.glob(os.path.join(self.ExtraDir,"level_??/level_classes%s.xmd"%fnSubset))
        if levelFiles:
            levelFiles.sort()
            lastLevelFile=levelFiles[-1]
            if self.DoShowLast:
                runShowJ("classes@"+lastLevelFile)
            else:
                listOfLevels = getListFromRangeString(self.LevelsToShow)
                lastLevel=int(re.search('level_(\d\d)',lastLevelFile).group(1))
                if max(listOfLevels)<=lastLevel:
                    files = "";
                    for level in listOfLevels:
                        fn=os.path.join(self.ExtraDir,"level_%02d/level_classes%s.xmd"%(level,fnSubset))
                        if os.path.exists(fn):
                            files+="classes_sorted@"+fn+" "
                    if files!="":
                        runShowJ(files)
        if self.DoShowHierarchy:
            fnHierarchy=self.getFilename("hierarchy",subset=fnSubset)
            if os.path.exists(fnHierarchy):
                from protlib_gui_ext import showTextfileViewer
                showTextfileViewer(fnHierarchy,[fnHierarchy])
                
    def insertCl2dStep(self):
        params= '-i %(InSelFile)s --odir %(ExtraDir)s --oroot level '+\
                ' --nref %(NumberOfReferences)d'+\
                ' --iter %(NumberOfIterations)d'+\
                ' %(AdditionalParameters)s'
        if self.ComparisonMethod=='correlation':
            params+= ' --distance correlation '
        if self.ClusteringMethod=='classical':
            params+= ' --classicalMultiref '
        if self.AdditionalParameters.find("--ref0")==-1:
            params+=' --nref0 %(NumberOfInitialReferences)d'
    
        self.insertRunJobStep("xmipp_classify_CL2D", params % self.ParamsDict, 
                              [self.workingDirPath("extra/images.xmd")])
        self.insertStep('postCl2d', WorkingDir=self.WorkingDir, 
                        NumberOfReferences=self.NumberOfReferences)
    
def postCl2d(log, WorkingDir, NumberOfReferences):
    levelDirs=glob.glob(os.path.join(WorkingDir,"extra/level_??"))
    createLink(log,os.path.join(WorkingDir,"extra/images.xmd"),os.path.join(WorkingDir,"images.xmd"))
    if levelDirs:
        levelDirs.sort()
        lastLevelDir=levelDirs[-1]
        fnLastLevelFile=os.path.join(lastLevelDir,"level_classes.xmd")
        from protlib_xmipp import getMdSize
        size = getMdSize("classes@"+fnLastLevelFile)
        if size == NumberOfReferences:
            createLink(log, fnLastLevelFile, os.path.join(WorkingDir,"classes.xmd"))

def sortClasses(log,ExtraDir,Nproc,suffix):
    if Nproc==1:
        Nproc=2
    for filename in glob.glob(os.path.join(ExtraDir,"level_??/level_classes%s.xmd"%suffix)):
        level=int(re.search('level_(\d\d)',filename).group(1))
        fnRoot=os.path.join(ExtraDir,"level_%02d/level_classes%s_sorted"%(level,suffix))
        params= "-i classes@"+filename+" --oroot "+fnRoot
        runJob(log,"xmipp_image_sort",params,Nproc)
        mD=MetaData(fnRoot+".xmd")
        mD.write("classes_sorted@"+filename,MD_APPEND)
        deleteFile(log,fnRoot+".xmd")

def evaluateClasses(log,WorkingDir,ExtraDir,subset):
    levelFiles=glob.glob(os.path.join(ExtraDir,"level_??/level_classes%s.xmd"%subset))
    levelFiles.sort()
    for filename in levelFiles:
        runJob(log,"xmipp_classify_evaluate_classes","-i "+filename)
        level=int(re.search('level_(\d\d)',filename).group(1))
        if level>0:
            previousFile=os.path.join(ExtraDir,"level_%02d/level_classes%s.xmd"%(level-1,subset))
            if os.path.exists(previousFile):
                fnOut="%s/classes%s_hierarchy.txt"%(ExtraDir,subset)
                args="--i1 %s --i2 %s -o %s"%(previousFile,filename,fnOut)
                if os.path.exists(fnOut):
                    args+=" --append"
                runJob(log,"xmipp_classify_compare_classes",args)

def postEvaluation(log,WorkingDir):
        levelFiles=glob.glob(os.path.join(WorkingDir,"extra/level_??/level_classes_core.xmd"))
        if levelFiles:
            levelFiles.sort()
            lastLevelFile=levelFiles[-1]
            createLink(log,lastLevelFile,os.path.join(WorkingDir,'classes_core.xmd'))
        levelFiles=glob.glob(os.path.join(WorkingDir,"extra/level_??/level_classes_stable_core.xmd"))
        if levelFiles:
            levelFiles.sort()
            lastLevelFile=levelFiles[-1]
            createLink(log,lastLevelFile,os.path.join(WorkingDir,'classes_stable_core.xmd'))
