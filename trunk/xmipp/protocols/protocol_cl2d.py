#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using hierarchical clustering principles
# Author: Carlos Oscar Sanchez Sorzano, August 2011
#

import glob,os,sys,shutil,time
from protlib_base import *
from config_protocols import protDict
from protlib_utils import runJob,getRangeValuesFromString
from protlib_filesystem import createLink
from xmipp import MetaData
from protocol_particle_pick_auto import createLinkToMicrographs

class ProtCL2D(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.cl2d.name, scriptname, project)
        self.Import = 'from protocol_cl2d import *'    

    def defineSteps(self):
        fnClasses=os.path.join(self.WorkingDir,"results_classes.xmd")
        fnImages=os.path.join(self.WorkingDir,"results_images.xmd")
        self.Db.insertStep('cl2d',verifyfiles=[fnClasses,fnImages],Selfile=self.InSelFile,WorkingDir=self.WorkingDir,
                           NumberOfReferences=self.NumberOfReferences,NumberOfInitialReferences=self.NumberOfInitialReferences,
                           NumberOfIterations=self.NumberOfIterations,ComparisonMethod=self.ComparisonMethod,
                           ClusteringMethod=self.ClusteringMethod,AdditionalParameters=self.AdditionalParameters,
                           Nproc=self.NumberOfMpi)
        if self.NumberOfReferences>self.NumberOfInitialReferences:
            self.Db.insertStep('core_analysis',verifyfiles=[],WorkingDir=self.WorkingDir,
                               thGoodClass=self.thGoodClass,thZscore=self.thZscore,thPCAZscore=self.thPCAZscore,
                               Nproc=self.NumberOfMpi)
    
    def summary(self):
        message=[]
        levelFiles=glob.glob(self.WorkingDir+"/results_level_*.xmd")
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
            mD=MetaData("classes@"+lastLevelFile)
            message.append("Last iteration is %s from level %d with %d classes (at %s)"%(iteration,lastLevel,mD.size(),date))
        return message
    
    def validate(self):
        errors = []
        if self.NumberOfInitialReferences>self.NumberOfReferences:
            errors.append("The number of initial classes cannot be larger than the number of final classes")
        if self.thGoodClass<0 or self.thGoodClass>100:
            errors.append("The good class threshold must be between 0 and 100")
        return errors
    
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
        import re
        levelFiles.sort()
        lastLevelFile=levelFiles[-1]
        mD=xmipp.MetaData("classes@"+lastLevelFile)
        if mD.size()==NumberOfReferences:
            createLink(log, lastLevelFile, WorkingDir+"results_classes.sel")

def core_analysis(log,WorkingDir,thGoodClass,thZscore,thPCAZscore,Nproc):
    return
    params= WorkingDir+'/results '+\
            str(thGoodClass)+' '+\
            str(thZscore)+' '+\
            str(thPCAZscore)
    runJob(log,"xmipp_classify_CL2D_core_analysis",params,Nproc)
    