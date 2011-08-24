#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using hierarchical clustering principles
# Author: Carlos Oscar Sanchez Sorzano, August 2011
#

import os,sys,shutil,time
from protlib_base import *
from config_protocols import protDict

class ProtCL2D(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.cl2d.name, scriptname, project)
        self.Import = 'from protocol_cl2d import *'    

    def defineSteps(self):
        self.Db.insertStep('cl2d',verifyfiles=[],Selfile=self.InSelFile,WorkingDir=self.WorkingDir,
                           NumberOfReferences=self.NumberOfReferences,NumberOfInitialReferences=self.NumberOfInitialReferences,
                           NumberOfIterations=self.NumberOfIterations,ComparisonMethod=self.ComparisonMethod,
                           ClusteringMethod=self.ClusteringMethod,AdditionalParameters=self.AdditionalParameters,
                           Nproc=self.NumberOfMpi)
        if self.NumberOfReferences>self.NumberOfReferences0:
            self.Db.insertStep('core_analysis',verifyfiles=[],WorkingDir=self.WorkingDir,
                               thGoodClass=self.thGoodClass,thZscore=self.thZscore,thPCAZscore=self.thPCAZscore,
                               Nproc=self.NumberOfMpi)
    
    def summary(self):
        pass
    
    def validate(self):
        errors = []
        if self.NumberOfReferences0>self.NumberOfReferences:
            errors.append("The number of initial classes cannot be larger than the number of final classes")
        if self.thGoodClass<0 or self.thGoodClass>100:
            errors.append("The good class threshold must be between 0 and 100")
        return errors
    
    def visualize(self):
        pass
    
def cl2d(log,Selfile,WorkingDir,NumberOfReferences,NumberOfInitialReferences,NumberOfIterations,
         ComparisonMethod,ClusteringMethod,AdditionalParameters,Nproc):
    params= '-i '+str(selFileToUse)+' -o '+WorkingDir+'/class '+\
            ' --codes '+str(NumberOfReferences)+\
            ' --codes0 '+str(NumberOfInitialReferences)+\
            ' --iter '+str(NumberOfIterations)+\
            ' '+AdditionalParameters
    if ComparisonMethod=='correlation':
        params+= ' --useCorrelation '
    if ClusteringMethod=='classical':
        params+= ' --classicalMultiref '

    runJob(log,"xmipp_classify_CL2D",params,Nproc)

def core_analysis(log,WorkingDir,thGoodClass,thZscore,thPCAZscore,Nproc):
    return
    params= WorkingDir+'/class '+\
            str(thGoodClass)+' '+\
            str(thZscore)+' '+\
            str(thPCAZscore)
    runJob(log,"xmipp_classify_CL2D_core_analysis",params,Nproc)
    