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
from protlib_filesystem import renameFile
from protlib_utils import runJob
from xmipp import MetaData, Image

class ProtDenoiseParticles(XmippProtocol):
    
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.denoise_particles.name, scriptname, project)
        self.Import = 'from protocol_denoise_particles import *'    

    def defineSteps(self):
        self.insertStep('createDir',path=self.ExtraDir)
        self.Db.insertStep("linkAcquisitionInfo",InputFile=self.InStack,dirDest=self.WorkingDir)
        fnRoot=self.extraPath("pca")
        args="-i %s --oroot %s --eigenvectors %d --maxImages %d"%(self.InClasses,fnRoot,self.NumberOfBases,self.MaxNumberOfClasses)
        self.insertStep('runJob',programname="xmipp_image_rotational_pca",params=args,verifyfiles=[fnRoot+".stk"],
                        NumberOfMpi=self.NumberOfMpi,NumberOfThreads=1)
        fnRootDenoised=self.workingDirPath("imagesDenoised")
        N=min(self.NumberOfBases,self.NumberOfBasesDenoising)
        args="-i %s -o %s.stk --save_metadata_stack %s.xmd --basis %s.stk %d"\
             %(self.InStack,fnRootDenoised,fnRootDenoised,fnRoot,N)
        self.insertStep('runJob', programname="xmipp_transform_filter",params=args,verifyfiles=[fnRootDenoised+".stk"],
                        NumberOfMpi=1,NumberOfThreads=1)
    
    def summary(self):
        message=["Input: [%s]"%self.InStack]
        fnRootDenoised=self.workingDirPath("imagesDenoised")
        message.append("Output: [%s.xmd]"%fnRootDenoised)
        return message
    
    def papers(self):
        papers=[]
        papers.append('Ponce, IEEE TIP (2011) [http://www.ncbi.nlm.nih.gov/pubmed/21536533]')
        return papers

    def validate(self):
        errors = []
        return errors
    
    def visualize(self):
        from protlib_utils import runShowJ
        fnOut=self.workingDirPath("imagesDenoised")+".xmd"
        if exists(fnOut):
            runShowJ(fnOut)
