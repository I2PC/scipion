#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for flexible angular alignment
#
# Author: Carlos Oscar Sanchez Sorzano, May 2013
#         Qiyu Jin
#         Slavica Jonic
#

import glob,os,re,sys,shutil,time
from protlib_base import *
from config_protocols import protDict
from protlib_utils import runJob
from protlib_filesystem import changeDir, createLink, getExt, moveFile, createDir, deleteFile
from xmipp import MetaData, MDL_NMA

class ProtNMAAlignment(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.nma_alignment.name, scriptname, project)
        self.Import = 'from protocol_nma_alignment import *'    

    def defineSteps(self):
        fnImgs=self.workingDirPath("images.xmd")
        self.insertStep("copyFile",verifyfiles=[fnImgs],source=self.InSelFile,dest=fnImgs)
        self.insertStep("performNMA",WorkingDir=self.WorkingDir, InSelFile=fnImgs, 
                        PDBfile=self.PDBfile, Modesfile=self.Modesfile,
                        SamplingRate=self.SamplingRate, TrustRegionScale=self.TrustRegionScale,
                        ProjMatch=self.ProjMatch,
                        MinAngularSampling=self.MinAngularSampling,
                        NProc=self.NumberOfMpi)
        self.insertStep("deleteFile",filename=self.workingDirPath("nmaTodo.xmd"))
        self.insertStep("projectOntoLowerDim",WorkingDir=self.WorkingDir,OutputDim=self.OutputDim)
    
    def summary(self):
        message=[];
        message.append("Input images: ["+self.InSelFile+"]")
        message.append("PDB volume:   ["+self.PDBfile+"]")
        message.append("Normal modes: ["+self.Modesfile+"]")
        return message
    
    def validate(self):
        errors = []
        return errors
    
    def visualize(self):
        pass
    
def performNMA(log, WorkingDir, InSelFile, PDBfile, Modesfile, SamplingRate,
               TrustRegionScale,ProjMatch,MinAngularSampling,NProc):
    arguments="-i "+InSelFile+" --pdb "+PDBfile+" --modes "+Modesfile+" --sampling_rate "+\
              str(SamplingRate) +" --discrAngStep " +str(MinAngularSampling) +" --odir "+WorkingDir+" --centerPDB "+\
              " --trustradius_scale "+str(TrustRegionScale)
    if PDBfile.find("pseudoatoms.pdb")!=-1:
        arguments+=" --fixed_Gaussian"
    if ProjMatch:
        arguments+=" --projMatch"

    runJob(log,"xmipp_nma_alignment",arguments,NProc)

def projectOntoLowerDim(log,WorkingDir,OutputDim):
    MD=MetaData(os.path.join(WorkingDir,"images.xmd"))
    deformations=MD.getColumnValues(MDL_NMA)
    fnDef=os.path.join(WorkingDir,"deformations.txt")
    fhDef=open(fnDef,'w')
    Ydim=len(deformations)
    Xdim=-1
    for deformation in deformations:
        if Xdim==-1:
            Xdim=len(deformation)
        for coef in deformation:
            fhDef.write("%f "%coef)
        fhDef.write("\n")
    fhDef.close()
    fnDefLow=os.path.join(WorkingDir,"deformationsProjected.txt")
    runJob(log,"xmipp_matrix_dimred","-i %s -o %s --din %d --dout %d --samples %d"%(fnDef,fnDefLow,Xdim,OutputDim,Ydim))
