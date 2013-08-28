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
from xmipp import MetaData, MDL_NMA, MDL_NMA_MODEFILE

class ProtNMAAlignment(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.nma_alignment.name, scriptname, project)
        self.Import = 'from protocol_nma_alignment import *'    

    def defineSteps(self):
        self.insertStep('createDir',path=self.ExtraDir)
        fnImgs=self.workingDirPath("images.xmd")
        self.insertStep("copyFile",verifyfiles=[fnImgs],source=self.InSelFile,dest=fnImgs)
        fnModes=self.workingDirPath("modes.xmd")
        self.insertStep("copyFile",verifyfiles=[fnModes],source=self.Modesfile,dest=fnModes)
        fnAtoms=self.workingDirPath("atoms.pdb")
        self.insertStep("copyFile",verifyfiles=[fnAtoms],source=self.PDBfile,dest=fnAtoms)
        self.insertStep("performNMA",WorkingDir=self.WorkingDir, InSelFile=fnImgs, 
                        PDBfile=self.PDBfile, Modesfile=self.Modesfile,
                        SamplingRate=self.SamplingRate, TrustRegionScale=self.TrustRegionScale,
                        ProjMatch=self.ProjMatch,
                        DiscreteAngularSampling=self.DiscreteAngularSampling,
                        NProc=self.NumberOfMpi)
        self.insertStep("deleteFile",filename=self.workingDirPath("nmaTodo.xmd"))
	self.insertStep("extractDeformations",WorkingDir=self.WorkingDir)
    
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
        from protlib_gui_figure import XmippArrayPlotter1D, XmippArrayPlotter2D, XmippArrayPlotter3D
        components=self.DisplayRawDeformation.split()
        dim=len(components)
        if dim>0:
            modeList=[]
            modeNameList=[]
            # Get modes
            MD=MetaData(self.Modesfile)
            MD.removeDisabled()
            for modeComponent in components:
                mode=int(modeComponent)
                if mode>MD.size():
                    from protlib_gui_ext import showWarning
                    showWarning('Warning', "You don't have so many modes",parent=self.master)
                else:
                    mode-=1
                    currentMode=0;
                    modeName=""
                    for id in MD:
                        modeName=MD.getValue(MDL_NMA_MODEFILE,id)
                        currentMode+=1
                        if currentMode>mode:
                            break
                    modeNameList.append(modeName)
                    modeList.append(mode)
            
            # Actually plot
            if dim==1:
                XmippArrayPlotter1D(self.extraPath("deformations.txt"),modeList[0],"Histogram for mode %s"%modeNameList[0],"Deformation value","Number of images")
            elif dim==2:
                XmippArrayPlotter2D(self.extraPath("deformations.txt"),modeList[0],modeList[1],"",
                                    modeNameList[0],modeNameList[1])
            elif dim==3:
                XmippArrayPlotter3D(self.extraPath("deformations.txt"),modeList[0],modeList[1],modeList[2],"",
                                    modeNameList[0],modeNameList[1],modeNameList[2])
  
    def visualizeVar(self, varName):
        if varName=="AnalyzeMATLAB" and self.AnalyzeMATLAB:
            os.system("matlab -r \"xmipp_nma_selection_tool(\'"+self.WorkingDir+"\')\"")

def performNMA(log, WorkingDir, InSelFile, PDBfile, Modesfile, SamplingRate,
               TrustRegionScale,ProjMatch,DiscreteAngularSampling,NProc):
    arguments="-i "+InSelFile+" --pdb "+PDBfile+" --modes "+Modesfile+" --sampling_rate "+\
              str(SamplingRate) +" --discrAngStep " +str(DiscreteAngularSampling) +" --odir "+WorkingDir+"/tmp --centerPDB "+\
              " --trustradius_scale "+str(TrustRegionScale)
    if PDBfile.find("pseudoatoms.pdb")!=-1:
        arguments+=" --fixed_Gaussian"
    if ProjMatch:
        arguments+=" --projMatch"

    runJob(log,"xmipp_nma_alignment",arguments,NProc)

def extractDeformations(log,WorkingDir):
    MD=MetaData(os.path.join(WorkingDir,"images.xmd"))
    deformations=MD.getColumnValues(MDL_NMA)
    fnDef=os.path.join(WorkingDir,"extra/deformations.txt")
    fhDef=open(fnDef,'w')
    for deformation in deformations:
        for coef in deformation:
            fhDef.write("%f "%coef)
        fhDef.write("\n")
    fhDef.close()

