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
        from protlib_gui_figure import XmippArrayPlotter1D, XmippArrayPlotter2D, XmippArrayPlotter3D
        components=self.DisplayRawDeformation.split()
        dim=len(components)
        if dim>0:
            modeList=[]
            modeNameList=[]
            # Get modes
            MD=MetaData(self.Modesfile)
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
                XmippArrayPlotter1D(self.workingDirPath("deformations.txt"),modeList[0],"Histogram for mode %s"%modeNameList[0],"Deformation value","Number of images")
            elif dim==2:
                XmippArrayPlotter2D(self.workingDirPath("deformations.txt"),modeList[0],modeList[1],"",
                                    modeNameList[0],modeNameList[1])
            elif dim==3:
                XmippArrayPlotter3D(self.workingDirPath("deformations.txt"),modeList[0],modeList[1],modeList[2],"",
                                    modeNameList[0],modeNameList[1],modeNameList[2])
        components=self.DisplayCombinedDeformation.split()
        dim=len(components)
        if dim>0:
            modeList=[]
            fnDeformationsProjected=self.workingDirPath("deformationsProjected.txt")
            if not os.path.exists(fnDeformationsProjected):
                from protlib_gui_ext import showWarning
                showWarning('Warning', "You don't have enough images to compute combined modes",parent=self.master)
            else:
                # Get modes
                for modeComponent in components:
                    mode=int(modeComponent)
                    if mode>self.OutputDim:
                        from protlib_gui_ext import showWarning
                        showWarning('Warning', "You don't have so many combined modes",parent=self.master)
                    else:
                        mode-=1
                        modeList.append(mode)
                if dim==1:
                    XmippArrayPlotter1D(fnDeformationsProjected,modeList[0],"Histogram for combined mode %d"%(modeList[0]+1),
                                        "Deformation value","Number of images")
                elif dim==2:
                    XmippArrayPlotter2D(fnDeformationsProjected,modeList[0],modeList[1],
                                        "",
                                        "Combined mode %d"%(modeList[0]+1),"Combined mode %d"%(modeList[1]+1))
                elif dim==3:
                    XmippArrayPlotter3D(fnDeformationsProjected,modeList[0],modeList[1],modeList[2],
                                        "",
                                        "Combined mode %d"%(modeList[0]+1),"Combined mode %d"%(modeList[1]+1),
                                        "Combined mode %d"%(modeList[2]+1))
    
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
    Nimgs=MD.size()
    if Nimgs>=10:
        fnDefLow=os.path.join(WorkingDir,"deformationsProjected.txt")
        runJob(log,"xmipp_matrix_dimred","-i %s -o %s --din %d --dout %d --samples %d"%(fnDef,fnDefLow,Xdim,OutputDim,Ydim))
    else:
        print("The number of images ("+str(Nimgs)+") is smaller than 10, combined modes are not calculated")
