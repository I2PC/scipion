#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based classification with KerDenSOM
# Author: Carlos Oscar Sanchez Sorzano, September 2011
#

import glob,os,sys,shutil,time
from protlib_base import *
from config_protocols import protDict
from protlib_utils import runJob
from protlib_filesystem import deleteFiles
from protlib_sql import SqliteDb
from xmipp import MetaData, MDL_IMAGE, MD_APPEND

class ProtKerdensom(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.kerdensom.name, scriptname, project)
        self.Import = 'from protocol_kerdensom import *'    

    def defineSteps(self):
        self.Db.insertStep('img2vector',[self.workingDirPath("vectors.xmd")],
                           Selfile=self.InSelFile,Mask=self.Mask,WorkingDir=self.WorkingDir)
        self.Db.insertStep('kerdensom',[self.workingDirPath("results_vectors.xmd"),
                                        self.workingDirPath("results_classes.xmd"),
                                        self.workingDirPath("results_images.xmd")],
                           WorkingDir=self.WorkingDir,SomXdim=self.SomXdim,SomYdim=self.SomYdim,
                           SomReg0=self.SomReg0,SomReg1=self.SomReg1,SomSteps=self.SomSteps,
                           KerdensomExtraCommand=self.KerdensomExtraCommand)
        self.Db.insertStep('vector2img',[self.workingDirPath("results_classes.stk")],Mask=self.Mask,WorkingDir=self.WorkingDir)
        self.Db.insertStep('rewriteClassBlock',WorkingDir=self.WorkingDir)

    def summary(self):
        message=[]
        message.append("Classification of "+self.InSelFile+" into a map of size "+str(self.SomYdim)+"x"+str(self.SomXdim))
        if self.getRunState()==SqliteDb.RUN_STARTED:
            lines=[]
            for line in open(self.LogPrefix+".err").readlines():
                if "Training Deterministic Annealing" in line:
                    lines.append(line)
            message.append("Currently at iteration "+str(len(lines))+" out of "+str(self.SomSteps))
        return message
    
    def validate(self):
        errors = []
        if self.SomReg0<self.SomReg1:
            errors.append("Regularization must decrease over iterations: Initial regularization must be larger than final")
        return errors
    
    def visualize(self):
        if self.getRunState()==SqliteDb.RUN_FINISHED:
            os.system("xmipp_showj -i "+self.workingDirPath("results_classes.stk")+" --columns "+str(self.SomXdim)+" &")
            os.system("xmipp_metadata_showj -i "+self.workingDirPath("results_classes.xmd")+" &")
        else:
            tkMessageBox.showwarning("Warning", "The algorithm has not finished yet", parent=self.master)

def img2vector(log,Selfile,Mask,WorkingDir):
     args=' -i '+ Selfile + ' -o ' + os.path.join(WorkingDir,"vectors.xmd")
     if Mask!='':
         args+=' --mask binary_file '+Mask
     runJob(log,"xmipp_image_vectorize", args)

def kerdensom(log,WorkingDir,SomXdim,SomYdim,SomReg0,SomReg1,SomSteps,KerdensomExtraCommand):
    args='-i '+os.path.join(WorkingDir,"vectors.xmd")+\
         ' --oroot '+os.path.join(WorkingDir,"results")+\
         ' --xdim ' + str(SomXdim) + \
         ' --ydim ' + str(SomYdim) + \
         ' --deterministic_annealing %f %f %f'%(SomSteps,SomReg0,SomReg1) + \
         ' '+ str(KerdensomExtraCommand)
    runJob(log,"xmipp_classify_kerdensom",args)
    deleteFiles(log, [os.path.join(WorkingDir,"vectors.xmd"),os.path.join(WorkingDir,"vectors.xmd.raw")], True)
   
def vector2img(log,Mask,WorkingDir):
    args=' -i '+os.path.join(WorkingDir,"results_vectors.xmd")+\
         ' -o '+os.path.join(WorkingDir,"results_classes.stk")
    if Mask!='':
         args+=' --mask binary_file '+Mask
    runJob(log,"xmipp_image_vectorize", args)
    deleteFiles(log, [os.path.join(WorkingDir,"results_vectors.xmd"),os.path.join(WorkingDir,"results_vectors.xmd.raw")], True)

def rewriteClassBlock(log,WorkingDir):
    fnClass="classes@%s"%os.path.join(WorkingDir,"results_classes.xmd")
    fnClassStack=os.path.join(WorkingDir,"results_classes.stk")
    mD=MetaData(fnClass)
    counter=1
    for id in mD:
        mD.setValue(MDL_IMAGE,"%06d@%s"%(counter,fnClassStack),id)
        counter+=1
    mD.write(fnClass,MD_APPEND)