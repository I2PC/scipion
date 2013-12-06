#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based classification with KerDenSOM
# Author: Carlos Oscar Sanchez Sorzano, September 2011
#

from protlib_base import *
from protlib_utils import runJob
from xmipp import MetaData, MDL_IMAGE, MD_APPEND
from protlib_filesystem import createLink

class ProtKerdensom(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.kerdensom.name, scriptname, project)
        self.Import = 'from protocol_kerdensom import *'    

    def defineSteps(self):
        self.insertStep('createDir',path=self.ExtraDir)
        self.Db.insertStep("linkAcquisitionInfo",InputFile=self.InSelFile,dirDest=self.WorkingDir)
        self.Db.insertStep('img2vector',[self.extraPath("vectors.xmd")],
                           Selfile=self.InSelFile,Mask=self.Mask,ExtraDir=self.ExtraDir)
        self.Db.insertStep('kerdensom',[self.extraPath("kerdensom_vectors.xmd"),
                                        self.extraPath("kerdensom_classes.xmd"),
                                        self.extraPath("kerdensom_images.xmd")],
                           ExtraDir=self.ExtraDir,SomXdim=self.SomXdim,SomYdim=self.SomYdim,
                           SomReg0=self.SomReg0,SomReg1=self.SomReg1,SomSteps=self.SomSteps,
                           KerdensomExtraCommand=self.KerdensomExtraCommand)
        self.Db.insertStep('vector2img',[self.extraPath("classes.stk")],Mask=self.Mask,ExtraDir=self.ExtraDir)
        self.Db.insertStep('rewriteClassBlock',WorkingDir=self.WorkingDir,ExtraDir=self.ExtraDir)

    def summary(self):
        message = ["KerdenSOM classification"]
        message.append("  Input classes: [%s]" % self.InSelFile)
        message.append("  Map size: <%(SomYdim)d> x <%(SomXdim)d>" % self.ParamsDict)

        if self.getRunState()==SqliteDb.RUN_STARTED:
            lines=[]
            for line in open(self.LogPrefix+".err").readlines():
                if "Training Deterministic Annealing" in line:
                    lines.append(line)
            message.append("Currently at iteration "+str(len(lines))+" out of "+str(self.SomSteps))
        return message
    
    def papers(self):
        papers=[]
        papers.append('Pascual-Montano, JSB (2001) [http://www.ncbi.nlm.nih.gov/pubmed/11472094]')
        papers.append('Pascual-Montano, JSB (2002) [http://www.ncbi.nlm.nih.gov/pubmed/12160707]')
        return papers

    def validate(self):
        errors = []
        if self.SomReg0 < self.SomReg1:
            errors.append("Regularization must decrease over iterations:")
            errors.append("    Initial regularization must be larger than final")
        return errors
    
    def visualize(self):
        if self.getRunState() == SqliteDb.RUN_FINISHED:
            from protlib_utils import runShowJ
            runShowJ("classes@"+self.workingDirPath("classes.xmd"), extraParams=" --columns %d" % self.SomXdim)
        else:
            from protlib_gui_ext import showWarning
            showWarning("Warning", "The algorithm has not finished yet", parent=self.master)

def img2vector(log,Selfile,Mask,ExtraDir):
    args=' -i '+ Selfile + ' -o ' + os.path.join(ExtraDir,"vectors.xmd")
    if Mask!='':
        args+=' --mask binary_file '+Mask
    runJob(log,"xmipp_image_vectorize", args)

def kerdensom(log,ExtraDir,SomXdim,SomYdim,SomReg0,SomReg1,SomSteps,KerdensomExtraCommand):
    args='-i ' + os.path.join(ExtraDir,"vectors.xmd")+\
         ' --oroot ' + os.path.join(ExtraDir,"kerdensom")+\
         ' --xdim ' + str(SomXdim) + \
         ' --ydim ' + str(SomYdim) + \
         ' --deterministic_annealing %f %f %f'%(SomSteps,SomReg0,SomReg1) + \
         ' '+ str(KerdensomExtraCommand)
    runJob(log,"xmipp_classify_kerdensom",args)
    deleteFiles(log, [os.path.join(ExtraDir,"vectors.xmd"),os.path.join(ExtraDir,"vectors.vec")], True)
   
def vector2img(log, Mask, ExtraDir):
    args=' -i ' + os.path.join(ExtraDir,"kerdensom_vectors.xmd")+\
         ' -o ' + os.path.join(ExtraDir,"classes.stk")
    if Mask != '':
        args += ' --mask binary_file ' + Mask
    runJob(log,"xmipp_image_vectorize", args)
    deleteFiles(log, [os.path.join(ExtraDir,"kerdensom_vectors.xmd"),os.path.join(ExtraDir,"kerdensom_vectors.vec")], True)

def rewriteClassBlock(log,WorkingDir,ExtraDir):
    fnClassMetadata=os.path.join(ExtraDir,"kerdensom_classes.xmd")
    fnClass="classes@%s"%fnClassMetadata
    fnClassStack=os.path.join(ExtraDir,"classes.stk")
    mD = MetaData(fnClass)
    counter = 1
    for id in mD:
        mD.setValue(MDL_IMAGE,"%06d@%s"%(counter,fnClassStack),id)
        counter += 1
    mD.write(fnClass,MD_APPEND)
    createLink(log,fnClassMetadata,os.path.join(WorkingDir,"classes.xmd"))
    createLink(log,os.path.join(ExtraDir,"kerdensom_images.xmd"),os.path.join(WorkingDir,"images.xmd"))
