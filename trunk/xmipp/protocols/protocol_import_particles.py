#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based importing single-particles: 
#  - normalization
#  - sort_junk

# Author: Carlos Oscar, August 2011
#
from protlib_base import *
import xmipp
import glob
import os
from protlib_utils import runJob
from protlib_filesystem import deleteFile, createDir, copyFile, fixPath, findProjectInPathTree
from httplib import CREATED

class ProtImportParticles(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.import_particles.name, scriptname, project)
        self.Import = 'from protocol_import_particles import *'

    def defineSteps(self):
        fnOut=self.workingDirPath("micrographs.xmd")
        self.Db.insertStep('createEmptyMicrographSel',verifyfiles=[fnOut],fnOut=fnOut)
        selfileRoot=self.workingDirPath(self.Family)
        fnOut=selfileRoot+".xmd"
        self.Db.insertStep('linkOrCopy', verifyfiles=[fnOut],
                           Family=self.Family,InputFile=self.InputFile, WorkingDir=self.WorkingDir, DoCopy=self.DoCopy,
                           ImportAll=self.ImportAll, SubsetMode=self.SubsetMode, Nsubset=self.Nsubset, TmpDir=self.TmpDir)
        if self.DoInvert and self.ImportAll:
            self.Db.insertStep('invert', FamilySel=fnOut,Nproc=self.NumberOfMpi)
        if self.DoRemoveDust and self.ImportAll:
            self.Db.insertStep('removeDust', FamilySel=fnOut,threshold=self.DustRemovalThreshold,Nproc=self.NumberOfMpi)
        if self.DoNorm and self.ImportAll:
            self.Db.insertStep('normalize', FamilySel=fnOut,bgRadius=self.BackGroundRadius,Nproc=self.NumberOfMpi)
        self.Db.insertStep('sortImageInFamily', verifyfiles=[selfileRoot+".xmd"],selfileRoot=selfileRoot)
            
    def validate(self):
        errors = []
        inputExt=os.path.splitext(self.InputFile)[1]
        if not inputExt in ['.mrc','.stk','.sel','.xmd','.hed','.img', '.ctfdat']:
            errors.append("Input file must be stack or metadata")
        else:
            if inputExt in ['.sel', '.xmd', '.ctfdat']:
                mD=xmipp.MetaData(self.InputFile)
                if not mD.containsLabel(xmipp.MDL_IMAGE):
                    errors.append("Cannot find label for images in the input file")
        return errors

    def summary(self):
        message=[]
        message.append("Stack imported from: "+self.InputFile)
        if self.DoCopy:
            message.append("Copied into "+self.WorkingDir)
        steps=[]
        if self.DoInvert:
            steps.append("Constrast inversion")
        if self.DoRemoveDust:
            steps.append("Dust removal")
        if self.DoNorm:
            steps.append("Ramp normalization")
        if len(steps)>0:
            message.append("Steps applied: "+",".join(steps))
            
        return message

    def visualize(self):
        from protlib_utils import runShowJ
        runShowJ(self.workingDirPath(self.Family+".xmd"), memory='1024m')        

def createEmptyMicrographSel(log,fnOut):
    mD = xmipp.MetaData()
    id = mD.addObject()
    mD.setValue(xmipp.MDL_IMAGE,"ImportedImages",id)
    mD.write(fnOut)

def linkOrCopy(log,Family,InputFile,WorkingDir,DoCopy,ImportAll,SubsetMode,Nsubset,TmpDir):
    familySel = os.path.join(WorkingDir,Family+".xmd")
    if DoCopy:
        familyDir=os.path.join(WorkingDir,Family)
        createDir(log,familyDir)
    if xmipp.FileName.isMetaData(xmipp.FileName(InputFile)):
        mD=xmipp.MetaData(InputFile)
        inputRelativePath=os.path.split(os.path.relpath(InputFile,'.'))[0]
        projectDir=findProjectInPathTree(InputFile)
        for id in mD:
            file=mD.getValue(xmipp.MDL_IMAGE,id)
            if '@' in file:
                (num,file)=file.split('@')
                file=os.path.relpath(fixPath(file,'.',inputRelativePath),'.')
                file='%(num)s@%(file)s'%locals()
            else:
                file=os.path.relpath(fixPath(file,'.',inputRelativePath),'.')
            mD.setValue(xmipp.MDL_IMAGE,file,id)
        if not ImportAll:
            if SubsetMode=="Random subset":
                mDaux=xmipp.MetaData()
                mDaux.randomize(mD)
            else:
                mDaux=xmipp.MetaData(mD)
            mD.selectPart(mDaux,0,Nsubset)
        mD.write(familySel)
        if DoCopy:
            newStack=os.path.join(familyDir,Family+".stk")
            runJob(log,"xmipp_image_convert","-i "+familySel+" -o "+newStack)
            mD=xmipp.MetaData(newStack)
            mD.write(familySel)
    else:
        mD=xmipp.MetaData(InputFile)
        if not ImportAll:
            if SubsetMode=="Random subset":
                mDaux=xmipp.MetaData()
                mDaux.randomize(mD)
            else:
                mDaux=xmipp.MetaData(mD)
            mD.selectPart(mDaux,0,Nsubset)
        mD.write(familySel)
        if DoCopy:
            if ImportAll:
                fileName=os.path.split(InputFile)[1]
                destination=os.path.join(familyDir,fileName)
                copyFile(log,InputFile,destination)
                (inputRoot,inputExt)=os.path.splitext(fileName)
                if inputExt=='hed':
                    copyFile(log,InputFile.replace(".hed",".img"),destination.replace(".hed",".img"))
                elif inputExt=="img":
                    copyFile(log,InputFile.replace(".img",".hed"),destination.replace(".img",".hed"))
            else:
                destination=os.path.join(familyDir,Family+".stk")
                runJob(log,"xmipp_image_convert","-i "+familySel+" -o "+destination)
            mD=xmipp.MetaData(destination)
            mD.write(familySel)

def invert(log,FamilySel,Nproc):
    runJob(log,'xmipp_image_operate','-i '+FamilySel+" --mult -1",Nproc)

def removeDust(log,FamilySel,threshold,Nproc):
    runJob(log,'xmipp_transform_filter','-i '+FamilySel+" --bad_pixels outliers %f"%threshold,Nproc)

def normalize(log,FamilySel,bgRadius,Nproc):
    if bgRadius==0:
        particleSize=xmipp.ImgSize(FamilySel)[0]
        bgRadius=int(particleSize/2)
    runJob(log,"xmipp_transform_normalize","-i "+FamilySel+' --method Ramp --background circle '+str(bgRadius),Nproc)

def sortImageInFamily(log,selfileRoot):
    runJob(log,"xmipp_image_sort_by_statistics","-i "+selfileRoot+".xmd --multivariate --addToInput -o "+selfileRoot+"_sorted")
