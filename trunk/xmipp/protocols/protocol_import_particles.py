#!/usr/bin/env python
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
from protlib_filesystem import deleteFile, createDir, copyFile
from httplib import CREATED

class ProtImportParticles(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.import_particles.name, scriptname, project)
        self.Import = 'from protocol_import_particles import *'

    def defineSteps(self):
        fnOut=os.path.join(self.WorkingDir,"micrographs.sel")
        self.Db.insertStep('createEmptyMicrographSel',verifyfiles=[fnOut],fnOut=fnOut)
        fnOut=os.path.join(self.WorkingDir,self.Family+".sel")
        self.Db.insertStep('linkOrCopy', verifyfiles=[fnOut],
                           Family=self.Family,InputFile=self.InputFile, WorkingDir=self.WorkingDir, DoCopy=self.DoCopy,
                           TmpDir=self.TmpDir)
        if self.DoInvert:
            self.Db.insertStep('invert', FamilySel=fnOut)
            
    def validate(self):
        errors = []
        inputExt=os.path.splitext(self.InputFile)[1]
        if not inputExt in ['.mrc','.stk','.sel','.xmd']:
            error.append("Input file must be stack or metadata")
        else:
            if inputExt in ['.sel','.xmd']:
                mD=xmipp.MetaData(self.InputFile)
                if not mD.containsLabel(xmipp.MDL_IMAGE):
                    error.append("Cannot find label for images in the input file")
        return errors

    def summary(self):
        message=[]
        return message

    def visualize(self):
        pass

def createEmptyMicrographSel(log,fnOut):
    mD=xmipp.MetaData()
    id=mD.addObject()
    mD.setValue(xmipp.MDL_IMAGE,"ImportedImages",id)
    mD.write(fnOut)

def linkOrCopy(log,Family,InputFile,WorkingDir,DoCopy,TmpDir):
    familySel=os.path.join(WorkingDir,Family+".sel")
    if DoCopy:
        familyDir=os.path.join(WorkingDir,Family)
        createDir(log,familyDir)
    if xmipp.FileName.isMetaData(xmipp.FileName(InputFile)):
        mD=xmipp.MetaData(InputFile)
        imageName=mD.getValue(xmipp.MDL_IMAGE,mD.firstObject())
        if imageName[0]!='/':
            # They are relative paths
            relativeDir=os.path.split(os.path.relpath(InputFile,'.'))[0]
            mD.operate("image='%s'||image"%(relativeDir+"/"))
        mD.write(familySel)
        if DoCopy:
            newStack=os.path.join(familyDir,Family+".stk")
            runJob(log,"xmipp_image_convert","-i "+familySel+" -o "+newStack)
            mD=xmipp.MetaData(newStack)
            mD.write(familySel)
    else:
        if DoCopy:
            fileName=os.path.split(InputFile)[1]
            destination=os.path.join(familyDir,fileName)
            copyFile(log,InputFile,destination)
        else:
            destination=InputFile
        mD=xmipp.MetaData(destination)
        mD.write(familySel)

def invert(log,FamilySel):
    runJob()