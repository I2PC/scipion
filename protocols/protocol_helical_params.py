#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using hierarchical clustering principles
# Author: Carlos Oscar Sanchez Sorzano, August 2011
#

import glob,os,re,sys,shutil,time
from protlib_base import *
from config_protocols import protDict
from protlib_xmipp import getMdSize
from protlib_utils import runJob, getListFromRangeString
from protlib_filesystem import createLink, deleteFile, linkAcquisitionInfo, getXmippPath
from xmipp import MetaData, getImageSize, MDL_ANGLE_ROT, MDL_SHIFT_Z

class ProtHelicalParams(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.helical_params.name, scriptname, project)
        self.Import = 'from protocol_helical_params import *'
        self.fnSym=self.workingDirPath("volume_symmetrized.vol")
        self.fnCoarse=os.path.join(self.WorkingDir,'extra','coarseParams.xmd')
        self.fnFine=os.path.join(self.WorkingDir,'extra','fineParams.xmd')

    def defineSteps(self):
        self.insertStep('createDir',path=self.ExtraDir)
        self.insertStep('coarseSearch',verifyfiles=[self.fnCoarse],
                        WorkingDir=self.WorkingDir,InputVolume=self.InputVolume,CylinderRadius=self.CylinderRadius,
                        Rot0=self.Rot0, RotF=self.RotF, RotStep=self.RotStep, Z0=self.Z0, ZF=self.ZF, ZStep=self.ZStep,
                        NumberOfThreads=self.NumberOfThreads, fnCoarse=self.fnCoarse)
        self.insertStep('fineSearch',verifyfiles=[self.fnFine],WorkingDir=self.WorkingDir,InputVolume=self.InputVolume,
                        CylinderRadius=self.CylinderRadius,fnCoarse=self.fnCoarse,fnFine=self.fnFine)
        self.insertStep('symmetrize',verifyfiles=[self.fnSym],WorkingDir=self.WorkingDir,InputVolume=self.InputVolume,
                        OutputVolume=self.fnSym,CylinderRadius=self.CylinderRadius, fnFine=self.fnFine)
        if self.Dihedrical:
            self.insertStep('findDihedrical',WorkingDir=self.WorkingDir,InputVolume=self.fnSym,CylinderRadius=self.CylinderRadius)
    
    def summary(self):
        message=[]
        message.append("Input volume: [%s] "%self.InputVolume)
        if os.path.exists(self.fnSym):
            message.append("Symmetrized volume: [%s] "%self.fnSym)
        if os.path.exists(self.fnFine):
            md=MetaData(self.fnFine)
            id=md.firstObject()
            rot0=md.getValue(MDL_ANGLE_ROT,id)
            z0=md.getValue(MDL_SHIFT_Z,id)
            message.append("DeltaRot=%f"%rot0)
            message.append("DeltaZ=%f"%z0)
        return message
    
    def validate(self):
        errors = []
        return errors
    
    def visualize(self):
        from protlib_utils import runShowJ
        if os.path.exists(self.fnSym):
            runShowJ(self.fnSym)

def coarseSearch(log,WorkingDir,InputVolume,CylinderRadius,Rot0, RotF, RotStep, Z0, ZF, ZStep, NumberOfThreads,fnCoarse):
    args="-i %s --sym helical -z %f %f %f --rotHelical %f %f %f --thr %d -o %s"%(InputVolume,float(Z0),float(ZF),float(ZStep),
                                                                  float(Rot0),float(RotF),float(RotStep),NumberOfThreads,fnCoarse)
    if CylinderRadius>0:
        [xdim,ydim,zdim,ndim]=getImageSize(InputVolume)
        args+=" --mask cylinder %d %d"%(int(-CylinderRadius),int(-xdim))
    runJob(log,'xmipp_volume_find_symmetry',args)

def fineSearch(log,WorkingDir,InputVolume,CylinderRadius,fnCoarse,fnFine):
    md=MetaData(fnCoarse)
    id=md.firstObject()
    rot0=md.getValue(MDL_ANGLE_ROT,id)
    z0=md.getValue(MDL_SHIFT_Z,id)
    args="-i %s --sym helical --localHelical %f %f -o %s"%(InputVolume,z0,rot0,fnFine)
    if CylinderRadius>0:
        [xdim,ydim,zdim,ndim]=getImageSize(InputVolume)
        args+=" --mask cylinder %d %d"%(int(-CylinderRadius),int(-xdim))
    runJob(log,'xmipp_volume_find_symmetry',args)

def symmetrize(log,WorkingDir,InputVolume,OutputVolume,CylinderRadius,fnFine):
    md=MetaData(fnFine)
    id=md.firstObject()
    rot0=md.getValue(MDL_ANGLE_ROT,id)
    z0=md.getValue(MDL_SHIFT_Z,id)
    args="-i %s --sym helical --helixParams %f %f -o %s"%(InputVolume,z0,rot0,OutputVolume)
    runJob(log,'xmipp_transform_symmetrize',args)
    if CylinderRadius>0:
        [xdim,ydim,zdim,ndim]=getImageSize(InputVolume)
        args="-i %s --mask cylinder %d %d"%(OutputVolume,int(-CylinderRadius),int(-xdim))
        runJob(log,'xmipp_transform_mask',args)

def findDihedrical(log,WorkingDir,InputVolume,CylinderRadius):
    fnRotated=os.path.join(WorkingDir,'tmp','volume_rotated.vol')
    runJob(log,"xmipp_transform_geometry","-i %s -o %s --rotate_volume axis 180 1 0 0"%(InputVolume,fnRotated))
    if CylinderRadius>0:
        [xdim,ydim,zdim,ndim]=getImageSize(InputVolume)
        maskArgs=" --mask cylinder %d %d"%(int(-CylinderRadius),int(-xdim))
    else:
        maskArgs=""
    runJob(log,"xmipp_volume_align","--i1 %s --i2 %s --rot 0 360 3 -z -2 2 0.5 --apply"%(InputVolume,fnRotated)+maskArgs)
    runJob(log,"xmipp_volume_align","--i1 %s --i2 %s --local --apply"%(InputVolume,fnRotated)+maskArgs)
    runJob(log,"xmipp_image_operate","-i %s --plus %s -o %s"%(InputVolume,fnRotated,InputVolume))
    runJob(log,"xmipp_image_operate","-i %s --divide 2"%(InputVolume))
    runJob(log,"xmipp_transform_mask","-i %s "%(InputVolume)+maskArgs)
    deleteFile(log,fnRotated)
