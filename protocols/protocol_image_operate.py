#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Author: Carlos Oscar Sanchez Sorzano, September 2013
#

import glob,os,re,sys,shutil,time
from protlib_base import *
from config_protocols import protDict
from protlib_xmipp import getMdSize
from protlib_utils import runJob
from protlib_filesystem import findAcquisitionInfo
from xmipp import MetaData, MetaDataInfo

class ProtImageOperate(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.image_operate.name, scriptname, project)
        self.Import = 'from protocol_image_operate import *'    

    def defineSteps(self):
        acquisitionInfo = findAcquisitionInfo(self.Operand1)
        if acquisitionInfo is not None:
            self.Db.insertStep("linkAcquisitionInfo",InputFile=self.Operand1,dirDest=self.WorkingDir)
        self.Db.insertStep("runImageOperate",WorkingDir=self.WorkingDir,O1=self.Operand1,Operation=self.Operation,O2=self.Operand2,
                           Nproc=self.NumberOfMpi)
    
    def summary(self):
        message=[]
        message.append("Operand 1: [%s]"%self.Operand1)
        message.append("Operation: %s"%self.Operation)
        Operation=self.Operation
        if Operation=='plus' or Operation=='minus' or Operation=='multiply' or Operation=='divide' or Operation=='minimum' or \
           Operation=='maximum' or Operation=='dot product' or Operation=='column' or Operation=='slice' or \
           Operation=='row':
           if not self.Operand2.isdigit():
               message.append("Operand 2: [%s]"%self.Operand2)
           else:
               message.append("Operand 2: %s"%self.Operand2)
        return message
    
    def validate(self):
        errors = []
        N1=getMdSize(self.Operand1)
        Operation=self.Operation
        checkDimensions=False
        if Operation=='column' or Operation=='slice' or Operation=='row':
            if not self.Operand2.isdigit():
                errors.append('You should give a number for the column, slice or row')
        elif Operation=='dot product':
            if self.Operand2.isdigit():
                errors.append('Second operand cannot be a number')
            else:
                checkDimensions=True
        elif Operation=='plus' or Operation=='minus' or Operation=='multiply' or Operation=='divide' or Operation=='minimum' or \
           Operation=='maximum':
            if not self.Operand2.isdigit():
                checkDimensions=True
        if checkDimensions:
            md1=MetaData(self.Operand1)
            md2=MetaData(self.Operand2)
            x1, y1, z1, _, _ = MetaDataInfo(md1)    
            x2, y2, z2, _, _ = MetaDataInfo(md2)
            if x1!=x2 or y1!=y2 or z1!=z2:
                errors.append("Image/Volume sizes in the two operands are not the same")    
            if md2.size()>1:
                if md2.size()!=md1.size():
                    errors.append("The number of images/volumes in the two operands are not the same")
        return errors
    
    def visualize(self):
        from protlib_utils import runShowJ
        fn=self.workingDirPath('results.xmd')
        if os.path.exists(fn):
            runShowJ(fn)

def runImageOperate(log,WorkingDir,O1,Operation,O2,Nproc):
    fnStkOut=os.path.join(WorkingDir,"results.stk")
    fnXmdOut=os.path.join(WorkingDir,"results.xmd")
    args="-i %s"%O1
    if Operation!='dot product':
        args+=" -o %s --save_metadata_stack %s --track_origin"%(fnStkOut,fnXmdOut)
    if Operation=='multiply':
        args+=" --mult"
    elif Operation=="minimum":
        args+=" --min"
    elif Operation=="maximum":
        args+=" --max"
    elif Operation=="dot product":
        args+=" --dot_product"
    elif Operation=="radial average":
        args+=" --radial_avg"
    else:
        args+=" --"+Operation
        
    if Operation=='plus' or Operation=='minus' or Operation=='multiply' or Operation=='divide' or Operation=='minimum' or \
       Operation=='maximum' or Operation=='dot product' or Operation=='column' or Operation=='slice' or \
       Operation=='row':
       args+=" "+O2
    runJob(log,"xmipp_image_operate",args,Nproc)
