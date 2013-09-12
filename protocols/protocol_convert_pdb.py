#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based pre-processing of volumes 

# Author: Carlos Oscar, August 2013
#
from protlib_base import *
import os
from xmipp import MetaData, MDL_SAMPLINGRATE
from os.path import exists, split, splitext
from protlib_utils import runJob, runShowJ
from protlib_filesystem import copyFile
import glob
from protlib_gui_ext import showError

class ProtConvertPDB(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.convert_pdb.name, scriptname, project)
        self.Import = 'from protocol_convert_pdb import *'
        self.OutModel=self.workingDirPath("volume.vol")

    def defineSteps(self):
        self.insertStep('createAcquisition',WorkingDir=self.WorkingDir,Ts=self.FinalTs)
        self.insertStep("convertFromPDB",InModel=self.InModel,WorkingDir=self.WorkingDir,Ts=self.FinalTs,Size=self.FinalSize)
        
    def validate(self):
        errors = []
        return errors
        
    def summary(self):
        messages = []      
        messages.append("Input model: [%s]" % self.InModel)
        messages.append("Output: [%s]" % self.OutModel)
        return messages

    def visualize(self):
        from protlib_utils import runShowJ
        fnVolume=self.workingDirPath('volume.vol')
        if os.path.exists(fnVolume):
            runShowJ(fnVolume)

def createAcquisition(log,WorkingDir,Ts):
    md=MetaData()
    id=md.addObject()
    md.setValue(MDL_SAMPLINGRATE,float(Ts),id)
    md.write(os.path.join(WorkingDir,"acquisition_info.xmd"))

def convertFromPDB(log,InModel,WorkingDir,Ts,Size):
    args="-i %s -o %s/volume --centerPDB"%(InModel,WorkingDir)
    if Size>0:
        args+=" --size %d"%(int(Size))
    if Ts>4:
        args+=" --poor_Gaussian"
    runJob(log,"xmipp_volume_from_pdb",args)
