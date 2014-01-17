#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for exploring classes with a volume
#
# Example use:
# ./xmipp_protocol_rct.py
#
# Author: Carlos Oscar Sorzano, March 2013 
#

from os.path import join, exists
from xmipp import MetaData, MetaDataInfo, MDL_IMAGE, MDL_IMAGE1, MDL_IMAGE_REF, MDL_ANGLE_ROT, MDL_ANGLE_TILT, MDL_ANGLE_PSI, MDL_REF, \
        MDL_SHIFT_X, MDL_SHIFT_Y, MDL_FLIP, MD_APPEND, MDL_MAXCC, MDL_ENABLED, Euler_angles2matrix, Image, FileName, getBlocksInMetaDataFile, \
        getImageSize
from protlib_base import *
from protlib_utils import getListFromRangeString, runJob, runShowJ
from protlib_filesystem import copyFile, deleteFile, removeFilenamePrefix
from protlib_projmatch import projMatch,produceAlignedImages

class ProtScreenClasses(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.screen_classes.name, scriptname, project)
        self.Import = 'from protocol_screen_classes import *'
        self.trueClass=False
        if self.Classes.find('@')==-1:
            blocks=getBlocksInMetaDataFile(self.Classes)
            if 'classes' in blocks:
                self.fnImages='classes@'+self.Classes
                self.trueClass=True
            else:
                self.fnImages=self.Classes
        else:
            self.fnImages=self.Classes
        self.Xdim= MetaDataInfo(self.fnImages)[0]
        
    def defineSteps(self):
        fnOutputClass=self.workingDirPath('classes.xmd')
        self.insertStep('createDir',path=self.ExtraDir)
        self.insertStep("linkAcquisitionInfo",InputFile=self.Classes,dirDest=self.WorkingDir)
        if self.trueClass:
            self.insertStep("copyFile",source=removeFilenamePrefix(self.Classes),dest=fnOutputClass)
        else:
            self.insertRunJobStep("xmipp_metadata_utilities","-i %s -o classes@%s"%(self.Classes,fnOutputClass),NumberOfMpi=1,NumberOfThreads=1)

        # Projection Matching        
        fnAngles=self.extraPath('angles.xmd')
        self.insertStep("projMatch",Volume=self.Volume,AngularSampling=self.AngularSampling,SymmetryGroup=self.SymmetryGroup,
                        Images="classes@%s"%fnOutputClass,ExtraDir=self.ExtraDir,fnAngles=fnAngles,NumberOfMpi=self.NumberOfMpi)        

        # Reorganize output and produce difference images 
        self.insertStep("runJob",programname="xmipp_metadata_utilities", params="-i classes@%s --set join %s --mode append"%(fnOutputClass,fnAngles),NumberOfMpi=1)  
        self.insertStep("produceAlignedImages",fnIn='classes@'+fnOutputClass, fnOut='classes_aligned@'+fnOutputClass, fnDiff=self.extraPath("diff.stk"),
                        volumeIsCTFCorrected=False)
        self.insertStep("runJob",programname="xmipp_metadata_utilities", params="-i classes_aligned@%s --operate sort maxCC desc --mode append"%(fnOutputClass),NumberOfMpi=1)  
   
    def summary(self):
        message = []
        message.append("Set of classes: [%s] " % self.Classes)
        message.append("Volume: [%s] " % self.Volume)
        message.append("Symmetry: %s " % self.SymmetryGroup)
        return message
    
    def validate(self):
        errors = []
        (Xdim,Ydim,Zdim,Ndim)=getImageSize(self.Volume)
        if Xdim!=self.Xdim:
            errors.append("Make sure that the volume and the classes have the same size")
        return errors    

    def visualize(self):
        fnAligned = 'classes_aligned@' + self.workingDirPath('classes.xmd')
        runShowJ(fnAligned, extraParams="--mode metadata --render first")
        os.system('xmipp_chimera_client -i '+self.Volume+' --mode projector 256 &')
