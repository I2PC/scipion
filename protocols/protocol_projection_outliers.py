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

class ProtProjectionOutliers(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.projection_outliers.name, scriptname, project)
        self.Import = 'from protocol_projection_outliers import *'
        self.Xdim= MetaDataInfo(self.Images)[0]
        
    def defineSteps(self):
        self.insertStep('createDir',path=self.ExtraDir)
        self.insertStep("linkAcquisitionInfo",InputFile=self.Images,dirDest=self.WorkingDir)

        # Projection matching
        fnAngles=self.tmpPath('angles.xmd')
        self.insertStep("projMatch",Volume=self.Volume,AngularSampling=self.AngularSampling,SymmetryGroup=self.SymmetryGroup,
                        Images=self.Images,ExtraDir=self.ExtraDir,fnAngles=fnAngles,NumberOfMpi=self.NumberOfMpi)        

        # Prepare output
        fnOutputImages=self.workingDirPath('images.xmd')
        self.insertRunJobStep("xmipp_metadata_utilities","-i %s --set join %s -o %s"%(fnAngles,self.Images,fnOutputImages),NumberOfMpi=1)
        
        # Produce difference images
        fnDiff=self.extraPath("diff.stk")
        fnAligned=self.extraPath("images_aligned.xmd")
        self.insertStep("produceAlignedImages",fnIn=fnOutputImages, fnOut=fnAligned, fnDiff=self.extraPath("diff.stk"),
                        volumeIsCTFCorrected=self.VolumeIsCTFCorrected)
        
        # Evaluate each image
        fnAutoCorrelations=self.extraPath("autocorrelations.xmd")
        self.insertRunJobStep("xmipp_image_residuals", "-i %s -o %s --save_metadata_stack %s"%
                              (fnDiff,self.extraPath("autocorrelations.stk"),fnAutoCorrelations),
                              [fnAutoCorrelations],NumberOfMpi=1)
        self.insertRunJobStep("xmipp_metadata_utilities", "-i %s --set merge %s"%(fnAligned,fnAutoCorrelations),NumberOfMpi=1)
        self.insertStep("deleteFile",filename=fnAutoCorrelations)

        # Prepare output
        self.insertRunJobStep("xmipp_metadata_utilities","-i %s --set join %s"%(fnOutputImages,fnAligned),NumberOfMpi=1)
        self.insertRunJobStep("xmipp_metadata_utilities","-i %s --operate sort zScoreResCov desc"%fnOutputImages,NumberOfMpi=1)  

    def summary(self):
        message = []
        message.append("Set of images: [%s] " % self.Images)
        message.append("Volume: [%s] " % self.Volume)
        message.append("Symmetry: %s " % self.SymmetryGroup)
        return message
    
    def validate(self):
        errors = []
        (Xdim,Ydim,Zdim,Ndim)=getImageSize(self.Volume)
        if Xdim!=self.Xdim:
            errors.append("Make sure that the volume and the images have the same size")
        return errors    

    def visualize(self):
        fnImages = self.workingDirPath("images.xmd")
        runShowJ(fnImages, extraParams="--mode metadata --render first")
