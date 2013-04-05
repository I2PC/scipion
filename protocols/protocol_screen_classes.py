#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for creaing a symmetric initial volume
#
# Example use:
# ./xmipp_protocol_rct.py
#
# Author: Carlos Oscar Sorzano, March 2013 
#

from os.path import join, exists
from xmipp import MetaData, MetaDataInfo, MDL_IMAGE, MDL_IMAGE_REF, MDL_ANGLE_ROT, MDL_ANGLE_TILT, MDL_ANGLE_PSI, MDL_REF, MDL_SHIFT_X, MDL_SHIFT_Y, \
        MDL_FLIP, MD_APPEND, MDL_MAXCC, MDL_ENABLED, Euler_angles2matrix

from protlib_base import *
from protlib_utils import getListFromRangeString, runJob, runShowJ
from protlib_filesystem import copyFile, deleteFile, removeFilenamePrefix

class ProtScreenClasses(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.screen_classes.name, scriptname, project)
        self.Import = 'from protocol_screen_classes import *'
        if self.Classes.find('@')==-1:
            self.fnImages='classes@'+self.Classes
        else:
            self.fnImages=self.Classes
        self.Xdim= MetaDataInfo(self.fnImages)[0]
        
    def defineSteps(self):
        fnOutputClass=self.workingDirPath('classes.xmd')
        self.insertStep("linkAcquisitionInfo",InputFile=self.Classes,dirDest=self.WorkingDir)
        self.insertStep('copyFile',source=removeFilenamePrefix(self.Classes),dest=fnOutputClass)

        # Generate gallery of projections        
        fnGallery=self.workingDirPath('gallery.stk')
        self.insertRunJobStep("xmipp_angular_project_library", "-i %s -o %s --sampling_rate %f --sym %s --method fourier 1 0.25 bspline --compute_neighbors --angular_distance -1 --experimental_images %s"\
                              %(self.Volume,fnGallery,float(self.AngularSampling),self.SymmetryGroup,self.fnImages),[fnGallery])

        # Assign angles
        fnAngles=self.workingDirPath('angles.xmd')        
        self.insertRunJobStep("xmipp_angular_projection_matching", "-i %s -o %s --ref %s --Ri 0 --Ro %s --max_shift 1000 --search5d_shift %s --search5d_step  1 --append"\
                              %(self.fnImages,fnAngles,fnGallery,str(self.Xdim/2),str(self.Xdim/10)),[fnAngles])  
        self.insertStep("deleteFile",filename=self.workingDirPath('gallery_sampling.xmd') )
        self.insertStep("deleteFile",filename=self.workingDirPath('gallery_angles.doc') )
        self.insertStep("deleteFile",filename=self.workingDirPath('gallery.doc') )

        # Write angles in the original file and sort
        self.insertStep("substituteReferenceImages",fnAngles=fnAngles,fnGallery=fnGallery)
        self.insertStep("runJob",programname="xmipp_metadata_utilities", params="-i classes@%s --set join %s --mode append"%(fnOutputClass,fnAngles),NumberOfMpi=1)  
        self.insertStep("produceAlignedImages",fnOutputClass=fnOutputClass)
        self.insertStep("runJob",programname="xmipp_metadata_utilities", params="-i classes_aligned@%s --operate sort maxCC --mode append"%(fnOutputClass),NumberOfMpi=1)  
        self.insertStep("deleteFile",filename=fnAngles)
   
    def summary(self):
        message = []
        message.append("Set of classes: [%s] " % self.Classes)
        message.append("Volume: [%s] " % self.Volume)
        message.append("Symmetry: %s " % self.SymmetryGroup)
        return message
    
    def validate(self):
        errors = []
        return errors    

    def visualize(self):
        fnAligned='classes_aligned@'+self.workingDirPath('classes.xmd')
        runShowJ(fnAligned)

def substituteReferenceImages(log,fnAngles,fnGallery):
    MD=MetaData(fnAngles)
    for id in MD:
        galleryReference=MD.getValue(MDL_REF,id)
        MD.setValue(MDL_IMAGE_REF,"%05d@%s"%(galleryReference+1,fnGallery),id)
    MD.write(fnAngles)

def produceAlignedImages(log,fnOutputClass):
    from numpy import array, dot
    MDin=MetaData('classes@'+fnOutputClass)
    MDout=MetaData()
    for i in MDin:
        fnImg=MDin.getValue(MDL_IMAGE,i)
        fnImgRef=MDin.getValue(MDL_IMAGE_REF,i)
        maxCC=MDin.getValue(MDL_MAXCC,i)
        psi =-1.*MDin.getValue(MDL_ANGLE_PSI,i)
        flip = MDin.getValue(MDL_FLIP,i)
        if(flip):
            psi =-psi
        eulerMatrix = Euler_angles2matrix(0.,0.,psi)
        x = MDin.getValue(MDL_SHIFT_X,i)
        y = MDin.getValue(MDL_SHIFT_Y,i)
        shift = array([x, y, 0])
        shiftOut = dot(eulerMatrix, shift)
        [x,y,z]= shiftOut
        if flip:
            x = -x
        id=MDout.addObject()
        MDout.setValue(MDL_IMAGE, fnImg, id)
        MDout.setValue(MDL_IMAGE_REF, fnImgRef, id)
        MDout.setValue(MDL_MAXCC, maxCC, id)
        MDout.setValue(MDL_ANGLE_PSI, psi, id)
        MDout.setValue(MDL_SHIFT_X, x,id)
        MDout.setValue(MDL_SHIFT_Y, y,id)
        MDout.setValue(MDL_FLIP,flip,id)
        MDout.setValue(MDL_ENABLED,1,id)
    MDout.write('classes_aligned@'+fnOutputClass,MD_APPEND)
