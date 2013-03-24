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
from xmipp import MetaData, ImgSize, MDL_IMAGE, MDL_ANGLE_ROT, MDL_ANGLE_TILT, MDL_ANGLE_PSI, MDL_COST

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
        (self.Xdim, Ydim, Zdim, Ndim) = ImgSize(self.fnImages)
        
    def defineSteps(self):
        fnOutputClass=self.workingDirPath('classes.xmd')
        self.insertStep("linkAcquisitionInfo",InputFile=self.Classes,dirDest=self.WorkingDir)
        self.insertStep('copyFile',source=removeFilenamePrefix(self.Classes),dest=fnOutputClass)

        # Generate gallery of projections        
        fnGallery=self.workingDirPath('gallery.stk')
        self.insertRunJobStep("xmipp_angular_project_library", "-i %s -o %s --sampling_rate 10 --sym %s --method fourier 1 0.25 bspline --compute_neighbors --angular_distance -1 --experimental_images %s"\
                              %(self.Volume,fnGallery,self.SymmetryGroup,self.fnImages),[fnGallery])

        # Assign angles
        fnAngles=self.workingDirPath('angles.xmd')        
        self.insertRunJobStep("xmipp_angular_projection_matching", "-i %s -o %s --ref %s --Ri 0 --Ro %s --max_shift 1000 --search5d_shift %s --search5d_step  1 --append"\
                              %(self.fnImages,fnAngles,fnGallery,str(self.Xdim/2),str(self.Xdim/10)),[fnAngles])  

        # Write angles in the original file and sort
        self.insertRunJobStep("xmipp_metadata_utilities -i %s --fill ref lineal 0 1"%self.workingDirPath('gallery.doc'))
        self.insertRunJobStep("xmipp_metadata_utilities", "-i classes@%s --set join %s --mode append"%(fnOutputClass,fnAngles))  
        self.insertRunJobStep("xmipp_metadata_utilities", "-i classes@%s --operate sort maxCC --mode append"%(fnOutputClass,fnAngles))  
   
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
        return
