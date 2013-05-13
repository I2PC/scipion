#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for creating an initial volume using a dimensionality reduction method and RANSAC
#
# Example use:
# ./xmipp_protocol_screen_classes.py
#
# Author: Javier Vargas and Carlos Oscar May 2013 
#

from xmipp import MetaData
from protlib_base import *

from os.path import join, exists
from xmipp import MetaData, MetaDataInfo, MDL_IMAGE, MDL_IMAGE_REF, MDL_ANGLE_ROT, MDL_ANGLE_TILT, MDL_ANGLE_PSI, MDL_REF, MDL_SHIFT_X, MDL_SHIFT_Y, \
        MDL_FLIP, MD_APPEND, MDL_MAXCC, MDL_ENABLED, Euler_angles2matrix

from protlib_base import *
from protlib_utils import getListFromRangeString, runJob, runShowJ
from protlib_filesystem import copyFile, deleteFile, removeFilenamePrefix

class ProtInitVolRANSAC(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.initvolume_ransac.name, scriptname, project)
        self.Import = 'from protocol_initvolume_ransac import *'
        self.numInliers = 0
        if self.Classes.find('@')==-1:
            self.fnImages='classes@'+self.Classes
        else:
            self.fnImages=self.Classes        
        self.Xdim= MetaDataInfo(self.fnImages)[0]
        
    def defineSteps(self):
        fnOutputClass=self.workingDirPath(self.Classes)
        print fnOutputClass
        self.insertStep("linkAcquisitionInfo",InputFile=self.Classes,dirDest=self.WorkingDir)
        self.insertStep('copyFile',source=removeFilenamePrefix(self.Classes),dest=fnOutputClass)
        #fnRootMaxInliers = ""
        #maxInliers = 0
        for n in range(self.NRansac):

            fnRoot = "ransac%05d"%n

            #self.insertStep("runJob",programname="xmipp_metadata_utilities", params="-i classes@%s --set join %s --mode append"%(fnOutputClass,fnAngles),NumberOfMpi=1)
            print fnOutputClass
            print fnRoot
            parent_id = self.insertParallelRunJobStep(programname="xmipp_transform_dimred", 
                                                params="-i %s -o %s.xmd -m LTSA -d 2"%(fnOutputClass,fnRoot),
                                                verifyfiles = [fnRoot+".xmd"],
                                                NumberOfMpi = 1)            
            #self.insertStep('runJob', 
            #         programname="xmipp_transform_dimred", 
            #         params="-i %s -o %s.xmd -m LTSA -d 2"%(fnOutputClass,fnRoot),
            #         verifyfiles = [fnRoot+".xmd"],
            #         NumberOfMpi = 1)
            parent_id=self.insertReconstruction(fnRoot, parent_id)
            numInliers = self.insertValidateVol(fnRoot, parent_id)
            
            if (numInliers > maxInliers) :
                fnRootMaxInliers = fnRoot
            
        for it in range(5):    
            self.insertProjMatch(fnRoot);
            self.insertReconstruction(fnRoot);
            

    def insertReconstruction(self,fnRoot, parent_id):
        id=self.insertParallelRunJobStep("xmipp_reconstruct_fourier","-i %s.xmd -o %s.vol"%(fnRoot,fnRoot),
                                      parent_step_id=parent_id)
        id=self.insertParallelRunJobStep("xmipp_transform_mask","-i %s.vol --mask circular -%d "%(fnRoot,self.Xdim/2),
                                         parent_step_id=id)
        return id

    def insertProjMatch(self,fnRoot):
        # Generate gallery of projections        
        fnGallery=self.workingDirPath('gallery_'+fnRoot+'.stk')
        self.insertRunJobStep("xmipp_angular_project_library", "-i %s.vol -o %s --sampling_rate %f --sym %s --method fourier 1 0.25 bspline --compute_neighbors --angular_distance -1 --experimental_images %s"\
                              %(fnRoot,fnGallery,float(self.AngularSampling),self.SymmetryGroup,self.fnImages),[fnGallery])

        # Assign angles
        fnAngles=self.workingDirPath('angles_'+fnRoot+'.xmd')        
        self.insertRunJobStep("xmipp_angular_projection_matching", "-i %s -o %s --ref %s --Ri 0 --Ro %s --max_shift 10 --search5d_shift %s --search5d_step  1 --append"\
                              %(fnRoot+".xmd",fnAngles,fnGallery,str(self.Xdim/2),str(self.Xdim/20)),[fnAngles])  
        self.insertStep("deleteFile",filename=self.workingDirPath('gallery_sampling.xmd') )
        self.insertStep("deleteFile",filename=self.workingDirPath('gallery_angles.doc') )
        self.insertStep("deleteFile",filename=self.workingDirPath('gallery.doc') )
        
    def insertValidateVol(self,fnRoot):
        
        fnGallery=self.workingDirPath('gallery_'+fnRoot+'.stk')
        self.insertRunJobStep("xmipp_angular_project_library", "-i %s.vol -o %s --sampling_rate %f --sym %s --method fourier 1 0.25 bspline --compute_neighbors --angular_distance -1 --experimental_images %s"\
                              %(fnRoot,fnGallery,float(self.AngularSampling),self.SymmetryGroup,self.fnImages),[fnGallery])
        # Assign angles
        fnAngles=self.workingDirPath('angles_'+fnRoot+'.xmd')      
        self.insertRunJobStep("xmipp_angular_projection_matching", "-i %s -o %s --ref %s --Ri 0 --Ro %s --max_shift 10 --search5d_shift %s --search5d_step  1 --append"\
                              %(self.fnImages,fnAngles,fnGallery,str(self.Xdim/2),str(self.Xdim/20)),[fnAngles])  
        self.insertStep("deleteFile",filename=self.workingDirPath('gallery_sampling.xmd') )
        self.insertStep("deleteFile",filename=self.workingDirPath('gallery_angles.doc') )
        self.insertStep("deleteFile",filename=self.workingDirPath('gallery.doc') )
        
        numInliers = self.validateVolumen()            
        fnInliers = 'inliers_'+fnRoot+'.stk'
        
        return numInliers

        
    def validateVolumen(self):
        
        md = MetaData('angles.xmd')
        numInliers = 0        
        for objId in md:
            corr = md.getValue(MDL_MAXCC, objId)
        
            if (corr >= CorrThresh) :
                numInliers = numInliers+1
            
        return numInliers                            
   
    def summary(self):
        message = []
        message.append("Set of classes: [%s] " % self.Classes)
        message.append("Symmetry: %s " % self.SymmetryGroup)
        return message
    
    def validate(self):
        errors = []
        return errors    

    def visualize(self):
        fnAligned='classes_aligned@'+self.workingDirPath('classes.xmd')
        runShowJ(fnAligned)
