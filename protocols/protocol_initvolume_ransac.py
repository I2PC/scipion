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
from xmipp import MetaData, MetaDataInfo, MD_APPEND, MDL_MAXCC, MDL_WEIGHT

from protlib_base import *
from protlib_utils import getListFromRangeString, runJob, runShowJ
from protlib_filesystem import copyFile, deleteFile, removeFilenamePrefix

class ProtInitVolRANSAC(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.initvolume_ransac.name, scriptname, project)
        self.Import = 'from protocol_initvolume_ransac import *'
        if self.Classes.find('@')==-1:
            self.fnImages='classes@'+self.Classes
        else:
            self.fnImages=self.Classes        
        self.Xdim= MetaDataInfo(self.fnImages)[0]
        
    def defineSteps(self):
        fnOutputClass=self.workingDirPath(self.Classes)

        self.insertStep("linkAcquisitionInfo",InputFile=self.Classes,dirDest=self.WorkingDir)
        self.insertStep('copyFile',source=removeFilenamePrefix(self.Classes),dest=fnOutputClass)
        fnOutputReducedClass = self.workingDirPath("reducedClasses") 
        
        self.MaxFreq = float(self.MaxFreq)
        self.Ts = float(self.Ts)
       
        K = 0.25*(self.MaxFreq/self.Ts)
        Xdim2 = self.Xdim/K
        
        if (Xdim2 < 32):
            self.Xdim2 = 32
            K = self.Xdim/Xdim2
        else:
            self.Xdim2 = Xdim2
        
        self.Ts = K*self.Ts
        freq = self.Ts/self.MaxFreq

        if (self.InitialVolume != '') :
            self.fnOutputInitVolume=self.workingDirPath(self.InitialVolume)
            self.insertStep('copyFile',source=removeFilenamePrefix(self.InitialVolume),dest=self.fnOutputInitVolume)
            self.insertRunJobStep("xmipp_image_resize","-i %s -o %s --dim %d %d" 
                                  %(self.fnOutputInitVolume,self.fnOutputInitVolume,self.Xdim2,self.Xdim2 ))
                        
        self.insertRunJobStep("xmipp_image_resize","-i %s -o %s.xmd --dim %d %d --oroot %s" 
                                                %(fnOutputClass,fnOutputReducedClass,self.Xdim2,self.Xdim2,
                                                  self.workingDirPath("tmp/reduced")))
        
        self.insertRunJobStep("xmipp_transform_filter","-i %s.xmd --fourier low_pass %f" 
                                                %(fnOutputReducedClass,freq))
                
        for n in range(self.NRansac):

            fnRoot=self.workingDirPath("ransac%05d"%n)                        
            parent_id = self.insertParallelRunJobStep("xmipp_transform_dimred",
                                                      "-i %s.xmd -o %s.xmd -m LTSA -d 2 -n %d"%(fnOutputReducedClass,fnRoot,self.NumGrids),
                                                      verifyfiles = [fnRoot+".xmd"],parent_step_id=XmippProjectDb.FIRST_STEP)
            
            if (self.InitialVolume != '') :
                parent_id = self.insertProjMatchInitialVolume(fnRoot,parent_id)
            
            parent_id=self.insertReconstruction(fnRoot, parent_id)
            
            fnRoot = "ransac%05d"%n
            parent_id = self.insertValidateVol(fnRoot, parent_id)
            
        self.insertStep("getBestVolume",WorkingDir=self.WorkingDir, NRansac=self.NRansac, NumVolumes=self.NumVolumes)        
        
        for n in range(self.NumVolumes):
            
            fnRoot=self.workingDirPath('bestAssignment%05d'%n)
            self.insertStep('runJob',
                     programname="xmipp_metadata_utilities", 
                     params="-i %s.xmd -o %s.xmd --query select \"maxCC>%f \"" %(fnRoot,fnRoot,self.CorrThresh),
                     verifyfiles = [fnRoot+".xmd"],NumberOfMpi = 1)
            
            parent_id=XmippProjectDb.FIRST_STEP
            
            for it in range(self.NumIter):    
                parent_id = self.insertReconstruction(fnRoot, parent_id)
                parent_id = self.insertProjMatch(fnRoot,parent_id)
            
            self.insertRunJobStep("xmipp_image_resize","-i %s.vol -o %s.vol --dim %d %d" 
                                              %(fnRoot,fnRoot,self.Xdim,self.Xdim))
        
    def insertReconstruction(self,fnRoot, parent_id):

        id=self.insertParallelRunJobStep("xmipp_reconstruct_fourier","-i %s.xmd -o %s.vol --sym %s " %(fnRoot,fnRoot,self.SymmetryGroup),
                                      parent_step_id=parent_id)
        id=self.insertParallelRunJobStep("xmipp_transform_mask","-i %s.vol --mask circular -%d "%(fnRoot,self.Xdim2/2),
                                         parent_step_id=id)
        return id

    def insertProjMatch(self,fnRoot, parent_id):

        fnGallery=self.workingDirPath('gallery_BestVolume'+'.stk')
        fnOutputReducedClass = self.workingDirPath("reducedClasses") 
        
        id=self.insertParallelRunJobStep("xmipp_angular_project_library", "-i %s.vol -o %s --sampling_rate %f --sym %s --method fourier 1 0.25 bspline --compute_neighbors --angular_distance -1 --experimental_images %s.xmd"\
                              %(fnRoot,fnGallery,float(self.AngularSampling),self.SymmetryGroup,fnOutputReducedClass),[fnGallery],parent_step_id=parent_id)

        #fnAngles=self.workingDirPath('angles_'+fnRoot+'.xmd')
        fnAngles=self.workingDirPath('angles_BestVolume'+'.xmd')
        id=self.insertParallelRunJobStep("xmipp_angular_projection_matching", "-i %s.xmd -o %s --ref %s --Ri 0 --Ro %s --max_shift %s --append"\
                              %(fnRoot,fnAngles,fnGallery,str(self.Xdim/2),str(self.Xdim/20)),[fnAngles],parent_step_id=id)
                
        id = self.insertParallelStep("deleteFile",filename=self.workingDirPath('gallery_sampling.xmd'),parent_step_id=id)
        id = self.insertParallelStep("deleteFile",filename=self.workingDirPath('gallery_angles.doc'),parent_step_id=id)
        id = self.insertParallelStep("deleteFile",filename=self.workingDirPath('gallery.doc'),parent_step_id=id)
        
        return id
    
    def insertProjMatchInitialVolume(self,fnRoot, parent_id):

        fnGallery=self.workingDirPath('gallery_InitialVolume'+'.stk')
        fnOutputReducedClass = self.workingDirPath("reducedClasses") 
        
        id=self.insertParallelRunJobStep("xmipp_angular_project_library", "-i %s -o %s --sampling_rate %f --sym %s --method fourier 1 0.25 bspline --compute_neighbors --angular_distance -1 --experimental_images %s.xmd"\
                              %(self.fnOutputInitVolume,fnGallery,float(self.AngularSampling),self.SymmetryGroup,fnOutputReducedClass),[fnGallery],parent_step_id=parent_id)

        id=self.insertParallelRunJobStep("xmipp_angular_projection_matching", "-i %s.xmd -o %s.xmd --ref %s --Ri 0 --Ro %s --max_shift %s --append"\
                              %(fnRoot,fnRoot,fnGallery,str(self.Xdim/2),str(self.Xdim/20)),parent_step_id=id)
                
        id = self.insertParallelStep("deleteFile",filename=self.workingDirPath('gallery_sampling.xmd'),parent_step_id=id)
        id = self.insertParallelStep("deleteFile",filename=self.workingDirPath('gallery_angles.doc'),parent_step_id=id)
        id = self.insertParallelStep("deleteFile",filename=self.workingDirPath('gallery.doc'),parent_step_id=id)
        
        return id

        
    def insertValidateVol(self,fnRoot, parent_id):
        
        fnGallery=self.workingDirPath('gallery_'+fnRoot+'.stk')
        fnVol = self.workingDirPath(fnRoot+'.vol')
        fnOutputReducedClass = self.workingDirPath("reducedClasses")
        
        id=self.insertParallelRunJobStep("xmipp_angular_project_library", "-i %s -o %s --sampling_rate %f --sym %s --method fourier 1 0.25 bspline --compute_neighbors --angular_distance -1 --experimental_images %s.xmd"\
                                  %(fnVol,fnGallery,float(self.AngularSampling),self.SymmetryGroup,fnOutputReducedClass),[fnGallery],parent_step_id=parent_id)
            
        # Assign angles
        fnAngles=self.workingDirPath('angles_'+fnRoot+'.xmd')
        
        id=self.insertParallelRunJobStep("xmipp_angular_projection_matching", "-i %s.xmd -o %s --ref %s --Ri 0 --Ro %s --max_shift %s --append"\
                              %(fnOutputReducedClass,fnAngles,fnGallery,str(self.Xdim/2),str(self.Xdim/20)),[fnAngles],parent_step_id=id)
       
        id = self.insertParallelStep("deleteFile",filename=self.workingDirPath('gallery_sampling.xmd'),parent_step_id=id)
        id = self.insertParallelStep("deleteFile",filename=self.workingDirPath('gallery_angles.doc'),parent_step_id=id)
        id = self.insertParallelStep("deleteFile",filename=self.workingDirPath('gallery.doc'),parent_step_id=id)
        id = self.insertParallelStep("validateVolume",WorkingDir=self.WorkingDir,fnRoot=fnRoot,CorrThresh=self.CorrThresh,
                                     parent_step_id=id)
        return id
        
    def summary(self):
        message = []
        return message
    
    def validate(self):
        errors = []
        return errors    

    def visualize(self):
        runShowJ(self.fnBestVolume)

def validateVolume(log,WorkingDir,fnRoot,CorrThresh):
    
    fnAngles=os.path.join(WorkingDir,"angles_"+fnRoot+".xmd")    
    md=MetaData(fnAngles)
    
    numInliers=0
    for objId in md:
        corr = md.getValue(MDL_MAXCC, objId)
       
        if (corr >= CorrThresh) :
            numInliers = numInliers+corr

    md= MetaData()
    objId = md.addObject()
    md.setValue(MDL_WEIGHT,float(numInliers),objId)
    md.write("inliers@"+fnAngles,MD_APPEND)
    
def getCorrThresh(WorkingDir,CorrThresh):
    
    corrVector = []
    for n in range(self.NRansac):
        
        fnRoot=self.workingDirPath("ransac%05d"%n)
        fnAngles=os.path.join(WorkingDir,"angles_"+fnRoot+".xmd")                            

        md=MetaData(fnAngles)

        for objId in md:
            corr = md.getValue(MDL_MAXCC, objId)
            corrVector.append(corr)
        
    sortedCorrVector = sorted(corrVector)
    indx = int(sum(sortedCorrVector)*self.CorrThresh)
    percentil = sortedCorrVector[idnx]
           
def getBestVolume(log,WorkingDir,NRansac,NumVolumes):
       
    volumes = []
    inliers = []
    
    for n in range(NRansac):
        fnAngles = os.path.join(WorkingDir,"angles_ransac%05d"%n+".xmd")
        md=MetaData("inliers@"+fnAngles)
        numInliers=md.getValue(MDL_WEIGHT,md.firstObject())
        volumes.append(fnAngles)
        inliers.append(numInliers)
    
    index = sorted(range(inliers.__len__()), key=lambda k: inliers[k])
    indx = 0
    fnBestAngles = ''
    
    for i in index[-NumVolumes:]:
        fnBestAngles = volumes[i]
        copyFile(log,fnBestAngles,os.path.join(WorkingDir,("bestAssignment%05d"%indx+".xmd")))       
        indx += 1