#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for creating an initial volume using a dimensionality reduction method and RANSAC
#
# Example use:
# ./xmipp_protocol_screen_classes.py
#
# Author: Javier Vargas and Carlos Oscar May 2013 
#

from protlib_base import *
from os.path import join, exists, split, basename
from xmipp import MetaData, MetaDataInfo, MD_APPEND, MDL_MAXCC, MDL_WEIGHT, \
    MDL_IMAGE, \
    MDL_VOLUME_SCORE_SUM, MDL_VOLUME_SCORE_SUM_TH, MDL_VOLUME_SCORE_MEAN, MDL_VOLUME_SCORE_MIN

from math import floor
from numpy import array, savetxt, sum, zeros
from protlib_xmipp import getMdSize
from protlib_utils import getListFromRangeString, runJob, runShowJ
from protlib_filesystem import copyFile, deleteFile, moveFile, removeFilenamePrefix
from protlib_initvolume import *

class ProtHG3D(ProtHG3DBase):
    def __init__(self, scriptname, project):
        ProtInitVolumeBase.__init__(self, protDict.initvolume_ransac.name, scriptname, project)
        self.Import += 'from protocol_hg3d import *'
        
    def defineSteps(self):
        ProtInitVolumeBase.defineSteps(self)

        self.NumHg = 1;
        # Generate projection gallery from the initial volume
        if (self.InitialVolume != ''):
            self.insertStep("projectInitialVolume",WorkingDir=self.WorkingDir,InitialVolume=self.InitialVolume,Xdim2=self.Xdim2,
                            AngularSampling=self.AngularSampling,SymmetryGroup=self.SymmetryGroup)

        # RANSAC iterations
        for n in range(self.NRansac):
            self.insertParallelStep('ransacIteration',WorkingDir=self.WorkingDir,n=n,SymmetryGroup=self.SymmetryGroup,Xdim=self.Xdim,
                                    Xdim2=self.Xdim2,NumGrids=self.NumGrids,NumSamples=self.NumSamples,DimRed=self.DimRed,
                                    InitialVolume=self.InitialVolume,AngularSampling=self.AngularSampling,
                                    UseSA=self.UseSA, NIterRandom=self.NIterRandom, Rejection=self.Rejection,
                                    parent_step_id=XmippProjectDb.FIRST_STEP)
        
        # Look for threshold, evaluate volumes and get the best
        if (self.InitialVolume != ''):
            self.insertStep("runJob",programname="rm", params=self.tmpPath("gallery_InitialVolume*"), NumberOfMpi=1)
        
        self.insertStep("getCorrThresh",WorkingDir=self.WorkingDir, NRansac=self.NRansac, CorrThresh=self.CorrThresh)
        self.insertStep("evaluateVolumes",WorkingDir=self.WorkingDir,NRansac=self.NRansac)
        self.insertStep("getBestVolumes",WorkingDir=self.WorkingDir, NRansac=self.NRansac, NumVolumes=self.NumVolumes, UseAll=self.UseAll,NumHg=self.NumHg)        
        
        # Refine the best volumes
        for n in range(self.NumVolumes):
            fnBase='proposedVolume%05d'%n
            fnRoot=self.workingDirPath(fnBase)
            parent_id=XmippProjectDb.FIRST_STEP

            # Simulated annealing
            parent_id = self.insertParallelStep('reconstruct',fnRoot=fnRoot,symmetryGroup=self.SymmetryGroup,maskRadius=self.Xdim2/2,
                                                parent_step_id=parent_id)
            if self.UseSA:
                parent_id = self.insertParallelRunJobStep("xmipp_volume_initial_simulated_annealing","-i %s.xmd --initial %s.vol --oroot %s_sa --sym %s --randomIter %d --rejection %f --dontApplyPositive"\
                          %(fnRoot,fnRoot,fnRoot,self.SymmetryGroup,self.NIterRandom,self.Rejection),parent_step_id=parent_id)
                parent_id = self.insertParallelStep('moveFile', source=fnRoot+"_sa.vol", dest=fnRoot+".vol",parent_step_id=parent_id)
                parent_id = self.insertParallelStep('deleteFile', filename=fnRoot+"_sa.xmd",parent_step_id=parent_id)
        
            for it in range(self.NumIter):    
                parent_id = self.insertParallelStep('reconstruct',fnRoot=fnRoot,symmetryGroup=self.SymmetryGroup,maskRadius=self.Xdim2/2,
                                                    parent_step_id=parent_id)
                parent_id = self.insertParallelStep('projMatch',WorkingDir=self.WorkingDir,fnBase=fnBase,AngularSampling=self.AngularSampling,
                                                    SymmetryGroup=self.SymmetryGroup, Xdim=self.Xdim2, parent_step_id=parent_id)
            self.insertParallelRunJobStep("xmipp_image_resize","-i %s.vol -o %s.vol --dim %d %d" 
                                          %(fnRoot,fnRoot,self.Xdim,self.Xdim),parent_step_id=parent_id)
        
        # Score each of the final volumes
        self.insertStep("scoreFinalVolumes",WorkingDir=self.WorkingDir,NumVolumes=self.NumVolumes)
        
    def validate(self):
        errors = []
        if float(self.MaxFreq)<2*float(self.Ts):
            errors.append("Maximum frequency cannot be smaller than twice the sampling rate")
        return errors

    def summary(self):
        message=ProtInitVolumeBase.summary(self)
        message.append("RANSAC iterations: %d"%self.NRansac)
        
        for n in range(self.NumVolumes):

            fnBase='proposedVolume%05d'%n
            fnRoot=self.workingDirPath(fnBase+".xmd")
                           
            if os.path.isfile(fnRoot):
                md=MetaData(fnRoot)
                if (md.size()< 5) :
                    message.append("Num of inliers for %s too small and equal to %d"%(fnRoot,md.size()))
                    message.append("Decrease the value of Inlier Threshold parameter and run again")
                                
        fnBase="ransac00000.xmd"
        fnRoot=self.workingDirPath("tmp/"+fnBase)
        
        if os.path.isfile(fnRoot):
            md=MetaData(fnRoot)
        
            if (md.size()< 5) :
                message.append("Num of random samples too small and equal to %d"%(md.size()))
                message.append("If the option Dimensionality reduction is on, increase the number of grids per dimension")
                message.append("If the option Dimensionality reduction is off, increase the number of random samples")
            
        if self.UseSA:
            message.append("Simulated annealing used")
        return message
    
    def visualize(self):
        if self.DoShowList:
            fnVolumes = self.workingDirPath('proposedVolumes.xmd')
            runShowJ(fnVolumes)
        if self.VolumesToShow!="":
            listOfVolumes = getListFromRangeString(self.VolumesToShow)
            for volume in listOfVolumes:
                fnRoot=self.workingDirPath('proposedVolume%05d'%int(volume))
                os.system("xmipp_chimera_client -i %s.vol --mode projector 256 --angulardist %s.xmd"%(fnRoot,fnRoot))

def evaluateVolumes(log,WorkingDir,NRansac):
    fnCorr=os.path.join(WorkingDir,"extra/correlations.xmd")
    fnCorr = 'corrThreshold@'+fnCorr
    mdCorr= MetaData(fnCorr)
    objId = mdCorr.firstObject()    
    CorrThresh = mdCorr.getValue(MDL_WEIGHT,objId)
    for n in range(NRansac):        
        fnRoot=os.path.join("ransac%05d"%n)              
        fnAngles=os.path.join(WorkingDir,"tmp/angles_"+fnRoot+".xmd")    
        md=MetaData(fnAngles)
        numInliers=0
        for objId in md:
            corr = md.getValue(MDL_MAXCC, objId)
           
            if (corr >= CorrThresh) :
#JV: QUIZA cambiar esta linea a: numInliers += numInliers
                numInliers = numInliers+corr

        md= MetaData()
        objId = md.addObject()
        md.setValue(MDL_WEIGHT,float(numInliers),objId)
        md.write("inliers@"+fnAngles,MD_APPEND)
    
def getCorrThresh(log,WorkingDir,NRansac,CorrThresh):
    corrVector = []
    fnCorr=os.path.join(WorkingDir,"extra/correlations.xmd")               
    mdCorr= MetaData()

    for n in range(NRansac):
        fnRoot=os.path.join("ransac%05d"%n)
        fnAngles=os.path.join(WorkingDir,"tmp/angles_"+fnRoot+".xmd")
        md=MetaData(fnAngles)
        
        for objId in md:
            corr = md.getValue(MDL_MAXCC, objId)
            corrVector.append(corr)
            objIdCorr = mdCorr.addObject()
            mdCorr.setValue(MDL_MAXCC,float(corr),objIdCorr)

    mdCorr.write("correlations@"+fnCorr,MD_APPEND)                            
    mdCorr= MetaData()
    sortedCorrVector = sorted(corrVector)
    indx = int(floor(CorrThresh*(len(sortedCorrVector)-1)))    
    
    #With the line below commented the percentil is not used for the threshold and is used the value introduced in the form
    
    #CorrThresh = sortedCorrVector[indx]#
        
    objId = mdCorr.addObject()
    mdCorr.setValue(MDL_WEIGHT,float(CorrThresh),objId)
    mdCorr.write("corrThreshold@"+fnCorr,MD_APPEND)
    print "Correlation threshold: "+str(CorrThresh)

def getCCThreshold(WorkingDir):
    fnCorr=os.path.join(WorkingDir,"extra/correlations.xmd")               
    mdCorr=MetaData("corrThreshold@"+fnCorr)
    return mdCorr.getValue(MDL_WEIGHT, mdCorr.firstObject())
    
def getBestVolumes(log,WorkingDir,NRansac,NumVolumes,UseAll,NumHg):
    
    StructureDir = os.path.join(WorkingDir,"Structure%05d"%NumHg) 
    createDir(log,StructureDir)
    
    volumes = []
    inliers = []
    
    for n in range(NRansac):
        fnAngles = os.path.join(WorkingDir,"tmp/angles_ransac%05d"%n+".xmd")
        md=MetaData("inliers@"+fnAngles)
        numInliers=md.getValue(MDL_WEIGHT,md.firstObject())
        volumes.append(fnAngles)
        inliers.append(numInliers)
    
    index = sorted(range(inliers.__len__()), key=lambda k: inliers[k])
    fnBestAngles = ''
    threshold=getCCThreshold(WorkingDir)
 
    i=NRansac-1
    indx = 0
    while i>=0 and indx<NumVolumes:
        fnBestAngles = volumes[index[i]]
        #fnBestAnglesOut=os.path.join(WorkingDir,"proposedVolume%05d"%indx+".xmd")
        fnBestAnglesOut=os.path.join(StructureDir,"proposedVolume%05d"%indx+".xmd")
        copyFile(log,fnBestAngles,fnBestAnglesOut)
        print("Best volume "+str(indx)+" = "+fnBestAngles)
        if not UseAll:
            runJob(log,"xmipp_metadata_utilities","-i %s -o %s --query select \"maxCC>%f \" --mode append" %(fnBestAnglesOut,fnBestAnglesOut,threshold))
            if getMdSize(fnBestAnglesOut) > 0:
                indx += 1
        else:
            indx += 1
        i -= 1
        
    # Remove unnecessary files
    for n in range(NRansac):
        fnAngles = os.path.join(WorkingDir,"tmp/angles_ransac%05d"%n+".xmd")
        deleteFile(log, fnAngles)
   
def projectInitialVolume(log,WorkingDir,InitialVolume,Xdim2,AngularSampling,SymmetryGroup):
    fnOutputInitVolume=os.path.join(WorkingDir,"tmp/initialVolume.vol")
    runJob(log,'xmipp_image_convert',"-i %s -o %s"%(removeFilenamePrefix(InitialVolume),fnOutputInitVolume))
    runJob(log,"xmipp_image_resize","-i %s --dim %d %d"%(fnOutputInitVolume,Xdim2,Xdim2))
    fnGallery=os.path.join(WorkingDir,'tmp/gallery_InitialVolume.stk')
    fnOutputReducedClass = os.path.join(WorkingDir,"extra/reducedClasses.xmd") 
    runJob(log,"xmipp_angular_project_library", "-i %s -o %s --sampling_rate %f --sym %s --method fourier 1 0.25 bspline --compute_neighbors --angular_distance -1 --experimental_images %s"\
                          %(fnOutputInitVolume,fnGallery,float(AngularSampling),SymmetryGroup,fnOutputReducedClass))

def reconstruct(log,fnRoot,symmetryGroup,maskRadius):
    runJob(log,"xmipp_reconstruct_fourier","-i %s.xmd -o %s.vol --sym %s " %(fnRoot,fnRoot,symmetryGroup))
    runJob(log,"xmipp_transform_mask","-i %s.vol --mask circular -%d "%(fnRoot,maskRadius))

def ransacIteration(log,WorkingDir,n,SymmetryGroup,Xdim,Xdim2,NumGrids,NumSamples,DimRed,InitialVolume,AngularSampling,UseSA,
                    NIterRandom,Rejection):
    
    NumSameFile = 20
    fnBase="ransac%05d"%n
    TmpDir=os.path.join(WorkingDir,"tmp")
    fnRoot=os.path.join(TmpDir,fnBase)
    fnOutputReducedClass = os.path.join(WorkingDir,"extra/reducedClasses.xmd")

    if ((n % NumSameFile)==0):
        runJob(log,"xmipp_metadata_utilities","-i %s -o %s.xmd  --operate random_subset %d --mode overwrite "%(fnOutputReducedClass,fnRoot,NumSamples))        
    else:
        fnBaseBef="ransac%05d"%(n-1)
        TmpDir=os.path.join(WorkingDir,"tmp")
        fnRoot=os.path.join(TmpDir,fnBaseBef)
        runJob(log,"cp","%s.xmd  %s.xmd "%(fnRoot,fnRoot))            
    
    runJob(log,"xmipp_metadata_utilities","-i %s.xmd --fill angleRot rand_uniform -180 180 "%(fnRoot))
    runJob(log,"xmipp_metadata_utilities","-i %s.xmd --fill angleTilt rand_uniform 0 180 "%(fnRoot))
    runJob(log,"xmipp_metadata_utilities","-i %s.xmd --fill anglePsi  rand_uniform 0 360 "%(fnRoot)) 
 
    # If there is an initial volume, assign angles        
    if (InitialVolume != ''):
        fnGallery=os.path.join(TmpDir,'gallery_InitialVolume.stk')
        runJob(log,"xmipp_angular_projection_matching", "-i %s.xmd -o %s.xmd --ref %s --Ri 0 --Ro %s --max_shift %s --append"\
               %(fnRoot,fnRoot,fnGallery,str(Xdim/2),str(Xdim/20)))

    # Reconstruct with the small sample
    reconstruct(log,fnRoot,SymmetryGroup,Xdim2/2)
    fnVol = fnRoot+'.vol'
    
    # Simulated annealing
    if UseSA:
        smallIter=int(min(floor(NIterRandom/5.0),0));
        runJob(log,"xmipp_volume_initial_simulated_annealing","-i %s --initial %s --oroot %s_sa --sym %s --randomIter %d --rejection %f --dontApplyPositive"
                  %(fnRoot+".xmd",fnVol,fnRoot,SymmetryGroup,smallIter,Rejection))
        moveFile(log, fnRoot+"_sa.vol", fnVol)
        deleteFile(log, fnRoot+"_sa.xmd")

    # Generate projections from this reconstruction
    fnGallery=os.path.join(TmpDir,'gallery_'+fnBase+'.stk')
    runJob(log,"xmipp_angular_project_library", "-i %s -o %s --sampling_rate %f --sym %s --method fourier 1 0.25 bspline --compute_neighbors --angular_distance -1 --experimental_images %s"\
                %(fnVol,fnGallery,float(AngularSampling),SymmetryGroup,fnOutputReducedClass))
        
    # Assign angles to the rest of images
    fnAngles=os.path.join(TmpDir,'angles_'+fnBase+'.xmd')
    runJob(log,"xmipp_angular_projection_matching", "-i %s -o %s --ref %s --Ri 0 --Ro %s --max_shift %s --append"\
                          %(fnOutputReducedClass,fnAngles,fnGallery,str(Xdim/2),str(Xdim/20)))
   
    # Delete intermediate files 
    deleteFile(log,fnGallery)
    deleteFile(log,os.path.join(TmpDir,'gallery_'+fnBase+'_sampling.xmd'))
    deleteFile(log,os.path.join(TmpDir,'gallery_'+fnBase+'.doc'))
    deleteFile(log,fnVol)
    deleteFile(log,os.path.join(TmpDir,fnBase+'.xmd'))

def projMatch(log, WorkingDir, fnBase, AngularSampling, SymmetryGroup, Xdim):
    fnRoot=os.path.join(WorkingDir,fnBase)
    fnGallery=os.path.join(WorkingDir,'tmp/gallery_'+fnBase+'.stk')
    fnOutputReducedClass = os.path.join(WorkingDir,"extra/reducedClasses.xmd") 
    
    AngularSampling=int(max(floor(AngularSampling/2.0),2));
    runJob(log,"xmipp_angular_project_library", "-i %s.vol -o %s --sampling_rate %f --sym %s --method fourier 1 0.25 bspline --compute_neighbors --angular_distance -1 --experimental_images %s"\
                          %(fnRoot,fnGallery,float(AngularSampling),SymmetryGroup,fnOutputReducedClass))

    runJob(log,"xmipp_angular_projection_matching", "-i %s.xmd -o %s.xmd --ref %s --Ri 0 --Ro %s --max_shift %s --append"\
           %(fnRoot,fnRoot,fnGallery,str(Xdim/2),str(Xdim/20)))
            
    deleteFile(log,os.path.join(WorkingDir,'tmp/gallery_'+fnBase+'_sampling.xmd'))
    deleteFile(log,os.path.join(WorkingDir,'tmp/gallery_'+fnBase+'.doc'))
    deleteFile(log,os.path.join(WorkingDir,'tmp/gallery_'+fnBase+'.stk'))

def scoreFinalVolumes(log,WorkingDir,NumVolumes):
    threshold=getCCThreshold(WorkingDir)
    mdOut=MetaData()
    for n in range(NumVolumes):
        fnRoot=os.path.join(WorkingDir,'proposedVolume%05d'%n)
        fnAssignment=fnRoot+".xmd"
        if exists(fnAssignment):
            runJob(log,"xmipp_metadata_utilities","-i %s --fill weight constant 1"%fnAssignment)
            MDassignment=MetaData(fnAssignment)
            sum=0
            thresholdedSum=0
            N=0
            minCC=2
            for id in MDassignment:
                cc=MDassignment.getValue(MDL_MAXCC,id)
                sum+=cc
                thresholdedSum+=cc-threshold
                if cc<minCC:
                    minCC=cc
                N+=1
            avg=sum/N
            id=mdOut.addObject()
            mdOut.setValue(MDL_IMAGE,fnRoot+".vol",id)
            mdOut.setValue(MDL_VOLUME_SCORE_SUM,float(sum),id)
            mdOut.setValue(MDL_VOLUME_SCORE_SUM_TH,float(thresholdedSum),id)
            mdOut.setValue(MDL_VOLUME_SCORE_MEAN,float(avg),id)
            mdOut.setValue(MDL_VOLUME_SCORE_MIN,float(minCC),id)
    mdOut.write(os.path.join(WorkingDir,"proposedVolumes.xmd"))

