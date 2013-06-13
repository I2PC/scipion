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

from os.path import join, exists, split, basename
from xmipp import MetaData, MetaDataInfo, MD_APPEND, MDL_MAXCC, MDL_WEIGHT, \
    MDL_IMAGE, MDL_VOLUME_SCORE1, MDL_VOLUME_SCORE2, MDL_VOLUME_SCORE3, MDL_VOLUME_SCORE4

from protlib_base import *
from math import floor
from numpy import array, savetxt, sum, zeros
from protlib_utils import getListFromRangeString, runJob, runShowJ
from protlib_filesystem import copyFile, deleteFile, removeFilenamePrefix

class ProtInitVolRANSAC(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.initvolume_ransac.name, scriptname, project)
        self.Import = 'from protocol_initvolume_ransac import *'
        self.Xdim=MetaDataInfo(self.Classes)[0]
        
    def defineSteps(self):
        self.insertStep('createDir',path=self.ExtraDir)
        self.insertStep("linkAcquisitionInfo",InputFile=self.Classes,dirDest=self.WorkingDir)
        fnOutputReducedClass = self.tmpPath("reducedClasses.stk") 

        # Low pass filter and resize        
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
        
        self.insertRunJobStep("xmipp_transform_filter","-i %s -o %s --fourier low_pass %f"
                                                %(self.Classes,fnOutputReducedClass,freq*K))
        self.insertRunJobStep("xmipp_image_resize","-i %s --dim %d %d " %(fnOutputReducedClass,self.Xdim2,self.Xdim2))

        if (self.InitialVolume != ''):
            self.fnOutputInitVolume=self.extraPath("initialVolume.vol")
            self.insertStep('runJob',programname='xmipp_image_convert',
                            params="-i %s -o %s"%(removeFilenamePrefix(self.InitialVolume),self.fnOutputInitVolume),
                            NumberOfMpi=1)
            self.insertRunJobStep("xmipp_image_resize","-i %s --dim %d %d" 
                                  %(self.fnOutputInitVolume,self.Xdim2,self.Xdim2 ))
            fnGallery=self.tmpPath('gallery_InitialVolume.stk')
            self.insertRunJobStep("xmipp_angular_project_library", "-i %s -o %s --sampling_rate %f --sym %s --method fourier 1 0.25 bspline --compute_neighbors --angular_distance -1 --experimental_images %s"\
                                  %(self.fnOutputInitVolume,fnGallery,float(self.AngularSampling),self.SymmetryGroup,fnOutputReducedClass),[fnGallery])

        for n in range(self.NRansac):
            fnRoot=self.tmpPath("ransac%05d"%n)                        
            parent_id = self.insertParallelRunJobStep("xmipp_transform_dimred",
                                                      "-i %s --randomSample %s.xmd  %d -m LTSA "%(fnOutputReducedClass,fnRoot,self.NumGrids),
                                                      verifyfiles = [fnRoot+".xmd"],parent_step_id=XmippProjectDb.FIRST_STEP)
            
            if (self.InitialVolume != ''):
               parent_id=self.insertParallelRunJobStep("xmipp_angular_projection_matching", "-i %s.xmd -o %s.xmd --ref %s --Ri 0 --Ro %s --max_shift %s --append"\
                                                       %(fnRoot,fnRoot,fnGallery,str(self.Xdim/2),str(self.Xdim/20)),parent_step_id=parent_id)
                
            parent_id = self.insertParallelStep('reconstruct',fnRoot=fnRoot,symmetryGroup=self.SymmetryGroup,maskRadius=Xdim2/2, parent_step_id=parent_id)
            parent_id = self.insertParallelStep('evaluateVolume',fnRoot="ransac%05d"%n, TmpDir=self.TmpDir, AngularSampling=self.AngularSampling,
                                                SymmetryGroup=self.SymmetryGroup, Xdim=self.Xdim, parent_step_id=parent_id)
        
        if (self.InitialVolume != ''):
            parent_id=self.insertStep("runJob",programname="rm", params=self.tmpPath("gallery_InitialVolume*"), NumberOfMpi=1)
        self.insertStep("getCorrThresh",WorkingDir=self.WorkingDir, NRansac=self.NRansac, CorrThresh=self.CorrThresh)
        self.insertStep("validateVolumes",WorkingDir=self.WorkingDir,NRansac=self.NRansac)
        self.insertStep("getBestVolumes",WorkingDir=self.WorkingDir, NRansac=self.NRansac, NumVolumes=self.NumVolumes, UseAll=self.UseAll)        
        
        for n in range(self.NumVolumes):
            
            fnRoot=self.workingDirPath('proposedVolume%05d'%n)
            parent_id=XmippProjectDb.FIRST_STEP
            
            for it in range(self.NumIter):    
                parent_id = self.insertParallelStep('reconstruct',fnRoot=fnRoot,symmetryGroup=self.SymmetryGroup,maskRadius=Xdim2/2)
                parent_id = self.insertProjMatch(fnRoot,parent_id)
            
            self.insertParallelRunJobStep("xmipp_image_resize","-i %s.vol -o %s.vol --dim %d %d" 
                                          %(fnRoot,fnRoot,self.Xdim,self.Xdim),parent_step_id=parent_id)
        self.insertStep("scoreFinalVolumes",WorkingDir=self.WorkingDir,NumVolumes=self.NumVolumes)
        self.insertStep("deleteFile",filename=fnOutputReducedClass)
        
    def insertReconstruction(self,fnRoot, parent_id):

        id=self.insertParallelRunJobStep("xmipp_reconstruct_fourier","-i %s.xmd -o %s.vol --sym %s " %(fnRoot,fnRoot,self.SymmetryGroup),
                                      parent_step_id=parent_id)
        id=self.insertParallelRunJobStep("xmipp_transform_mask","-i %s.vol --mask circular -%d "%(fnRoot,self.Xdim2/2),
                                         parent_step_id=id)
        return id

    def insertProjMatch(self,fnRoot, parent_id):
        fnGallery=self.tmpPath('gallery_'+fnRoot+'.stk')
        fnOutputReducedClass = self.tmpPath("reducedClasses.stk") 
        
        id=self.insertParallelRunJobStep("xmipp_angular_project_library", "-i %s.vol -o %s --sampling_rate %f --sym %s --method fourier 1 0.25 bspline --compute_neighbors --angular_distance -1 --experimental_images %s"\
                              %(fnRoot,fnGallery,float(self.AngularSampling),self.SymmetryGroup,fnOutputReducedClass),[fnGallery],parent_step_id=parent_id)

        fnAngles=self.tmpPath('angles_BestVolume'+'.xmd')
        id=self.insertParallelRunJobStep("xmipp_angular_projection_matching", "-i %s.xmd -o %s --ref %s --Ri 0 --Ro %s --max_shift %s --append"\
                              %(fnRoot,fnAngles,fnGallery,str(self.Xdim/2),str(self.Xdim/20)),[fnAngles],parent_step_id=id)
                
        id = self.insertParallelStep("deleteFile",filename=self.tmpPath('gallery_'+fnRoot+'_sampling.xmd'),parent_step_id=id)
        id = self.insertParallelStep("deleteFile",filename=self.tmpPath('gallery_'+fnRoot+'.doc'),parent_step_id=id)
        id = self.insertParallelStep("deleteFile",filename=self.tmpPath('gallery_'+fnRoot+'.stk'),parent_step_id=id)
        
        return id
    
    def summary(self):
        message = ['Initial volume reconstruction by RANSAC using %d iterations and symmetry %s'%(self.NRansac,self.SymmetryGroup) ]
        return message
    
    def validate(self):
        errors = []
        return errors

    def visualize(self):
        n=0    
        fnVolumes = self.workingDirPath('proposedVolumes.xmd')
        runShowJ(fnVolumes)

def validateVolumes(log,WorkingDir,NRansac):
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
    CorrThresh = sortedCorrVector[indx]    
    objId = mdCorr.addObject()
    mdCorr.setValue(MDL_WEIGHT,float(CorrThresh),objId)
    mdCorr.write("corrThreshold@"+fnCorr,MD_APPEND)
    print "Correlation threshold: "+str(CorrThresh)

def getCCThreshold(WorkingDir):
    fnCorr=os.path.join(WorkingDir,"extra/correlations.xmd")               
    mdCorr=MetaData("corrThreshold@"+fnCorr)
    return mdCorr.getValue(MDL_WEIGHT, mdCorr.firstObject())
    
def getBestVolumes(log,WorkingDir,NRansac,NumVolumes,UseAll):
    volumes = []
    inliers = []
    
    Nimgs=MetaData(os.path.join(WorkingDir,"tmp/angles_ransac00000.xmd")).size()
    allCC=zeros((Nimgs,NRansac))
    for n in range(NRansac):
        fnAngles = os.path.join(WorkingDir,"tmp/angles_ransac%05d"%n+".xmd")
        md=MetaData("inliers@"+fnAngles)
        numInliers=md.getValue(MDL_WEIGHT,md.firstObject())
        volumes.append(fnAngles)
        inliers.append(numInliers)

        # Get the column of CCs 
        allCC[:,n]=MetaData(fnAngles).getColumnValues(MDL_MAXCC)
    
    savetxt(os.path.join(WorkingDir,'allCC.txt'),allCC)
    
    # Now each image score each volume
    for i in range(Nimgs):
        ccList=allCC[i,:]
        idx=sorted(range(len(ccList)),key=lambda x:ccList[x])
        allCC[i,:]=1.0-array(idx)/float(NRansac)
    savetxt(os.path.join(WorkingDir,'allCCScore.txt'),allCC)
    
    # Final volume score
    volumeRankScores=[]
    for n in range(NRansac):
        volumeRankScores.append(sum(allCC[:,n]))

    # index = sorted(range(volumeRankScores.__len__()), key=lambda k: volumeRankScores[k])
    index = sorted(range(inliers.__len__()), key=lambda k: inliers[k])
    fnBestAngles = ''
    threshold=getCCThreshold(WorkingDir)
    
    i=NRansac-1
    indx = 0
    while i>=0 and indx<NumVolumes:
        fnBestAngles = volumes[index[i]]
        fnBestAnglesOut=os.path.join(WorkingDir,"proposedVolume%05d"%indx+".xmd")
        copyFile(log,fnBestAngles,fnBestAnglesOut)
        print("Best volume "+str(indx)+" = "+fnBestAngles)
        if not UseAll:
            runJob(log,"xmipp_metadata_utilities","-i %s -o %s --query select \"maxCC>%f \" --mode append" %(fnBestAnglesOut,fnBestAnglesOut,threshold))
            md=MetaData(fnBestAnglesOut)
            if md.size()>0:
                indx += 1
        else:
            indx+=1
        i-=1
        
    # Remove unnecessary files
    for n in range(NRansac):
        fnAngles = os.path.join(WorkingDir,"tmp/angles_ransac%05d"%n+".xmd")
        deleteFile(log, fnAngles)
   
def scoreFinalVolumes(log,WorkingDir,NumVolumes):
    threshold=getCCThreshold(WorkingDir)
    mdOut=MetaData()
    for n in range(NumVolumes):
        fnRoot=os.path.join(WorkingDir,'proposedVolume%05d'%n)
        if exists(fnRoot+".xmd"):
            MDassignment=MetaData(fnRoot+".xmd")
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
            mdOut.setValue(MDL_VOLUME_SCORE1,float(sum),id)
            mdOut.setValue(MDL_VOLUME_SCORE2,float(thresholdedSum),id)
            mdOut.setValue(MDL_VOLUME_SCORE3,float(avg),id)
            mdOut.setValue(MDL_VOLUME_SCORE4,float(minCC),id)
    mdOut.write(os.path.join(WorkingDir,"proposedVolumes.xmd"))

def reconstruct(log,fnRoot,symmetryGroup,maskRadius):
    runJob(log,"xmipp_reconstruct_fourier","-i %s.xmd -o %s.vol --sym %s " %(fnRoot,fnRoot,symmetryGroup))
    runJob(log,"xmipp_transform_mask","-i %s.vol --mask circular -%d "%(fnRoot,maskRadius))

def evaluateVolume(log, fnRoot, TmpDir, AngularSampling, SymmetryGroup, Xdim):
    fnGallery=os.path.join(TmpDir,'gallery_'+fnRoot+'.stk')
    fnVol = os.path.join(TmpDir,fnRoot+'.vol')
    fnOutputReducedClass = os.path.join(TmpDir,"reducedClasses.stk")
    
    runJob(log,"xmipp_angular_project_library", "-i %s -o %s --sampling_rate %f --sym %s --method fourier 1 0.25 bspline --compute_neighbors --angular_distance -1 --experimental_images %s"\
                %(fnVol,fnGallery,float(AngularSampling),SymmetryGroup,fnOutputReducedClass))
        
    # Assign angles
    fnAngles=os.path.join(TmpDir,'angles_'+fnRoot+'.xmd')
    
    runJob(log,"xmipp_angular_projection_matching", "-i %s -o %s --ref %s --Ri 0 --Ro %s --max_shift %s --append"\
                          %(fnOutputReducedClass,fnAngles,fnGallery,str(Xdim/2),str(Xdim/20)))
   
    deleteFile(log,fnGallery)
    deleteFile(log,os.path.join(TmpDir,'gallery_'+fnRoot+'_sampling.xmd'))
    deleteFile(log,os.path.join(TmpDir,'gallery_'+fnRoot+'.doc'))
    deleteFile(log,fnVol)
    deleteFile(log,os.path.join(TmpDir,fnRoot+'.xmd'))
