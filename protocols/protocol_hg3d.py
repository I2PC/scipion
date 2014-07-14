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
from xmipp import MetaData, MetaDataInfo, MD_APPEND, MDL_MAXCC, MDL_WEIGHT, MDL_IMAGE, MDL_VOLUME_SCORE_SUM, MDL_VOLUME_SCORE_SUM_TH, \
                  MDL_VOLUME_SCORE_MEAN, MDL_VOLUME_SCORE_MIN
from math import floor
from numpy import array, savetxt, sum, zeros,  empty, loadtxt, ones, where, percentile
from protlib_xmipp import getMdSize
from protlib_utils import getListFromRangeString, runJob, runShowJ
from protlib_filesystem import copyFile, deleteFile, moveFile, removeFilenamePrefix, createDir
from protlib_hg3d import *


class ProtHG3D(ProtHG3DBase):
    def __init__(self, scriptname, project):
        ProtHG3DBase.__init__(self, protDict.hg3d.name, scriptname, project)
        self.Import += 'from protocol_hg3d import *'
        
    def defineSteps(self):
        ProtHG3DBase.defineSteps(self)
        
        md=MetaData(self.RemainingClasses)
        nI=0
        
        # Generate projection gallery from the initial volume
        if (self.InitialVolume != ''):
            self.insertStep("projectInitialVolume",WorkingDir=self.WorkingDir,InitialVolume=self.InitialVolume,
                            RemainingClasses=self.RemainingClasses,Xdim2=self.Xdim2,
                            AngularSampling=self.AngularSampling,SymmetryGroup=self.SymmetryGroup)
        
        # Compute RANSAC with all remaining images
        WorkingDirStructure = os.path.join(self.WorkingDir,"Structure%05d"%nI)
        self.insertStep('createDir',path=WorkingDirStructure)
        fnClasses=self.RemainingClasses
        self.defineRANSACSteps(fnClasses,WorkingDirStructure,False,self.NRansacInitial,self.NumVolumesInitial,True)
        self.defineRefinementSteps(fnClasses,WorkingDirStructure,self.NumVolumesInitial)
        
        if (self.InitialVolume != ''):
            self.insertStep("runJob",programname="rm", params=self.tmpPath("gallery_InitialVolume*"), NumberOfMpi=1)
        
        # Compute imagesCore for this structure
        WorkingDirStructureCore = os.path.join(self.WorkingDir,"Structure%05d_core"%nI)
        self.insertStep('createDir',path=WorkingDirStructureCore)
        self.insertStep("coocurenceMatrix",RemainingClasses=self.RemainingClasses,WorkingDirStructure=WorkingDirStructure,
                            NumVolumes=self.NumVolumesInitial,nI=nI,CorePercentile=self.CorePercentile,
                            CorrThresh=self.CorrThresh)

        # Compute RANSAC with the core
        fnClasses=os.path.join(WorkingDirStructureCore,'imagesCore.xmd')
        self.defineRANSACSteps(fnClasses,WorkingDirStructureCore,True,self.NRansacCore,self.NumVolumesFinal,False)
        self.defineRefinementSteps(fnClasses,WorkingDirStructureCore,self.NumVolumesFinal)
                               
        # Match the rest of images and extend the core
        WorkingDirStructureExtended = os.path.join(self.WorkingDir,"Structure%05d_core_extended"%nI)
        self.insertStep('createDir',path=WorkingDirStructureExtended)
        self.insertStep('extendCore',WorkingDirStructureCore=WorkingDirStructureCore,WorkingDirStructureExtended=WorkingDirStructureExtended,
                        NumVolumes=self.NumVolumesFinal,RemainingClasses=self.RemainingClasses,ComplementaryClasses=self.Classes,
                        AngularSampling=self.AngularSampling,SymmetryGroup=self.SymmetryGroup, CorrThresh=self.CorrThresh,
                        NumberOfMpi=self.NumberOfMpi)
        
        # Refine the extended cores
        fnClasses=self.RemainingClasses
        self.defineRefinementSteps(fnClasses,WorkingDirStructureExtended,self.NumVolumesFinal)
        
        # Resize the output volumes
        for i in range(self.NumVolumesFinal):
            fnRoot=os.path.join(WorkingDirStructureExtended,"proposedVolume%05d"%i)
            self.insertParallelRunJobStep("xmipp_image_resize","-i %s.vol -o %s.vol --dim %d %d" 
                                          %(fnRoot,fnRoot,self.Xdim,self.Xdim))
        
        self.insertStep("scoreFinalVolumes",WorkingDir=self.WorkingDir,NumVolumes=self.NumVolumesFinal)

    def defineRANSACSteps(self,fnClasses,WorkingDirStructure,useAll,NRansac,NumVolumes,UseInitial):
        # RANSAC iterations
        if UseInitial:
            InitialVolume=self.InitialVolume
        else:
            InitialVolume=""
        for n in range(NRansac):
            self.insertParallelStep('ransacIteration',WorkingDir=self.WorkingDir,fnClasses=fnClasses,
                                    n=n,SymmetryGroup=self.SymmetryGroup,Xdim=self.Xdim,
                                    Xdim2=self.Xdim2,NumSamples=self.NumSamples,InitialVolume=InitialVolume,
                                    AngularSampling=self.AngularSampling,
                                    parent_step_id=XmippProjectDb.FIRST_STEP)
        
        self.insertStep("getCorrThresh",WorkingDir=self.WorkingDir, WorkingDirStructure=WorkingDirStructure,
                        NRansac=NRansac, CorrThresh=self.CorrThresh)
        self.insertStep("evaluateVolumes",WorkingDir=self.WorkingDir,NRansac=NRansac,WorkingDirStructure=WorkingDirStructure)        
        self.insertStep("getBestVolumes",WorkingDir=self.WorkingDir, WorkingDirStructure=WorkingDirStructure,
                        NRansac=NRansac, NumVolumes=NumVolumes, UseAll=useAll)        
    
    def defineRefinementSteps(self,fnClasses,WorkingDirStructure,NumVolumes):
        # Refine the best volumes
        for n in range(NumVolumes):
            fnBase='proposedVolume%05d'%n
            fnRoot=os.path.join(WorkingDirStructure,fnBase)
            parent_id=XmippProjectDb.FIRST_STEP
            parent_id = self.insertParallelStep('reconstruct',fnRoot=fnRoot,symmetryGroup=self.SymmetryGroup,maskRadius=self.Xdim2/2,
                                                parent_step_id=parent_id)
            for it in range(self.NumIter):    
                parent_id = self.insertParallelStep('projMatch',WorkingDir=self.WorkingDir,WorkingDirStructure=WorkingDirStructure,
                                                    fnClasses=fnClasses, fnBase=fnBase, AngularSampling=self.AngularSampling,
                                                    SymmetryGroup=self.SymmetryGroup, Xdim=self.Xdim2, parent_step_id=parent_id)
                parent_id = self.insertParallelStep('reconstruct',fnRoot=fnRoot,symmetryGroup=self.SymmetryGroup,maskRadius=self.Xdim2/2,
                                                    parent_step_id=parent_id)

    def validate(self):
        errors = []
        return errors

    def summary(self):
        message=ProtHG3DBase.summary(self)
        message.append("RANSAC iterations: %d"%self.NRansacInitial)
        for n in range(self.NumVolumesFinal):
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
        
        return message
    
    def visualize(self):
        if self.DoShowList:
            fnVolumes = self.workingDirPath('proposedVolumes.xmd')
            runShowJ(fnVolumes)
        if self.VolumesToShow!="":
            listOfVolumes = getListFromRangeString(self.VolumesToShow)
            for volume in listOfVolumes:
                fnRoot=os.path.join(WorkingDir,'Structure00000_core_extended','proposedVolume%05d'%volume)
                os.system("xmipp_chimera_client -i %s.vol --mode projector 256 --angulardist %s.xmd"%(fnRoot,fnRoot))

def evaluateVolumes(log,WorkingDir,NRansac,WorkingDirStructure):
    fnCorr=os.path.join(WorkingDirStructure,"correlations.xmd")
    fnCorr = 'corrThreshold@'+fnCorr
    mdCorr= MetaData(fnCorr)
    objId = mdCorr.firstObject()    
    CorrThresh = mdCorr.getValue(MDL_WEIGHT,objId)
    for n in range(NRansac):        
        fnRoot="ransac%05d"%n              
        fnAngles=os.path.join(WorkingDir,"tmp/angles_"+fnRoot+".xmd")
        if os.path.exists(fnAngles):    
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
    
def getCorrThresh(log,WorkingDir,WorkingDirStructure,NRansac,CorrThresh):
    corrVector = []
    fnCorr=os.path.join(WorkingDirStructure,"correlations.xmd")               
    mdCorr= MetaData()

    for n in range(NRansac):
        fnRoot=os.path.join("ransac%05d"%n)
        fnAngles=os.path.join(WorkingDir,"tmp/angles_"+fnRoot+".xmd")
        if os.path.exists(fnAngles):
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
    
    objId = mdCorr.addObject()
    mdCorr.setValue(MDL_WEIGHT,float(CorrThresh),objId)
    mdCorr.write("corrThreshold@"+fnCorr,MD_APPEND)
    print "Correlation threshold: "+str(CorrThresh)

def getBestVolumes(log,WorkingDir,WorkingDirStructure,NRansac,NumVolumes,UseAll):    
    volumes = []
    inliers = []
    
    for n in range(NRansac):
        fnAngles = os.path.join(WorkingDir,"tmp/angles_ransac%05d"%n+".xmd")
        if os.path.exists(fnAngles):
            md=MetaData("inliers@"+fnAngles)
            numInliers=md.getValue(MDL_WEIGHT,md.firstObject())
            volumes.append(fnAngles)
            inliers.append(numInliers)
        else:
            inliers.append(0)
    
    index = sorted(range(inliers.__len__()), key=lambda k: inliers[k])
    fnBestAngles = ''

    fnCorr=os.path.join(WorkingDirStructure,"correlations.xmd")               
    mdCorr=MetaData("corrThreshold@"+fnCorr)
    threshold=mdCorr.getValue(MDL_WEIGHT, mdCorr.firstObject())
 
    i=NRansac-1
    indx = 0
    while i>=0 and indx<NumVolumes:
        fnBestAngles = volumes[index[i]]
        #fnBestAnglesOut=os.path.join(WorkingDir,"proposedVolume%05d"%indx+".xmd")
        fnBestAnglesOut=os.path.join(WorkingDirStructure,"proposedVolume%05d"%indx+".xmd")
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
   
def projectInitialVolume(log,WorkingDir,InitialVolume,RemainingClasses,Xdim2,AngularSampling,SymmetryGroup):
    fnOutputInitVolume=os.path.join(WorkingDir,"tmp/initialVolume.vol")
    runJob(log,'xmipp_image_convert',"-i %s -o %s"%(removeFilenamePrefix(InitialVolume),fnOutputInitVolume))
    runJob(log,"xmipp_image_resize","-i %s --dim %d %d"%(fnOutputInitVolume,Xdim2,Xdim2))
    fnGallery=os.path.join(WorkingDir,'tmp/gallery_InitialVolume.stk')
    runJob(log,"xmipp_angular_project_library", "-i %s -o %s --sampling_rate %f --sym %s --method fourier 1 0.25 bspline --compute_neighbors --angular_distance -1 --experimental_images %s"\
                          %(fnOutputInitVolume,fnGallery,float(AngularSampling),SymmetryGroup,RemainingClasses))

def reconstruct(log,fnRoot,symmetryGroup,maskRadius):
    runJob(log,"xmipp_reconstruct_fourier","-i %s.xmd -o %s.vol --sym %s -v 0" %(fnRoot,fnRoot,symmetryGroup))
    runJob(log,"xmipp_transform_mask","-i %s.vol --mask circular -%d -v 0"%(fnRoot,maskRadius))

def ransacIteration(log,WorkingDir,n,fnClasses,SymmetryGroup,Xdim,Xdim2,NumSamples,InitialVolume,AngularSampling):
    fnBase="ransac%05d"%n
    TmpDir=os.path.join(WorkingDir,"tmp")
    fnRoot=os.path.join(TmpDir,fnBase)

    runJob(log,"xmipp_metadata_utilities","-i %s -o %s.xmd  --operate random_subset %d --mode overwrite "%(fnClasses,fnRoot,NumSamples))        
          
    
    runJob(log,"xmipp_metadata_utilities","-i %s.xmd --fill angleRot  rand_uniform -180 180 "%(fnRoot))
    runJob(log,"xmipp_metadata_utilities","-i %s.xmd --fill angleTilt rand_uniform 0 180 "%(fnRoot))
    runJob(log,"xmipp_metadata_utilities","-i %s.xmd --fill anglePsi  rand_uniform 0 360 "%(fnRoot)) 
 
    # If there is an initial volume, assign angles        
    if (InitialVolume != ''):
        fnGallery=os.path.join(TmpDir,'gallery_InitialVolume.stk')
        runJob(log,"xmipp_angular_projection_matching", "-i %s.xmd -o %s.xmd --ref %s --Ri 0 --Ro %s --max_shift %s --append -v 0"\
               %(fnRoot,fnRoot,fnGallery,str(Xdim/2),str(Xdim/20)))

    # Reconstruct with the small sample
    reconstruct(log,fnRoot,SymmetryGroup,Xdim2/2)
    fnVol = fnRoot+'.vol'
    
    # Generate projections from this reconstruction
    fnGallery=os.path.join(TmpDir,'gallery_'+fnBase+'.stk')
    runJob(log,"xmipp_angular_project_library", "-i %s -o %s --sampling_rate %f --sym %s --method fourier 1 0.25 bspline --compute_neighbors --angular_distance -1 --experimental_images %s -v 0"\
                %(fnVol,fnGallery,float(AngularSampling),SymmetryGroup,fnClasses))
        
    # Assign angles to the rest of images
    fnAngles=os.path.join(TmpDir,'angles_'+fnBase+'.xmd')
    runJob(log,"xmipp_angular_projection_matching", "-i %s -o %s --ref %s --Ri 0 --Ro %s --max_shift %s --append -v 0"\
                          %(fnClasses,fnAngles,fnGallery,str(Xdim/2),str(Xdim/20)))
   
    # Delete intermediate files 
    deleteFile(log,fnGallery)
    deleteFile(log,os.path.join(TmpDir,'gallery_'+fnBase+'_sampling.xmd'))
    deleteFile(log,os.path.join(TmpDir,'gallery_'+fnBase+'.doc'))
    deleteFile(log,fnVol)
    deleteFile(log,os.path.join(TmpDir,fnBase+'.xmd'))

def projMatch(log, WorkingDir,WorkingDirStructure, fnClasses, fnBase, AngularSampling, SymmetryGroup, Xdim):
    fnRoot=os.path.join(WorkingDirStructure,fnBase)
    fnGallery=os.path.join(WorkingDir,'tmp/gallery_'+fnBase+'.stk')
    
    AngularSampling=int(max(floor(AngularSampling/2.0),2));
    runJob(log,"xmipp_angular_project_library", "-i %s.vol -o %s --sampling_rate %f --sym %s --method fourier 1 0.25 bspline --compute_neighbors --angular_distance -1 --experimental_images %s -v 0"\
                          %(fnRoot,fnGallery,float(AngularSampling),SymmetryGroup,fnClasses))

    runJob(log,"xmipp_angular_projection_matching", "-i %s.xmd -o %s.xmd --ref %s --Ri 0 --Ro %s --max_shift %s --append -v 0"\
           %(fnRoot,fnRoot,fnGallery,str(Xdim/2),str(Xdim/20)))
            
    deleteFile(log,os.path.join(WorkingDir,'tmp/gallery_'+fnBase+'_sampling.xmd'))
    deleteFile(log,os.path.join(WorkingDir,'tmp/gallery_'+fnBase+'.doc'))
    deleteFile(log,os.path.join(WorkingDir,'tmp/gallery_'+fnBase+'.stk'))

def coocurenceMatrix(log,RemainingClasses,WorkingDirStructure,NumVolumes,nI,CorePercentile,CorrThresh):
    import numpy
    mdRemaining = MetaData(RemainingClasses)
    Nimgs=mdRemaining.size()
    allNames=mdRemaining.getColumnValues(MDL_IMAGE)
    matrixTotal = numpy.zeros([Nimgs,Nimgs])
    for n in range(NumVolumes):
        fnBase='proposedVolume%05d'%n
        fnRoot=os.path.join(WorkingDirStructure,fnBase)
        
        md = MetaData(fnRoot+".xmd")
        size = md.size()
        
        num=[]
        corr=[]
        for objId in md:
            name = md.getValue(MDL_IMAGE, objId)
            if name in allNames:
                num.append(allNames.index(name))
                corr.append(md.getValue(MDL_MAXCC, objId))
            else:
                print "Cannot find ",name      
        if size!=len(num):
            print "Error when processing: ",fnRoot+".xmd"
            aaa
        
        matrix = numpy.zeros([Nimgs,Nimgs])
        for i in range(size):
            for j in range(size):
                matrix[num[i],num[j]]=((corr[i]+corr[j])/2)
        
        #numpy.savetxt(os.path.join(WorkingDirStructure,'coocurrenceMatrix_%05d.txt'%n), matrix) 
        matrixTotal=matrixTotal+matrix
    matrixTotal=matrixTotal/NumVolumes
    numpy.savetxt(os.path.join(WorkingDirStructure,'coocurrenceMatrix.txt'),matrixTotal)
    largestComponent=procCoocurenceMatrix(matrixTotal,CorePercentile,CorrThresh)

    md = MetaData()
    for idx in largestComponent:
        id=md.addObject()
        md.setValue(MDL_IMAGE,allNames[idx],id)
    md.write(os.path.join(WorkingDirStructure+"_core","imagesCore.xmd"))
    if md.size()==0:
        print "There are no images in the core"
        aaa

#http://wiki.scipy.org/NumPy_for_Matlab_Users        
def procCoocurenceMatrix(sumMatrix,CorePercentile,CorrThresh):
    import numpy
    import Queue as q
    
    # Binarize sum matrix to the 80%
    y=sumMatrix[numpy.where(sumMatrix!=0)]
    threshold=numpy.percentile(y,CorePercentile)
    print "Correlation threshold due to CorePercentile=",threshold
    intSumMatrix = sumMatrix
    intSumMatrix[sumMatrix<threshold]=0
    intSumMatrix[sumMatrix>=threshold]=1
    
    #Compute connected components
    num = intSumMatrix.shape[0]
    rowSum=sum(sumMatrix,0)
    remaining = numpy.ones(num)
    queue = q.Queue()
    components = numpy.empty((num,num))
    components[:] = -1
    nComp = -1
    idxComp = 0
    
    for i in range(num):
       if rowSum[i]==0:
          remaining[i]=0
    
    while (sum(remaining) != 0):
        if (queue.empty()):
            idy = 0
            while ((remaining[idy])==0):
                idy += 1
            nComp += 1
            idxComp = 0
            remaining[idy] = 0
        else:
            idy=int(queue.get())
        for idx in range(idy,num):
            if (intSumMatrix[idy,idx]!=0 and remaining[idx]==1):
                queue.put(idx)
                remaining[idx] = 0
                components[nComp,idxComp]=idx
                idxComp += 1
    
    # Compute largest component
    nBest=-1
    nBestSize=-1
    for i in range(nComp+1):
        currentSize=0
        j=0
        while components[i,j]!=-1:
           currentSize+=1
           j+=1
           if j==num:
                 break
        if currentSize>nBestSize:
            nBestSize=currentSize
            nBest=i
    
    # Isolate largest component
    largestComponent=[]
    for j in range(nBestSize):
        largestComponent.append(int(components[nBest,j]))
    largestComponent.sort()
    return largestComponent

def extendCore(log,WorkingDirStructureCore,WorkingDirStructureExtended,NumVolumes,RemainingClasses,ComplementaryClasses,AngularSampling,\
               SymmetryGroup,CorrThresh,NumberOfMpi):
    from protlib_projmatch import projMatch
    
    fnCore=os.path.join(WorkingDirStructureCore,"imagesCore.xmd")
    fnComplementary=os.path.join(WorkingDirStructureExtended,"imagesComplementary.xmd")
    runJob(log,"xmipp_metadata_utilities","-i %s --set subtraction %s image image -o %s"%(RemainingClasses,fnCore,fnComplementary))
    runJob(log,"xmipp_metadata_utilities","-i %s --operate keep_column image"%fnComplementary)
    if ComplementaryClasses!="":
        fnAux=os.path.join(WorkingDirStructureExtended,"aux.xmd")
        runJob(log,"xmipp_metadata_utilities","-i %s -o %s --operate keep_column image"%(ComplementaryClasses,fnAux))
        runJob(log,"xmipp_metadata_utilities", "-i %s  --set union %s image image"%(fnComplementary,fnAux))
        deleteFile(log,fnAux)
        runJob(log,"xmipp_metadata_utilities", "-i %s  --set subtraction %s image image"%(fnComplementary,fnCore))
        
    for i in range(NumVolumes):
        # Projection matching with the complementary images
        fnVolume=os.path.join(WorkingDirStructureCore,"proposedVolume%05d.vol"%i)
        fnAngles=os.path.join(WorkingDirStructureExtended,"anglesComplementary%05d.xmd"%i)
        projMatch(log,fnVolume,AngularSampling,SymmetryGroup,fnComplementary,WorkingDirStructureExtended,fnAngles,NumberOfMpi)
        
        # Calculate minimum correlation
        fnInliers=fnVolume=os.path.join(WorkingDirStructureCore,"proposedVolume%05d.xmd"%i)
        md=MetaData(fnInliers)
        cc=md.getColumnValues(MDL_MAXCC)
        cc.sort()
        minCC=min(cc[0],CorrThresh)
        
        # Choose those complementary images whose correlation is larger than the minimum
        fnAnglesChosen=os.path.join(WorkingDirStructureExtended,"anglesComplementaryChosen%05d.xmd"%i)
        runJob(log,"xmipp_metadata_utilities",'-i %s --query select "maxCC>%f" -o %s'%(fnAngles,minCC,fnAnglesChosen))
        fnExtended=os.path.join(WorkingDirStructureExtended,"proposedVolume%05d.xmd"%i)
        runJob(log,"xmipp_metadata_utilities",'-i %s --set union %s -o %s'%(fnInliers,fnAnglesChosen,fnExtended))
        deleteFile(log,fnAngles)
        deleteFile(log,fnAnglesChosen)
        fnOutOfCore=os.path.join(WorkingDirStructureExtended,"imagesOutOfCore%05d.xmd"%i)
        runJob(log,"xmipp_metadata_utilities",'-i %s --set subtraction %s -o %s'%(RemainingClasses,fnExtended,fnOutOfCore))
        runJob(log,"xmipp_metadata_utilities",'-i %s --operate sort'%fnOutOfCore)
        runJob(log,"xmipp_metadata_utilities",'-i %s --operate drop_column scale'%fnOutOfCore)
        runJob(log,"xmipp_metadata_utilities",'-i %s --operate sort'%fnExtended)
    
    deleteFile(log,os.path.join(WorkingDirStructureExtended,"gallery.stk"))

def scoreFinalVolumes(log,WorkingDir,NumVolumes):
    mdOut=MetaData()
    for n in range(NumVolumes):
        fnRoot=os.path.join(WorkingDir,'Structure00000_core_extended','proposedVolume%05d'%n)
        fnAssignment=fnRoot+".xmd"
        if exists(fnAssignment):
            MDassignment=MetaData(fnAssignment)
            sum=0
            N=0
            minCC=2
            for id in MDassignment:
                cc=MDassignment.getValue(MDL_MAXCC,id)
                sum+=cc
                if cc<minCC:
                    minCC=cc
                N+=1
            avg=sum/N
            id=mdOut.addObject()
            mdOut.setValue(MDL_IMAGE,fnRoot+".vol",id)
            mdOut.setValue(MDL_VOLUME_SCORE_SUM,float(sum),id)
            mdOut.setValue(MDL_VOLUME_SCORE_MEAN,float(avg),id)
            mdOut.setValue(MDL_VOLUME_SCORE_MIN,float(minCC),id)
    mdOut.write(os.path.join(WorkingDir,"proposedVolumes.xmd"))
