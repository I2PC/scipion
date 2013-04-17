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
from xmipp import MetaData, MetaDataInfo, MDL_IMAGE, MDL_ANGLE_ROT, MDL_ANGLE_TILT, MDL_ANGLE_PSI, MDL_COST

from protlib_base import *
from protlib_utils import getListFromRangeString, runJob, runShowJ
from protlib_filesystem import copyFile, deleteFile

class ProtSymmetric(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.symmetric_initial.name, scriptname, project)
        self.Import = 'from protocol_symmetric_initial import *'
        (self.Xdim,_,_,_,_) = MetaDataInfo(self.SideViews)
        
    def defineSteps(self):
        self.insertStep('createDir',path=self.ExtraDir)
        self.insertStep("linkAcquisitionInfo",InputFile=self.SideViews,dirDest=self.WorkingDir)
        
        # Align all side views
        dirSide=self.extraPath('sideViews')
        self.insertStep('createDir',path=dirSide)
        fnAlignedSide=os.path.join(dirSide,'level_00/class_classes.stk')
        self.insertRunJobStep("xmipp_classify_CL2D","-i %s --odir %s --nref0 1 --nref 1 --iter 10"%(self.SideViews,dirSide),[fnAlignedSide])
        self.insertStep("deleteFile",filename=os.path.join(dirSide,'images.xmd'))
        
        fnAlignedTop=""
        if self.TopViews!="":
            # Align all top views
            dirTop=self.extraPath('topViews')
            self.insertStep('createDir',path=dirTop)
            fnAlignedTop=os.path.join(dirTop,'level_00/class_classes.stk')
            self.insertRunJobStep("xmipp_classify_CL2D","-i %s --odir %s --nref0 1 --nref 1 --iter 10"%(self.TopViews,dirTop),[fnAlignedTop])
            self.insertStep("deleteFile",filename=os.path.join(dirTop,'images.xmd'))

        # Try to find the in-plane orientation of the average side view
        fnPsiProfile=self.extraPath('psiProfile.xmd')
        self.insertStep('calculatePsiProfile',verifyfiles=[fnPsiProfile],
                        avgSideView=fnAlignedSide,avgTopView=fnAlignedTop,ExtraDir=self.ExtraDir,Symmetry=self.SymmetryGroup,Nproc=self.NumberOfMpi)

        # Construct mask
        fnMask=self.extraPath('cylindricalMask.vol')
        self.insertStep('constructCylindricalMask',verifyfiles=[fnMask],
                        ExtraDir=self.ExtraDir,avgSideView=fnAlignedSide,avgTopView=fnAlignedTop,radius=self.Xdim/2)

        # Try to find the in-plane orientation of the average side view
        # Construct bestFirstVolume.vol
        fnFirstVolume=self.extraPath('volume.stk')
        self.insertStep('constructFirstVolume',verifyfiles=[fnFirstVolume],
                        avgSideView=fnAlignedSide,avgTopView=fnAlignedTop,fnMask=fnMask,
                        ExtraDir=self.ExtraDir,Symmetry=self.SymmetryGroup,Nproc=self.NumberOfMpi)
                
        # Now align all side views with respect to this volume
        Niter=5
        for i in range(2,Niter+1):
            self.insertStepsAlignment(i)
        
        # Take the last volume and make it the final result
        self.insertStep("runJob",programname="xmipp_image_convert",
                        params="-i %d@%s -o %s"%(Niter,self.extraPath('volume.stk'),self.workingDirPath('volume.vol')))

    def insertStepsAlignment(self,i):
        fnVolumeStack=self.extraPath('volume.stk')
        fnAngles=self.extraPath('angles.xmd')
        fnSurface=self.extraPath('cylindricalMaskInverted.vol')
        
        # Generate projection gallery for side views
        fnPreviousVolume='%d@%s'%(i-1,fnVolumeStack)
        fnGallery=self.extraPath('gallery.stk')
        self.insertRunJobStep("xmipp_angular_project_library", "-i %s -o %s --sampling_rate 2 --sym %s --method fourier 1 0.25 bspline --min_tilt_angle 88 --max_tilt_angle 92 --compute_neighbors --angular_distance -1 --experimental_images %s"\
                              %(fnPreviousVolume,fnGallery,self.SymmetryGroup,self.SideViews),[fnGallery])

        # Assign angles for side views
        fnAnglesSide='angles_side_%02d@'%i+fnAngles
        self.insertRunJobStep("xmipp_angular_projection_matching", "-i %s -o %s --ref %s --Ri 0 --Ro %s --max_shift 1000 --search5d_shift %s --search5d_step  1 --append"\
                              %(self.SideViews,fnAnglesSide,fnGallery,str(self.Xdim/2),str(self.Xdim/10)),[fnAnglesSide])  
        fnAnglesToReconstruct=fnAnglesSide

        if self.TopViews!="":
            # Generate projection gallery for top views
            self.insertRunJobStep("xmipp_angular_project_library", "-i %s -o %s --sampling_rate 2 --sym %s --method fourier 1 0.25 bspline --min_tilt_angle 0 --max_tilt_angle 2 --compute_neighbors --angular_distance -1 --experimental_images %s"\
                                  %(fnPreviousVolume,fnGallery,self.SymmetryGroup,self.SideViews),[fnGallery])
    
            # Assign angles for top views
            fnAnglesTop='angles_top_%02d@'%i+fnAngles
            self.insertRunJobStep("xmipp_angular_projection_matching", "-i %s -o %s --ref %s --Ri 0 --Ro %s --max_shift 1000 --search5d_shift %s --search5d_step  1 --append"\
                                  %(self.TopViews,fnAnglesTop,fnGallery,str(self.Xdim/2),str(self.Xdim/10)),[fnAnglesTop])
            fnAnglesToReconstruct='angles_%02d@'%i+fnAngles
            self.insertStep("runJob",programname="xmipp_metadata_utilities",params="-i %s --set union %s -o %s --mode append"%(fnAnglesSide,fnAnglesTop,fnAnglesToReconstruct))

        # Reconstruct
        fnCurrentVolume=self.extraPath('volume%02d.vol'%i)
        self.insertRunJobStep("xmipp_reconstruct_art","-i %s -o %s --sym %s -n 15 --surface %s"%(fnAnglesToReconstruct,fnCurrentVolume,self.SymmetryGroup,fnSurface))
        # Fourier: self.insertRunJobStep("xmipp_reconstruct_fourier","-i %s -o %s --sym %s"%(fnAnglesToReconstruct,fnCurrentVolume,self.SymmetryGroup))
        self.insertStep("runJob",programname="xmipp_transform_symmetrize",params="-i %s --sym %s"%(fnCurrentVolume,self.SymmetryGroup))
        self.insertStep("runJob",programname="xmipp_transform_mask",params="-i %s --mask circular -%s"%(fnCurrentVolume,str(self.Xdim/2)))
        fnMask=self.extraPath('cylindricalMask.vol')
        self.insertStep("runJob",programname="xmipp_transform_mask",params="-i %s --mask binary_file %s"%(fnCurrentVolume,fnMask))

        # Delete intermediate files
        self.insertStep("runJob",programname="rm",params="%s/gallery*"%(self.ExtraDir))
        self.insertStep("runJob",programname="rm",params="%s/volume*.hist"%(self.ExtraDir))
        
        # Stack reconstruction
        self.insertStep("runJob",programname="xmipp_image_convert",params="-i %s -o %02d@%s"%(fnCurrentVolume,i,fnVolumeStack))
        self.insertStep('deleteFile',filename=fnCurrentVolume)
        
    def summary(self):
        message = []
        message.append("Top  views: [%s] " % self.TopViews)
        message.append("Side views: [%s] " % self.SideViews)
        message.append("Symmetry: %s " % self.SymmetryGroup)
        return message
    
    def validate(self):
        errors = []
        return errors    

    def visualize(self):
        fnVol=self.workingDirPath('volume.vol')
        if os.path.exists():
            runShowJ(fnVol)
        return

def getReconstructionError(fnHist):
    error=1e38
    for line in open(fnHist):
        if "Global mean squared error:" in line:
            error=float(line.split(':')[1])
    return error

def calculatePsiProfile(log,avgSideView,avgTopView,ExtraDir,Symmetry,Nproc):
    MD=MetaData()
    if avgTopView!="":
        objId=MD.addObject()
        MD.setValue(MDL_IMAGE,avgTopView,objId)
        MD.setValue(MDL_ANGLE_ROT,0.,objId)
        MD.setValue(MDL_ANGLE_TILT,0.,objId)
        MD.setValue(MDL_ANGLE_PSI,0.,objId)
    objId=MD.addObject()
    MD.setValue(MDL_IMAGE,avgSideView,objId)
    MD.setValue(MDL_ANGLE_ROT,0.,objId)
    MD.setValue(MDL_ANGLE_TILT,90.,objId)
    psi=0.
    fnInitial=os.path.join(ExtraDir,'testPsi.xmd')
    fnVolume=os.path.join(ExtraDir,'testReconstruction.vol')
    fnHist=os.path.join(ExtraDir,'testReconstruction.hist')
    MDpsi=MetaData()
    while psi<180.:
        MD.setValue(MDL_ANGLE_PSI,psi,objId)
        MD.write(fnInitial)
        runJob(log,"xmipp_reconstruct_art","-i %s -o %s -n 3 --sym %s"%(fnInitial,fnVolume,Symmetry),Nproc)
        error=getReconstructionError(fnHist)
        objIdPsi=MDpsi.addObject()
        MDpsi.setValue(MDL_ANGLE_PSI,psi,objIdPsi)
        MDpsi.setValue(MDL_COST,error,objIdPsi)
        MDpsi.write(os.path.join(ExtraDir,'psiProfile.xmd'))
        psi+=3.
    deleteFile(log,fnVolume)
    deleteFile(log,fnInitial)
    runJob(log,"rm","%s/*.hist"%ExtraDir)

def getBestPsi(ExtraDir):
    MDpsi=MetaData(os.path.join(ExtraDir,'psiProfile.xmd'))
    bestError=1e38
    bestPsi=0
    for objId in MDpsi:
        error=MDpsi.getValue(MDL_COST,objId)
        if error<bestError:
            bestError=error
            bestPsi=MDpsi.getValue(MDL_ANGLE_PSI,objId)
    return float(bestPsi)

def constructCylindricalMask(log,ExtraDir,avgSideView,avgTopView,radius):
    bestPsi=getBestPsi(ExtraDir)
    
    # Construct metadata with two images
    MD=MetaData()
    objId=MD.addObject()
    MD.setValue(MDL_IMAGE,avgSideView,objId)
    MD.setValue(MDL_ANGLE_ROT,0.,objId)
    MD.setValue(MDL_ANGLE_TILT,90.,objId)
    MD.setValue(MDL_ANGLE_PSI,bestPsi,objId)
    if avgTopView!="":
        objId=MD.addObject()
        MD.setValue(MDL_IMAGE,avgTopView,objId)
        MD.setValue(MDL_ANGLE_ROT,0.,objId)
        MD.setValue(MDL_ANGLE_TILT,0.,objId)
        MD.setValue(MDL_ANGLE_PSI,0.,objId)
    fnAngles=os.path.join(ExtraDir,'anglesMask.xmd')
    MD.write(fnAngles)

    # Construct symmetry file
    fnSym=os.path.join(ExtraDir,'cylindricalSymmetry.sym')
    fhSym = open(fnSym, "w")
    fhSym.write("rot_axis 120 0 0 1")
    fhSym.close()
    
    # Construct cylindrical mask
    fnMask=os.path.join(ExtraDir,'cylindricalMask.vol')
    runJob(log,"xmipp_reconstruct_fourier","-i %s -o %s --sym %s"%(fnAngles,fnMask,fnSym))
    deleteFile(log,fnAngles)
    deleteFile(log,fnSym)
    
    # Binarize mask
    runJob(log,"xmipp_transform_range_adjust","-i %s --range 0 1"%fnMask)
    runJob(log,"xmipp_transform_threshold","-i %s --select below 0.5 --substitute value 0"%fnMask)
    runJob(log,"xmipp_transform_threshold","-i %s --select above 0.5 --substitute value 1"%fnMask)

    # Process it
    runJob(log,"xmipp_transform_mask","-i %s --mask circular %s"%(fnMask,str(-(radius-3))))
    runJob(log,"xmipp_transform_morphology","-i %s --binaryOperation opening --size 2"%fnMask)
    runJob(log,"xmipp_transform_morphology","-i %s --binaryOperation closing --size 2"%fnMask)
    runJob(log,"xmipp_transform_morphology","-i %s --binaryOperation dilation"%fnMask)
    
    # Invert it
    fnMaskInverted=os.path.join(ExtraDir,'cylindricalMaskInverted.vol')
    runJob(log,"xmipp_image_operate","-i %s --mult -1 -o %s"%(fnMask,fnMaskInverted))
    runJob(log,"xmipp_image_operate","-i %s --plus 1 "%fnMaskInverted)

def constructFirstVolume(log,avgSideView,avgTopView,fnMask,ExtraDir,Symmetry,Nproc):
    bestPsi=getBestPsi(ExtraDir)
    fnInitial=os.path.join(ExtraDir,'initial.xmd')
    fnVolume=os.path.join(ExtraDir,'firstReconstruction.vol')
    fnHist=os.path.join(ExtraDir,'firstReconstruction.hist')
    fnSurface=os.path.join(ExtraDir,'cylindricalMaskInverted.vol')

    MD=MetaData()
    if avgTopView!="":
        objId=MD.addObject()
        MD.setValue(MDL_IMAGE,avgTopView,objId)
        MD.setValue(MDL_ANGLE_ROT,0.,objId)
        MD.setValue(MDL_ANGLE_TILT,0.,objId)
        MD.setValue(MDL_ANGLE_PSI,0.,objId)
    objId=MD.addObject()
    MD.setValue(MDL_IMAGE,avgSideView,objId)
    MD.setValue(MDL_ANGLE_ROT,0.,objId)
    MD.setValue(MDL_ANGLE_TILT,90.,objId)
    MD.setValue(MDL_ANGLE_PSI,bestPsi,objId)
    MD.write(fnInitial)

    runJob(log,"xmipp_reconstruct_art","-i %s -o %s -n 25 --sym %s --thr %d --surface %s"\
           %(fnInitial,fnVolume,Symmetry,Nproc,fnSurface ))
    # Fourier: runJob(log,"xmipp_reconstruct_fourier","-i %s -o %s --sym %s"%(fnInitial,fnVolume,Symmetry))
    deleteFile(log,fnHist)
    deleteFile(log,fnInitial)

    runJob(log,"xmipp_transform_mask","-i %s -o 1@%s/volume.stk --mask binary_file %s"%(fnVolume,ExtraDir,fnMask))
    deleteFile(log,fnVolume)

    runJob(log,"xmipp_transform_symmetrize","-i 1@%s/volume.stk --sym %s"%(ExtraDir,Symmetry))
