from xmipp import *
import os
from protlib_utils import runJob

#No defocus group
def joinImageCTFscale(_log, CTFgroupName,DocFileExp,inputSelfile):
    ctfMD = MetaData(CTFgroupName)
    MDaux = MetaData(inputSelfile)
    outMD = MetaData()
    outMD.join(MDaux, ctfMD, MDL_IMAGE_ORIGINAL, MDL_IMAGE, INNER_JOIN)
    outMD.write(DocFileExp, MD_APPEND)

#No defocus group
def joinImageCTF(_log, CTFgroupName,DocFileExp,inputSelfile):
    ctfMD = MetaData(CTFgroupName)
    MDaux = MetaData(inputSelfile)
    outMD = MetaData()
    outMD.join(MDaux, ctfMD, MDL_IMAGE, MDL_IMAGE, INNER_JOIN)
    outMD.write(DocFileExp, MD_APPEND)

#reconstruct
def reconstructVolume(_log 
                     ,DocFileExp
                     , DoParallel
                     , MpiJobSize
                     , NumberOfMpiProcesses
                     , NumberOfThreads
                     , reconstructedVolume
                     , SymmetryGroup
                     , SystemFlavour
                     ):
    #xmipp_reconstruct_fourier -i ctfgroup_1@ctfgroups_Iter_13_current_angles.doc -o rec_ctfg01.vol --sym i3 --weight
    parameters  = ' -i ' +  DocFileExp 
    parameters += ' -o ' +  reconstructedVolume 
    parameters += ' --sym ' + SymmetryGroup +'h'
    doParallel = DoParallel
    if (doParallel):
            parameters += ' --mpi_job_size ' + MpiJobSize
            parameters += ' --thr ' + str(NumberOfThreads)

    runJob(_log,'xmipp_reconstruct_fourier',
                             parameters,
                             doParallel,
                             NumberOfMpiProcesses,
                             NumberOfThreads,
                             SystemFlavour)

def maskVolume(_log
              ,dRradiusMax
              ,dRradiusMin
              ,maskReconstructedVolume
              ,reconstructedVolume
):
    parameters  = ' -i ' +  reconstructedVolume 
    parameters += ' -o ' +  maskReconstructedVolume 
    parameters += ' --mask raised_crown -%d -%d 2' %(dRradiusMin, dRradiusMax)

    runJob( _log,"xmipp_transform_mask",
                             parameters,
                             False,1,1,'')
    
    #project
def createProjections(_log
                      ,AngSamplingRateDeg
                      ,DocFileExp
                      ,DoParallel
                      ,maskReconstructedVolume
                      ,MaxChangeInAngles
                      , MpiJobSize
                      ,NumberOfMpiProcesses
                      ,NumberOfThreads
                      ,referenceStack
                      ,SymmetryGroup
                      ,SystemFlavour
):
    doParallel = DoParallel

    parameters  = ' -i ' +  maskReconstructedVolume 
    parameters += ' --experimental_images ' +  DocFileExp
    
    parameters += ' -o ' +  referenceStack
    parameters += ' --sampling_rate ' + AngSamplingRateDeg
    parameters += ' --sym ' + SymmetryGroup +'h'
    parameters += ' --compute_neighbors --near_exp_data ' 
    parameters += ' --angular_distance ' + str(MaxChangeInAngles)
    if (doParallel):
            parameters += ' --mpi_job_size ' + MpiJobSize

    runJob(_log,'xmipp_angular_project_library',
                             parameters,
                             doParallel,
                             NumberOfMpiProcesses *NumberOfThreads,
                             1,
                             SystemFlavour)
                


def subtractionScript(_log
                      ,DocFileExp
                      ,referenceStackDoc
                      ,subtractedStack
                      ):
    #printLog("subtractionScript, docfile %s referenceStackDoc %s substracted stack %s"%
    #                   (DocFileExp
    #                  ,referenceStackDoc
    #                  ,subtractedStack))

    md = MetaData(DocFileExp)#experimental images
    mdRotations = MetaData(md) #rotations
    
    # Save Metadata with just rotations (shifts will be applied when reading)
    mdRotations.setValueCol(MDL_SHIFTX, 0.)
    mdRotations.setValueCol(MDL_SHIFTY, 0.)
    
    #reference projection for a given defocus group
    mdRef = MetaData(referenceStackDoc) #reference library
    subtractedStackName = subtractedStack #output
    
    imgExp = Image()
    imgRef = Image()
    imgSub = Image()
    #apply ctf to a temporary reference    
    for id in md:
        #refNum = md.getValue(MDL_REF, id)
        angRot = md.getValue(MDL_ANGLEROT, id)
        angTilt = md.getValue(MDL_ANGLETILT, id)
        psi = md.getValue(MDL_ANGLEPSI, id)
        
        # Search for the closest idRef
        dist = -1.
        distMin = 999.
        for idRef in mdRef:
            angRotRef  = mdRef.getValue(MDL_ANGLEROT, idRef)
            angTiltRef = mdRef.getValue(MDL_ANGLETILT, idRef)
            
            dist = abs(float(angRotRef) - float(angRot)) +  abs(float(angTiltRef) - float(angTilt))
            if(dist < distMin or dist == -1):
                refNum = idRef
                distMin = dist
            
        #print id,str(refNum)+'@' + referenceStack
        #refImgName = '%06d@%s'%(refNum,referenceStackName)
        
        # Apply alignment as follows: shifts firt (while reading), then rotations
        print "step 1"
        imgExp.readApplyGeo(md, id, True, DATA, ALL_IMAGES, False) # only_apply_shifts = true
        print "step 2"
        imgExp.write("image_tmp.spi")
        #imgExp = Image("image_tmp.spi")
        imgExp.applyGeo(mdRotations, id, False, False)  
        print "step 3"
        #void applyGeo(const MetaData &md, size_t objId, bool only_apply_shifts = false, bool wrap = WRAP)
        
        imgRef.readApplyGeo(mdRef,refNum, False, DATA, ALL_IMAGES,False)
        imgSub = imgExp - imgRef
        imgSub.write('%06d@%s'%(id,subtractedStackName))
        imgExp.write('%06d@%s'%(id,subtractedStackName+'exp'))
        imgRef.write('%06d@%s'%(id,subtractedStackName+'ref'))

