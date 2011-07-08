
from xmipp import *
import launch_job
import os
from protlib_utils import runJob

#'([A-z]*)'[ ]*:(.*)\n
#$1\n

def scaleImages(_log
                    , dimX
                    , dimY
                    , DoParallel
                    , filename_currentAngles
                    , MpiJobSize
                    , NumberOfMpiProcesses
                    , NumberOfThreads
                    , scaledImages
                    , SystemFlavour
                    ):

    outFileName=scaledImages + ".stk"
    if os.path.exists(outFileName):
        os.remove(outFileName)

    if (dimY<0):
        dimY = dimX
    
    parameters  = ' -i ' +  filename_currentAngles 
    parameters += ' -o ' + outFileName 
    parameters += ' --scale fourier ' + str(dimX) + ' ' + str(dimY) + ' ' + str(NumberOfThreads)
    parameters += ' --disable_metadata' 

    runJob(_log,'xmipp_transform_geometry',
                     parameters,
                     DoParallel,
                     NumberOfMpiProcesses,
                     1, # Threads go in --scale option
                     SystemFlavour)

    #merge new scaled metadata with original metadata
    fnScaledImages = scaledImages+".xmd"
    mdScaled = MetaData(fnScaledImages)
    md = MetaData(filename_currentAngles)
    
    # copy images to original images column
    md.addLabel(MDL_IMAGE_ORIGINAL)
    for id in md:
        imageScale=mdScaled.getValue(MDL_IMAGE, id)
        imageOriginal=md.getValue(MDL_IMAGE, id)
        md.setValue(MDL_IMAGE, imageScale, id)
        md.setValue(MDL_IMAGE_ORIGINAL, imageOriginal, id)
    
    # Calculate scale factor
    (x,y,z,n) =ImgSize(filename_currentAngles)
    factorX = float(x) / dimX
    factorY = float(y) / dimY

    scaleXstr = 'shiftX=(shiftX /  ' + str(factorX) + ')'
    scaleYstr = 'shiftY=(shiftY /  ' + str(factorY) + ')'
    md.operate(scaleXstr)
    md.operate(scaleYstr)
    
    md.write(fnScaledImages)
    

def joinImageCTF(_log, CTFgroupName,DocFileExp,inputSelfile):
    ctfMD = MetaData(CTFgroupName)
    MDaux = MetaData(inputSelfile)
    outMD = MetaData()
    outMD.join(MDaux, ctfMD, MDL_IMAGE_ORIGINAL, MDL_IMAGE, INNER_JOIN)
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
    printLog("subtractionScript, docfile %s referenceStackDoc %s substracted stack %s"%
                       (DocFileExp
                      ,referenceStackDoc
                      ,subtractedStack))
    md = MetaData(DocFileExp)#experimental images
    #referenceStackName = referenceStack#reference projection for a given defocus group
    mdRef = MetaData(referenceStackDoc)#refrence ibrary
    subtractedStackName = subtractedStack#output
    
    #Set shifts to 0
    md.operate("shiftX=(shiftX*-1)")
    md.operate("shiftY=(shiftY*-1)")
    md.operate("anglePsi=(anglePsi*-1)")
    
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
        imgExp.readApplyGeo(md,       id, False, DATA, ALL_IMAGES,False)
        #read only
        imgRef.readApplyGeo(mdRef,refNum, False, DATA, ALL_IMAGES,False)
        imgSub = imgExp - imgRef
        imgSub.write('%06d@%s'%(id,subtractedStackName))
        imgExp.write('%06d@%s'%(id,subtractedStackName+'exp'))
        imgRef.write('%06d@%s'%(id,subtractedStackName+'ref'))

