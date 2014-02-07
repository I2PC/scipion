#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Xmipp protocol for subtraction
#
# Example use:
# ./protocol_subtraction_header.py
#
# Authors: Roberto Marabini,
#          Alejandro Echeverria Rey.
#
from xmipp import *
import os
from protlib_utils import runJob

#No defocus group
def joinImageCTFscale(_log, CTFgroupName,DocFileExp,inputSelfile):
    ctfMD = MetaData(CTFgroupName)
    MDaux = MetaData(inputSelfile)
    outMD = MetaData()
    outMD.join2(MDaux, ctfMD, MDL_IMAGE_ORIGINAL, MDL_IMAGE, INNER_JOIN)
    outMD.write(DocFileExp, MD_APPEND)

#No defocus group
def joinImageCTF(_log, CTFgroupName,DocFileExp,inputSelfile):
    ctfMD = MetaData(CTFgroupName)
    MDaux = MetaData(inputSelfile)
    outMD = MetaData()
    outMD.join1(MDaux, ctfMD, MDL_IMAGE, INNER_JOIN)
    outMD.write(DocFileExp, MD_APPEND)

#reconstruct
def reconstructVolume(_log 
                     ,DocFileExp
                     , MpiJobSize
                     , NumberOfMpi
                     , NumberOfThreads
                     , reconstructedVolume
                     , SymmetryGroup
                     ):
    #xmipp_reconstruct_fourier -i ctfgroup_1@ctfgroups_Iter_13_current_angles.doc -o rec_ctfg01.vol --sym i3 --weight
    parameters  = ' -i ' +  DocFileExp 
    parameters += ' -o ' +  reconstructedVolume 
    parameters += ' --sym ' + SymmetryGroup +'h'
    if (NumberOfMpi>1):
            parameters += ' --mpi_job_size ' + MpiJobSize
            parameters += ' --thr ' + str(NumberOfThreads)

    runJob(_log,'xmipp_reconstruct_fourier',
                             parameters,
                             NumberOfMpi,
                             NumberOfThreads)

def maskVolume(_log
              ,dRradiusMax
              ,dRradiusMin
              ,maskReconstructedVolume
              ,reconstructedVolume
              ,NumberOfMpi
              ,NumberOfThreads
):
    parameters  = ' -i ' +  reconstructedVolume 
    parameters += ' -o ' +  maskReconstructedVolume 
    parameters += ' --mask raised_crown -%d -%d 2' %(dRradiusMin, dRradiusMax)

    runJob( _log,"xmipp_transform_mask",
                             parameters,
                             NumberOfMpi *NumberOfThreads)
    
    #project
def createProjections(_log
                      ,AngSamplingRateDeg
                      ,DocFileExp
                      ,maskReconstructedVolume
                      ,MaxChangeInAngles
                      , MpiJobSize
                      ,NumberOfMpi
                      ,NumberOfThreads
                      ,referenceStack
                      ,SymmetryGroup
):

    parameters  = ' -i ' +  maskReconstructedVolume 
    parameters += ' --experimental_images ' +  DocFileExp
    
    parameters += ' -o ' +  referenceStack
    parameters += ' --sampling_rate ' + AngSamplingRateDeg
    parameters += ' --sym ' + SymmetryGroup +'h'
    parameters += ' --compute_neighbors --near_exp_data ' 
    parameters += ' --angular_distance ' + str(MaxChangeInAngles)
    if ((NumberOfMpi *NumberOfThreads)>1):
            parameters += ' --mpi_job_size ' + MpiJobSize

    runJob(_log,'xmipp_angular_project_library',
                             parameters,
                             NumberOfMpi *NumberOfThreads)
                


def subtractionScript(_log
                      ,DocFileExp
                      ,referenceStackDoc
                      ,subtractedStack
                      ,resultsImagesName
                      ):
    md = MetaData(DocFileExp)#experimental images
    mdRotations = MetaData(md) #rotations
    
    # Save Metadata with just rotations (shifts will be applied when reading)
    mdRotations.setValueCol(MDL_SHIFT_X, 0.)
    mdRotations.setValueCol(MDL_SHIFT_Y, 0.)
    mdRotations.operate('anglePsi=-anglePsi')
    
    mdResults = MetaData(mdRotations)
    mdResults.operate('anglePsi=0')
    
    #reference projection for a given defocus group
    mdRef = MetaData(referenceStackDoc) #reference library
    
    imgExp = Image()
    imgRef = Image()
    imgSub = Image()
    
    stackResults = resultsImagesName+'.stk'
    
    if (os.path.exists(stackResults)):
        from protlib_xmipp import getMdSize
        idResults = getMdSize(stackResults)
    else:
        idResults = 0
        
    xmdResults = resultsImagesName+'.xmd'
    if (os.path.exists(xmdResults)):
        xmdOldResults = MetaData(xmdResults)
        xmdOldResults.unionAll(mdResults)
        mdResults = MetaData(xmdOldResults)
    
    for id in md:
        angRot = md.getValue(MDL_ANGLE_ROT, id)
        angTilt = md.getValue(MDL_ANGLE_TILT, id)
        psi = md.getValue(MDL_ANGLE_PSI, id)
        
        # Search for the closest idRef
        dist = -1.
        distMin = 999.
        for idRef in mdRef:
            angRotRef  = mdRef.getValue(MDL_ANGLE_ROT, idRef)
            angTiltRef = mdRef.getValue(MDL_ANGLE_TILT, idRef)
            
            dist = abs(float(angRotRef) - float(angRot)) +  abs(float(angTiltRef) - float(angTilt))
            if(dist < distMin or dist == -1):
                refNum = idRef
                distMin = dist
                    
        # Apply alignment as follows: shifts firt (while reading), then rotations
        imgExp.readApplyGeo(md, id, True, DATA, ALL_IMAGES, False) # only_apply_shifts = true
        imgExp.applyGeo(mdRotations, id, False, False)  
        
        imgRef.readApplyGeo(mdRef,refNum, False, DATA, ALL_IMAGES,False)
        imgSub = imgExp - imgRef
        imgSub.write('%06d@%s'%(id,subtractedStack))
        imgExp.write('%06d@%s'%(id,subtractedStack+'exp'))
        imgRef.write('%06d@%s'%(id,subtractedStack+'ref'))
        
        mdRotations.setValue(MDL_IMAGE, '%06d@%s'%(id,subtractedStack), id)
        
        # write protocols results
        imgSub.write('%06d@%s'%(id + idResults, stackResults))
        mdResults.setValue(MDL_IMAGE, '%06d@%s'%(id + idResults, stackResults), id + idResults)
        
    
    mdRotations.operate('anglePsi=0')
    mdRotations.write(subtractedStack.replace('.stk','.xmd'))
    
    mdResults.write(xmdResults)
        

