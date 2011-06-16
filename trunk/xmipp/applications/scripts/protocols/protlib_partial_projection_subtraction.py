
from xmipp import *
import launch_job
import os

def joinImageCTF(_log, CTFgroupName, DocFileRef):
    tmpMD = MetaData(CTFgroupName)
    MDaux = MetaData(filename_currentAngles)
    outMD = MetaData()
    outMD.join(MDaux, tmpMD, MDL_IMAGE, NATURAL_JOIN)
    outMD.write(DocFileRef, MD_APPEND)


#reconstruct
def reconstructVolume(_log, dict):
    #xmipp_reconstruct_fourier -i ctfgroup_1@ctfgroups_Iter_13_current_angles.doc -o rec_ctfg01.vol --sym i3 --weight
    parameters  = ' -i ' +  DocFileRef 
    parameters += ' -o ' +  reconstructedVolume 
    parameters += ' --sym ' + SymmetryGroup +'h'
    parameters += ' --weight'
    doParallel = DoParallel
    if (doParallel):
            parameters += ' --mpi_job_size ' + MpiJobSize
            parameters += ' --thr ' + str(NumberOfThreads)

    launchJob('xmipp_reconstruct_fourier',
                             parameters,
                             _log,
                             doParallel,
                             NumberOfMpiProcesses,
                             NumberOfThreads,
                             SystemFlavour)


def maskVolume(_log, dict):
    parameters  = ' -i ' +  reconstructedVolume 
    parameters += ' -o ' +  maskReconstructedVolume 
    
    parameters += ' --mask raised_crown -%d -%d 2' %(dRradiusMin, dRradiusMax)
    launchJob("xmipp_transform_mask",
                             parameters,
                             _log,
                             False,1,1,'')
    
    #project
def createProjections(_log, dict):
    #dirname = 'ReferenceLibrary'

    doParallel = DoParallel

    parameters  = ' -i ' +  maskReconstructedVolume 
    parameters += ' --experimental_images ' +  DocFileRef
    
    parameters += ' -o ' +  referenceStack
    parameters += ' --sampling_rate ' + AngSamplingRateDeg
    parameters += ' --sym ' + SymmetryGroup +'h'
    parameters += ' --compute_neighbors --near_exp_data ' 
    parameters += ' --angular_distance ' + str(MaxChangeInAngles)
    if (doParallel):
            parameters += ' --mpi_job_size ' + MpiJobSize

    launchJob('xmipp_angular_project_library',
                             parameters,
                             _log,
                             doParallel,
                             NumberOfMpiProcesses *NumberOfThreads,
                             1,
                             SystemFlavour)
                


def subtractionScript(_log, dict):
    md = MetaData(DocFileRef)#experimental images
    #referenceStackName = referenceStack#reference projection for a given defocus group
    mdRef = MetaData(referenceStackDoc)#experimental images
    subtractedStackName = subtractedStack
    
    
    imgExp = Image()
    imgRef = Image()
    imgSub = Image()
    #apply ctf to a temporary reference    
    for id in md:
        #refNum = md.getValue(MDL_REF, id)
        angRot = md.getValue(MDL_ANGLEROT, id)
        angTilt = md.getValue(MDL_ANGLETILT, id)
        psi = md.getValue(MDL_ANGLEPSI, id)
        md.setValue(MDL_ANGLEPSI, psi * -1.,id)
        
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
        imgRef.readApplyGeo(mdRef,refNum, False, DATA, ALL_IMAGES,False)
        imgSub = imgExp - imgRef
        imgSub.write('%06d@%s'%(id,subtractedStackName))
        imgExp.write('%06d@%s'%(id,subtractedStackName+'exp'))
        imgRef.write('%06d@%s'%(id,subtractedStackName+'ref'))


        
