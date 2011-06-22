
from xmipp import *
import launch_job
import os

def joinImageCTF(_log, dict):
    tmpMD = MetaData(dict['CTFgroupName'])
    MDaux = MetaData(dict['filename_currentAngles'])
    outMD = MetaData()
    outMD.join(MDaux, tmpMD, MDL_IMAGE, NATURAL_JOIN)
    outMD.write(dict['DocFileRef'], MD_APPEND)

#reconstruct
def reconstructVolume(_log, dict):
    #xmipp_reconstruct_fourier -i ctfgroup_1@ctfgroups_Iter_13_current_angles.doc -o rec_ctfg01.vol --sym i3 --weight
    parameters  = ' -i ' +  dict['DocFileRef'] 
    parameters += ' -o ' +  dict['reconstructedVolume'] 
    parameters += ' --sym ' + dict['SymmetryGroup'] +'h'
    parameters += ' --weight'
    doParallel = dict['DoParallel']
    if (doParallel):
            parameters += ' --mpi_job_size ' + dict['MpiJobSize']
            parameters += ' --thr ' + str(dict['NumberOfThreads'])

    launch_job.launch_job('xmipp_reconstruct_fourier',
                             parameters,
                             _log,
                             doParallel,
                             dict['NumberOfMpiProcesses'],
                             dict['NumberOfThreads'],
                             dict['SystemFlavour'])


def maskVolume(_log, dict):
    parameters  = ' -i ' +  dict['reconstructedVolume'] 
    parameters += ' -o ' +  dict['maskReconstructedVolume'] 
    
    parameters += ' --mask raised_crown -%d -%d 2' %(dict['dRradiusMin'], dict['dRradiusMax'])
    launch_job.launch_job("xmipp_transform_mask",
                             parameters,
                             _log,
                             False,1,1,'')
    
    #project
def createProjections(_log, dict):
    #dirname = 'ReferenceLibrary'

    doParallel = dict['DoParallel']

    parameters  = ' -i ' +  dict['maskReconstructedVolume'] 
    parameters += ' --experimental_images ' +  dict['DocFileRef']
    
    parameters += ' -o ' +  dict['referenceStack']
    parameters += ' --sampling_rate ' + dict['AngSamplingRateDeg']
    parameters += ' --sym ' + dict['SymmetryGroup'] +'h'
    parameters += ' --compute_neighbors --near_exp_data ' 
    parameters += ' --angular_distance ' + str(dict['MaxChangeInAngles'])
    if (doParallel):
            parameters += ' --mpi_job_size ' + dict['MpiJobSize']

    launch_job.launch_job('xmipp_angular_project_library',
                             parameters,
                             _log,
                             doParallel,
                             dict['NumberOfMpiProcesses'] *dict['NumberOfThreads'],
                             1,
                             dict['SystemFlavour'])
                


def subtractionScript(_log, dict):
    md = MetaData(dict['DocFileRef'])#experimental images
    #referenceStackName = dict['referenceStack']#reference projection for a given defocus group
    mdRef = MetaData(dict['referenceStackDoc'])#experimental images
    subtractedStackName = dict['subtractedStack']
    
    
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
            
        #print id,str(refNum)+'@' + dict['referenceStack']
        #refImgName = '%06d@%s'%(refNum,referenceStackName)
        imgExp.readApplyGeo(md,       id, False, DATA, ALL_IMAGES,False)
        imgRef.readApplyGeo(mdRef,refNum, False, DATA, ALL_IMAGES,False)
        imgSub = imgExp - imgRef
        imgSub.write('%06d@%s'%(id,subtractedStackName))
        imgExp.write('%06d@%s'%(id,subtractedStackName+'exp'))
        imgRef.write('%06d@%s'%(id,subtractedStackName+'ref'))


        