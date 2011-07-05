
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
    #create blanck file
    outFileName=scaledImages + ".stk"
    if os.path.exists(outFileName):
        os.remove(outFileName)
    dimN=MetaData(filename_currentAngles).size()
    if (dimY<0):
        dimY = dimX
    createEmptyFile(outFileName,dimX,dimY,1,dimN)

    parameters  = ' -i ' +  filename_currentAngles 
    parameters += ' -o ' + outFileName 
    parameters += ' --scale fourier ' + str(dimX) + ' ' + str(dimY) 
    parameters += ' ' +  str(NumberOfThreads) + ' 1 ' 

#    launch_job.launch_job('xmipp_transform_geometry',
#                             parameters,
#                             _log,
#                             False,1,1,'')
    runJob(_log,'xmipp_transform_geometry',
                     parameters,
                     DoParallel,
                     NumberOfMpiProcesses,
                     1, # Threads go in --scale option
                     SystemFlavour)
    #create selfile, the one provided by defualt is buggy when mpi is used
    
    
    scaledImages = scaledImages+".xmd"
    parameters  = ' -p ' +  outFileName + ' -o' + scaledImages + ' -s'

    runJob(_log,'xmipp_metadata_selfile_create',
                     parameters,
                     False,
                     1,
                     1,
                     None)
    
    md = MetaData(filename_currentAngles)
    mdScaled = MetaData(scaledImages)
    
    #md.addLabel("original_image")
    mdScaled.addLabel(MDL_IMAGE_ORIGINAL)
    mdScaled.setValueCol(MDL_IMAGE_ORIGINAL,".")
    
    # Calculate scale factor
    img = Image()
    img.readApplyGeo(md,       1, False, DATA, ALL_IMAGES,False)
    x=y=z=n=0
    (x,y,z,n) = img.getDimensions()
    print "X = ",x 
    print "Y = ",y 
    print "Z = ",z 
    print "N = ",n 
    factorX = x / dimX
    factorY = y / dimY
    # Copy images to original images column
    for id in md:
        image = mdScaled.getValue(MDL_IMAGE, id)
        #print image, id
        mdScaled.setValue(MDL_IMAGE_ORIGINAL, image, id)
        
    mdScaled.write(scaledImages)
    
    print "scaledImages file written"
    #mdScaled.deleteRow(MDL_IMAGE_ORIGINAL)
    
    #parameters  = ' -i ' +  filename_currentAngles 
    #parameters += ' --set merge ' +  scaledImages 
    parameters  = ' -i ' +  scaledImages 
    parameters += ' --set merge ' +  filename_currentAngles 
    parameters += ' -o ' + scaledImages
    
    launch_job.launch_job('xmipp_metadata_utilities',
                             parameters,
                             _log,
                             False,1,1,'')
    
    mdScaled = MetaData(scaledImages)
    scaleXstr = 'shiftX=(shiftX /  ' + str(factorX) + ')'
    scaleYstr = 'shiftY=(shiftY /  ' + str(factorY) + ')'
    print "scaleXstr: ",scaleXstr
    print "scaleYstr: ",scaleYstr
    mdScaled.operate(scaleXstr)
    mdScaled.operate(scaleYstr)
    
    mdScaled.write(scaledImages)
    
    # xmipp_metadata_utilities  -i Iter_6_current_angles.doc --operate modify_values "shiftX=(shiftX/4)" -o Iter_6_current_angles_scaled.doc
    # xmipp_metadata_utilities  -i Iter_6_current_angles_scaled.doc --operate modify_values "shiftY=(shiftY/4)" -o Iter_6_current_angles_scaled.doc

    # ctfgroups sels and current angles, images names, need to be changed to the new scaled names.
    

def joinImageCTF(_log, CTFgroupName,DocFileExp,inputSelfile):
    tmpMD = MetaData(CTFgroupName)
    MDaux = MetaData(inputSelfile)
    outMD = MetaData()
    outMD.join(MDaux, tmpMD, MDL_IMAGE, NATURAL_JOIN)
    outMD.write(DocFileExp, MD_APPEND)


def updateImageLabel(_log, inputSelfile):
    md = MetaData(inputSelfile)
    
    for id in md:
        image = md.getValue(MDL_IMAGE_ORIGINAL, id)
        md.setValue(MDL_IMAGE, image, id)

    md.write(inputSelfile)

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
    parameters += ' --weight'
    doParallel = DoParallel
    if (doParallel):
            parameters += ' --mpi_job_size ' + MpiJobSize
            parameters += ' --thr ' + str(NumberOfThreads)

    launch_job.launch_job('xmipp_reconstruct_fourier',
                             parameters,
                             _log,
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
    launch_job.launch_job("xmipp_transform_mask",
                             parameters,
                             _log,
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
    #dirname = 'ReferenceLibrary'

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

    launch_job.launch_job('xmipp_angular_project_library',
                             parameters,
                             _log,
                             doParallel,
                             NumberOfMpiProcesses *NumberOfThreads,
                             1,
                             SystemFlavour)
                


def subtractionScript(_log
                      ,DocFileExp
                      ,referenceStackDoc
                      ,subtractedStack
                      ):
    md = MetaData(DocFileExp)#experimental images
    #referenceStackName = referenceStack#reference projection for a given defocus group
    mdRef = MetaData(referenceStackDoc)#experimental images
    subtractedStackName = subtractedStack
    
    #Set shifts to 0
    a = md.getValue(MDL_SHIFTX, 1)
    print "antes  : ", a
    md.operate("shiftX=(shiftX*-1)")
    b = md.getValue(MDL_SHIFTX, 1)
    print "despues: ", b
    
    md.operate("shiftY=(shiftY*-1)")
    
    mdRef.setValueCol(MDL_SHIFTX,0.)
    mdRef.setValueCol(MDL_SHIFTY,0.)
    
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
        #read only
        imgRef.readApplyGeo(mdRef,refNum, False, DATA, ALL_IMAGES,False)
        imgSub = imgExp - imgRef
        imgSub.write('%06d@%s'%(id,subtractedStackName))
        imgExp.write('%06d@%s'%(id,subtractedStackName+'exp'))
        imgRef.write('%06d@%s'%(id,subtractedStackName+'ref'))

