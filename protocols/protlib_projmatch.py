from xmipp import MetaData, MetaDataInfo, MDL_IMAGE, MDL_IMAGE1, MDL_IMAGE_REF, MDL_ANGLE_ROT, MDL_ANGLE_TILT, MDL_ANGLE_PSI, MDL_REF, \
        MDL_SHIFT_X, MDL_SHIFT_Y, MDL_FLIP, MD_APPEND, MDL_MAXCC, MDL_ENABLED, MDL_CTF_MODEL, MDL_SAMPLINGRATE, DT_DOUBLE, \
        Euler_angles2matrix, Image, FileName, getBlocksInMetaDataFile
from protlib_utils import runJob
from protlib_filesystem import deleteFile, findAcquisitionInfo
import os

def projMatch(log,Volume,AngularSampling,SymmetryGroup,Images,ExtraDir,fnAngles,NumberOfMpi):
    Xdim=MetaDataInfo(Images)[0]

    # Generate gallery of projections        
    fnGallery=os.path.join(ExtraDir,'gallery.stk')
    runJob(log,"xmipp_angular_project_library", "-i %s -o %s --sampling_rate %f --sym %s --method fourier 1 0.25 bspline --compute_neighbors --angular_distance -1 --experimental_images %s"\
               %(Volume,fnGallery,float(AngularSampling),SymmetryGroup,Images),NumberOfMpi)

    # Assign angles
    runJob(log,"xmipp_angular_projection_matching", "-i %s -o %s --ref %s --Ri 0 --Ro %s --max_shift 1000 --search5d_shift %s --search5d_step  %s --append"\
               %(Images,fnAngles,fnGallery,str(Xdim/2),str(int(Xdim/10)),str(int(Xdim/25))),NumberOfMpi)
    deleteFile(log, os.path.join(ExtraDir,'gallery_sampling.xmd'))
    deleteFile(log, os.path.join(ExtraDir,'gallery_angles.doc'))
    deleteFile(log, os.path.join(ExtraDir,'gallery.doc'))

    # Write angles in the original file and sort
    MD=MetaData(fnAngles)
    for id in MD:
        galleryReference=MD.getValue(MDL_REF,id)
        MD.setValue(MDL_IMAGE_REF,"%05d@%s"%(galleryReference+1,fnGallery),id)
    MD.write(fnAngles)

def produceAlignedImages(log,fnIn,fnOut,fnDiff,volumeIsCTFCorrected):
    from numpy import array, dot
    MDin=MetaData(fnIn)
    MDout=MetaData()
    n=1
    hasCTF=MDin.containsLabel(MDL_CTF_MODEL)
    for i in MDin:
        fnImg=MDin.getValue(MDL_IMAGE,i)
        fnImgRef=MDin.getValue(MDL_IMAGE_REF,i)
        maxCC=MDin.getValue(MDL_MAXCC,i)
        rot =  MDin.getValue(MDL_ANGLE_ROT,i)
        tilt = MDin.getValue(MDL_ANGLE_TILT,i)
        psi =-1.*MDin.getValue(MDL_ANGLE_PSI,i)
        flip = MDin.getValue(MDL_FLIP,i)
        if(flip):
            psi =-psi
        eulerMatrix = Euler_angles2matrix(0.,0.,psi)
        x = MDin.getValue(MDL_SHIFT_X,i)
        y = MDin.getValue(MDL_SHIFT_Y,i)
        shift = array([x, y, 0])
        shiftOut = dot(eulerMatrix, shift)
        [x,y,z]= shiftOut
        if flip:
            x = -x
        id=MDout.addObject()
        MDout.setValue(MDL_IMAGE, fnImg, id)
        MDout.setValue(MDL_IMAGE_REF, fnImgRef, id)
        MDout.setValue(MDL_IMAGE1, "%05d@%s"%(n,fnDiff), id)
        if hasCTF:
            fnCTF=MDin.getValue(MDL_CTF_MODEL,i)
            MDout.setValue(MDL_CTF_MODEL,fnCTF,id)
        MDout.setValue(MDL_MAXCC, maxCC, id)
        MDout.setValue(MDL_ANGLE_ROT, rot, id)
        MDout.setValue(MDL_ANGLE_TILT, tilt, id)
        MDout.setValue(MDL_ANGLE_PSI, psi, id)
        MDout.setValue(MDL_SHIFT_X, x,id)
        MDout.setValue(MDL_SHIFT_Y, y,id)
        MDout.setValue(MDL_FLIP,flip,id)
        MDout.setValue(MDL_ENABLED,1,id)
        n+=1
    MDout.write(fnOut,MD_APPEND)
    
    # Actually create the differences
    img=Image()
    imgRef=Image()
    if hasCTF and volumeIsCTFCorrected:
        fnAcquisition=findAcquisitionInfo(fnOut)
        if not fnAcquisition:
            hasCTF=False
        else:
            mdAcquisition=MetaData(fnAcquisition)
            Ts=mdAcquisition.getValue(MDL_SAMPLINGRATE,mdAcquisition.firstObject())

    for i in MDout:
        img.readApplyGeo(MDout,i)
        imgRef.read(MDout.getValue(MDL_IMAGE_REF,i))
        if hasCTF and volumeIsCTFCorrected:
            fnCTF=MDout.getValue(MDL_CTF_MODEL,i)
            imgRef.applyCTF(fnCTF,Ts)
            img.convert2DataType(DT_DOUBLE)
        imgDiff=img-imgRef
        imgDiff.write(MDout.getValue(MDL_IMAGE1,i))
