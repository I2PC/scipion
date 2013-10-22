import os, shutil, string, glob, math
from os.path import join, exists

#import launch_job, utils_xmipp
from distutils.dir_util import mkpath
from xmipp import *
from protlib_utils import runJob, getMemoryAvailable
from protlib_filesystem import copyFile
from protlib_xmipp import emptyMd

from math import floor

CtfBlockName = 'ctfGroup'
RefBlockName = 'refGroup'

def executeMask(_log,
                  DoMask,
                  DoSphericalMask,
                  maskedFileName,
                  maskRadius,
                  ReconstructedFilteredVolume,
                  userSuppliedMask
                  ):
    _log.debug("executeMask")
    print "executeMask", maskRadius
    if DoMask:
        command = ' -i ' + ReconstructedFilteredVolume + \
                  ' -o ' + maskedFileName
        if DoSphericalMask:
            command += ' --mask circular -' + str(maskRadius)
        else:
            command += ' --mask binary_file ' + userSuppliedMask
        runJob(_log, "xmipp_transform_mask", command)
    else:
        copyFile(_log, ReconstructedFilteredVolume, maskedFileName)


def angular_project_library(_log
                                , AngSamplingRateDeg
                                , BlockWithAllExpImages
                                , ConstantToAddToFiltration
                                , CtfGroupSubsetFileName
                                , DoCtfCorrection
                                , DocFileInputAngles
                                , DoParallel
                                , DoRestricSearchbyTiltAngle
                                , FourierMaxFrequencyOfInterest
                                , KernelAngularProjection
                                , MaxChangeInAngles
                                , maskedFileNamesIter
                                , MpiJobSize
                                , NumberOfMpi
                                , NumberOfThreads
                                , OnlyWinner
                                , PaddingAngularProjection
                                , PerturbProjectionDirections
                                , ProjectLibraryRootName
                                , ProjectionMethod
                                , ResolSam
                                , ResolutionXmdPrevIterMax
                                , SymmetryGroup
                                , SymmetryGroupNeighbourhood
                                , Tilt0
                                , TiltF):
    _log.debug("execute_projection_matching")
    ###need one block per reference
    # Project all references
    print '* Create projection library'
    (Xdim, Ydim, Zdim, Ndim, _) = MetaDataInfo(maskedFileNamesIter)
    memoryUsed=(Xdim*Xdim*Xdim*8.0)/pow(2,20)
    parameters = ' -i ' + maskedFileNamesIter + \
              ' --experimental_images ' + BlockWithAllExpImages + '@' + DocFileInputAngles + \
              ' -o ' + ProjectLibraryRootName + \
              ' --sampling_rate ' + AngSamplingRateDeg + \
              ' --sym ' + SymmetryGroup + 'h' + \
              ' --compute_neighbors' + \
              ' --method ' + ProjectionMethod 
    if ProjectionMethod == 'fourier':
        memoryUsed=memoryUsed*6
        if FourierMaxFrequencyOfInterest == -1:
                md = MetaData(ResolutionXmdPrevIterMax)
                id = md.firstObject()
                FourierMaxFrequencyOfInterest = md.getValue(MDL_RESOLUTION_FREQREAL, id)
                FourierMaxFrequencyOfInterest = ResolSam / FourierMaxFrequencyOfInterest + float(ConstantToAddToFiltration)
                if FourierMaxFrequencyOfInterest > 0.5:
                    FourierMaxFrequencyOfInterest = 0.5
                elif FourierMaxFrequencyOfInterest < 0.:
                    FourierMaxFrequencyOfInterest = 0.001

        parameters += " " + str(PaddingAngularProjection)
        parameters += " " + str(FourierMaxFrequencyOfInterest)
        parameters += " " + str(KernelAngularProjection)

    if (string.atof(MaxChangeInAngles) < 181.):
        parameters += \
              ' --near_exp_data --angular_distance ' + str(MaxChangeInAngles)
    else:
        parameters += \
              ' --angular_distance -1'

    if (PerturbProjectionDirections):
        perturb = math.sin(math.radians(float(AngSamplingRateDeg))) / 4.
        parameters += \
           ' --perturb ' + str(perturb)

    if (DoRestricSearchbyTiltAngle):
        parameters += \
              ' --min_tilt_angle ' + str(Tilt0) + \
              ' --max_tilt_angle ' + str(TiltF)

    if (DoCtfCorrection):
        parameters += \
              ' --groups ' + CtfGroupSubsetFileName
    processorsToUse=NumberOfMpi * NumberOfThreads
    if processorsToUse>1:
        memoryAvailable=getMemoryAvailable()
        processorsToUse=min(processorsToUse,floor(memoryAvailable/memoryUsed))
    if (DoParallel and processorsToUse>1):
        parameters = parameters + ' --mpi_job_size ' + str(MpiJobSize)
    if (len(SymmetryGroupNeighbourhood) > 1):
        parameters += \
          ' --sym_neigh ' + SymmetryGroupNeighbourhood + 'h'
    if (OnlyWinner):
        parameters += \
              ' --only_winner '

    runJob(_log, 'xmipp_angular_project_library',
                         parameters,
                         processorsToUse)
    if (not DoCtfCorrection):
        src = ProjectLibraryRootName.replace(".stk", '_sampling.xmd')
        dst = src.replace('sampling.xmd', 'group%06d_sampling.xmd' % 1)
        copyFile(_log, src, dst)


def projection_matching(_log
                            , AvailableMemory
                            , CtfGroupRootName
                            , CtfGroupDirectory
                            , DocFileInputAngles
                            , DoComputeResolution
                            , DoCtfCorrection
                            , DoScale
                            , DoParallel
                            , InnerRadius
                            , MaxChangeOffset
                            , MpiJobSize
                            , NumberOfMpi
                            , NumberOfThreads
                            , OuterRadius
                            , PaddingFactor
                            , ProjectLibraryRootName
                            , ProjMatchRootName
                            , ReferenceIsCtfCorrected
                            , ScaleStep
                            , ScaleNumberOfSteps
                            , Search5DShift
                            , Search5DStep
                            ):
    # Loop over all CTF groups
    # Use reverse order to have same order in add_to docfiles from angular_class_average
    # get all ctf groups
    _DoCtfCorrection = DoCtfCorrection
    _ProjMatchRootName = ProjMatchRootName
    refname = str(ProjectLibraryRootName)
    file_name = join(CtfGroupDirectory, CtfGroupRootName) + 'Info.xmd'
    if exists(file_name):
        auxMD = MetaData("numberGroups@" + file_name)
        NumberOfCtfGroups = auxMD.getValue(MDL_COUNT, auxMD.firstObject())
    else:
        NumberOfCtfGroups = 1

    
    CtfGroupName = CtfGroupRootName #,ictf+1,'')
    CtfGroupName = CtfGroupDirectory + '/' + CtfGroupName
    #remove output metadata
    if os.path.exists(_ProjMatchRootName):
        os.remove(_ProjMatchRootName)
    
    for ii in range(NumberOfCtfGroups):
        if NumberOfCtfGroups > 1 :
            print 'Focus Group: ', ii + 1, '/', NumberOfCtfGroups
        ictf = NumberOfCtfGroups - ii 
        
        inputdocfile = CtfBlockName + str(ictf).zfill(FILENAMENUMBERLENGTH) + '@' + DocFileInputAngles
        outputname = CtfBlockName + str(ictf).zfill(FILENAMENUMBERLENGTH) + '@' + _ProjMatchRootName
        baseTxtFile = refname[:-len('.stk')] 
        neighbFile = baseTxtFile + '_sampling.xmd'
        if (os.path.exists(neighbFile)):
            os.remove(neighbFile)
        neighbFileb = baseTxtFile + '_group' + str(ictf).zfill(FILENAMENUMBERLENGTH) + '_sampling.xmd'
        print 'neighbFileb: ', neighbFileb
        copyFile(_log, neighbFileb, neighbFile)

        parameters = ' -i ' + inputdocfile + \
                    ' -o ' + outputname + \
                    ' --ref ' + refname + \
                    ' --Ri ' + str(InnerRadius) + \
                    ' --Ro ' + str(OuterRadius) + \
                    ' --max_shift ' + str(MaxChangeOffset) + \
                    ' --search5d_shift ' + str(Search5DShift) + \
                    ' --search5d_step  ' + str(Search5DStep) + \
                    ' --mem ' + str(AvailableMemory * NumberOfThreads) + \
                    ' --thr ' + str(NumberOfThreads) + \
                    ' --append '

        
        if (DoScale):
            parameters += \
                    ' --scale ' + str(ScaleStep) + ' ' + str(ScaleNumberOfSteps) 
        
        if (_DoCtfCorrection and ReferenceIsCtfCorrected):
            ctffile = str(ictf).zfill(FILENAMENUMBERLENGTH) + '@' + CtfGroupName + '_ctf.stk'
            parameters += \
                      ' --pad ' + str(PaddingFactor) + \
                      ' --ctf ' + ctffile
        
        if (DoParallel):
            parameters = parameters + ' --mpi_job_size ' + str(MpiJobSize)
        
        runJob(_log, 'xmipp_angular_projection_matching',
                            parameters,
                            NumberOfMpi,
                            NumberOfThreads
                            )
        
def assign_images_to_references(_log
                         , BlockWithAllExpImages
                         , CtfGroupDirectory
                         , CtfGroupRootName
                         , DocFileInputAngles #outputfile
                         , ProjMatchRootName #input file
                         , NumberOfReferences
                         ):
    ''' assign the images to the different references based on the crosscorrelation coeficient
        #if only one reference it just copy the docfile generated in the previous step
        '''
    file_name = join(CtfGroupDirectory, CtfGroupRootName) + 'Info.xmd'
    if exists(file_name):
        auxMD = MetaData("numberGroups@" + file_name)
        NumberOfCtfGroups = auxMD.getValue(MDL_COUNT, auxMD.firstObject())
    else:
        NumberOfCtfGroups = 1

    #!a
    #DocFileInputAngles  = DocFileInputAngles
    #ProjMatchRootName   = ProjMatchRootName#
    #NumberOfCtfGroups   = NumberOfCtfGroups
    #NumberOfReferences  = NumberOfReferences
    
    #first we need a list with the references used. That is,
    #read all docfiles and map referecendes to a mdl_order
    MDaux = MetaData()
    MDSort = MetaData()
    MD = MetaData()
    MD1 = MetaData()
    MDout = MetaData()
    MDout.setComment("metadata with  images, the winner reference as well as the ctf group")
    print "one"
    """ compute auxiliary index order, it may become handy to match projections and
    projection directions
    """
    
    
    
    mycounter = 1L
    for iCTFGroup in range(1, NumberOfCtfGroups + 1):
        auxInputdocfile = CtfBlockName + str(iCTFGroup).zfill(FILENAMENUMBERLENGTH) + '@'
        for iRef3D in range(1, NumberOfReferences + 1):
            inputFileName = ProjMatchRootName[iRef3D]
            inputdocfile = auxInputdocfile + inputFileName
            MD.read(inputdocfile)
            for id in MD:
                t = MD.getValue(MDL_REF, id)
                i = MDSort.addObject()
                MDSort.setValue(MDL_REF, t, i)
    
    MDSort.removeDuplicates()
    
    for id in MDSort:
        MDSort.setValue(MDL_ORDER, mycounter, id)
        mycounter += 1
    ####################
    outputdocfile = DocFileInputAngles
    if os.path.exists(outputdocfile):
        os.remove(outputdocfile)
        
        
    MDout2 = MetaData()
    for iCTFGroup in range(1, NumberOfCtfGroups + 1):
        MDaux.clear()
        auxInputdocfile = CtfBlockName + str(iCTFGroup).zfill(FILENAMENUMBERLENGTH) + '@'
        for iRef3D in range(1, NumberOfReferences + 1):
            inputFileName = ProjMatchRootName[iRef3D]
            inputdocfile = auxInputdocfile + inputFileName
            MD.clear()
            MD.read(inputdocfile)
            #In practice you should not get duplicates
            MD.removeDuplicates()
            MD.setValueCol(MDL_REF3D, iRef3D)
            MD.setValueCol(MDL_DEFGROUP, iCTFGroup)
            #MD.setValueCol(MDL_CTF_MODEL,auxInputdocfile[:-1])
            MDaux.unionAll(MD)
        MDaux.sort()
        MD.aggregate(MDaux, AGGR_MAX, MDL_IMAGE, MDL_MAXCC, MDL_MAXCC)
        #if a single image is assigned to two references with the same 
        #CC use it in both reconstruction
        #recover atribbutes after aggregate function
        
        MD1.join  (MD, MDaux, MDL_UNDEFINED, MDL_UNDEFINED, NATURAL)
        MDout.join(MD1, MDSort, MDL_UNDEFINED, MDL_UNDEFINED, NATURAL)
        print 'write file: ', auxInputdocfile + outputdocfile
        MDout.write(auxInputdocfile + outputdocfile, MD_APPEND)
        MDout2.unionAll(MDout)
        
    MDout2.write(BlockWithAllExpImages + '@' + outputdocfile, MD_APPEND)
    #!a original_angles too
    
    #we are done but for the future it is convenient to create more blocks
    #with the pairs ctf_group reference    
    for iCTFGroup in range(1, NumberOfCtfGroups + 1):
        auxInputdocfile = CtfBlockName + str(iCTFGroup).zfill(FILENAMENUMBERLENGTH) + '@'
        #print 'read file: ', auxInputdocfile+outputdocfile
        MDaux.read(auxInputdocfile + outputdocfile)
        for iRef3D in range(1, NumberOfReferences + 1):
            auxOutputdocfile = CtfBlockName + \
                                str(iCTFGroup).zfill(FILENAMENUMBERLENGTH)
            auxOutputdocfile += '_' + RefBlockName + \
                                      str(iRef3D).zfill(FILENAMENUMBERLENGTH) + '@'
            #select images with ref3d=iRef3D
            MDout.importObjects(MDaux, MDValueEQ(MDL_REF3D, iRef3D))
            MDout.write(auxOutputdocfile + outputdocfile, MD_APPEND)

def angular_class_average(_log
                         , Align2DIterNr
                         , Align2dMaxChangeOffset
                         , Align2dMaxChangeRot
                         , CtfGroupDirectory
                         , CtfGroupRootName
                         , DiscardImages
                         , DiscardPercentage
                         , DiscardPercentagePerClass
                         , DoAlign2D
                         , DoComputeResolution
                         , DoCtfCorrection
                         , DocFileInputAngles
                         , DoParallel
                         , DoSaveImagesAssignedToClasses
                         , DoSplitReferenceImages
                         , InnerRadius
                         , MaxChangeOffset
                         , MinimumCrossCorrelation
                         , MpiJobSize
                         , NumberOfMpi
                         , NumberOfThreads
                         , OutClasses
                         , PaddingFactor
                         , ProjectLibraryRootName
                         ):
                             

    CtfGroupName = CtfGroupDirectory + '/' + CtfGroupRootName
    refname = str(ProjectLibraryRootName)
    baseTxtFile = refname[:-len('.stk')] 
    neighbFile = baseTxtFile + '.xmd'
	
    if emptyMd(DocFileInputAngles):
        print "Empty metadata file: %s" % DocFileInputAngles
        return

    parameters = ' -i ctfGroup[0-9][0-9][0-9][0-9][0-9][0-9]\$@' + DocFileInputAngles + \
                 ' --lib ' + refname.replace(".stk", ".doc") + \
                 ' -o ' + OutClasses
    if(DoSaveImagesAssignedToClasses):
        parameters += ' --save_images_assigned_to_classes'  
                  
    if(DiscardImages == 'maxCC'):
        parameters += ' --limit0 ' + MinimumCrossCorrelation
    elif(DiscardImages == 'percentage'):
        parameters += ' --limitRper ' + DiscardPercentage
    elif(DiscardImages == 'classPercentage'):
        parameters += ' --limitRclass ' + DiscardPercentagePerClass
    #else 'none'
        
    # On-the fly apply Wiener-filter correction and add all CTF groups together
    if (DoCtfCorrection):
        parameters += \
                   ' --wien ' + CtfGroupName + '_wien.stk' + \
                   ' --pad ' + str(PaddingFactor)
                   
    if (DoAlign2D == '1'):
        parameters += \
                  ' --iter ' + Align2DIterNr + \
                  ' --Ri ' + str(InnerRadius) + \
                  ' --Ro ' + str(OuterRadius)
                  
    if (DoComputeResolution and DoSplitReferenceImages):
        parameters += ' --split '
                  
    if (DoParallel):
        parameters = parameters + ' --mpi_job_size ' + str(MpiJobSize)

    runJob(_log,
           'xmipp_angular_class_average',
           parameters,
           NumberOfMpi * NumberOfThreads
           )
        
def reconstruction(_log
                   , ARTReconstructionExtraCommand
                   , WBPReconstructionExtraCommand
                   , FourierReconstructionExtraCommand
                   , Iteration_number
                   , DoParallel
                   , maskedFileNamesIter
                   , MpiJobSize
                   , NumberOfMpi
                   , NumberOfThreads
                   , ReconstructionMethod
                   , FourierMaxFrequencyOfInterest
                   , ARTLambda
                   , SymmetryGroup
                   , ReconstructionXmd
                   , ReconstructedVolume
                   , ResolSam
                   , ResolutionXmdPrevIterMax
                   , PaddingFactor
                   , ConstantToAddToFiltration
                   ):
    
    #if inout metadata is empty create a Blanck image
    if emptyMd(ReconstructionXmd):
        from protlib_utils import printLog
        img = Image()
        img.read(maskedFileNamesIter, DATA)
        #(x,y,z,n) = img.getDimensions()
        printLog("Metadata %s is empty. Creating a Black file named %s" % (ReconstructionXmd, ReconstructedVolume))
        #createEmptyFile(ReconstructedVolume,x,y,z,n)
        img.initRandom()
        img.write(ReconstructedVolume)
        return

    
    print '*********************************************************************'
    print '* Reconstruct volume using '
    if ReconstructionMethod == 'wbp':
        program = 'xmipp_reconstruct_wbp'
        parameters = ' -i ' + ReconstructionXmd + \
                    ' --doc ' + ReconstructionXmd + \
                    ' -o ' + ReconstructedVolume + \
                    ' --sym ' + SymmetryGroup + \
                    ' --weight --use_each_image '
        parameters = parameters + WBPReconstructionExtraCommand
                  
    elif ReconstructionMethod == 'art':
        program = 'xmipp_reconstruct_art'

        parameters = ' -i ' + ReconstructionXmd + \
                   ' -o ' + ReconstructedVolume + ' ' + \
                   ' --sym ' + SymmetryGroup + \
                   ' --thr ' + str(NumberOfThreads) + \
                   ' --WLS '
        if len(ARTLambda) > 1:
           parameters = parameters + ' -l ' + ARTLambda + ' '
        parameters = parameters + ARTReconstructionExtraCommand
        
        NumberOfMpi = 1
        NumberOfThreads = 1
        DoParallel = False
                
    elif ReconstructionMethod == 'fourier':
        
        if FourierMaxFrequencyOfInterest == -1:
                md = MetaData(ResolutionXmdPrevIterMax)
                id = md.firstObject()
                FourierMaxFrequencyOfInterest = md.getValue(MDL_RESOLUTION_FREQREAL, id)
                FourierMaxFrequencyOfInterest = ResolSam / FourierMaxFrequencyOfInterest + float(ConstantToAddToFiltration)
                if FourierMaxFrequencyOfInterest > 0.5:
                    FourierMaxFrequencyOfInterest = 0.5
                elif FourierMaxFrequencyOfInterest < 0.:
                    FourierMaxFrequencyOfInterest = 0.001

        program = 'xmipp_reconstruct_fourier'
        parameters = ' -i ' + ReconstructionXmd + \
                   ' -o ' + ReconstructedVolume + \
                   ' --sym ' + SymmetryGroup + \
                   ' --thr ' + str(NumberOfThreads) + \
                   ' --weight ' + \
                   ' --max_resolution ' + str(FourierMaxFrequencyOfInterest) + \
                   ' --padding ' + str(PaddingFactor) + ' ' + str(PaddingFactor)
 
    if (DoParallel):
        parameters = parameters + ' --mpi_job_size ' + str(MpiJobSize)
            
    runJob(_log
           , program
           , parameters
           , NumberOfMpi
           , NumberOfThreads
           )


def  compute_resolution(_log
                         , ConstantToAddToFiltration
                         , FourierMaxFrequencyOfInterest
                         , ResolutionXmdCurrIter
                         , ResolutionXmdCurrIterMax
                         , ResolutionXmdPrevIterMax
                         , OuterRadius
                         , ReconstructedVolumeSplit1
                         , ReconstructedVolumeSplit2
                         , ReconstructionMethod
                         , ResolSam
                         ):

        # Prevent high-resolution correlation because of discrete mask from wbp
        innerrad = int(OuterRadius) - 2
        Outputvolumes = [ReconstructedVolumeSplit1, ReconstructedVolumeSplit2]
        for i in range(len(Outputvolumes)):
           print '*********************************************************************'
           print '* Applying a soft mask'
           command = " -i " + Outputvolumes[i] + \
                     " --mask  raised_cosine -" + str(innerrad) + \
                     " -" + str(OuterRadius)
                     
           runJob(_log, "xmipp_transform_mask", command)
      
        print '**************************************************************'
        print '* Compute resolution ' 
        command = " --ref " + ReconstructedVolumeSplit1 + \
                  " -i " + ReconstructedVolumeSplit2 + ' --sampling_rate ' + str(ResolSam) + \
                  " -o " + ResolutionXmdCurrIter
        if ReconstructionMethod == 'fourier':
#            img = Image()
#            img.read(ReconstructedVolumeSplit1, HEADER)
#            (x,y,z,n) = img.getDimensions()
            #2.5 is the blob radius 
            #aux4 = 2.5 * 0.5 / x
            if FourierMaxFrequencyOfInterest == -1:
                md = MetaData(ResolutionXmdPrevIterMax)
                id = md.firstObject()
                FourierMaxFrequencyOfInterest = md.getValue(MDL_RESOLUTION_FREQREAL, id)
                print "1FourierMaxFrequencyOfInterest", FourierMaxFrequencyOfInterest
                normalizedFreq = ResolSam / FourierMaxFrequencyOfInterest
                print "2normalizedFreq", normalizedFreq
                normalizedFreq += float(ConstantToAddToFiltration)
                print "3normalizedFreq", normalizedFreq
                FourierMaxFrequencyOfInterest = ResolSam / normalizedFreq
                print "4FourierMaxFrequencyOfInterest", FourierMaxFrequencyOfInterest
            else:
                FourierMaxFrequencyOfInterest = ResolSam / FourierMaxFrequencyOfInterest
            command += " --max_sam " + str (FourierMaxFrequencyOfInterest)
        
        runJob(_log, "xmipp_resolution_fsc", command)
        print "compute resolution1"
        #compute resolution
        mdRsol = MetaData(ResolutionXmdCurrIter)
        mdResolOut = MetaData()
        mdResolOut.importObjects(mdRsol, MDValueLT(MDL_RESOLUTION_FRC, 0.5))
        print "compute resolution2"
        if mdResolOut.size()==0:
            mdResolOut.clear()
            mdResolOut.addObject()
            id=mdResolOut.firstObject()
            mdResolOut.setValue(MDL_RESOLUTION_FREQREAL,ResolSam*2.,id)
            mdResolOut.setValue(MDL_RESOLUTION_FRC,0.5,id)
        else:
            mdResolOut.sort()

#        filter_frequence = mdResolOut.aggregateSingle(AGGR_MAX,MDL_RESOLUTION_FREQREAL)
#        mdResolOut.clear()
#        mdResolOut.importObjects(mdRsol,MDValueEQ(MDL_RESOLUTION_FREQREAL, filter_frequence))
        id = mdResolOut.firstObject()
        filter_frequence = mdResolOut.getValue(MDL_RESOLUTION_FREQREAL, id)
        frc = mdResolOut.getValue(MDL_RESOLUTION_FRC, id)

        md = MetaData()
        id = md.addObject()
        md.setColumnFormat(False)

        md.setValue(MDL_RESOLUTION_FREQREAL, filter_frequence, id)
        md.setValue(MDL_RESOLUTION_FRC, frc, id)
        md.setValue(MDL_SAMPLINGRATE, ResolSam, id)
        md.write(ResolutionXmdCurrIterMax, MD_APPEND)
        

def  filter_volume(_log
                   , FourierMaxFrequencyOfInterest
                   , ReconstructedVolume
                   , ReconstructedFilteredVolume
                   , DoComputeResolution
                   , OuterRadius
                   , DoLowPassFilter
                   , UseFscForFilter
                   , ConstantToAddToFiltration
                   , ResolutionXmdPrevIterMax
                   , ResolSam
                   ):

    if (not DoLowPassFilter):
        copyFile(_log, ReconstructedVolume, ReconstructedFilteredVolume)
    else:   
        if (UseFscForFilter):
           if (FourierMaxFrequencyOfInterest == -1):
               md = MetaData(ResolutionXmdPrevIterMax)
               id = md.firstObject()
               FourierMaxFrequencyOfInterest = md.getValue(MDL_RESOLUTION_FREQREAL, id)
               FourierMaxFrequencyOfInterest = ResolSam / FourierMaxFrequencyOfInterest
            
           filter_in_pixels_at = float(FourierMaxFrequencyOfInterest) + float(ConstantToAddToFiltration)
        else:
           filter_in_pixels_at = float(ConstantToAddToFiltration)

        if (filter_in_pixels_at > 0.5):
           copyFile(_log, ReconstructedVolume, ReconstructedFilteredVolume)
        else:
           command = " -i " + ReconstructedVolume + \
                     " -o " + ReconstructedFilteredVolume + ' --fourier low_pass ' + \
                     str (filter_in_pixels_at)
           runJob(_log, "xmipp_transform_filter", command)

