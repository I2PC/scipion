import os, shutil, string, glob, math
#import launch_job, utils_xmipp
from distutils.dir_util import mkpath
from xmipp import *
from protlib_utils import runJob

CtfBlockName = 'ctfGroup'
RefBlockName = 'refGroup'

def executeMask(_log, 
                  DoMask, 
                  DoSphericalMask, 
                  maskedFileName, 
                  maskRadius,
                  reconstructedFileName,
                  userSuppliedMask
                  ):
    _log.debug("executeMask")
    if DoMask:
#        print '*********************************************************************'
#        print '* Mask the reference volume'
        command = ' -i ' + reconstructedFileName + \
                  ' -o ' + maskedFileName

        if DoSphericalMask:
            command += ' --mask circular -' + str(maskRadius)
        else:
            command += ' --mask ' + userSuppliedMask

        runJob(_log,"xmipp_transform_mask", command)

    else:
        shutil.copy(reconstructedFileName, maskedFileName)
        _log.info("Skipped Mask")
        _log.info("cp " + reconstructedFileName + " " + maskedFileName)
#        print '*********************************************************************'
        print '* Skipped Mask'


def angular_project_library(_log
                                ,AngSamplingRateDeg
                                ,CtfGroupSubsetFileName
                                ,DoCtfCorrection
                                ,DocFileInputAngles
                                ,DoParallel
                                ,DoRestricSearchbyTiltAngle
                                ,MaxChangeInAngles
                                ,maskedFileNamesIter
                                ,MpiJobSize
                                ,NumberOfMpi
                                ,NumberOfThreads
                                ,OnlyWinner
                                ,PerturbProjectionDirections
                                ,ProjectLibraryRootName
                                ,SymmetryGroup
                                ,SymmetryGroupNeighbourhood
                                ,Tilt0
                                ,TiltF):
    _log.debug("execute_projection_matching")
    ###need one block per reference
    # Project all references
    print '* Create projection library'
    parameters=' -i '                   + maskedFileNamesIter + \
              ' --experimental_images ' + DocFileInputAngles + \
              ' -o '                    + ProjectLibraryRootName + \
              ' --sampling_rate '       + AngSamplingRateDeg  + \
              ' --sym '                 + SymmetryGroup + 'h' + \
              ' --compute_neighbors'

    if ( string.atof(MaxChangeInAngles) < 181.):
        parameters+= \
              ' --near_exp_data --angular_distance '    + str(MaxChangeInAngles)
    else:
        parameters+= \
              ' --angular_distance -1'

    if (PerturbProjectionDirections):
        perturb=math.sin(math.radians(float(AngSamplingRateDeg)))/4.
        parameters+= \
           ' --perturb ' + str(perturb)

    if (DoRestricSearchbyTiltAngle):
        parameters+=  \
              ' --min_tilt_angle '      + str(Tilt0) + \
              ' --max_tilt_angle '      + str(TiltF)

    if (DoCtfCorrection):
        parameters+=  \
              ' --groups '              + CtfGroupSubsetFileName
    _DoParallel=DoParallel
    if (DoParallel):
        parameters = parameters + ' --mpi_job_size ' + str(MpiJobSize)
    if (len(SymmetryGroupNeighbourhood)>1):
        parameters+= \
          ' --sym_neigh ' + SymmetryGroupNeighbourhood + 'h'
    if (OnlyWinner):
        parameters+= \
              ' --only_winner '

    runJob(_log,'xmipp_angular_project_library',
                         parameters,
                         NumberOfMpi*NumberOfThreads)
    if (not DoCtfCorrection):
        print "a1"
        src=ProjectLibraryRootName.replace(".stk",'_sampling.xmd')
        dst = src.replace('sampling.xmd','group%06d_sampling.xmd' % 1)
        shutil.copy(src, dst)


def projection_matching(_log
                            , AvailableMemory
                            , CtfGroupRootName
                            , CtfGroupDirectory
                            , DoComputeResolution
                            , DoCtfCorrection
                            , DoScale
                            , DoParallel
                            , InnerRadius
                            , MaxChangeOffset
                            , MpiJobSize
                            , NumberOfCtfGroups
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
    _DoCtfCorrection    = DoCtfCorrection
    _ProjMatchRootName  = ProjMatchRootName
    refname = str(ProjectLibraryRootName)
    NumberOfCtfGroups=NumberOfCtfGroups
    
    CtfGroupName = CtfGroupRootName #,ictf+1,'')
    CtfGroupName = CtfGroupDirectory + '/' + CtfGroupName
    #remove output metadata
    if os.path.exists(_ProjMatchRootName):
        os.remove(_ProjMatchRootName)
    
    for ii in range(NumberOfCtfGroups):
        if NumberOfCtfGroups>1 :
            print 'Focus Group: ', ii+1,'/',NumberOfCtfGroups
        ictf    = NumberOfCtfGroups - ii 
#        if (_DoCtfCorrection):
        #outputname   = _ProjMatchRootName + '_' + CtfGroupName 
        inputdocfile    = CtfBlockName+str(ictf).zfill(FILENAMENUMBERLENGTH) + '@' + CtfGroupName + '_images.sel'
        outputname   = CtfBlockName+str(ictf).zfill(FILENAMENUMBERLENGTH) + '@'+ _ProjMatchRootName
        #inputdocfile = (os.path.basename(inselfile)).replace('.sel','.doc')
        baseTxtFile  = refname[:-len('.stk')] 
        neighbFile      = baseTxtFile + '_sampling.xmd'
        if (os.path.exists(neighbFile)):
            os.remove(neighbFile)
        neighbFileb     = baseTxtFile + '_group'+str(ictf).zfill(FILENAMENUMBERLENGTH) + '_sampling.xmd'
        shutil.copy(neighbFileb, neighbFile)
#        else:
#            inputdocfile    = 'ctfGroup'+str(1).zfill(utils_xmipp.FILENAMENUMBERLENTGH) + '@' + CtfGroupName + '_images.sel'
#            outputname   = 'ctfGroup'+str(1).zfill(utils_xmipp.FILENAMENUMBERLENTGH) + '@'+ _ProjMatchRootName

        parameters= ' -i '               + inputdocfile + \
                    ' -o '               + outputname + \
                    ' --ref '            + refname + \
                    ' --Ri '             + str(InnerRadius)           + \
                    ' --Ro '             + str(OuterRadius)           + \
                    ' --max_shift '      + str(MaxChangeOffset) + \
                    ' --search5d_shift ' + str(Search5DShift) + \
                    ' --search5d_step  ' + str(Search5DStep) + \
                    ' --mem '            + str(AvailableMemory * NumberOfThreads) + \
                    ' --thr '            + str(NumberOfThreads) +\
                    ' --append '

        
        if (DoScale):
            parameters += \
                    ' --scale '          + str(ScaleStep) + ' ' + str(ScaleNumberOfSteps) 
        
        if (_DoCtfCorrection and ReferenceIsCtfCorrected):
            ctffile = str(ictf).zfill(FILENAMENUMBERLENGTH) + '@' + CtfGroupName + '_ctf.stk'
            parameters += \
                      ' --pad '            + str(PaddingFactor) + \
                      ' --ctf '            + ctffile
        
        if (DoParallel):
            parameters = parameters + ' --mpi_job_size ' + str(MpiJobSize)
        
        runJob(_log,'xmipp_angular_projection_matching',
                            parameters,
                            NumberOfMpi,
                            NumberOfThreads
                            )
        
def assign_images_to_references(_log
                         , DocFileInputAngles #outputfile
                         , NumberOfCtfGroups
                         , ProjMatchRootName #input file
                         , NumberOfReferences
                         ):
    ''' assign the images to the different references based on the crosscorrelation coeficient
        #if only one reference it just copy the docfile generated in the previous step
        '''
    #!a
    #DocFileInputAngles  = DocFileInputAngles
    #ProjMatchRootName   = ProjMatchRootName#
    #NumberOfCtfGroups   = NumberOfCtfGroups
    #NumberOfReferences  = NumberOfReferences
    
    #first we need a list with the references used. That is,
    #read all docfiles and map referecendes to a mdl_order
    MDaux  = MetaData()
    MDSort = MetaData()
    MD     = MetaData()
    MD1    = MetaData()
    MDout  = MetaData()
    MDout.setComment("metadata with  images, the winner reference as well as the ctf group")
    print "one"
    """ compute auxiliary index order, it may become handy to match projections and
    projection directions
    """
    mycounter=1L
    for iCTFGroup in range(1,NumberOfCtfGroups+1):
        auxInputdocfile = CtfBlockName + str(iCTFGroup).zfill(FILENAMENUMBERLENGTH)+'@'
        for iRef3D in range(1,NumberOfReferences+1):
            inputFileName = ProjMatchRootName[iRef3D]
            print "[",iCTFGroup," ", iRef3D, "] inputFileName:", inputFileName
            inputdocfile    = auxInputdocfile+ inputFileName
            print "inputdocfile: ", inputdocfile 
            MD.read(inputdocfile)
            for id in MD:
                t=MD.getValue(MDL_REF,id)
                i=MDSort.addObject()
                MDSort.setValue(MDL_REF,t,i)
    print "two"
    MDSort.write("mdsort_wt_duplicates.xmd")
    
    MDSort.removeDuplicates()
    
    MDSort.write("mdsort_wo_duplicates.xmd")
    
    for id in MDSort:
        MDSort.setValue(MDL_ORDER,mycounter,id)
        mycounter += 1
    MDSort.write("mdsort_wo_duplicates_sorted.xmd")
    ####################
    print "three"
    print "bbb",ProjMatchRootName[1], DocFileInputAngles
    outputdocfile =  DocFileInputAngles
    if os.path.exists(outputdocfile):
        os.remove(outputdocfile)
    print "four"
    for iCTFGroup in range(1,NumberOfCtfGroups+1):
        print "five",iCTFGroup
        MDaux.clear()
        auxInputdocfile = CtfBlockName + str(iCTFGroup).zfill(FILENAMENUMBERLENGTH)+'@'
        for iRef3D in range(1,NumberOfReferences+1):
            inputFileName = ProjMatchRootName[iRef3D]
            inputdocfile    = auxInputdocfile+ inputFileName
            MD.clear()
            MD.read(inputdocfile)
            #In practice you should not get duplicates
            MD.removeDuplicates()
            MD.setValueCol(MDL_REF3D,iRef3D)
            #MD.setValueCol(MDL_CTFMODEL,auxInputdocfile[:-1])
            MDaux.unionAll(MD)
        MDaux.sort()
        MD.write("md_test.xmd")
        MDaux.write("mdaux_test.xmd")
        MD.aggregate(MDaux,AGGR_MAX,MDL_IMAGE,MDL_MAXCC,MDL_MAXCC)
        MD.write("md_test_afteraggregate.xmd")
        MDaux.write("mdaux_test_afteraggregate.xmd")
        #if a single image is assigned to two references with the same 
        #CC use it in both reconstruction
        #recover atribbutes after aggregate function
        
        MD1.join  (MD,  MDaux,  MDL_UNDEFINED, MDL_UNDEFINED, NATURAL)
        MD1.write("md1.xmd")
        MDout.join(MD1, MDSort, MDL_UNDEFINED, MDL_UNDEFINED, NATURAL)
        MDout.write("mdout.xmd")
        MDout.write(auxInputdocfile+outputdocfile,MD_APPEND)
        
    #we are done but for the future it is convenient to create more blocks
    #with the pairs ctf_group reference    
    for iCTFGroup in range(1,NumberOfCtfGroups+1):
        auxInputdocfile  = CtfBlockName + str(iCTFGroup).zfill(FILENAMENUMBERLENGTH)+'@'
        MDaux.read(auxInputdocfile+outputdocfile)
        for iRef3D in range(1,NumberOfReferences+1):
            auxOutputdocfile  = CtfBlockName + \
                                str(iCTFGroup).zfill(FILENAMENUMBERLENGTH)
            auxOutputdocfile += '_' + RefBlockName +\
                                      str(iRef3D).zfill(FILENAMENUMBERLENGTH)+'@'
            #select images with ref3d=iRef3D
            MDout.importObjects(MDaux,MDValueEQ(MDL_REF3D, iRef3D))
            MDout.write(auxOutputdocfile+outputdocfile,MD_APPEND)

def angular_class_average(_log
                         , Action
                         , Align2DIterNr
                         , Align2dMaxChangeRot
                         , Align2dMaxChangeOffset
                         , CtfGroupDirectory
                         , CtfGroupRootName
                         , DiscardPercentage
                         , DoAlign2D
                         , DoComputeResolution
                         , DoCtfCorrection
                         , DocFileInputAngles
                         , DoParallel
                         , DoSplitReferenceImages
                         , InnerRadius
                         , MaxChangeOffset
                         , MinimumCrossCorrelation
                         , NumberOfReferences
                         , NumberOfCtfGroups
                         , NumberOfMpi
                         , NumberOfThreads
                         , PaddingFactor
                         , ProjectLibraryRootName
                         , ProjMatchRootName
                         , refN
                         ):
    
    if Action == "preprocessing" or Action == "postprocessing" :
        NumberOfCtfGroups = 0
        
                         
    # Now make the class averages
    CtfGroupName        = CtfGroupDirectory + '/' +\
                          CtfGroupRootName
    refname             = str(ProjectLibraryRootName)
    MD = MetaData()
    ProjMatchRootName = ProjMatchRootName
    for iCTFGroup in range(1,NumberOfCtfGroups+1):
#        for iRef3D in range(1,NumberOfReferences+1):
        auxInputdocfile  = CtfBlockName + \
                            str(iCTFGroup).zfill(FILENAMENUMBERLENGTH)
        auxInputdocfile += '_' + RefBlockName +\
                            str(refN).zfill(FILENAMENUMBERLENGTH)+'@'
        MD.read(auxInputdocfile+DocFileInputAngles)
        if MD.size()==0:
            print "Empty metadata, remember to copy the reference",iCTFGroup,refN
            continue;
        #Md.write("test.xmd" + str(iCTFGroup).zfill(2) +'_'+str(refN).zfill(2))
        parameters =  ' -i '       + auxInputdocfile  + DocFileInputAngles +\
                      ' --lib '    + refname.replace(".stk",".doc") + \
                      ' --write_selfiles ' + \
                      ' --limit0 ' + MinimumCrossCorrelation + \
                      ' --limitR ' + DiscardPercentage
        if (DoCtfCorrection):
            # On-the fly apply Wiener-filter correction and add all CTF groups together
            parameters += \
                       ' --wien '   + str(iCTFGroup).zfill(FILENAMENUMBERLENGTH)+'@' + CtfGroupName + '_wien.stk' + \
                       ' --pad '    + str(PaddingFactor) + \
                       ' --add_to ' + ProjMatchRootName.replace('.doc','__')
        else:
            parameters += \
                      ' -o '                + ProjMatchRootName

                       
        if Action == "preprocessing":
            parameters += \
                       ' --do_preprocess '
             
        if Action == "postprocessing" :
            parameters += \
                       ' --do_postprocess --number_3dreferences ' + str(NumberOfReferences) + ' '

        if (DoAlign2D == '1'):
            parameters += \
                      ' --iter '             + Align2DIterNr  + \
                      ' --Ri '               + str(InnerRadius)           + \
                      ' --Ro '               + str(OuterRadius)           + \
                      ' --max_shift '        + MaxChangeOffset + \
                      ' --max_shift_change ' + Align2dMaxChangeOffset + \
                      ' --max_psi_change '   + Align2dMaxChangeRot 
        if (DoComputeResolution and DoSplitReferenceImages):
            parameters += \
                      ' --split '
        
        runJob(_log,
               'xmipp_angular_class_average',
               parameters,
               NumberOfMpi,
               NumberOfThreads
               )

