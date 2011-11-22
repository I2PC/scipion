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
    
    MDSort.removeDuplicates()
    
    
    for id in MDSort:
        MDSort.setValue(MDL_ORDER,mycounter,id)
        mycounter += 1
    ####################
    outputdocfile =  DocFileInputAngles
    if os.path.exists(outputdocfile):
        os.remove(outputdocfile)
    for iCTFGroup in range(1,NumberOfCtfGroups+1):
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
            MD.setValueCol(MDL_DEFGROUP,iCTFGroup)
            #MD.setValueCol(MDL_CTFMODEL,auxInputdocfile[:-1])
            MDaux.unionAll(MD)
        MDaux.sort()
        MD.aggregate(MDaux,AGGR_MAX,MDL_IMAGE,MDL_MAXCC,MDL_MAXCC)
        #if a single image is assigned to two references with the same 
        #CC use it in both reconstruction
        #recover atribbutes after aggregate function
        
        MD1.join  (MD,  MDaux,  MDL_UNDEFINED, MDL_UNDEFINED, NATURAL)
        MDout.join(MD1, MDSort, MDL_UNDEFINED, MDL_UNDEFINED, NATURAL)
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
                         , Align2DIterNr
                         , Align2dMaxChangeOffset
                         , Align2dMaxChangeRot
                         , CtfGroupDirectory
                         , CtfGroupRootName
                         , DiscardPercentage
                         , DoAlign2D
                         , DoComputeResolution
                         , DoCtfCorrection
                         , DocFileInputAngles
                         , DoSplitReferenceImages
                         , InnerRadius
                         , MaxChangeOffset
                         , MinimumCrossCorrelation
                         , NumberOfCtfGroups
                         , NumberOfMpi
                         , NumberOfThreads
                         , OutClasses
                         , PaddingFactor
                         , ProjectLibraryRootName
                         , Ref3dNum
                         ):
                             

    CtfGroupName = CtfGroupDirectory + '/' + CtfGroupRootName
    refname      = str(ProjectLibraryRootName)

    MD = MetaData()
    MD.read(DocFileInputAngles)
    if MD.size()==0:
        print "Empty metadata, remember to copy the reference ",NumberOfCtfGroups,Ref3dNum
        return

    parameters =  ' -i '       + DocFileInputAngles +\
                  ' --lib '    + refname.replace(".stk",".doc") + \
                  ' --write_selfiles ' + \
                  ' --limit0 ' + MinimumCrossCorrelation + \
                  ' --limitR ' + DiscardPercentage + \
                  ' --ctfNum ' + str(NumberOfCtfGroups) + \
                  ' --ref3dNum ' + str(Ref3dNum) + \
                  ' -o '        + OutClasses
                  
        # On-the fly apply Wiener-filter correction and add all CTF groups together
    if (DoCtfCorrection):
        parameters += \
                   ' --wien '   + str(NumberOfCtfGroups).zfill(FILENAMENUMBERLENGTH)+'@' + CtfGroupName + '_wien.stk' + \
                   ' --pad '    + str(PaddingFactor)
                   
    if (DoAlign2D == '1'):
        parameters += \
                  ' --iter '             + Align2DIterNr  + \
                  ' --Ri '               + str(InnerRadius)           + \
                  ' --Ro '               + str(OuterRadius)           + \
                  ' --max_shift '        + MaxChangeOffset + \
                  ' --max_shift_change ' + Align2dMaxChangeOffset + \
                  ' --max_psi_change '   + Align2dMaxChangeRot 
                  
    if (DoComputeResolution and DoSplitReferenceImages):
        parameters += ' --split '
                  
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
                   #, DisplayReconstruction
                   , DoParallel
                   , NumberOfMpi
                   , NumberOfThreads
                   , ReconstructionMethod
                   #, FourierMaxFrequencyOfInterest
                   , ARTLambda
                   , SymmetryGroup
                   #, ReconstructedandfilteredVolume
                   , ReconstructionXmd
                   , ReconstructedVolume
                   , DoComputeResolution
                   , DoSplitReferenceImages
                   , PaddingFactor
                   ):

    InputVolume = ReconstructionXmd
    OutputVolume = ReconstructedVolume
    
    print '*********************************************************************'
    print '* Reconstruct volume using '
    if ReconstructionMethod=='wbp':
        OutputVolume = OutputVolume+".vol"
        program = 'xmipp_reconstruct_wbp'
        parameters= ' -i '    + ForReconstructionSel + \
                    ' --doc '  + ForReconstructionDoc + \
                    ' -o '    + OutputVolume + \
                    ' --sym '  + SymmetryGroup + \
                    ' --weight -use_each_image '
        parameters = parameters + WBPReconstructionExtraCommand
        #MyNumberOfThreads = 1
                  
    elif ReconstructionMethod=='art':
        program = 'xmipp_reconstruct_art'
        #DoParallel=False
        parameters=' -i '    + ForReconstructionSel + \
                   ' -o '    + OutputVolume + ' ' + \
                   ' --sym '  + SymmetryGroup + \
                   ' --thr '  + str(NumberOfThreads) + \
                   ' --WLS '
        if len(_ARTLambda)>1:
           parameters = parameters + ' -l '   + _ARTLambda + ' '
        parameters = parameters + _ARTReconstructionExtraCommand
    elif ReconstructionMethod=='fourier':
        #if ( _MyNumberOfMpiProcesses ==1):
            #DoParallel=False
        program = 'xmipp_reconstruct_fourier'
        parameters=' -i '    + ForReconstructionSel + \
                   ' -o '    + OutputVolume + '.vol ' + \
                   ' --sym '  + SymmetryGroup + \
                   ' --thr '  + str(NumberOfThreads) + \
                   ' --weight ' + \
                   ' --max_resolution ' + str(FourierMaxFrequencyOfInterest) +\
                   ' --pad_proj ' + str(PaddingFactor) +\
                   ' --pad_vol ' + str(PaddingFactor)
 
        if (DoParallel):
            parameters = parameters + ' --mpi_job_size ' + str(MpiJobSize)

        #if (_DoComputeResolution and not _DoSplitReferenceImages):
            #myFileName =  ProjMatchDir + '/' + ProjMatchName
            #parameters = parameters + ' -prepare_fsc ' + myFileName + ' '
            #rand_command  = ' xmipp_selfile_split -i '
            #rand_command += ForReconstructionSel + ' -dont_sort -n 1 ' 
            #os.system(rand_command)
            #_mylog.info(rand_command)
            #parameters = parameters + FourierReconstructionExtraCommand
        #else:
            #_mylog.error("Reconstruction method unknown. Quiting")
            #print "Reconstruction method unknown. Quiting"
            #exit(1)
            
    runJob(_log
           , program
           , parameters
           , NumberOfMpi
           , NumberOfThreads
           )

#    if DisplayReconstruction==True:
#        command='xmipp_show -vol '+ Outputvolume + '&'
#        print '*********************************************************************'
#        print '* ',command
#        _mylog.info(command)
#        os.system(command)

#------------------------------------------------------------------------
#def  compute_resolution_and_filter(_log
#                                   , DoComputeResolution
#                                   , ARTReconstructionExtraCommand
#                                   , WBPReconstructionExtraCommand
#                                   , FourierReconstructionExtraCommand
#                                   , ReconstructionMethod
#                                   , FourierMaxFrequencyOfInterest
#                                   , iteration_number
#                                   #, DisplayReconstruction
#                                   , ResolSam
#                                   , DoParallel
#                                   , NumberOfMpi
#                                   , NumberOfThreads
#                                   , MyMpiJobSize
#                                   , SymmetryGroup
#                                   #, DisplayResolution
#                                   , ReconstructedVolume
#                                   , ARTLambda
#                                   , OuterRadius
#                                   , DoSplitReferenceImages
#                                   , PaddingFactor
#                                   , DoLowPassFilter
#                                   , UseFscForFilter
#                                   , ConstantToAddToFiltration
#                                   #, filter_frequence
#                                   #, ReconstructedVolume
#                                   #, ReconstructedandfilteredVolume
#                                   ):
#
#    ##
#    if (DoComputeResolution):
#        
#        #import os,shutil,math
#        PerformReconstruction=True
#        split_sel_root_name=ProjMatchRootName+'_split'
#        Outputvolumes=[]
#        Outputvolumes.append(split_sel_root_name+'_1')
#        Outputvolumes.append(split_sel_root_name+'_2')
#        
#        Selfiles=[]
#        Selfiles.append(split_sel_root_name+'_1_classes.sel')
#        Selfiles.append(split_sel_root_name+'_2_classes.sel')
#        Docfiles=[]
#        Docfiles.append(split_sel_root_name+'_1_classes.doc')
#        Docfiles.append(split_sel_root_name+'_2_classes.doc')
#        for i in range(len(Outputvolumes)):
#           print '*********************************************************************'
#           print '* Reconstruct volume'
#           if _ReconstructionMethod=='wbp':
#              program = 'xmipp_reconstruct_wbp'
#              parameters= ' -i '    + Selfiles[i] + \
#                          ' -doc '  + Docfiles[i] + \
#                          ' -o '    + Outputvolumes[i] + ".vol" + \
#                          ' --sym '  + _SymmetryGroup + \
#                          ' --weight --use_each_image '
#              parameters = parameters + _WBPReconstructionExtraCommand
#              _MyNumberOfThreads = 1
#           elif _ReconstructionMethod=='art':
#              program = 'xmipp_reconstruct_art'
#              _DoParallel=False
#              parameters=' -i '    + Selfiles[i] + \
#                         ' -o '    + Outputvolumes[i] + \
#                         ' -sym '  + _SymmetryGroup + \
#                 ' -thr '  + str(_MyNumberOfThreads) + \
#                         ' -WLS '
#              if len(_ARTLambda)>1:
#                 parameters = parameters + ' -l '   + _ARTLambda + ' '
#              parameters = parameters + _ARTReconstructionExtraCommand
#           elif _ReconstructionMethod=='fourier':
#              if ( _MyNumberOfMpiProcesses ==1):
#                  _DoParallel=False
#              program = 'xmipp_reconstruct_fourier'
#              parameters=' -i '    +  Selfiles[i] + \
#                         ' -o '    +  Outputvolumes[i] + '.vol ' + \
#                         ' --sym '  + _SymmetryGroup + \
#                 ' --thr '  + str(_MyNumberOfThreads) + \
#                         ' --weight ' + \
#                         ' --max_resolution ' + str(_FourierMaxFrequencyOfInterest) +\
#                 ' --pad_proj ' + str(_PaddingFactor) +\
#                 ' --pad_vol ' + str(_PaddingFactor)
#              if (_DoParallel):
#                 parameters = parameters + ' --mpi_job_size ' + str(_MyMpiJobSize)
#              if ( not _DoSplitReferenceImages):
#                  PerformReconstruction=False
#              parameters = parameters + _FourierReconstructionExtraCommand
#           else:
#              _mylog.error("Reconstruction method unknown. Quiting")
#              print "Reconstruction method unknown. Quiting"
#              exit(1)
#    
#           import launch_job
#           if(PerformReconstruction):
#               launch_job.launch_job(program,
#                                 parameters,
#                                 _mylog,
#                                 _DoParallel,
#                                 _MyNumberOfMpiProcesses,
#                                 _MyNumberOfThreads,
#                                 _MySystemFlavour)
#    
#        # Prevent high-resolution correlation because of discrete mask from wbp
#        innerrad = _OuterRadius - 2
#        for i in range(len(Outputvolumes)):
#           Outputvolumes[i]+=".vol"
#           print '*********************************************************************'
#           print '* Applying a soft mask'
#           command = " -i " + Outputvolumes[i] + \
#                     " --mask  raised_cosine -" + str(innerrad) + \
#                     " -" + str(_OuterRadius)
#           launch_job.launch_job("xmipp_mask",
#                                 command,
#                                 _mylog,
#                                 False,1,1,_MySystemFlavour)
#      
#        print '**************************************************************'
#        print '* Compute resolution ' 
#        command = " --ref " + Outputvolumes[0] +\
#                  " -i " +Outputvolumes[1]  + ' --sam ' + str(_ResolSam)
#        if ReconstructionMethod=='fourier':
#            import spider_header
#            myheader=spider_header.spiderheader(Outputvolumes[i] )
#            ncolumns=myheader.nx
#            #2.5 is the blob radius 
#            aux4 = 2.5 * 0.5 / ncolumns
#            command += " --max_sam " + str (_ResolSam/(aux4+_FourierMaxFrequencyOfInterest))
#        
#        launch_job.launch_job("xmipp_resolution_fsc",
#                              command,
#                              _mylog,
#                              False,1,1,_MySystemFlavour)
#        import visualization
#        if _DisplayResolution==True:
#          plot=visualization.gnuplot()
#          plot.plot_xy_file(Outputvolumes[1]+'.frc',
#                              Title="Resolution",
#                              X_Label="Armstrong^-1",
#                              Y_Label="y",
#                              X_col=1,
#                              Y_col=2)
#          print '*********************************************************************'
#          print '* plot resolution'
#          _mylog.info(" plot resolution")
#    
#        # Copy FSC to standard name file
#        outputfsc=_ReconstructedVolume.replace(ReconstructedVolume,OutputFsc)
#        shutil.copy(Outputvolumes[1]+'.frc',outputfsc) 
#    
#        #compute resolution
#        resolution_fsc_file = Outputvolumes[1]+'.frc'
#        f = open(resolution_fsc_file, 'r')
#        #skip first line
#        fi=f.readline()
#          
#        filter_frequence=0. 
#        for line in f:
#            line = line.strip()
#            if not line.startswith('#'):
#                mylist = (line.split())
#                if( float(mylist[1]) < 0.5):
#                   break
#                else:
#                  filter_frequence=float(mylist[0])
#    
#    
#        f.close()
#        print '* maximum resolution (A^-1): ', filter_frequence
#        filter_frequence *= _ResolSam
#        print '* maximum resolution (px^-1): ', filter_frequence
#        return filter_frequence
#    
#    else:
#        filter_frequence = 0
#    
#
##------------------------------------------------------------------------
##           filter_at_given_resolution
##------------------------------------------------------------------------
#    import os,shutil
#    import launch_job
#    Inputvolume   =_ReconstructedVolume+'.vol'
#    Outputvolume  =_ReconstructedandfilteredVolume+'.vol'
#
#    filter_in_pixels_at=0.5
#    if (not _DoLowPassFilter):
#        shutil.copy(Inputvolume,Outputvolume) 
#        command ="shutilcopy" + Inputvolume + ' ' + Outputvolume
#        _mylog.info(command)
#    else:   
#        print '**************************************************************'
#        print '* Filter reconstruction ' 
#        if (_UseFscForFilter):
#           filter_in_pixels_at = float(_filter_frequence) +\
#                                 float(_ConstantToAddToFiltration)
#        else:
#           filter_in_pixels_at = float(_ConstantToAddToFiltration)
#
#        if (filter_in_pixels_at>0.5):
#           shutil.copy(Inputvolume,Outputvolume) 
#           command ="shutilcopy" + Inputvolume + ' ' + Outputvolume
#           _mylog.info(command)
#           filter_in_pixels_at=0.5
#        else:
#           command = " -i " + Inputvolume +\
#                     " -o " + Outputvolume + ' --low_pass ' +\
#                     str (filter_in_pixels_at)
#           launch_job.launch_job("xmipp_fourier_filter",
#                                 command,
#                                 _mylog,
#                                 False,1,1,_MySystemFlavour)
#
#    #return filter_in_pixels_at


