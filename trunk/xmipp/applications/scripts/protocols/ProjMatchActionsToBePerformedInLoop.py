import os, shutil, string, glob, math
import launch_job, utils_xmipp
from distutils.dir_util import mkpath
from xmipp import *

def execute_mask(_log, dict):
    _mylog = _log
    _mylog.debug("execute_mask")
    if (dict['DoMask']):
#        print '*********************************************************************'
#        print '* Mask the reference volume'
        command = ' -i ' + dict['reconstructedFileName'] + \
                  ' -o ' + dict['maskedFileName']

        if (dict['DoSphericalMask']):
            command += ' --mask circular -' + str(dict['maskRadius'])
        else:
            command += ' --mask ' + dict['userSuppliedMask']

        launch_job.launch_job("xmipp_transform_mask",
                              command,
                              _mylog,
                              False, 1, 1, '')

    else:
        shutil.copy(dict['reconstructedFileName'], dict['maskedFileName'])
        _mylog.info("Skipped Mask")
        _mylog.info("cp " + dict['reconstructedFileName'] + " " + dict['maskedFileName'])
#        print '*********************************************************************'
        print '* Skipped Mask'


def angular_project_library(_log,dict):
    _log.debug("execute_projection_matching")
    ###need one block per reference
    # Project all references
    print '* Create projection library'
    parameters=' -i '                   + dict['maskedFileNamesIter'] + \
              ' --experimental_images ' + dict['DocFileInputAngles'] + \
              ' -o '                    + dict['ProjectLibraryRootName'] + \
              ' --sampling_rate '       + dict['AngSamplingRateDeg']  + \
              ' --sym '                 + dict['SymmetryGroup'] + 'h' + \
              ' --compute_neighbors'

    if ( string.atof(dict['MaxChangeInAngles']) < 181.):
        parameters+= \
              ' --near_exp_data --angular_distance '    + str(dict['MaxChangeInAngles'])
    else:
        parameters+= \
              ' --angular_distance -1'

    if (dict['PerturbProjectionDirections']):
        perturb=math.sin(math.radians(float(dict['AngSamplingRateDeg'])))/4.
        parameters+= \
           ' --perturb ' + str(perturb)

    if (dict['DoRestricSearchbyTiltAngle']):
        parameters+=  \
              ' --min_tilt_angle '      + str(dict['Tilt0']) + \
              ' --max_tilt_angle '      + str(dict['TiltF'])

    if (dict['DoCtfCorrection']):
        parameters+=  \
              ' --groups '              + dict['CtfGroupSubsetFileName']
    _DoParallel=dict['DoParallel']
    if (dict['DoParallel']):
        parameters = parameters + ' --mpi_job_size ' + str(dict['MpiJobSize'])
    if (len(dict['SymmetryGroupNeighbourhood'])>1):
        parameters+= \
          ' --sym_neigh ' + dict['SymmetryGroupNeighbourhood'] + 'h'
    if (dict['OnlyWinner']):
        parameters+= \
              ' --only_winner '

    launch_job.launch_job('xmipp_angular_project_library',
                         parameters,
                         _log,
                         _DoParallel,
                         dict['NumberOfMpiProcesses']*dict['NumberOfThreads'],
                         1,
                         dict['SystemFlavour'])
    if (not dict['DoCtfCorrection']):
        src=dict['ProjectLibraryRootName'].replace(".stk",'_sampling.txt')
        dst = src.replace('sampling.txt','group'+
                              str(1).zfill(utils_xmipp.FILENAMENUMBERLENTGH)+
                              '_sampling.txt')
        print "aaa", src,dst
        shutil.copy(src, dst)


def projection_matching(_log,dict):
    # Loop over all CTF groups
    # Use reverse order to have same order in add_to docfiles from angular_class_average
    # get all ctf groups
    _DoCtfCorrection    = dict['DoCtfCorrection']
    _ProjMatchRootName  = dict['ProjMatchRootName']
    refname = str(dict['ProjectLibraryRootName'])
    NumberOfCtfGroups=dict['NumberOfCtfGroups']
    
    CtfGroupName = dict['CtfGroupRootName'] #,ictf+1,'')
    CtfGroupName = dict['CtfGroupDirectory'] + '/' + CtfGroupName
    #remove output metadata
    if os.path.exists(_ProjMatchRootName):
        os.remove(_ProjMatchRootName)
    
    for ii in range(NumberOfCtfGroups):
        if NumberOfCtfGroups>1 :
            print 'Focus Group: ', ii+1,'/',NumberOfCtfGroups
        ictf    = NumberOfCtfGroups - ii 
#        if (_DoCtfCorrection):
        #outputname   = _ProjMatchRootName + '_' + CtfGroupName 
        inputdocfile    = 'ctfGroup'+str(ictf).zfill(utils_xmipp.FILENAMENUMBERLENTGH) + '@' + CtfGroupName + '_images.sel'
        outputname   = 'ctfGroup'+str(ictf).zfill(utils_xmipp.FILENAMENUMBERLENTGH) + '@'+ _ProjMatchRootName
        #inputdocfile = (os.path.basename(inselfile)).replace('.sel','.doc')
        baseTxtFile  = refname[:-len('.stk')] 
        txtfile      = baseTxtFile + '_sampling.txt'
        if (os.path.exists(txtfile)):
            os.remove(txtfile)
        txtfileb     = utils_xmipp.composeFileName(baseTxtFile + '_group',ictf,'')
        txtfileb     += '_sampling.txt'
        shutil.copy(txtfileb, txtfile)
#        else:
#            inputdocfile    = 'ctfGroup'+str(1).zfill(utils_xmipp.FILENAMENUMBERLENTGH) + '@' + CtfGroupName + '_images.sel'
#            outputname   = 'ctfGroup'+str(1).zfill(utils_xmipp.FILENAMENUMBERLENTGH) + '@'+ _ProjMatchRootName

        parameters= ' -i '               + inputdocfile + \
                    ' -o '               + outputname + \
                    ' --ref '            + refname + \
                    ' --Ri '             + str(dict['InnerRadius'])           + \
                    ' --Ro '             + str(dict['OuterRadius'])           + \
                    ' --max_shift '      + str(dict['MaxChangeOffset']) + \
                    ' --search5d_shift ' + str(dict['Search5DShift']) + \
                    ' --search5d_step  ' + str(dict['Search5DStep']) + \
                    ' --mem '            + str(dict['AvailableMemory'] * dict['NumberOfThreads']) + \
                    ' --thr '            + str(dict['NumberOfThreads'])
        
        if (dict['DoScale']):
            parameters += \
                    ' --scale '          + str(dict['ScaleStep']) + ' ' + str(dict['ScaleNumberOfSteps']) 
        
        if (_DoCtfCorrection and dict['ReferenceIsCtfCorrected']):
            ctffile = str(ictf).zfill(utils_xmipp.FILENAMENUMBERLENTGH) + '@' + CtfGroupName + '_ctf.stk'
            parameters += \
                      ' --pad '            + str(dict['PaddingFactor']) + \
                      ' --ctf '            + ctffile
        
        if (dict['DoParallel']):
            parameters = parameters + ' --mpi_job_size ' + str(dict['MpiJobSize'])
        
        launch_job.launch_job('xmipp_angular_projection_matching',
                            parameters,
                            _log,
                            dict['DoParallel'],
                            dict['NumberOfMpiProcesses'],
                            dict['NumberOfThreads'],
                            dict['SystemFlavour'])
        
def assign_images_to_references(_log,dict):
    ''' assign the images to the different references based on the crosscorrelation coeficient
        #if only one reference it just copy the docfile generated in the previous step
        '''
    DocFileInputAngles  = dict['DocFileInputAngles']
    ProjMatchRootName   = dict['ProjMatchRootName']#
    NumberOfCtfGroups   = dict['NumberOfCtfGroups']
    NumberOfReferences  = dict['NumberOfReferences']

    #print "bbb",ProjMatchRootName[1], DocFileInputAngles
    MDaux = MetaData()
    MD    = MetaData()
    MD1   = MetaData()

    for iRef3D in range(1,NumberOfReferences+1):
        inputFileName = ProjMatchRootName[iRef3D]#skip first Null element
        for iCTFGroup in range(1,NumberOfCtfGroups+1):
###                MDaux.clear()
###                MD1.clear()
                inputdocfile    = 'ctfGroup' + \
                                  str(iCTFGroup).zfill(utils_xmipp.FILENAMENUMBERLENTGH) + \
                                  '@' + inputFileName
                MD.clear()
                MD.read(inputdocfile)
                MD.setValueCol(MDL_CTF_GROUP,iCTFGroup)
                MD.setValueCol(MDL_REF3D,iRef3D)
                MDaux.unionAll(MD)
    MD.aggregate(MDaux,AGGR_MAX,MDL_IMAGE,MDL_MAXCC,MDL_MAXCC)
    MD1.join(MD,MDaux,MDL_UNDEFINED,NATURAL)
        
    outputdocfile =  DocFileInputAngles
    MD1.setComment("metadata with  images, the winner reference as well as the ctf group")
    MD1.write(outputdocfile)
    print "www",outputdocfile
            
def angular_class_average(_log,dict):
    # Now make the class averages
    DocFileInputAngles  = dict['DocFileInputAngles']#number of references
    NumberOfCtfGroups   = dict['NumberOfCtfGroups']
    CtfGroupName        = dict['CtfGroupDirectory'] + '/' + dict['CtfGroupRootName']
    ProjMatchRootName   = dict['ProjMatchRootName']
    NumberOfReferences  = len(ProjMatchRootName)
    refname = str(dict['ProjectLibraryRootName'])

    
    Md=MetaData()
    MdSelect=MetaData()
    for iRef3D in range(1,NumberOfReferences):#already has plus 1
        ProjMatchRootName = dict['ProjMatchRootName'][iRef3D]
        for iCTFGroup in range(1,NumberOfCtfGroups+1):
            #extract from metadata relevant images
            Md.read(DocFileInputAngles)
            MdSelect.importObjects(Md, MDValueEQ(MDL_REF3D, iRef3D))
            Md.clear()
            Md.importObjects(MdSelect,MDValueEQ(MDL_CTF_GROUP, iCTFGroup))
            tmpFileName=DocFileInputAngles+"_tmp"
            Md.write(tmpFileName)
            #Md.write("test.xmd" + str(iCTFGroup).zfill(2) +'_'+str(iRef3D).zfill(2))
            parameters =  ' -i '      + tmpFileName  + \
                          ' --lib '    + refname + '.doc' + \
                          ' --dont_write_selfiles ' + \
                          ' --limit0 ' + dict['MinimumCrossCorrelation'] + \
                          ' --limitR ' + dict['DiscardPercentage']
            if (dict['DoCtfCorrection']):
                # On-the fly apply Wiener-filter correction and add all CTF groups together
                parameters += \
                           ' --wien '   + CtfGroupName + '.wien' + \
                           ' --pad '    + str(dict['PaddingFactor']) + \
                           ' --add_to ' + ProjMatchRootName
            else:
                parameters += \
                          ' -o '                + ProjMatchRootName
            if (dict['Align2DIterNr'] == '1'):
                parameters += \
                          ' --iter '             + dict['Align2DIterNr']  + \
                          ' --Ri '               + str(dict['InnerRadius'])           + \
                          ' --Ro '               + str(dict['OuterRadius'])           + \
                          ' --max_shift '        + dict['MaxChangeOffset'] + \
                          ' --max_shift_change ' + dict['Align2dMaxChangeOffset'] + \
                          ' --max_psi_change '   + dict['Align2dMaxChangeRot'] 
    if (dict['DoComputeResolution'] and dict['DoSplitReferenceImages']):
        parameters += \
                  ' --split '
    
    launch_job.launch_job('xmipp_angular_class_average',
                          parameters,
                          _log,
                          dict['DoParallel'],
                          dict['NumberOfMpiProcesses'] * dict['NumberOfThreads'],
                          1,
                          dict['SystemFlavour'])
######    
######    if (_DoAlign2D == '1'):
######       outputdocfile =  ProjMatchRootName + '_realigned.doc'
######    else:
######       outputdocfile =  ProjMatchRootName + '.doc'
######    
######    if (_DoCtfCorrection):
######       os.remove(outputname + '.doc')
######       os.remove(inputdocfile)

#            shutil.copy(ProjMatchRootName[1], DocFileInputAngles)

#
#        # Now make the class averages
#        parameters =  ' -i '      + outputname + '.doc'  + \
#                    ' --lib '    + refname[:-len('.stk')] + '.doc' + \
#                    ' --dont_write_selfiles ' + \
#                    ' --limit0 ' + str(dict['MinimumCrossCorrelation']) + \
#                    ' --limitR ' + str(dict['DiscardPercentage'])
#        if (_DoCtfCorrection):
#            # On-the fly apply Wiener-filter correction and add all CTF groups together
#            parameters += \
#                       ' --wien '             + CtfGroupName + '.wien' + \
#                       ' --pad '              + str(dict['PaddingFactor']) + \
#                       ' --add_to '           + _ProjMatchRootName
#        else:
#            parameters += \
#                        ' -o '                + _ProjMatchRootName
#        if (dict['DoAlign2D']):
#            parameters += \
#                        ' --iter '             + str(_Align2DIterNr) + \
#                        ' --Ri '               + str(_Ri)           + \
#                        ' --Ro '               + str(_Ro)           + \
#                        ' --max_shift '        + str(_MaxChangeOffset) + \
#                        ' --max_shift_change ' + str(_Align2dMaxChangeOffset) + \
#                        ' --max_psi_change '   + str(_Align2dMaxChangeRot) 
#        if (dict['DoComputeResolution'] and dict['DoSplitReferenceImages']):
#            parameters += ' --split '
#
#        launch_job.launch_job('xmipp_angular_class_average',
#                            parameters,
#                            _log,
#                            dict['DoParallel'],
#                            dict['NumberOfMpiProcesses']*dict['NumberOfThreads'],
#                            1,
#                            dict['SystemFlavour'])
#
#      if (_DoAlign2D == '1'):
#         outputdocfile =  ProjMatchRootName + '_realigned.doc'
#      else:
#         outputdocfile =  ProjMatchRootName + '.doc'
#
#      if (_DoCtfCorrection):
#         os.remove(outputname + '.doc')
#         os.remove(inputdocfile)
