import os, shutil, string, glob, math
import launch_job, utils_xmipp
from distutils.dir_util import mkpath

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
              ' --sampling_rate '       + str(dict['AngSamplingRateDeg'])  + \
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
    if (dict['OnlyWinner']!='0'):
        parameters+= \
              ' --only_winner '

    launch_job.launch_job('xmipp_angular_project_library',
                         parameters,
                         _log,
                         _DoParallel,
                         dict['NumberOfMpiProcesses']*dict['NumberOfThreads'],
                         1,
                         dict['SystemFlavour'])


def projection_matching(_log,dict):
    # Loop over all CTF groups
    # Use reverse order to have same order in add_to docfiles from angular_class_average
    # get all ctf groups
   
    print dict['NumberOfCtfGroups']
    for ii in range(dict['NumberOfCtfGroups']):
        ictf = dict['NumberOfCtfGroups'] - ii + 1

#      refname          = ProjectLibraryRootName
#      if (_DoCtfCorrection):
#         CtfGroupName = utils_xmipp.composeFileName(CtfGroupRootName + '_group',ictf+1,'')
#         outputname   = ProjMatchRootName + '_' + CtfGroupName 
#         CtfGroupName = '../' + CtfGroupDirectory + '/' + CtfGroupName
#         inselfile    = CtfGroupName + '.sel'
#         inputdocfile = (os.path.basename(inselfile)).replace('.sel','.doc')
#         txtfile      = ProjectLibraryRootName + '_sampling.txt'
#         if (os.path.exists(txtfile)):
#            os.remove(txtfile)
#         txtfileb     = utils_xmipp.composeFileName(ProjectLibraryRootName + '_group',ictf+1,'')
#         txtfileb     += '_sampling.txt'
#         shutil.copy(txtfileb, txtfile)
#      else:
#         outputname   = ProjMatchRootName
#         inputdocfile = _InputDocFileName
#
#      print '*********************************************************************'
#      print '* Perform projection matching'
#      parameters= ' -i '              + inputdocfile + \
#                  ' -o '              + outputname + \
#                  ' --ref '            + refname + \
#                  ' --Ri '             + str(_Ri)           + \
#                  ' --Ro '             + str(_Ro)           + \
#                  ' --max_shift '      + str(_MaxChangeOffset) + \
#                  ' --search5d_shift ' + str(_Search5DShift) + \
#                  ' --search5d_step  ' + str(_Search5DStep) + \
#                  ' --mem '            + str(_AvailableMemory * _MyNumberOfThreads) + \
#                  ' --thr '            + str(_MyNumberOfThreads)
#      
#      if (_DoScale):
#         parameters += \
#                  ' --scale '          + str(_ScaleStep) + ' ' + str(_ScaleNumberOfSteps) 
#
#      if (_DoCtfCorrection and _ReferenceIsCtfCorrected):
#         ctffile = CtfGroupName + '.ctf'
#         parameters += \
#                  ' --pad '            + str(_PaddingFactor) + \
#                  ' --ctf '            + ctffile
#
#      if (_DoParallel):
#         parameters = parameters + ' --mpi_job_size ' + str(_MyMpiJobSize)
#
#      launch_job.launch_job('xmipp_angular_projection_matching',
#                            parameters,
#                            _mylog,
#                            _DoParallel,
#                            _MyNumberOfMpiProcesses,
#                            _MyNumberOfThreads,
#                            _MySystemFlavour)
#
#      # Now make the class averages
#      parameters =  ' -i '      + outputname + '.doc'  + \
#                    ' --lib '    + ProjectLibraryRootName + '.doc' + \
#                    ' --dont_write_selfiles ' + \
#                    ' --limit0 ' + str(_MinimumCrossCorrelation) + \
#                    ' --limitR ' + str(_DiscardPercentage)
#      if (_DoCtfCorrection):
#         # On-the fly apply Wiener-filter correction and add all CTF groups together
#         parameters += \
#                    ' --wien '             + CtfGroupName + '.wien' + \
#                    ' --pad '              + str(_PaddingFactor) + \
#                    ' --add_to '           + ProjMatchRootName
#      else:
#         parameters += \
#                    ' -o '                + ProjMatchRootName
#      if (_DoAlign2D == '1'):
#         parameters += \
#                    ' --iter '             + str(_Align2DIterNr) + \
#                    ' --Ri '               + str(_Ri)           + \
#                    ' --Ro '               + str(_Ro)           + \
#                    ' --max_shift '        + str(_MaxChangeOffset) + \
#                    ' --max_shift_change ' + str(_Align2dMaxChangeOffset) + \
#                    ' --max_psi_change '   + str(_Align2dMaxChangeRot) 
#      if (_DoComputeResolution and _DoSplitReferenceImages):
#         parameters += \
#                    ' --split '
#
#      launch_job.launch_job('xmipp_angular_class_average',
#                            parameters,
#                            _mylog,
#                            _DoParallel,
#                            _MyNumberOfMpiProcesses*_MyNumberOfThreads,
#                            1,
#                            _MySystemFlavour)
#
#      if (_DoAlign2D == '1'):
#         outputdocfile =  ProjMatchRootName + '_realigned.doc'
#      else:
#         outputdocfile =  ProjMatchRootName + '.doc'
#
#      if (_DoCtfCorrection):
#         os.remove(outputname + '.doc')
#         os.remove(inputdocfile)
