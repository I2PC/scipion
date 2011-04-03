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

##                                , 'ProjMatchDir' : ProjMatchDir
##                                , 'DocFileOutAngles':DocFileInputAngles[iterN]
#                                  'AngSamplingRateDeg':AngSamplingRateDeg[iterN]
#                                , 'CtfGroupRootName': CtfGroupRootName
#                                , 'CtfGroupDirectory': CtfGroupDirectory
#                                , 'CtfGroupSubsetFileName':CtfGroupSubsetFileName
#                                , 'DoCtfCorrection': DoCtfCorrection
#                                , 'DocFileInputAngles':DocFileInputAngles[iterN-1][refN]
#                                , 'NumberOfCtfGroups':NumberOfCtfGroups
#                                , 'ProjectLibraryRootName':ProjectLibraryRootNames[iterN][refN]
#                                , 'reconstructedFileName':reconstructedFileNamesIter[iterN-1][refN]
#                                , 'SymmetryGroup':SymmetryGroup
#                                }

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

    if (dict['DoRetricSearchbyTiltAngle']):
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
    a=0