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

#            _Parameters = {
#                          **       'DoCtfCorrection': DoCtfCorrection
#                                , 'ProjMatchDir' : ProjMatchDir
#                         **       , 'NumberOfCtfGroups':NumberOfCtfGroups
#                         **       , 'DocFileInputAngles':DocFileInputAngles[iterN]
#                                , 'DocFileOutAngles':DocFileInputAngles[iterN+1]
#                         **       , 'CtfGroupRootName': CtfGroupRootName
#                         **       , 'CtfGroupDirectory': CtfGroupDirectory
#                         **       , 'CtfGroupSubsetFileName':CtfGroupSubsetFileName
#                                }

def angular_project_library(_log,dict):
    _log.debug("execute_projection_matching")
#    ProjMatchDir=dict['ProjMatchDir']
    if dict['DoCtfCorrection']:
#        # To use -add_to in angular_class_average correctly, 
#        # make sure there are no proj_match_class* files from previous runs. 
#        print ' * CleanUp: deleting directory '+ ProjMatchDir
#        shutil.rmtree(ProjMatchDir, True)#true-> ignore error
#        mkpath(ProjMatchDir)
#        print ProjMatchDir
        # Create docfiles for each defocus group and corresponding selfile containing all of them      
        make_subset_docfiles(dict)

#def make_subset_docfiles(dict['mylog'],
#                         dict['InputDocFileName'],
#                         dict['NumberOfCtfGroups'],
#                         dict['CtfGroupRootName'],
#                         dict['CtfGroupDirectory'],
#                         dict['CtfGroupSubsetFileName']):

    # Project all references
    print '* Create projection library'
    refN=dict['refN']
    ###need one block per reference
    parameters=' -i '                    + dict['reconstructedFileName'] + \
               ' --experimental_images ' + dict['DocFileInputAngles'][refN] + \
               ' -o '                    + dict['ProjectLibraryRootName'] + '_'+str(refN)+\
               ' --sampling_rate '       + dict['AngSamplingRateDeg']  + \
               ' --sym '                 + dict['SymmetryGroup'] + 'h' + \
               ' --compute_neighbors'
    print parameters
#   if ( string.atof(_MaxChangeInAngles) < 181.):
#      parameters+= \
#              ' --near_exp_data --angular_distance '    + str(_MaxChangeInAngles)
#   else:
#      parameters+= \
#              ' --angular_distance -1'
#
#   if (_PerturbProjectionDirections):
#      # Just follow Roberto's suggestion
#      perturb=math.sin(math.radians(float(_AngSamplingRateDeg)))/4.
#      parameters+= \
#          ' --perturb ' + str(perturb)
#
#   if (_DoRetricSearchbyTiltAngle):
#     parameters+=  \
#              ' --min_tilt_angle '      + str(_Tilt0) + \
#              ' --max_tilt_angle '      + str(_TiltF)
#  
#   if (_DoCtfCorrection):
#     parameters+=  \
#              ' --groups '              + CtfGroupSubsetFileName
#
#   if (_DoParallel):
#      parameters = parameters + ' --mpi_job_size ' + str(_MyMpiJobSize)
#      
#   if (len(SymmetryGroupNeighbourhood)>1):
#      parameters+= \
#          ' --sym_neigh ' + SymmetryGroupNeighbourhood + 'h'
#   if (OnlyWinner):
#      parameters+= \
#              ' --only_winner '
#
#   launch_job.launch_job('xmipp_angular_project_library',
#                         parameters,
#                         _mylog,
#                         _DoParallel,
#                         _MyNumberOfMpiProcesses*_MyNumberOfThreads,
#                         1,
#                         _MySystemFlavour)

def make_subset_docfiles(_mylog,
                         _InputDocFileName,
                         _NumberOfCtfGroups,
                         _CtfGroupRootName,
                         _CtfGroupDirectory,
                         _CtfGroupSubsetFileName):
    
    # Loop over all CTF groups
    docselfile = []
    for ictf in range(_NumberOfCtfGroups):
       
        CtfGroupName=utils_xmipp.composeFileName(_CtfGroupRootName + '_group',ictf+1,'')
        CtfGroupName = '../' + _CtfGroupDirectory + '/' + CtfGroupName
        inselfile = CtfGroupName + '.sel'
        inputdocfile = (os.path.basename(inselfile)).replace('.sel','.doc')
        command=' --join ' + _InputDocFileName + ' ' + inselfile + ' --label image'
        command= command + ' -o ' + inputdocfile
        print '*********************************************************************'
        launch_job.launch_job("xmipp_metadata_utilities",
                      command,
                      _mylog,
                      False,1,1,'')
       
        docselfile.append(inputdocfile+' 1\n')
    
    # Write the selfile of all these docfiles
    fh = open(_CtfGroupSubsetFileName,'w')
    fh.writelines(docselfile)
    fh.close()
