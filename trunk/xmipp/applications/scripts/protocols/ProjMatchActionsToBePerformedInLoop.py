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


def projection_matching(_log,dict):
    # Loop over all CTF groups
    # Use reverse order to have same order in add_to docfiles from angular_class_average
    # get all ctf groups
    _DoCtfCorrection    = dict['DoCtfCorrection']
    _ProjMatchRootName  = dict['ProjMatchRootName']
    refname = str(dict['ProjectLibraryRootName'])
    NumberOfCtfGroups=dict['NumberOfCtfGroups']
    for ii in range(NumberOfCtfGroups):
        if NumberOfCtfGroups>1 :
            print 'Focus Group: ', ii+1,'/',NumberOfCtfGroups
        ictf    = NumberOfCtfGroups - ii 
        outputname   = _ProjMatchRootName
        if (_DoCtfCorrection):
            CtfGroupName = dict['CtfGroupRootName'] #,ictf+1,'')
            #outputname   = _ProjMatchRootName + '_' + CtfGroupName 
            CtfGroupName = dict['CtfGroupDirectory'] + '/' + CtfGroupName
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
        else:
            print "CORRECT THIS"
            exit(1)
            inputdocfile = dict['InputDocFileName']

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
            parameters = parameters + ' --mpi_job_size ' + str(dict['MyMpiJobSize'])
        
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
    DocFileInputAngles  = dict['DocFileInputAngles']#number of references
    ProjMatchRootName   = dict['ProjMatchRootName']#
    NumberOfCtfGroups   = dict['NumberOfCtfGroups']
    #print "cp",ProjMatchRootName[1],DocFileInputAngles
    #if number of references is one just copy file
    if(len(ProjMatchRootName)==2):#single reference
        shutil.copy(ProjMatchRootName[1], DocFileInputAngles[1])
    else:#multiple reference
        #add all ProjMatchRootName
        MDaux = MetaData()
        MD    = MetaData()
        MD1   = MetaData()
        print ProjMatchRootName

        for ii in range(1,NumberOfCtfGroups+1):
            ProjMatchRootNameIter = iter(ProjMatchRootName)
            ProjMatchRootNameIter.next()#skip first Null element
            element = ProjMatchRootNameIter.next()#skip first Null element
            MDaux.clear()
            MD1.clear()
            while True:
                fileNameIn    = str(ii).zfill(utils_xmipp.FILENAMENUMBERLENTGH) + '@' + element
                print ii, fileNameIn
                inputdocfile  = 'ctfGroup'+fileNameIn
                MD.read(inputdocfile)
                MDaux.unionAll(MD)
                try:
                    element = ProjMatchRootNameIter.next()
                except StopIteration:
                    break
            MD.aggregate(MDaux,AGGR_MAX,MDL_IMAGE,MDL_MAXCC,MDL_MAXCC)
            MD1.join(MD,MDaux,MDL_UNDEFINED,NATURAL)
        
            fileNameOut   = str(ii).zfill(utils_xmipp.FILENAMENUMBERLENTGH) + '@' + DocFileInputAngles[ii]
            outputdocfile = 'ref'+ fileNameOut
            MD1.write(outputdocfile,MD_APPEND)
            MDaux.write(outputdocfile+"aux",MD_APPEND)
            MD.write(outputdocfile+"MD",MD_APPEND)
        
        
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
