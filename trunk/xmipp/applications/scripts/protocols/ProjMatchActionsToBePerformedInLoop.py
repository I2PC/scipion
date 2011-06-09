import os, shutil, string, glob, math
import launch_job, utils_xmipp
from distutils.dir_util import mkpath
from xmipp import *

CtfBlockName = 'ctfGroup'

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
        print "a1"
        src=dict['ProjectLibraryRootName'].replace(".stk",'_sampling.xmd')
        dst = src.replace('sampling.xmd','group'+
                              str(1).zfill(utils_xmipp.FILENAMENUMBERLENTGH)+
                              '_sampling.xmd')
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
        inputdocfile    = CtfBlockName+str(ictf).zfill(utils_xmipp.FILENAMENUMBERLENTGH) + '@' + CtfGroupName + '_images.sel'
        outputname   = CtfBlockName+str(ictf).zfill(utils_xmipp.FILENAMENUMBERLENTGH) + '@'+ _ProjMatchRootName
        #inputdocfile = (os.path.basename(inselfile)).replace('.sel','.doc')
        baseTxtFile  = refname[:-len('.stk')] 
        neighbFile      = baseTxtFile + '_sampling.xmd'
        if (os.path.exists(neighbFile)):
            os.remove(neighbFile)
        neighbFileb     = utils_xmipp.composeFileName(baseTxtFile + '_group',ictf,'')
        neighbFileb     += '_sampling.xmd'
        shutil.copy(neighbFileb, neighbFile)
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
                    ' --thr '            + str(dict['NumberOfThreads']) +\
                    ' --append '

        
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
    #first we need a list with the references used. That is,
    #read all docfiles and map referecendes to a mdl_order
    MDaux  = MetaData()
    MDSort = MetaData()
    MD     = MetaData()
    MD1    = MetaData()
    MDout  = MetaData()
    MDout.setComment("metadata with  images, the winner reference as well as the ctf group")

    mycounter=0L
    for iCTFGroup in range(1,NumberOfCtfGroups+1):
        auxInputdocfile = CtfBlockName + str(iCTFGroup).zfill(utils_xmipp.FILENAMENUMBERLENTGH)+'@'
        for iRef3D in range(1,NumberOfReferences+1):
            inputFileName = ProjMatchRootName[iRef3D]
            inputdocfile    = auxInputdocfile+ inputFileName
            MD.read(inputdocfile)
            for id in MD:
                t=MD.getValue(MDL_REF,id)
                i=MDSort.addObject()
                MDSort.setValue(MDL_REF,t,i)
    MDSort.removeDuplicates()
    for id in MDSort:
        MDSort.setValue(MDL_ORDER,mycounter,id)
        mycounter += 1
    #print "bbb",ProjMatchRootName[1], DocFileInputAngles
    outputdocfile =  DocFileInputAngles
    if os.path.exists(outputdocfile):
        os.remove(outputdocfile)
    for iCTFGroup in range(1,NumberOfCtfGroups+1):
        MDaux.clear()
        auxInputdocfile = CtfBlockName + str(iCTFGroup).zfill(utils_xmipp.FILENAMENUMBERLENTGH)+'@'
        for iRef3D in range(1,NumberOfReferences+1):
            inputFileName = ProjMatchRootName[iRef3D]
            inputdocfile    = auxInputdocfile+ inputFileName
            MD.clear()
            MD.read(inputdocfile)
            #In practice you should not get duplicates
            MD.removeDuplicates()
            MD.setValueCol(MDL_REF3D,iRef3D)
            MDaux.unionAll(MD)
        MDaux.sort()
        MD.aggregate(MDaux,AGGR_MAX,MDL_IMAGE,MDL_MAXCC,MDL_MAXCC)
        #if a single image is assigned to two references with the same 
        #CC use it in both reconstruction
        #recover atributes after aggregate function
        MD1.join(MD,MDaux,MDL_UNDEFINED,NATURAL)
        MD1.write('MD1.xmd')
        MDSort.write('MDSort.xmd')
        #add a sorting number to make easier to create an stack of averaged classes
        MDout.join(MD1,MDSort,MDL_UNDEFINED,NATURAL)        
        MDout.write(auxInputdocfile+outputdocfile,MD_APPEND)

def angular_class_average(_log,dict):
    # Now make the class averages
    CtfGroupName        = dict['CtfGroupDirectory'] + '/' + dict['CtfGroupRootName']
    #DocFileInputAngles  = dict['DocFileInputAngles']#
    DoCtfCorrection     = dict['DoCtfCorrection']
    NumberOfCtfGroups   = dict['NumberOfCtfGroups']
    #NumberOfReferences  = dict['NumberOfReferences']
    #ProjMatchRootName   = dict['ProjMatchRootName']
    refname = str(dict['ProjectLibraryRootName'])

    
    #Md=MetaData()
    #MdSelect=MetaData()
    ProjMatchRootName = dict['ProjMatchRootName']
    for iCTFGroup in range(1,NumberOfCtfGroups+1):
        tmpFileName = CtfBlockName + str(iCTFGroup).zfill(utils_xmipp.FILENAMENUMBERLENTGH)+'@'
        #extract from metadata relevant images
        tmpFileName += ProjMatchRootName
        #Md.write("test.xmd" + str(iCTFGroup).zfill(2) +'_'+str(iRef3D).zfill(2))
        parameters =  ' -i '      + tmpFileName  + \
                      ' --lib '    + refname.replace(".stk",".doc") + \
                      ' --dont_write_selfiles ' + \
                      ' --limit0 ' + dict['MinimumCrossCorrelation'] + \
                      ' --limitR ' + dict['DiscardPercentage']
        if (DoCtfCorrection):
            # On-the fly apply Wiener-filter correction and add all CTF groups together
            parameters += \
                       ' --wien '   + str(iCTFGroup).zfill(utils_xmipp.FILENAMENUMBERLENTGH)+'@' + CtfGroupName + '_wien.stk' + \
                       ' --pad '    + str(dict['PaddingFactor']) + \
                       ' --add_to ' + ProjMatchRootName
        else:
            parameters += \
                      ' -o '                + ProjMatchRootName
        if (dict['DoAlign2D'] == '1'):
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
