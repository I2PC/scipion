from distutils.errors   import DistutilsInternalError
from types              import StringTypes
import os, sys
from xmipp import *

def createDirectoryTree(dirName):
    """ Create directory or directory branch"""
    from distutils.dir_util import mkpath
    #
    if not isinstance(dirName, StringTypes):
        raise DistutilsInternalError, "mkpath: 'name' must be a string (got %r)" % (dirName,)
    mkpath(dirName, 0777, True)

def createRequiredDirectories(_log, dict):
    """ Create all required directories, I do not want to bother about that 
        later """
    absPath = dict['ProjectDir'] + "/" + dict['WorkingDir']
    for i in range(1, dict['NumberofIterations'] + 1):
        createDirectoryTree(absPath + "/Iter_" + str(i).zfill(2))
    _log.info("Creating directory tree " + absPath)

def deleteWorkingDirectory(_mylog, dict):

    if (not dict['DoDeleteWorkingDir']):
        return
    if not dict['firstIteration'] :
        print "You can not delete the working directory"
        print " and start at a iteration different from 1 "
        exit(1)

    from distutils.dir_util import remove_tree
    dirName = dict['ProjectDir'] + "/" + dict['WorkingDir']
    _mylog.info("Delete working directory tree " + dirName)
    if not isinstance(dict['WorkingDir'], StringTypes):
        raise DistutilsInternalError, "mkpath: 'name' must be a string (got %r)" % (dirName,)
    if os.path.exists(dirName):
        remove_tree(dirName, True)

def checkVolumeProjSize(_log, dict):

    """ check references y projection size match"""
    #if this is not the firt iteration you may skip it
    if not dict['firstIteration']:
        return
    #5a check volumes have same size
    try:
        (xdim, ydim, zdim, ndim) = SingleImgSize(dict['ReferenceFileNames'][0])
        for reference in dict['ReferenceFileNames']:
            (xdim2, ydim2, zdim2, ndim2) = SingleImgSize(reference)
    except XmippError as e:
        print __name__, 'checkVolumeProjSize:', e
        exit(1)
    if (xdim2, ydim2, zdim2, ndim2) != (xdim, ydim, zdim, ndim):
        print "Reference %s and %s have not the same size" % \
              (dict['ReferenceFileNames'][0], reference)
        exit(1)

    #5b check volume and projections  have same size
    (xdim2, ydim2, zdim2, ndim2) = ImgSize(dict['SelFileName'])
    if (xdim2, ydim2) != (xdim, ydim):
            print "Volume and reference images have not the same size"
            exit(1)
    _log.debug("checkVolumeProjSize")


def initOuterRadius(_log, dict):
    """ init mask radius for volumes. If suggested value is negative set to half dim"""
    OuterRadius = dict['OuterRadius']

#    print "REMOVE -34"
#    OuterRadius = -34

    xdim = 0
    if (OuterRadius < 0):
        try:
            (xdim, ydim, zdim, ndim) = ImgSize(dict['SelFileName'])
        except XmippError as e:
            print __name__, 'initOuterRadius:', e
            exit(1)
        OuterRadius = (xdim / 2) - 1
        comment = "InitOuterRadius: Outer radius set to: " + str(OuterRadius)
        print '* ' + comment
        _log.info(comment)
    return OuterRadius

# wrapper for log.make_backup_of_script_file
def pm_make_backup_of_script_file(_log, dict):
    import log
    _log.info("make backup script file: " + dict['progName'])
    log.make_backup_of_script_file(dict['progName'], dict['ProjectDir'] + "/" + dict['WorkingDir'])

#------------------------------------------------------------------------
#make ctf groups
#------------------------------------------------------------------------
#    'CTFDatName': CTFDatName
#  , 'CtfGroupDirectory': CtfGroupDirectory
#  , 'CtfGroupMaxDiff': CtfGroupMaxDiff
#  , 'CtfGroupMaxResol': CtfGroupMaxResol
#  , 'CtfGroupRootName': CtfGroupRootName
#  , 'DataArePhaseFlipped': DataArePhaseFlipped
#  , 'DoAutoCtfGroup'      : DoAutoCtfGroup
#  , 'DoCtfCorrection'      : DoCtfCorrectio
#  , '_log'                :_log
#  , 'PaddingFactor'       : PaddingFactor
#  , 'SelFileName'         : SelFileName
#  , 'SplitDefocusDocFile' : SplitDefocusDocFile
#  , 'WienerConstant'      : WienerConstant
def execute_ctf_groups (_log, dict):

    if (not dict['DoCtfCorrection']):
        return 1
    import glob, sys
    import utils_xmipp
    import launch_job
    if not os.path.exists(dict['CtfGroupDirectory']):
        os.makedirs(dict['CtfGroupDirectory'])

    print '*********************************************************************'
    print '* Make CTF groups'
    command = \
              ' --ctfdat ' + dict['CTFDatName' ] + \
              ' -o ' + dict['CtfGroupDirectory'] + '/' + dict['CtfGroupRootName'] + \
              ' --wiener --wc ' + str(dict['WienerConstant']) + \
              ' --pad ' + str(dict['PaddingFactor'])

    if (dict['DataArePhaseFlipped']):
        command += ' --phase_flipped '

    if (dict['DoAutoCtfGroup']):
        command += ' --error ' + str(dict['CtfGroupMaxDiff']) + \
                   ' --resol ' + str(dict['CtfGroupMaxResol'])
    else:
        if (len(dict['SplitDefocusDocFile']) > 0):
            command += ' --split ' + dict['SplitDefocusDocFile']
        else:
            message = "Error: for non-automated ctf grouping, please provide a docfile!"
            print '* ', message
            _log.info(message)
            sys.exit()

    (_command, retcode) = launch_job.launch_job("xmipp_ctf_group",
                          command,
                          _log,
                          False, 1, 1, '')

    if(retcode):
        print "command", command, "failed with exit status", retcode
        exit(1)
    wildcardname = utils_xmipp.composeWildcardFileName(dict['CtfGroupDirectory'] + '/'
                                                     + dict['CtfGroupRootName']
                                                     + '_group', 'ctf')
    ctflist = glob.glob(wildcardname)
    return len(ctflist)

def checkOptionsCompatibility(_log, dict):
    import arg
    # Never allow DoAlign2D and DoCtfCorrection together
    if (int(arg.getComponentFromVector(dict['DoAlign2D'], 0)) and int (dict['DoCtfCorrection'])):
        error_message = "You cannot realign classes AND perform CTF-correction. Switch either of them off!"
        dict['_log'].error(error_message)
        print error_message
        exit(1)

def dummy(_log, dict):
    print dict['dummy']
