from distutils.errors   import DistutilsInternalError
from types              import StringTypes
import os
from xmipp import *

def createDir(_log, dict):
    """ Create directory """
    from distutils.dir_util import mkpath
    from distutils.errors import DistutilsFileError
    i = dict['Iter']
    if i == 0:
        _Path = dict['ProjectDir'] + "/" + dict['WorkingDir'] + "/"
    elif i > 0:
        _Path = dict['WorkingDir'] + "/" + "Iter_" + str(i).zfill(2)
    try:
        mkpath(_Path, 0777, True)
    except DistutilsFileError, e:
        print "could not create '%s': %s" % (os.path.abspath(_Path), e)
        exit(1)
    _log.info("Create directory " + _Path)

def createDir2(_log, dict):
    """ Create directory no add workingdir"""
    from distutils.dir_util import mkpath
    from distutils.errors import DistutilsFileError
    _path=dict['path']
    try:
        mkpath(_path, 0777, True)
    except DistutilsFileError, e:
        print "could not create '%s': %s" % (os.path.abspath(_path), e)
        exit(1)
    _log.info("Create directory " + _path)

def changeDir(_log, dict):
    """ Change to Directory """
    _Path = dict['ProjectDir'] + "/" + dict['WorkingDir'] + "/"
    _log.info("Create directory " + _Path)
    try:
        os.chdir(_Path)
    except os.error, (errno, errstr):
        print "could not change to directory '%s'" % _Path
        print "Error(%d): %s" % (errno, errstr)
        exit(1)
def deleteWorkingDirectory(_mylog, dict):

    if (not dict['DoDeleteWorkingDir']):
        return
    #in theory this cannot happen
#    if not dict['firstIteration'] :
#        print "You can not delete the working directory"
#        print " and start at a iteration different from 1 "
#        exit(1)

    from distutils.dir_util import remove_tree
    dirName = dict['ProjectDir'] + "/" + dict['WorkingDir']
    _mylog.info("Delete working directory tree " + dirName)
    if not isinstance(dict['WorkingDir'], StringTypes):
        raise DistutilsInternalError, "mkpath: 'name' must be a string (got %r)" % (dirName,)
    if os.path.exists(dirName):
        remove_tree(dirName, True)

def checkVolumeProjSize(_log, dict):
    """ check references and projection size match"""
    #5a check volumes have same size
    result = True
    message=""
    try:
        (xdim, ydim, zdim, ndim) = SingleImgSize(dict['ReferenceFileNames'][0])
        for reference in dict['ReferenceFileNames']:
            (xdim2, ydim2, zdim2, ndim2) = SingleImgSize(reference)
    except XmippError as e:
        print __name__, 'checkVolumeProjSize:', e
        exit(False,message)
    if (xdim2, ydim2, zdim2, ndim2) != (xdim, ydim, zdim, ndim):
        message = "Reference %s and %s have not the same size" % \
              (dict['ReferenceFileNames'][0], reference) 
        result = False
    
    if result:
        #5b check volume and projections  have same size
        (xdim2, ydim2, zdim2, ndim2) = ImgSize(dict['SelFileName'])
        if (xdim2, ydim2) != (xdim, ydim):
            message = "Volume and reference images have not the same size"
            result = False
    #exit(1)
    #print message
    if (_log):
        _log.debug("checkVolumeProjSize")
        if (not result):
            print message
            exit(1)
    return (result,message)


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
    _log.info("make backup script file: " + dict['ProgName'])
    log.make_backup_of_script_file(dict['ProgName'], dict['ProjectDir'] + "/" + dict['WorkingDir'])

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
    CtfGroupDirectory  = dict['CtfGroupDirectory']
    CtfGroupRootName   = dict['CtfGroupRootName']
    #Verify

    import glob, sys
    import utils_xmipp
    import launch_job
    if not os.path.exists(CtfGroupDirectory):
        os.makedirs(CtfGroupDirectory)

#    print '*********************************************************************'
#    print '* Make CTF groups'
    command = \
              ' --ctfdat ' + dict['CTFDatName' ] + \
              ' -o ' + CtfGroupDirectory + '/' + CtfGroupRootName + \
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
    fn = CtfGroupDirectory + '/'+\
                  CtfGroupRootName+\
                 'Info.xmd'
#    ctflist = glob.glob(wildcardname)
#    print ctflist, wildcardname
#    exit(1)
#    return len(ctflist)
    MD = MetaData(fn)
    return MD.size()

    

def checkOptionsCompatibility(_log, dict):
    # Never allow DoAlign2D and DoCtfCorrection together
    if (dict['DoAlign2D'] and int (dict['DoCtfCorrection'])):
        error_message = "You cannot realign classes AND perform CTF-correction. Switch either of them off!"
        _log.error(error_message)
        print error_message
        exit(1)

def initAngularReferenceFile(_log,dict):
    '''Create Initial angular file. Either fill it with zeros or copy input'''
    import shutil
    if len(dict['DocFileName'])>1:
        shutil.copy(dict['DocFileName'],dict['DocFileWithOriginalAngles'])
    else:
        MD=MetaData(dict['SelFileName'])
        MD.addLabel(MDL_ANGLEROT)
        MD.addLabel(MDL_ANGLETILT)
        MD.addLabel(MDL_ANGLEPSI)
        MD.write(dict['DocFileWithOriginalAngles'])

def dummy(_log, dict):
    print dict['dummy']
