from distutils.errors   import DistutilsInternalError
from types              import StringTypes
import os
from xmipp import *

def createDirectoryTree(dirName):
    """ Create directory or directory branch"""
    from distutils.dir_util import mkpath
    #
    if not isinstance(dirName, StringTypes):
        raise DistutilsInternalError, "mkpath: 'name' must be a string (got %r)" % (dirName,)
    mkpath(dirName, 0777, True)

def createRequiredDirectories(_log, ProjectDir, WorkingDir, NumberofIterations):
    """ Create all required directories, I do not want to bother about that 
        later """
    absPath = ProjectDir + "/" + WorkingDir
    for i in range(1, NumberofIterations + 1):
        createDirectoryTree(absPath + "/Iter_" + str(i).zfill(2))
    _log.info("Creating directory tree" + absPath)

def deleteWorkingDirectory(_mylog, _ProjectDir, _WorkingDir, ContinueAtIteration):

    print "delete"
    if ContinueAtIteration != 1 :
        print "You can not delete the working directory"
        print " and start at iteration: ", ContinueAtIteration
        exit(1)

    from distutils.dir_util import remove_tree
    dirName = _ProjectDir + "/" + _WorkingDir
    _mylog.info("Delete working directory tree" + dirName)
    if not isinstance(_WorkingDir, StringTypes):
        raise DistutilsInternalError, "mkpath: 'name' must be a string (got %r)" % (dirName,)
    if os.path.exists(dirName):
        remove_tree(dirName, True)

def checkVolumeProjSize(ReferenceFileNames, SelFileName):
    """ check references y projection size match"""
    #5a check volumes have same size
    (xdim, ydim, zdim, ndim) = SingleImgSize(ReferenceFileNames[0])
    for reference in ReferenceFileNames:
        (xdim2, ydim2, zdim2, ndim2) = SingleImgSize(reference)
        if (xdim2, ydim2, zdim2, ndim2) != (xdim, ydim, zdim, ndim):
            print "Reference %s and %s have not the same size" % \
                  (ReferenceFileNames[0], reference)
            exit(1)

    #5b check volume and projections  have same size
    (xdim2, ydim2, zdim2, ndim2) = ImgSize(SelFileName)
    if (xdim2, ydim2) != (xdim, ydim):
            print "Volume and reference images have not the same size"
            exit(1)

def initOuterRadius(OuterRadius, volName, selfFileName, _log):
    """ init mask radius for volumes. If suggested value is negative set to half dim"""
    if (OuterRadius < 0):
        (xdim, ydim, zdim, ndim) = ImgSize(SelFileName)
        OuterRadius = (xdim / 2) - 1
        comment = "InitOuterRadius: Outer radius set to: " + str(OuterRadius)
        print '* ' + comment
        _log.info(comment)
    return OuterRadius

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
#  , '_log'                :_log
#  , 'PaddingFactor'       : PaddingFactor
#  , 'SelFileName'         : SelFileName
#  , 'SplitDefocusDocFile' : SplitDefocusDocFile
#  , 'WienerConstant'      : WienerConstant
def execute_ctf_groups (dict):

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
            _mylog.info(message)
            sys.exit()

    launch_job.launch_job("xmipp_ctf_group",
                          command,
                          dict['_log'],
                          False, 1, 1, '')

    wildcardname = utils_xmipp.composeWildcardFileName(dict['CtfGroupDirectory'] + '/'
                                                     + dict['CtfGroupRootName']
                                                     + '_group', 'ctf')
    ctflist = glob.glob(wildcardname)
    return len(ctflist)
