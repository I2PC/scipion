#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     I. Foche Perez (ifoche@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'ifoche@cnb.csic.es'
# *
# **************************************************************************

#
# Brief summary: This file is intended to be the main skeleton which will guide all the installation process using a SCons python-based system.
#

import os
import sys
import platform
import SCons.Script


#############
# VARIABLES #
#############

# OS boolean vars
MACOSX = platform.system() == 'Darwin'
WINDOWS = platform.system() == 'Windows'
LINUX = platform.system() == 'Linux'
MANDATORY_PYVERSION = '2.7.7' #python version that is obligatory for executing Scipion
PYVERSION = platform.python_version() #python version present in the machine

# Big scipion structure dictionary and associated vars
# indexes
# folders & packages | libs  | index
SOFTWARE_FOLDER =      DEF =   0 # is built by default?               
CONFIG_FOLDER =        INCS =  1 # includes                           
INSTALL_FOLDER =       LIBS =  2 # libraries to create                
BIN_FOLDER =           SRC =   3 # source pattern                     
PACKAGES_FOLDER =      DIR =   4 # folder name in temporal directory  
LIB_FOLDER =           TAR =   5                                      
MAN_FOLDER =           DEPS =  6 # explicit dependencies              
TMP_FOLDER =           URL =   7 # URL to download from               
INCLUDE_FOLDER =       FLAGS = 8 # Other flags for the compiler
LOG_FOLDER =                   9 # 


# indexes for LIBS

SCIPION = {
    'FOLDERS': {SOFTWARE_FOLDER: 'software',
                CONFIG_FOLDER: os.path.join('software', 'cfg'),
                INSTALL_FOLDER: os.path.join('software', 'install'),
                BIN_FOLDER: os.path.join('software', 'bin'),
                PACKAGES_FOLDER: os.path.join('software', 'em'),
                LIB_FOLDER: os.path.join('software', 'lib'),
                MAN_FOLDER: os.path.join('software', 'man'),
                TMP_FOLDER: os.path.join('software', 'tmp'),
                INCLUDE_FOLDER: os.path.join('software', 'include'),
                LOG_FOLDER: os.path.join('software', 'log')},
    'LIBS': {},
    'PACKAGES': {'xmipp': {INSTALL_FOLDER: 'xmipp',
                           LIB_FOLDER: os.path.join('xmipp', 'lib'),
                           BIN_FOLDER: os.path.join('xmipp', 'bin'),}}}

######################
# AUXILIAR FUNCTIONS #
###################### 
def _addLibrary(env, name, dir=None, dft=True, src=None, incs=None, libs=None, tar=None, deps=None, url=None, flags=[]):
    """
    This method is for adding a library to the main dict
    """
    from os.path import splitext
    if dir is None: 
        dir = splitext(tar)[0] if tar is not None else name
    tar = '%s.tgz' % dir if tar is None else tar
    src = dir if src is None else dir
    url = 'http://scipionwiki.cnb.csic.es/files/scipion/software/external/%s' % tar if url is None else url
    incs = [] if incs is None else incs
    libs = ["lib%s.so" % name] if libs is None else libs
    deps = [] if deps is None else deps
    # First we add the library to the main dictionary
    SCIPION['LIBS'][name] = {DEF: dft,
                     SRC: src,
                     INCS: incs,
                     LIBS: libs,
                     DIR: dir,
                     TAR: tar,
                     DEPS: deps,
                     URL: url,
                     FLAGS: flags}
    # Then we add the option, so the user can call scons with this the library as option
    AddOption('--%s' % name,
              dest='%s' % name,
              action='store_true',
              help='Library %s option' % name)

def _delLibrary(env, name):
    """
    This method is for removing a library from the main dict
    """
    del SCIPION['LIBS'][name]

def _addPackage(env, name, installFolder, libFolder, binFolder):
    """
    This method is for adding a package to the main dict
    """
    SCIPION['PACKAGES'][name] = {INSTALL_FOLDER: installFolder,
                                 LIB_FOLDER: libFolder,
                                 BIN_FOLDER: binFolder}

def _delPackage(env, name):
    """
    This method is for removing a package from the main dict
    """
    del SCIPION['PACKAGES'][name]

def _downloadLibrary(env, name, verbose=False):
    """
    This method is for downloading a library and placing it in tmp dir
    It will download first the .md5
    """
    import urllib2, urlparse
    library = SCIPION['LIBS'].get(name)
    # Check whether the library must be downloaded
    if not GetOption('%s' % name):
        if not library[DEF]:
            return False
    tar = os.path.join(SCIPION['FOLDERS'][TMP_FOLDER], library[TAR])
    tarFile = File(tar)
    md5 = "%s.md5" % tar
    md5File = File(tar)
    folder = SCIPION['FOLDERS'][TMP_FOLDER]
    folderDir = Dir(folder)
    url = library[URL]
    urlMd5 = "%s.md5" % url
    go = True
    message = message2 = ''
    
    md5check = GetOption('update')
    if GetOption('update') is None:
        md5check = False
        go = False
    
    if not os.path.exists(md5):
        md5check = True
    
    if md5check:
        print "Downloading md5 file for %s library..." % name
        go, message = _downloadFile(urlMd5, md5)
        print "Checking md5 checksum..."
        if _checkMd5(tar, md5) == 1:
            go = False
    
    if go:
        while go:
            print "Downloading %s in folder %s" % (url, folder)
            down, message2 = _downloadFile(url, tar)
            if not down:
                print "\t ...Library %s not downloaded. Server says: \n %s" % (name, message2)
            if _checkMd5(tar, md5) == 0:
                if not _askContinue("Downloaded %s file doesn't match md5 its md5 checksum. Download it again?" % tar):
                    go = False
            else:
                go = False
        return tar
    elif verbose:
        if _askContinue("\t ...%s library was not downloaded. Proceed anyway?" % name):
            return True
        else:
            raise SCons.Errors.StopError("User defined stop")
    else:
        print "\t ...%s not downloaded" % (name)
        return tarFile

def _untarLibrary(env, name, tar=None, folder=None):
    """
    This method is for untar the downloaded library in the tmp directory
    """
    import tarfile
    # Add builders to deal with source code, donwloads, etc 
    libraryDict = SCIPION['LIBS'].get(name)
    # Check whether the library must be untar
    if not GetOption('%s' % name):
        if not libraryDict[DEF]:
            return False
    if tar is None:
        tar = os.path.join(SCIPION['FOLDERS'][TMP_FOLDER], libraryDict[TAR])
    if folder is None:
        folder = os.path.join(SCIPION['FOLDERS'][TMP_FOLDER],libraryDict[DIR])
    sourceTar = tarfile.open(tar,'r')
    tarContents = sourceTar.getmembers()
    tarFileContents = filter(lambda tarEntry: tarEntry.isfile(), tarContents)
    tarFileContentsNames = map(__tarInfoToNode, tarFileContents)
    for indx, item in enumerate(tarFileContentsNames): 
        tarFileContentsNames[indx] = os.path.join(SCIPION['FOLDERS'][TMP_FOLDER], item)
    sourceTar.close()
    env["UNPACK"]["EXTRACTDIR"] = SCIPION['FOLDERS'][TMP_FOLDER] 
    unpacked = env.Unpack(target=folder, 
                          source=tar, 
                          UNPACKLIST=tarFileContentsNames)
    #for dep in tarFileContentsNames:
    #    "%s depends on %s and %s" % (dep, folder, tar)
        #Depends(dep, folder)
        #Depends(dep, tar)
    #result = env.Unpack(target=tarFileContentsNames, 
    #                    source=tar, 
    #                    UNPACKLIST=tarFileContentsNames)
    return unpacked

def __tarInfoToNode(tarInfoObject):
    return tarInfoObject.name

def _compileLibrary(env, name, incs=None, libs=None, deps=None, flags=None, source=None, target=None, autoTarget=None, autoSource=None, makePath=None, makeTargets="all install", configStdOutLog=None, makeStdOutLog=None):
    """
    Function that implements pseudo-builder for executing AutoConfig and Make builders
    Args:
     * name -> name of the library as used for the key of the LIBS dictionary in SCIPION main dictionary
     * incs -> includes for compiling process. Must be a list of paths
     * libs -> list of libraries that will be built when compilation of this library is finished. Will be used as the target of the compilation process
     * deps -> list of dependencies. 2 ways may be used. You can give the main key of another library in the SCIPION main dictionary or you can use the path of the future library file
     * flags -> any other flag (in its complete way) that may be passed to the configure process, can be written as a list of flags
     * source -> 
     * target ->
     * autoTarget ->
     * autoSource ->
     * makePath ->
     * makeTargets ->
     * stdoutLog ->
     * stderrLog ->
      
    """
    from shutil import copyfile, rmtree
    
    # PREPARING ENVIRONMENT
    env['CROSS_BUILD'] = False
    libraryDict = SCIPION['LIBS'].get(name)
    # Check whether the library must be compiled
    if not GetOption('%s' % name):
        if not libraryDict[DEF]:
            return False
    tmp = SCIPION['FOLDERS'][TMP_FOLDER]
    incs = libraryDict[INCS] if incs is None else incs
    libs = libraryDict[LIBS] if libs is None else libs
    deps = libraryDict[DEPS] if deps is None else deps
    flags = libraryDict[FLAGS] if flags is None else flags
    autoSource = 'Makefile.in' if autoSource is None else autoSource
    autoTarget = 'Makefile' if autoTarget is None else autoTarget
    folder = os.path.join(tmp,
                          libraryDict[DIR])
    tar = os.path.join(tmp,
                       libraryDict[TAR])
    if source is None:
        source = Glob(os.path.join(tmp, 
                                   libraryDict[SRC]),
                      '*.c')
    target = 'lib%s.so' % name if target is None else target
    target = os.path.join(tmp,
                          libraryDict[DIR],
                          target)
    if makePath is None:
        makePath = Dir(os.path.join(tmp, 
                                    libraryDict[DIR]))
    else:
        makePath = Dir(os.path.join(tmp, 
                                    libraryDict[DIR], 
                                    makePath))
    incflags = []
    if not len(incs) == 0:
        incflagString = "CPPFLAGS=\'"
        for inc in incs:
            if os.path.exists(inc):
                incflagString += '-I%s ' % inc
        incflagString += "\'"
        incflags = [incflagString]
    flags += incflags
    
    # Output redirection to log file
    configStdOutLog = os.path.join(SCIPION['FOLDERS'][LOG_FOLDER], '%s_configure.log' % name) if configStdOutLog is None else configStdOutLog
    makeStdOutLog = os.path.join(SCIPION['FOLDERS'][LOG_FOLDER], '%s_make.log' % name) if makeStdOutLog is None else makeStdOutLog
    if os.path.exists(configStdOutLog):
        if os.path.exists("%s.old" % configStdOutLog):
            os.remove("%s.old" % configStdOutLog)
        copyfile(configStdOutLog, "%s.old" % configStdOutLog)
    if os.path.exists(makeStdOutLog):
        if os.path.exists("%s.old" % makeStdOutLog):
            os.remove("%s.old" % makeStdOutLog)
        copyfile(makeStdOutLog, "%s.old" % makeStdOutLog)
    env['AUTOCONFIGCOMSTR'] = "Configuring $TARGET from $SOURCES "
    env['MAKECOMSTR'] = "Compiling $TARGET from $SOURCES " 

    # DEPENDENCIES
    for dep in deps:
        # if we provided the key for the main dictionary
        if dep in SCIPION['LIBS']:
            librs = SCIPION['LIBS'][dep][LIBS]
            if len(librs) > 0:
                for dp in librs:
                    Depends(File(autoSource),
                            Dir(os.path.join(tmp,
                                             SCIPION['LIBS'][dep][DIR])))
                    Depends(File(target),
                            File(os.path.join(tmp,
                                              SCIPION['LIBS'][dep][DIR],
                                              dp)))
        # if we provided directly the path of a file
        else:
            Depends(File(target), File(dep))
    Depends(File(autoSource), Dir(folder))
    Depends(File(target), File(os.path.join(tmp,
                                        libraryDict[DIR],
                                        autoTarget)))

    # CONFIGURE
    configure = env.AutoConfig(source=Dir(folder),
                               target=target, 
                               AutoConfigTarget=autoTarget, 
                               AutoConfigSource=autoSource, 
                               AutoConfigParams=flags,
                               AutoConfigStdOut=configStdOutLog)

    # MAKE
    make = env.Make(source=os.path.join(tmp,
                                        libraryDict[DIR],
                                        autoTarget),
                    target=File(target), 
                    MakePath=makePath, 
                    MakeEnv=os.environ, 
                    MakeTargets=makeTargets,
                    MakeStdOut=makeStdOutLog)

    return make

def _compileWithSetupPy(env, name, deps=None, actions=['build','install'], setupPyFile=None, real_stdout=None):
    """
    This method enter in a folder where a setup.py file is placed and executes setup.py build and setup.py install with the given prefix
    Args:
     * name ->
     * source ->
     * target ->
     * actions ->
     * setupPyFile ->
    """
    from shutil import copyfile
    import subprocess
    libraryDict = SCIPION['LIBS'].get(name)
    # Check whether the module must be compiled 
    if not GetOption('%s' % name):
        if not libraryDict[DEF]:
            return False
    tmp = SCIPION['FOLDERS'][TMP_FOLDER]
    bin = SCIPION['FOLDERS'][BIN_FOLDER]

    deps = [] if deps is None else deps
    if setupPyFile is not None:
        copyFile(setupPyFile, os.path.join(tmp, libraryDict[DIR], 'setup.py'))
    setupPyFile = os.path.join(tmp, libraryDict[DIR], 'setup.py')
    
    setup = None
    command = 'cd %s && ' % os.path.join(tmp, libraryDict[DIR])
    command += os.path.join(File(os.path.join(bin, 'python')).abspath)
    command += ' '
    command += 'setup.py'
    for action in actions:
        command += ' '
        command += action
    command += ' && cd -'
    setup = env.Command(Dir(os.path.join(tmp, libraryDict[DIR], 'build')),
                    File(os.path.join(tmp, libraryDict[DIR], 'setup.py')),
                    command)
    
    for dep in deps:
        if dep in SCIPION['LIBS']:
            librs = SCIPION['LIBS'][dep][LIBS]
            if len(librs) > 0:
                for dp in librs:
                    Depends(Dir(os.path.join(tmp, libraryDict[DIR], 'build')), 
                            File(os.path.join(tmp,
                                              SCIPION['LIBS'][dep][DIR],
                                              dp)))
        else:
            Depends(Dir(os.path.join(tmp, libraryDict[DIR], 'build')),
                    File(dep))
    
    return setup

def _scipionLogo(env):
    print ""
    print "QQQQQQQQQT!^'::\"\"?$QQQQQQ" + "  S   S   S"
    print "QQQQQQQY`          ]4QQQQ"   + "  C   C   C"
    print "QQQQQD'              \"$QQ"  + "  I   I   I"
    print "QQQQP                 \"4Q"  + "  P   P   P"
    print "QQQP        :.,        -$"   + "  I   I   I"
    print "QQD       awQQQQwp      )"   + "  O   O   O"
    print "QQ'     qmQQQWQQQQg,   jm"   + "  N   N   N"
    print "Qf     QQQD^   -?$QQp jQQ"   + " ################################################"
    print "Q`    qQQ!        4WQmQQQ"   + " # Integrating image processing packages for EM #"
    print "F     QQ[          ~)WQQQ"   + " ################################################"
    print "[    ]QP             4WQQ"   + ""
    print "f    dQ(             -$QQ"   + " Installation SCons system"
    print "'    QQ              qQQQ"
    print ".   )QW            _jQQQQ"
    print "-   =QQ           jmQQQQQ"
    print "/   -QQ           QQQQQQQ"
    print "f    4Qr    jQk   )WQQQQQ"
    print "[    ]Qm    ]QW    \"QQQQQ"
    print "h     $Qc   jQQ     ]$QQQ"
    print "Q,  :aQQf qyQQQL    _yQQQ"
    print "QL jmQQQgmQQQQQQmaaaQWQQQ"
    print ""

def _downloadFile(url, file):
    """
    Function that downloads the content of a URL into a file
    Returns a boolean telling if the download succeed, and the htmllib.Message instance to the answer
    """
    import urllib
    import htmllib
    message=None
    try:
        response, message = urllib.urlretrieve(url, file)
    except:
        message = "Exception caught when downloading the URL %s. Are you connected to internet?"
        return False, message
    # If we get a html answer, then it is a server answer telling us that the file doesn't exist, but we return the message
    if message.dict['content-type'] == 'text/html':
        return False, message

    return True, message

def _installLibs(name, libs):
    """
    Function that copies the generated libs to the proper folder in the scipion architecture
    """
    print "Not implemented yet"
    return True

def _checkMd5(file, md5):
    """
    Function that checks if the md5sum from a given file match the md5sum from a md5 file
    md5 file should contain it as the first element.
    The function will return 1 if both sums match, 0 if they doesn't match and -1 if either md5 or file doesn't exist 
    """
    if not os.path.exists(file):
        print "Checksum error. File %s doesn't exist" % file
        return -1
    if not os.path.exists(md5):
        print "Checksum error. File %s doesn't exist" % md5
        return -1
    tarFileMd5 = _md5Sum(file)
    md5File = open(md5, 'r+')
    lines = md5File.readlines()
    md5FileMd5 = lines[0].replace("\n","").split(" ")[0]
    md5File.close()
    
    answer = 1 if tarFileMd5 == md5FileMd5 else 0

    if not answer:
        print "Checksum error. %s says %s, but %s says %s" % (file, tarFileMd5, md5, md5FileMd5)
    else:
        print "\t ...md5 OK"
    return answer

def _askContinue(msg="continue?"):
    """
    Function that ask the user in command line and returns the answer
    """
    answer = "-"
    while answer != "y" and answer != "n" and answer != "":
        answer = raw_input("%s (y/[n]): " % msg)
        if answer == "n" or answer == "N":
            return False
        elif answer == "y" or answer == "Y":
            return True
    return None

def _md5Sum(file):
    """
    Function that calculates the md5 sum of a given file
    """
    import hashlib
    md5sum = 0
    md5 = hashlib.md5()
    with open(file,'r+') as fileToCheck:
        for chunk in iter(lambda: fileToCheck.read(128*md5.block_size), b''):
            md5.update(chunk)
    md5sum = md5.hexdigest()
    return md5sum


#########################
# ENVIRONMENT AND TOOLS #
#########################

# We create the environment the whole build will use
env = None
env = Environment(ENV=os.environ,
                  tools=['Make',
                         'Unpack',
                         'AutoConfig',
#                         'ConfigureJNI',
#                         'install',
#                         'cuda'
                         ],
                  toolpath=[os.path.join('software', 'install', 'scons-tools')])
# To decide if a target must be rebuilt, both md5 and timestamp will be used together
env.Decider('MD5-timestamp')
# For certain files or folders which change could affect the compilation, here there could exist also a user-defined decider function. At this moment MD5-timestamp will be enough 

# To avoid the scanning of new dependencies for every file each time you rebuild, we set an implicit cache. So once a file has been scanned (its #includes pointers) then it won't be scanned again unless the file is changed
#SetOption('implicit_cache', 1)

#Depending on the system, we have to add to the environment, the path to where dynamic libraries are, so linker can find them 
if LINUX:
    env.AppendUnique(LIBPATH=os.environ.get('LD_LIBRARY_PATH'))
elif MACOSX:
    print "OS not tested yet"
    env.AppendUnique(LIBPATH=os.environ.get('DYLD_FALLBACK_LIBRARY_PATH'))
elif WINDOWS:
    print "OS not tested yet"


# Add methods to manage main dict to the environment to put them available from SConscript
env.AddMethod(_addLibrary, "AddLibrary")
env.AddMethod(_delLibrary, "DelLibrary")
env.AddMethod(_addPackage, "AddPackage")
env.AddMethod(_delPackage, "DelPackage")

# Add pseudo-builder methods in order to perfectly manage what we want, and not depend only on the builders methods
env.AddMethod(_downloadLibrary, "DownloadLibrary")
env.AddMethod(_untarLibrary, "UntarLibrary")
env.AddMethod(_compileLibrary, "CompileLibrary")
env.AddMethod(_compileWithSetupPy, "CompileWithSetupPy")

# Add other auxiliar functions to environment
env.AddMethod(_scipionLogo, "ScipionLogo")

# Add main dict to environment
env.AppendUnique(SCIPION)

########################
# Command-line options #
########################

AddOption('--update',
          dest='update',
          action='store_true',
          help='Check for packages or libraries updates')
AddOption('--purge',
          dest='purge',
          action='store_true',
          help='Completely clean the installation and its binaries')
AddOption('--binary',
          dest='binary',
          action='store_true',
          help='After doing the installation, create and package a binary for distribution')

env['MANDATORY_PYVERSION'] = MANDATORY_PYVERSION
env['PYVERSION'] = PYVERSION
Export('env', 'SCIPION')

# Only in case user didn't select help message, we run SConscript
#if not GetOption('help') and not GetOption('clean'):
if not GetOption('help'):
    env.SConscript('SConscript')