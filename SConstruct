#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     I. Foche Perez (ifoche@cnb.csic.es)
# *              J. Burguet Castell (jburguet@cnb.csic.es)
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
# Main skeleton which will guide all the installation process.
#

import os
from os.path import join
import sys
import platform
import SCons.Script
import shutil


#############
# VARIABLES #
#############

# OS boolean vars
MACOSX = (platform.system() == 'Darwin')
WINDOWS = (platform.system() == 'Windows')
LINUX = (platform.system() == 'Linux')
MANDATORY_PYVERSION = '2.7.7' #python version that is obligatory for executing Scipion
PYVERSION = platform.python_version() #python version present in the machine

# Big scipion structure dictionary and associated vars
# indexes
# folders             | libs | packages           | index
SOFTWARE_FOLDER =      DEF =   PKG_DEF =             0 
                    # is built by default?               
CONFIG_FOLDER =        INCS =  PKG_INSTALL_FOLDER =  1
                    # includes                           
INSTALL_FOLDER =       LIBS =  PKG_LIB_FOLDER =      2
                    # libraries to create                
BIN_FOLDER =           SRC =   PKG_BIN_FOLDER =      3 
                    # source pattern                     
PACKAGES_FOLDER =      DIR =   PKG_URL =             4 
                    # folder name in temporal directory  
LIB_FOLDER =           TAR =   PKG_TAR =             5
                    # tarfile name in temporal directory
MAN_FOLDER =           DEPS =                        6 
                    # explicit dependencies              
TMP_FOLDER =           URL =                         7
                    # URL to download from               
INCLUDE_FOLDER =       FLAGS =                       8
                    # Other flags for the compiler
LOG_FOLDER =                                         9 # 


# indexes for LIBS

SCIPION = {
    'FOLDERS': {SOFTWARE_FOLDER: 'software',
                CONFIG_FOLDER: join('software', 'cfg'),
                INSTALL_FOLDER: join('software', 'install'),
                BIN_FOLDER: join('software', 'bin'),
                PACKAGES_FOLDER: join('software', 'em'),
                LIB_FOLDER: join('software', 'lib'),
                MAN_FOLDER: join('software', 'man'),
                TMP_FOLDER: join('software', 'tmp'),
                INCLUDE_FOLDER: join('software', 'include'),
                LOG_FOLDER: join('software', 'log')},
    'LIBS': {},
    'EMPACKAGES': {}}


#######################
# AUXILIARY FUNCTIONS #
#######################

# We create the environment the whole build will use
env = Environment(ENV=os.environ,
                  tools=['Make',
                         'Unpack',
                         'AutoConfig',
#                         'ConfigureJNI',
#                         'install',
#                         'cuda'
                         ],
                  toolpath=[join('software', 'install', 'scons-tools')])
# To decide if a target must be rebuilt, both md5 and timestamp will be used together
env.Decider('MD5-timestamp')
# For certain files or folders which change could affect the
# compilation, here there could exist also a user-defined decider
# function. At this moment MD5-timestamp will be enough

env['CROSS_BUILD'] = False
# We need this so Make doesn't try to do a cross-building. Strange enough.

env['AUTOCONFIGCOMSTR'] = "Configuring $TARGET from $SOURCES"
env['MAKECOMSTR'] = "Compiling $TARGET from $SOURCES "
# Message from autoconf and make, so we don't see all its verbosity.

def _addLibrary(env, name, dirname=None, dft=True, src=None, incs=None,
                libs=None, tar=None, deps=None, url=None, flags=[]):

    tar = tar or '%s.tgz' % name
    url = url or 'http://scipionwiki.cnb.csic.es/files/scipion/software/external/%s' % tar
    dirname = dirname or tar.rsplit('.tar.gz', 1)[0].rsplit('.tgz', 1)[0]
    libs = libs or ['lib%s.so' % name]

    x = Download(env, 'software/tmp/%s' % tar, Value(url))
    y = Untar(env, 'software/tmp/%s/configure' % dirname, x)
    z = env.AutoConfig(source=Dir('software/tmp/%s' % dirname),
                       AutoConfigTarget='Makefile',
                       AutoConfigSource='configure',
                       AutoConfigParams=['CPPFLAGS=-w',
                                         'CFLAGS=-DSQLITE_ENABLE_UPDATE_DELETE_LIMIT=1',
                                         '--prefix=%s' % os.path.abspath('software')],
                       AutoConfigStdOut='software/log/%s_config.log' % name)
    # TODO: remove fixed flags (use the ones passed as an argument
    env.Make(source=z,
             target=['software/lib/%s' % x for x in libs],
             MakePath='software/tmp/%s' % dirname,
             MakeEnv=os.environ,
             MakeTargets='install',
             MakeStdOut='software/log/%s_make.log' % name)


# def _addLibrary(env, name, dir=None, dft=True, src=None, incs=None, libs=None, tar=None, deps=None, url=None, flags=[]):
#     """
#     Add a library to the main dict.
#     """
#     from os.path import splitext
#     if dir is None: 
#         dir = splitext(tar)[0] if tar is not None else name
#     tar = '%s.tgz' % dir if tar is None else tar
#     src = dir if src is None else src
#     url = 'http://scipionwiki.cnb.csic.es/files/scipion/software/external/%s' % tar if url is None else url
#     incs = [] if incs is None else incs
#     libs = ["lib%s.so" % name] if libs is None else libs
#     deps = [] if deps is None else deps
#     # First we add the library to the main dictionary
#     SCIPION['LIBS'][name] = {DEF: dft,
#                      SRC: src,
#                      INCS: incs,
#                      LIBS: libs,
#                      DIR: dir,
#                      TAR: tar,
#                      DEPS: deps,
#                      URL: url,
#                      FLAGS: flags}
#     # Then we add the option, so the user can call scons with this the library as option
#     AddOption('--%s' % name,
#               dest='%s' % name,
#               action='store_true',
#               help='Library %s option' % name)


# def _delLibrary(env, name):
#     """
#     This method is for removing a library from the main dict
#     """
#     del SCIPION['LIBS'][name]


# def _addPackage(env, name, dft=True, dir=None, tar=None, url=None, lib=None, bin=None):
#     """
#     This method is for adding a package to the main dict
#     """
#     from os.path import splitext
#     if dir is None:
#         dir = splitext(tar)[0] if tar is not None else name
#     tar = '%s.tgz' % dir if tar is None else tar
#     url = 'http://scipionwiki.cnb.csic.es/files/scipion/software/em/%s' % tar if url is None else url
    
#     SCIPION['EMPACKAGES'][name] = {PKG_DEF: dft,
#                                  PKG_TAR: tar,
#                                  PKG_URL: url,
#                                  PKG_INSTALL_FOLDER: dir,
#                                  PKG_LIB_FOLDER: lib,
#                                  PKG_BIN_FOLDER: bin}
#     AddOption('--%s' % name,
#               dest='%s' % name,
#               action='store_true',
#               help='EM Package %s option')


# def _delPackage(env, name):
#     """
#     This method is for removing a package from the main dict
#     """
#     del SCIPION['EMPACKAGES'][name]


# def _download(target, source, env, name, type='Library'):
#     """
#     Download either a package or a library.
#     """
#     # Used to construct a builder.
#     import urllib2, urlparse
#     targetDict = {}
#     library = (type == 'Library')
#     em = type == 'EMPackage' 
#     if library:
#         targetDict = SCIPION['LIBS'].get(name)
#         DEFAULT = DEF
#         TARFILE = TAR
#         URLADD = URL
#     elif em:
#         targetDict = SCIPION['EMPACKAGES'].get(name)
#         DEFAULT = PKG_DEF
#         TARFILE = PKG_TAR
#         URLADD = PKG_URL
#     else:
#         print "Error: type %s not allowed for downloading" % type
#         return []
        
#     # Check whether the library must be downloaded
#     if not GetOption('%s' % name):
#         if not targetDict[DEFAULT]:
#             return []
#     tar = join(SCIPION['FOLDERS'][TMP_FOLDER], targetDict[TARFILE])
#     tarFile = File(tar)
#     md5 = "%s.md5" % tar
#     md5File = File(tar)
#     folder = SCIPION['FOLDERS'][TMP_FOLDER]
#     folderDir = Dir(folder)
#     url = targetDict[URLADD]
#     urlMd5 = "%s.md5" % url
#     go = True
#     message = message2 = ''
    
#     md5check = GetOption('update')
#     if GetOption('update') is None:
#         md5check = False
#         go = False
    
#     if not os.path.exists(md5):
#         md5check = True
    
#     if md5check:
#         print "Downloading md5 file for %s %s..." % (name, type)
#         go, message = _downloadFile(urlMd5, md5)
#         print "Checking md5 checksum..."
#         if _checkMd5(tar, md5) == 1:
#             go = False
    
#     if go:
#         while go:
#             print "Downloading %s in folder %s" % (url, folder)
#             down, message2 = _downloadFile(url, tar)
#             if not down:
#                 print "\t ...%s %s not downloaded. Server says: \n %s" % (type, name, message2)
#             if _checkMd5(tar, md5) == 0:
#                 if _ask("Downloaded file %s doesn't match its md5 checksum. "
#                         "Download it again?" % tar) != 'y':
#                     go = False
#             else:
#                 go = False
#         return 0  # it worked

#         raise RuntimeError("\t ...%s %s was not downloaded" % (name, type))


# def _downloadEMPackage(target, source, env, name):
#     """
#     This method is for downloading an EM package and place it in the em directory, so it can be compiled and installed
#     """
#     return _download(target, source, env, name, type='EMPackage')


# def _downloadLibrary(target, source, env, name):
#     """
#     This method is for downloading a library and placing it in tmp dir
#     It will download first the .md5
#     """
#     return _download(target, source, env, name, type='Library')


# def _untarLibrary(target, source, env, name, tar=None, folder=None):
#     """
#     Untar the downloaded library in the tmp directory
#     """
#     import tarfile
#     # Add builders to deal with source code, donwloads, etc 
#     libraryDict = SCIPION['LIBS'].get(name)
#     # Check whether the library must be untar-ed
#     if not GetOption('%s' % name):
#         if not libraryDict[DEF]:
#             return []
#     print tar
#     if tar is None:
#         tar = join(SCIPION['FOLDERS'][TMP_FOLDER], libraryDict[TAR])
#     if folder is None:
#         folder = join(SCIPION['FOLDERS'][TMP_FOLDER],libraryDict[DIR])
#     sourceTar = tarfile.open(tar,'r')
#     tarContents = sourceTar.getmembers()
#     tarFileContents = filter(lambda tarEntry: tarEntry.isfile(), tarContents)
#     tarFileContentsNames = map(__tarInfoToNode, tarFileContents)
#     for indx, item in enumerate(tarFileContentsNames): 
#         tarFileContentsNames[indx] = join(SCIPION['FOLDERS'][TMP_FOLDER], item)
#     sourceTar.close()
#     env["UNPACK"]["EXTRACTDIR"] = SCIPION['FOLDERS'][TMP_FOLDER] 
#     unpacked = env.Unpack(target=folder, 
#                           source=tar, 
#                           UNPACKLIST=tarFileContentsNames)
#     print tarFileContentsNames
#     for dep in tarFileContentsNames:
# #       "%s depends on %s and %s" % (dep, folder, tar)
#         Depends(dep, folder)
#         Depends(dep, tar)
#     result = env.Unpack(target=tarFileContentsNames, 
#                        source=tar, 
#                        UNPACKLIST=tarFileContentsNames)
#     return unpacked


# def __tarInfoToNode(tarInfoObject):
#     return tarInfoObject.name


# def _compileLibrary(env, name, libs=None, deps=None, flags=None, source=None, target=None, autoTarget=None, autoSource=None, makePath=None, makeTargets="install", configStdOutLog=None, makeStdOutLog=None):
#     """
#     Pseudo-builder for executing AutoConfig and Make builders
#     Args:
#      * name -> name of the library as used for the key of the LIBS dictionary in SCIPION main dictionary
#      * incs -> includes for compiling process. Must be a list of paths
#      * libs -> list of libraries that will be built when compilation of this library is finished. Will be used as the target of the compilation process
#      * deps -> list of dependencies. 2 ways may be used. You can give the main key of another library in the SCIPION main dictionary or you can use the path of the future library file
#      * flags -> any other flag (in its complete way) that may be passed to the configure process, can be written as a list of flags
#      * source -> 
#      * target ->
#      * autoTarget ->
#      * autoSource ->
#      * makePath ->
#      * makeTargets ->
#      * stdoutLog ->
#      * stderrLog ->
      
#     """
#     from shutil import copyfile, rmtree
    
#     # PREPARING ENVIRONMENT
#     env['CROSS_BUILD'] = False
#     libraryDict = SCIPION['LIBS'].get(name)
#     # Check whether the library must be compiled
#     if not GetOption('%s' % name):
#         if not libraryDict[DEF]:
#             return []
#     tmp = SCIPION['FOLDERS'][TMP_FOLDER]
#     incs = libraryDict[INCS] if incs is None else incs
#     libs = libraryDict[LIBS] if libs is None else libs
#     deps = libraryDict[DEPS] if deps is None else deps
#     flags = libraryDict[FLAGS] if flags is None else flags
#     autoSource = 'Makefile.in' if autoSource is None else autoSource
#     autoTarget = 'Makefile' if autoTarget is None else autoTarget
#     folder = join(tmp, libraryDict[DIR])
#     tar = join(tmp, libraryDict[TAR])
#     if source is None:
#         source = Glob(join(tmp, libraryDict[SRC]), '*.c')
#     target = 'lib%s.so' % name if target is None else target
#     target = join(tmp,
#                           libraryDict[DIR],
#                           target)
#     if makePath is None:
#         makePath = Dir(join(tmp, 
#                                     libraryDict[DIR]))
#     else:
#         makePath = Dir(join(tmp, 
#                                     libraryDict[DIR], 
#                                     makePath))
#     incflags = []
#     if not len(incs) == 0:
#         incflagString = "CPPFLAGS=\'"
#         for inc in incs:
#             if os.path.exists(inc):
#                 incflagString += '-I%s ' % inc
#         incflagString += "\'"
#         incflags = [incflagString]
#     flags += incflags
    
#     # Output redirection to log file
#     configStdOutLog = join(SCIPION['FOLDERS'][LOG_FOLDER], '%s_configure.log' % name) if configStdOutLog is None else configStdOutLog
#     makeStdOutLog = join(SCIPION['FOLDERS'][LOG_FOLDER], '%s_make.log' % name) if makeStdOutLog is None else makeStdOutLog
#     if os.path.exists(configStdOutLog):
#         if os.path.exists("%s.old" % configStdOutLog):
#             os.remove("%s.old" % configStdOutLog)
#         copyfile(configStdOutLog, "%s.old" % configStdOutLog)
#     if os.path.exists(makeStdOutLog):
#         if os.path.exists("%s.old" % makeStdOutLog):
#             os.remove("%s.old" % makeStdOutLog)
#         copyfile(makeStdOutLog, "%s.old" % makeStdOutLog)
#     env['AUTOCONFIGCOMSTR'] = "Configuring $TARGET from $SOURCES "
#     env['MAKECOMSTR'] = "Compiling $TARGET from $SOURCES " 

#     # # DEPENDENCIES
#     # for dep in deps:
#     #     # if we provided the key for the main dictionary
#     #     if dep in SCIPION['LIBS']:
#     #         librs = SCIPION['LIBS'][dep][LIBS]
#     #         if len(librs) > 0:
#     #             for dp in librs:
#     #                 Depends(autoSource,
#     #                         Dir(join(tmp,
#     #                                          SCIPION['LIBS'][dep][DIR])))
#     #                 Depends(File(target),
#     #                         File(join(tmp,
#     #                                           SCIPION['LIBS'][dep][DIR],
#     #                                           dp)))
#     #     # if we provided directly the path of a file
#     #     else:
#     #         Depends(File(target), File(dep))
#     # Depends(autoSource, Dir(folder))
#     # Depends(File(target), File(join(tmp,
#     #                                     libraryDict[DIR],
#     #                                     autoTarget)))

#     # CONFIGURE
#     configure = env.AutoConfig(source=Dir(folder),
#                                target=target, 
#                                AutoConfigTarget=autoTarget, 
#                                AutoConfigSource=str(autoSource),
#                                AutoConfigParams=flags,
#                                AutoConfigStdOut=configStdOutLog)

#     # MAKE
#     make = env.Make(source=join(tmp,
#                                         libraryDict[DIR],
#                                         autoTarget),
#                     target=File(target), 
#                     MakePath=makePath, 
#                     MakeEnv=os.environ, 
#                     MakeTargets=makeTargets,
#                     MakeStdOut=makeStdOutLog)

#     return make

# def _compileWithSetupPy(env, name, deps=None, actions=['build','install'], setupPyFile=None, real_stdout=None):
#     """
#     This method enter in a folder where a setup.py file is placed and executes setup.py build and setup.py install with the given prefix
#     Args:
#      * name ->
#      * source ->
#      * target ->
#      * actions ->
#      * setupPyFile ->
#     """
#     from shutil import copyfile
#     import subprocess
#     libraryDict = SCIPION['LIBS'].get(name)
#     # Check whether the module must be compiled 
#     if not GetOption('%s' % name):
#         if not libraryDict[DEF]:
#             return []
#     tmp = SCIPION['FOLDERS'][TMP_FOLDER]
#     bin = SCIPION['FOLDERS'][BIN_FOLDER]
#     log = SCIPION['FOLDERS'][LOG_FOLDER]

#     deps = [] if deps is None else deps
#     if setupPyFile is not None:
#         copyFile(setupPyFile, join(tmp, libraryDict[DIR], 'setup.py'))
#     setupPyFile = join(tmp, libraryDict[DIR], 'setup.py')
    
#     setup = None
#     command = 'cd %s &> /dev/null && ' % join(tmp, libraryDict[DIR])
#     command += join(File(join(bin, 'python')).abspath)
#     command += ' '
#     command += 'setup.py'
#     for action in actions:
#         command += ' '
#         command += action
#     command += ' &> %s && cd - &> /dev/null' % (join(os.environ['SCIPION_HOME'], log, '%s.log' % name))
#     setup = env.Command(Dir(join(tmp, libraryDict[DIR], 'build')),
#                         File(join(tmp, libraryDict[DIR], 'setup.py')),
#                         Action(command, 'Compiling %s > %s' % (name, join(log, '%s.log' % name))))
    
#     for dep in deps:
#         if dep in SCIPION['LIBS']:
#             librs = SCIPION['LIBS'][dep][LIBS]
#             if len(librs) > 0:
#                 for dp in librs:
#                     Depends(Dir(join(tmp, libraryDict[DIR], 'build')), 
#                             File(join(tmp,
#                                               SCIPION['LIBS'][dep][DIR],
#                                               dp)))
#         else:
#             Depends(Dir(join(tmp, libraryDict[DIR], 'build')),
#                     File(dep))
    
#     return setup


def _scipionLogo(env):
    print """
QQQQQQQQQT!^'::""?$QQQQQQ  S   S   S
QQQQQQQY`          ]4QQQQ  C   C   C
QQQQQD'              "$QQ  I   I   I
QQQQP                 "4Q  P   P   P
QQQP        :.,        -$  I   I   I
QQD       awQQQQwp      )  O   O   O
QQ'     qmQQQWQQQQg,   jm  N   N   N
Qf     QQQD^   -?$QQp jQQ ################################################
Q`    qQQ!        4WQmQQQ # Integrating image processing packages for EM #
F     QQ[          ~)WQQQ ################################################
[    ]QP             4WQQ
f    dQ(             -$QQ Installation SCons system
'    QQ              qQQQ
.   )QW            _jQQQQ
-   =QQ           jmQQQQQ
/   -QQ           QQQQQQQ
f    4Qr    jQk   )WQQQQQ
[    ]Qm    ]QW    "QQQQQ
h     $Qc   jQQ     ]$QQQ
Q,  :aQQf qyQQQL    _yQQQ
QL jmQQQgmQQQQQQmaaaQWQQQ
"""


# def _downloadFile(url, file):
#     """
#     Function that downloads the content of a URL into a file
#     Returns a boolean telling if the download succeed, and the htmllib.Message instance to the answer
#     """
#     import urllib
#     import htmllib
#     message=None
#     try:
#         response, message = urllib.urlretrieve(url, file)
#     except:
#         message = "Exception caught when downloading the URL %s. Are you connected to internet?"
#         return False, message
#     # If we get a html answer, then it is a server answer telling us that the file doesn't exist, but we return the message
#     if message.dict['content-type'] == 'text/html':
#         return False, message

#     return True, message


# def _installLibs(name, libs):
#     """
#     Function that copies the generated libs to the proper folder in the scipion architecture
#     """
#     print "Not implemented yet"
#     return True


# def _removeInstallation(env):
#     """
#     Function that cleans the folders used by a scipion installation in order to completely remove everything related to that installation
#     """
#     # Dictionary to store the folder that need to be emptied (TOCLEAN) or deleted (TOREMOVE)
#     UNINSTALL = {'TOCLEAN': [join('software','lib'), 
#                              join('software', 'lib64'),
#                              join('software', 'bin'),
#                              join('software', 'man'),
#                              join('software', 'share'),
#                              join('software', 'tmp'),
#                              join('software', 'log')],
#                  'TOREMOVE': [join('software', 'install', 'scons-2.3.1')]}
# #    if _ask("Proceeding with Scipion purge process. Everything is going to be removed from the machine. Are you sure?") != 'y':
# #        return False
#     for dir in UNINSTALL.get('TOCLEAN'):
#         print "Cleaning %s" % dir
#         list = os.listdir(dir)
#         for thing in list:
#             path = join(dir, thing)
#             if thing == '.gitignore':
#                 continue
#             if os.path.isfile(path) or os.path.islink(path):
#                 os.unlink(path)
#             else:
#                 shutil.rmtree(path)
#     for dir in UNINSTALL.get('TOREMOVE'):
#         print "Deleting %s" % dir
#         shutil.rmtree(dir)
#     return True


# def _checkMd5(fname, fmd5):
#     """
#     Check if the md5sum of a given file matches the md5sum mentioned in fmd5
#     fmd5 is the name of a file with the md5 as the first element.
#     Return 1 if both sums match, 0 if they don't and -1 if either fmd5 or fname don't exist.
#     """
#     try:
#         md5Computed = _md5Sum(fname)
#         md5Read = open(fmd5).readline().split()[0]
#     except IOError as e:
#         print 'Error computing checksum: %s' % e
#         return -1

#     if md5Computed == md5Read:
#         print "\t...md5 OK"
#         return 1
#     else:
#         print ("Checksum error. %s has md5 %s, but %s says %s" %
#                (fname, md5Computed, fmd5, md5Read))
#         return 0


# def _ask(question='Continue? (y/n): ', allowed=None):
#     """ Ask the question until it returns one of the allowed responses """

#     while True:
#         ans = raw_input(question)
#         if ans.lower() in (allowed if allowed else ['y', 'n']):
#             return ans


# def _md5Sum(fname):
#     """ Return the md5 hash of file fname """

#     import hashlib
#     mhash = hashlib.md5()
#     with open(fname) as f:
#         for chunk in iter(lambda: f.read(128 * mhash.block_size), ''):
#             mhash.update(chunk)
#     return mhash.hexdigest()


#########################
# ENVIRONMENT AND TOOLS #
#########################

# We create the environment the whole build will use
# env = Environment(ENV=os.environ,
#                   tools=['Make',
#                          'Unpack',
#                          'AutoConfig',
# #                         'ConfigureJNI',
# #                         'install',
# #                         'cuda'
#                          ],
#                   toolpath=[join('software', 'install', 'scons-tools')])
# # To decide if a target must be rebuilt, both md5 and timestamp will be used together
# env.Decider('MD5-timestamp')
# # For certain files or folders which change could affect the compilation, here there could exist also a user-defined decider function. At this moment MD5-timestamp will be enough 

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


# # Add methods to manage main dict to the environment to put them available from SConscript
env.AddMethod(_addLibrary, "AddLibrary")
# env.AddMethod(_delLibrary, "DelLibrary")
# env.AddMethod(_addPackage, "AddPackage")
# env.AddMethod(_delPackage, "DelPackage")

# # # Add pseudo-builder methods in order to perfectly manage what we want, and not depend only on the builders methods
# # env.AddMethod(_download, "Download")


# #env.Append(BUILDERS={"DownloadBuilderLibrary": Builder(action=_downloadLibrary)})
# # env.Append(BUILDERS={"DownloadBuilderEMPackage": Builder(action=_downloadEMPackage)})

Download = Builder(action='wget $SOURCE -c -O $TARGET')
# #env.Append(BUILDERS={"Download": DOWNLOAD})

# #env.AddMethod(_downloadLibrary, "DownloadLibrary")
# #env.AddMethod(_downloadEMPackage, "DownloadEMPackage")
# #env.AddMethod(_untarLibrary, "UntarLibrary")
# #env.Append(BUILDERS={"UntarLibrary": Builder(action=_untarLibrary)})

Untar = Builder(action='tar -C software/tmp --skip-old-files -xzf $SOURCE')
# #env.Append(BUILDERS={"UntarLibrary": Untar})


# env.AddMethod(_compileLibrary, "CompileLibrary")
# env.AddMethod(_compileWithSetupPy, "CompileWithSetupPy")

# Add other auxiliar functions to environment
env.AddMethod(_scipionLogo, "ScipionLogo")

# # Add auxiliar function to completely clean the installation
# env.AddMethod(_removeInstallation, "RemoveInstallation")

# Add main dict to environment
env.AppendUnique(SCIPION)

# ########################
# # Command-line options #
# ########################

# AddOption('--update',
#           dest='update',
#           action='store_true',
#           help='Check for packages or libraries updates')
# AddOption('--purge',
#           dest='purge',
#           action='store_true',
#           help='Completely clean the installation and its binaries')
# AddOption('--binary',
#           dest='binary',
#           action='store_true',
#           help='After doing the installation, create and package a binary for distribution')

env['MANDATORY_PYVERSION'] = MANDATORY_PYVERSION
env['PYVERSION'] = PYVERSION
Export('env', 'SCIPION')

# # Purge option
# if GetOption('purge'):
#     print "Purge option implies clean. Activating clean..."
#     SetOption('clean', True)

# Only in case user didn't select help message, we run SConscript
#if not GetOption('help') and not GetOption('clean'):
if not GetOption('help'):
    env.SConscript('SConscript')
