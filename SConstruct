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
CONFIG_FOLDER =        DEF =   0 # is built by default?
INSTALL_FOLDER =       INCS =  1 # includes
BIN_FOLDER =           LIBS =  2 # libraries to create
PACKAGES_FOLDER =      SRC =   3 # source pattern
LIB_FOLDER =           DIR =   4 # folder name in temporal directory 
MAN_FOLDER =           TAR =   5
TMP_FOLDER =           DEPS =  6 # explicit dependencies
INCLUDE_FOLDER =       URL =   7 # URL to download from
FLAGS = 8 # Other flags for the compiler
# indexes for LIBS

SCIPION = {
    'FOLDERS': {CONFIG_FOLDER: os.path.join('software', 'cfg'),
                INSTALL_FOLDER: os.path.join('software', 'install'),
                BIN_FOLDER: os.path.join('software', 'bin'),
                PACKAGES_FOLDER: os.path.join('software', 'em'),
                LIB_FOLDER: os.path.join('software', 'lib'),
                MAN_FOLDER: os.path.join('software', 'man'),
                TMP_FOLDER: os.path.join('software', 'tmp'),
                INCLUDE_FOLDER: os.path.join('software', 'include')},
    'LIBS': {},
    'PACKAGES': {'xmipp': {INSTALL_FOLDER: 'xmipp',
                           LIB_FOLDER: os.path.join('xmipp', 'lib'),
                           BIN_FOLDER: os.path.join('xmipp', 'bin'),}}}

######################
# AUXILIAR FUNCTIONS #
######################
# TODO - AddExternalLibrary for auto-tools made compilation and addLibrary for ourself compiled libraries and 
def addLibrary(env, name, dir, dft=True, src=[], incs=[], libs=[], tar=[], deps=[], url='', flags=[]):
    """
    This method is for adding a library to the main dict
    """
    SCIPION['LIBS'][name] = {DEF: dft,
                     SRC: src,
                     INCS: incs,
                     LIBS: libs,
                     DIR: dir,
                     TAR: tar,
                     DEPS: deps,
                     URL: url,
                     FLAGS: flags}

def delLibrary(env, name):
    """
    This method is for removing a library from the main dict
    """
    del SCIPION['LIBS'][name]

def addPackage(env, name, installFolder, libFolder, binFolder):
    """
    This method is for adding a package to the main dict
    """
    SCIPION['PACKAGES'][name] = {INSTALL_FOLDER: installFolder,
                                 LIB_FOLDER: libFolder,
                                 BIN_FOLDER: binFolder}

def delPackage(env, name):
    """
    This method is for removing a package from the main dict
    """
    del SCIPION['PACKAGES'][name]

def downloadLibrary(env, name):
    """
    This method is for downloading a library and placing it in tmp dir
    """
    library = SCIPION['LIBS'].get(name)
    tar = File(os.path.join(SCIPION['FOLDERS'][TMP_FOLDER], library[TAR]))
    folder = Dir(SCIPION['FOLDERS'][TMP_FOLDER])
    print "Downloading %s in folder %s" % (library[URL], folder)
    return env.URLDownload(folder, library[URL])

def untarLibrary(env, name, tar=None, folder=None):
    """
    This method is for untar the downloaded library in the tmp directory
    """
    # Add builders to deal with source code, donwloads, etc 
    libraryDict = SCIPION['LIBS'].get(name)
    if tar is None:
        tar = File(os.path.join(SCIPION['FOLDERS'][TMP_FOLDER], libraryDict[TAR]))
    if folder is None:
        folder = Dir(os.path.join(SCIPION['FOLDERS'][TMP_FOLDER], libraryDict[DIR]))
#    folder = Dir(os.path.join(SCIPION['FOLDERS'][TMP_FOLDER],libraryDict[DIR]))
    print "Unpacking %s in folder %s" % (tar, folder)
    return env.UnTar(folder, tar)

def compileWithSetupPy(env, name, prefix=None):
    """
    This method enter in a folder where a setup.py file is placed and executes setup.py build and setup.py install with the given prefix
    """
    print "Not implemented yet"

def compileLibrary(env, name):
    """
    This method is for compiling the library
    """
    print "Not implemented yet"

def scipionLogo(env):
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



#########################
# ENVIRONMENT AND TOOLS #
#########################

# We create the environment the whole build will use
env = None
env = Environment(ENV=os.environ,
                  tools=['URLDownload',
#                         'disttar',
                         'Make',
                         'untar',
#                         'Unpack',
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
    env.AppendUnique(LIBPATH=os.environ['LD_LIBRARY_PATH'])
elif MACOSX:
    print "OS not tested yet"
    env.AppendUnique(LIBPATH=os.environ['DYLD_FALLBACK_LIBRARY_PATH'])
elif WINDOWS:
    print "OS not tested yet"


# Add methods to manage main dict to the environment to put them available from SConscript
env.AddMethod(addLibrary, "AddLibrary")
env.AddMethod(delLibrary, "DelLibrary")
env.AddMethod(addPackage, "AddPackage")
env.AddMethod(delPackage, "DelPackage")

# Add pseudo-builder methods in order to perfectly manage what we want, and not depend only on the builders methods
env.AddMethod(downloadLibrary, "DownloadLibrary")
env.AddMethod(untarLibrary, "UntarLibrary")
env.AddMethod(compileLibrary, "CompileLibrary")
env.AddMethod(compileWithSetupPy, "CompileWithSetupPy")

# Add other auxiliar functions to environment
env.AddMethod(scipionLogo, "ScipionLogo")

# Add main dict to environment
env.AppendUnique(SCIPION)
env['MANDATORY_PYVERSION'] = MANDATORY_PYVERSION
env['PYVERSION'] = PYVERSION
Export('env')
env.SConscript('SConscript')