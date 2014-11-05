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

Import('env')

# Required for custom functions
from glob import glob
import os
from os.path import join, exists
from types import *
import os
import sys

lastTarget = env.get('lastTarget')

# Read some flags
CYGWIN = env['PLATFORM'] == 'cygwin'
MACOSX = env['PLATFORM'] == 'darwin'
MINGW = env['PLATFORM'] == 'win32'
SCONSCRIPT_PATH = Dir('.').abspath

# XMIPP ADDITIONAL EXTERNAL LIBRARIES
# Tiff
tiff = env.AddLibrary(
     'tiff',
     tar='tiff-3.9.4.tgz',
#     targets=['lib/libtiff.so'],
     clean=[Dir('#software/tmp/tiff-3.9.4').abspath],
     url='http://scipionwiki.cnb.csic.es/files/xmipp/software/external/tiff-3.9.4.tgz',)
Depends(tiff, lastTarget)

# Hdf5
hdf5 = env.AddLibrary(
     'hdf5',
     tar='hdf5-1.8.10.tgz',
     buildDir=['hdf5-1.8.10', 'hdf5-1.8.10/c++'],
     configDir=['hdf5-1.8.10', 'hdf5-1.8.10'],
     flags=[['--enable-cxx'], ['--enable-cxx']],
     autoConfigTargets=['src/Makefile', 'c++/Makefile'],
     targets=[[File('#software/lib/libhdf5.so').abspath], [File('#software/lib/libhdf5_cpp.so').abspath]],
     clean=[Dir('#software/tmp/hdf5-1.8.10').abspath],
     url='http://scipionwiki.cnb.csic.es/files/xmipp/software/external/hdf5-1.8.10.tgz')
Depends(hdf5, lastTarget)

# OpenCV
opencv = env.AddLibrary(
       'opencv',
       tar='opencv.tgz',
       clean=[Dir('#software/tmp/opencv').abspath],
       url='http://scipionwiki.cnb.csic.es/files/xmipp/software/external/opencv.tgz',
       default=False)
Depends(opencv, lastTarget)

BASIC_LIBS = ['fftw3', 'tiff', 'jpeg', 'sqlite3', 'hdf5','hdf5_cpp', 'rt']

# XMIPP SHARED LIBRARIES
# Xmipp External Libraries (bilib, condor, alglib)
xmippExternal = env.AddPackageLibrary(
              'XmippExternal',
              tars=['external/bilib.tgz', 
                    'external/condor.tgz', 
                    'external/alglib-3.8.0.cpp.tgz'],
              untarTargets=['bilib/configs.h',
                            'condor/Matrix.h',
                            'alglib-3.8.0.cpp/gpl3.txt'],
              dirs=['external',
                    'external',
                    'external'],
              patterns=['bilib/sources/*.cc', 
                        'condor/*.cpp',
                        'alglib-3.8.0.cpp/src/*.cpp'],
              incs=['external/bilib' + s for s in ['', '/headers', '/types']],
              deps=tiff + hdf5 + lastTarget)

xmippSqliteExt = env.AddPackageLibrary(
               'XmippSqliteExt',
               dirs=['external/sqliteExt'],
               patterns=['extension-functions.c'],
               libs=['m'],
               deps=tiff + hdf5 + lastTarget)

# Data
xmippData = env.AddPackageLibrary(
          'XmippData',
          dirs=['libraries/data'],
          patterns=['*.cpp'],
          incs=[Dir('.').path, Dir('libraries').path],
          libs=['XmippExternal'] + BASIC_LIBS,
          deps=xmippExternal + xmippSqliteExt)

# Classification
xmippClassif = env.AddPackageLibrary(
             'XmippClassif',
             dirs=['libraries/classification'],
             patterns=['*.cpp'],
             incs=[Dir('.').path, Dir('libraries').path],
             libs=['XmippExternal', 'XmippData'] + BASIC_LIBS,
             deps=xmippExternal + xmippSqliteExt + xmippData)

# Dimred
xmippDimred = env.AddPackageLibrary(
            'XmippDimred',
            dirs=['libraries/dimred'],
            incs=[Dir('.').path, Dir('libraries').path],
            libs=['XmippExternal', 'XmippData'] + BASIC_LIBS,
            deps=xmippExternal + xmippSqliteExt + xmippData)

# Reconstruction
xmippRecons = env.AddPackageLibrary(
            'XmippRecons',
            dirs=['libraries/reconstruction'],
            patterns=['*.cpp'],
            incs=[Dir('.').path, Dir('libraries').path, Dir('libraries/reconstruction').path, Dir('external').path],
            libs=['XmippExternal', 'XmippData', 'XmippClassif'] + BASIC_LIBS,
            deps=xmippExternal + xmippSqliteExt + xmippData + xmippClassif)

# Interface
xmippInterface = env.AddPackageLibrary(
               'XmippInterface',
               dirs=['libraries/interface'],
               patterns=['*.cpp'],
               incs=[Dir('.').path, Dir('libraries').path, Dir('#software/include/python2.7').abspath, Dir('#software/lib/python2.7/site-packages').abspath, Dir('#software/lib/python2.7/site-packages/numpy/core/include').abspath],
               libs=['XmippExternal', 'XmippData', 'pthread'],
               deps=xmippExternal + xmippSqliteExt + xmippData)

# Python binding
xmippPyBinding = env.AddPackageLibrary(
               'xmipp',
               dirs=['libraries/bindings/python'],
               patterns=['*.cpp'],
               incs=[Dir('.').path, Dir('libraries').path, Dir('#software/include/python2.7').abspath, Dir('#software/lib/python2.7/site-packages/numpy/core/include').abspath],
               libs=['XmippExternal', 'XmippData', 'XmippRecons'] + BASIC_LIBS,
               deps=xmippExternal + xmippSqliteExt + xmippData + xmippRecons,
               prefix='',
               installDir=Dir('#software/lib').abspath)

# Java binding
xmippJNI = env.AddPackageLibrary(
         'XmippJNI',
         dirs=['libraries/bindings/java'],
         patterns=['*.cpp'],
         incs=[Dir('.').path, Dir('libraries').path],
         libs=['XmippData', 'pthread', 'XmippRecons', 'XmippClassif', 'XmippExternal'] + BASIC_LIBS,
         deps=xmippData + xmippSqliteExt + xmippRecons + xmippClassif + xmippExternal)

# Parallelization
xmippParallel = env.AddPackageLibrary(
              'XmippParallel',
              dirs=['libraries/parallel'],
              patterns=['*.cpp'],
              incs=[Dir('.').path, Dir('libraries').path],
              libs=['XmippExternal', 'XmippData', 'XmippClassif', 'XmippRecons', 'pthread', (env['MPI_LIB'] or '')] + BASIC_LIBS,
              deps=xmippData + xmippSqliteExt + xmippRecons + xmippClassif + xmippExternal,
              mpi=True)

# Java libraries

######################################
# VERY UGLY. CHANGE IT WHEN POSSIBLE #
######################################
# --- From here...
javaEnumDict = javaEnumDict = {
            'ImageWriteMode': [File('libraries/data/xmipp_image_base.h').abspath, 'WRITE_'],
            'CastWriteMode': [File('libraries/data/xmipp_image_base.h').abspath, 'CW_'],
            'MDLabel': [File('libraries/data/metadata_label.h').abspath, 'MDL_'],
            'XmippError': [File('libraries/data/xmipp_error.h').abspath, 'ERR_']
            }


def WriteJavaEnum(class_name, header_file, pattern, log):
    java_file = File(join(SCONSCRIPT_PATH, 'java/src/xmipp/jni/%s.java' % class_name)).abspath
    env.Depends(java_file, header_file)
    f = open(header_file)
    fOut = open(java_file, 'w+')
    counter = 0;
    last_label_pattern = pattern + 'LAST_LABEL'
    fOut.write('package xmipp.jni; \n')
    fOut.write('public class ' + class_name + ' {\n')

    for line in f:
        l = line.strip();
        if l.startswith(pattern):
            if '///' in l:
                l, comment = l.split('///')
            else:
                comment = ''
            if l.startswith(last_label_pattern):
                l = l.replace(last_label_pattern, last_label_pattern + " = " + str(counter) + ";")
            if (l.find("=") == -1):
                l = l.replace(",", " = %d;" % counter)
                counter = counter + 1;
            else:
                l = l.replace(",", ";")

            fOut.write("   public static final int %s ///%s\n" % (l, comment))
    fOut.write("}\n")
    fOut.close()
    f.close()
    # Write log file
    if log:
        from datetime import datetime
        d = str(datetime.now())
        #d = date.today();
        log.write("Java file '%s' successful generated at %s\n" % (java_file, d))


def ExtractEnumFromHeader(source, target, env):
    log = open(str(target[0]), 'w+')
    for (class_name, list) in javaEnumDict.iteritems():
        WriteJavaEnum(class_name, list[0], list[1], log)

    log.close()
    return None


env.Append(JAVACLASSPATH='java/lib')
env['JAVABUILDPATH'] = 'java/build'
env['JAVADIR'] = 'java'
env['ENV']['LANG'] = 'en_GB.UTF-8'
env['JARFLAGS'] = '-Mcf'    # Default "cf". "M" = Do not add a manifest file.
# Set -g debug options if debugging
if env['debug'] == True:
    env['JAVAC'] = 'javac -g'

javaBuild = Execute(Mkdir('java/build'))

# Update enums from c++ headers, if they don't exist, generate it
javaLog = open(File('java/build/javaLog').abspath, 'w+')
for (class_name, class_list) in javaEnumDict.iteritems():
    java_file = "java/src/xmipp/jni/%s.java" % class_name
    WriteJavaEnum(class_name, class_list[0], class_list[1], javaLog)

javaEnums = env.Alias('javaEnums', 
                      env.Command('libraries/bindings/java/src/xmipp/jni/enums.changelog',
                      [File('libraries/data/xmipp_image_base.h').abspath, File('libraries/data/metadata_label.h').abspath],
                      ExtractEnumFromHeader))
env.Default(javaEnums)

ijLink = env.SymLink('java/lib/ij.jar', 'external/imagej/ij.jar')
env.Default(ijLink)
# --- ...to here

xmippJavaJNI = env.AddJavaLibrary(
             'XmippJNI',
             #dirs=['java/src/xmipp/jni', 'java/src/xmipp/ij/commons', 'java/src/xmipp/ij/plugins/maskstoolbar'],
             dirs=['java/src/xmipp/jni', 'java/src/xmipp/ij/commons'],
             deps=['ij'],
             )

xmippUtils = env.AddJavaLibrary(
           'XmippUtils',
           dirs=['java/src/xmipp/utils', 'java/src/xmipp/jni/'],
           deps=['ij', 'commons-cli-1.1'],
           )

xmippIJ = env.AddJavaLibrary(
        'XmippIJ',
        dirs=['java/src/xmipp/ij/commons'],
        deps=['XmippUtils'],
        )

xmippViewer = env.AddJavaLibrary(
            'XmippViewer',
            dirs=['java/src/xmipp/viewer'],
            deps=['XmippIJ', 'XmippUtils', 'ij', 'jfreechart-1.0.13'],
            )

xmippTest = env.AddJavaLibrary(
          'XmippTest',
          dirs=['java/src/xmipp/test'],
          deps=['XmippViewer', 'junit4-4.8.2', 'core-1.1']
          )

xmippIJPlugin = env.AddJavaLibrary(
          'XmippIJPlugin_MasksToolbar',
          dirs=['java/src/xmipp/ij/plugins/maskstoolbar']
          )

pluginLink = env.SymLink('external/imagej/plugins/XmippIJPlugin_MasksToolbar.jar', str(xmippIJPlugin[0]))
env.Default(pluginLink)

env2 = Environment()
env2.AppendUnique(JAVACLASSPATH=":".join(glob(join(Dir('java/lib').abspath,'*.jar'))))
javaExtraFileTypes = env2.Java('java/build/HandleExtraFileTypes.class', 'java/src/HandleExtraFileTypes.java')
env2.Depends(javaExtraFileTypes, 'java/lib/XmippViewer.jar')
env2.Default(javaExtraFileTypes)

# Java tests
AddOption('--run-java-tests', dest='run_java_tests', action='store_true',
          help='Run all Java tests (not only default ones)')

env.AddJavaTest('FilenameTest', default=False)
env.AddJavaTest('ImageGenericTest', default=False)
env.AddJavaTest('MetadataTest', default=False)

# XMIPP PROGRAMS
# 
lastTarget = xmippParallel

Default(lastTarget)
Return('lastTarget')
