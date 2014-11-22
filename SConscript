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
import os
import sys

packageDeps = env.get('packageDeps')

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
Depends(tiff, packageDeps)

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
Depends(hdf5, packageDeps)

# OpenCV
opencv = env.AddLibrary(
       'opencv',
       tar='opencv.tgz',
       clean=[Dir('#software/tmp/opencv').abspath],
       url='http://scipionwiki.cnb.csic.es/files/xmipp/software/external/opencv.tgz',
       default=False)
Depends(opencv, packageDeps)

BASIC_LIBS = ['fftw3', 'fftw3_threads', 'tiff', 'jpeg', 'sqlite3', 'hdf5','hdf5_cpp', 'rt']

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
              deps=tiff + hdf5 + packageDeps)

xmippSqliteExt = env.AddPackageLibrary(
               'XmippSqliteExt',
               dirs=['external/sqliteExt'],
               patterns=['extension-functions.c'],
               libs=['m'],
               deps=tiff + hdf5 + packageDeps)

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
            patterns=['*.cpp'],
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

# Gtest
if int(env['gtest']):
    gtest = env.ManualInstall('gtest',
                              tar='gtest-1.6.0.tgz',
                              buildDir='gtest-1.6.0/fused-src/gtest',
                              url='http://scipionwiki.cnb.csic.es/files/xmipp/software/external/gtest-1.6.0.tgz',
                              extraActions=[(File('#/software/lib/libgtest.so'), 
                                             '%s -o %s -shared -L %s -L /usr/lib %s -fPIC %s' % (env['CC'],
                                                                                                 File('#/software/lib/libgtest.so').abspath,
                                                                                                 Dir('#software/lib').abspath,
                                                                                                 '-l '+' -l'.join(BASIC_LIBS), 
                                                                                                 File('#software/tmp/gtest-1.6.0/fused-src/gtest/gtest-all.cc').abspath))])


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
           deps=['ij', 'commons-cli-1.1', 'XmippJNI'],
           )

xmippIJ = env.AddJavaLibrary(
        'XmippIJ',
        dirs=['java/src/xmipp/ij/commons'],
        deps=['XmippUtils', 'XmippJNI'],
        )

xmippViewer = env.AddJavaLibrary(
            'XmippViewer',
            dirs=['java/src/xmipp/viewer'],
            deps=['XmippIJ', 'XmippUtils', 'ij', 'jfreechart-1.0.13', 'XmippJNI'],
            )

xmippTest = env.AddJavaLibrary(
          'XmippTest',
          dirs=['java/src/xmipp/test'],
          deps=['XmippViewer', 'junit4-4.8.2', 'core-1.1', 'XmippJNI']
          )

xmippIJPlugin = env.AddJavaLibrary(
              'XmippIJPlugin_MasksToolbar',
              dirs=['java/src/xmipp/ij/plugins/maskstoolbar'],
              deps=['XmippJNI']
              )

pluginLink = env.SymLink('external/imagej/plugins/XmippIJPlugin_MasksToolbar.jar', str(xmippIJPlugin[0]))
env.Default(pluginLink)

# For any yet unkown issue, using the environment used for the rest of the SCons is imposible to compile java code
# FIXME: Its needed to guess why and correct this. In the meanwhile we'll use an alternative environment
env2 = Environment()
env2.AppendUnique(JAVACLASSPATH=":".join(glob(join(Dir('java/lib').abspath,'*.jar'))))
javaExtraFileTypes = env2.Java('java/build/HandleExtraFileTypes.class', 'java/src/HandleExtraFileTypes.java')
env2.Depends(javaExtraFileTypes, 'java/lib/XmippViewer.jar')
env2.Default(javaExtraFileTypes)

# FIXME: For any yet unkown issue, java is being compiled putting in -d flag the class name, producing a folder with the same name as the class and putting the class file inside
fileTypesInstallation = env.Install('external/imagej/plugins/', 'java/build/HandleExtraFileTypes.class/HandleExtraFileTypes.class')
env.Default(fileTypesInstallation)

# Java tests
AddOption('--run-java-tests', dest='run_java_tests', action='store_true',
          help='Run all Java tests (not only default ones)')

env.AddJavaTest('FilenameTest', 'XmippTest.jar', default=False)
env.AddJavaTest('ImageGenericTest', 'XmippTest.jar', default=False)
env.AddJavaTest('MetadataTest', 'XmippTest.jar', default=False)

# XMIPP PROGRAMS
# 
# Idea about adding a set of programs defined in a dictionary
xmippPrograms = {
                 'xmipp_angular_commonline': ['XmippRecons'],
                 'xmipp_angular_continuous_assign':['XmippRecons'], 
                 }

_libsrecons = BASIC_LIBS+['XmippRecons', 'XmippClassif', 'XmippData', 'XmippExternal']
_depsrecons = ['lib/libXmippRecons.so', 'lib/libXmippClassif.so']
_libsinterf = BASIC_LIBS+['XmippInterface', 'python2.7', 'XmippRecons', 'XmippClassif', 'XmippData', 'XmippExternal']
_depsinterf = ['lib/libXmippInterface.so']
_libsclassif = BASIC_LIBS+['XmippClassif', 'XmippData', 'XmippExternal']
_depsclassif = ['lib/libXmippClassif.so']
_libsdimred = BASIC_LIBS+['XmippDimred', 'XmippClassif', 'XmippData', 'XmippExternal']
_depsdimred = ['lib/libXmippDimred.so', 'lib/libXmippClassif.so']
_libsnothing = BASIC_LIBS+['XmippData', 'XmippExternal']

env.AddProgram('xmipp_angular_commonline', 
               src=['applications/programs/angular_commonline/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_angular_continuous_assign',
               src=['applications/programs/angular_continuous_assign/'],
               incs=[Dir('.').path],
               libs=_libsrecons, 
               deps=_depsrecons)
env.AddProgram('xmipp_angular_discrete_assign',
               src=['applications/programs/angular_discrete_assign/'],
               incs=[Dir('.').path],
               libs=_libsrecons, 
               deps=_depsrecons)
env.AddProgram('xmipp_angular_distance',
               src=['applications/programs/angular_distance/'],
               incs=[Dir('.').path],
               libs=_libsrecons, 
               deps=_depsrecons)
env.AddProgram('xmipp_angular_distribution_show',
               src=['applications/programs/angular_distribution_show/'],
               incs=[Dir('.').path],
               libs=_libsinterf, 
               deps=_depsinterf)
env.AddProgram('xmipp_angular_neighbourhood',
               src=['applications/programs/angular_neighbourhood/'],
               incs=[Dir('.').path],
               libs=_libsrecons, 
               deps=_depsrecons)
env.AddProgram('xmipp_angular_projection_matching',
               src=['applications/programs/angular_projection_matching/'],
               incs=[Dir('.').path],
               libs=_libsrecons, 
               deps=_depsrecons)
env.AddProgram('xmipp_angular_project_library',
               src=['applications/programs/angular_project_library/'],
               incs=[Dir('.').path],
               libs=_libsrecons, 
               deps=_depsrecons)
env.AddProgram('xmipp_angular_rotate',
               src=['applications/programs/angular_rotate/'],
               incs=[Dir('.').path],
               libs=_libsnothing)
env.AddProgram('xmipp_classify_analyze_cluster',
               src=['applications/programs/classify_analyze_cluster/'],
               incs=[Dir('.').path],
               libs=_libsclassif, 
               deps=_depsclassif)
env.AddProgram('xmipp_classify_compare_classes',
               src=['applications/programs/classify_compare_classes/'],
               incs=[Dir('.').path],
               libs=_libsrecons, 
               deps=_depsrecons)
env.AddProgram('xmipp_classify_evaluate_classes',
               src=['applications/programs/classify_evaluate_classes/'],
               incs=[Dir('.').path],
               libs=_libsrecons, 
               deps=_depsrecons)
env.AddProgram('xmipp_classify_kerdensom',
               src=['applications/programs/classify_kerdensom/'],
               incs=[Dir('.').path],
               libs=_libsclassif, 
               deps=_depsclassif)
env.AddProgram('xmipp_ctf_correct_wiener3d',
               src=['applications/programs/ctf_correct_wiener3d/'],
               incs=[Dir('.').path],
               libs=_libsrecons, 
               deps=_depsrecons)
env.AddProgram('xmipp_ctf_correct_idr',
               src=['applications/programs/ctf_correct_idr/'],
               incs=[Dir('.').path],
               libs=_libsrecons, 
               deps=_depsrecons)
env.AddProgram('xmipp_ctf_create_ctfdat',
               src=['applications/programs/ctf_create_ctfdat/'],
               incs=[Dir('.').path],
               libs=_libsrecons, 
               deps=_depsrecons)
env.AddProgram('xmipp_ctf_enhance_psd',
               src=['applications/programs/ctf_enhance_psd/'],
               incs=[Dir('.').path],
               libs=_libsrecons, 
               deps=_depsrecons)
env.AddProgram('xmipp_ctf_estimate_from_micrograph',
               src=['applications/programs/ctf_estimate_from_micrograph/'],
               incs=[Dir('.').path],
               libs=_libsrecons, 
               deps=_depsrecons)
env.AddProgram('xmipp_ctf_estimate_from_psd',
               src=['applications/programs/ctf_estimate_from_psd/'],
               incs=[Dir('.').path],
               libs=_libsrecons, 
               deps=_depsrecons)
env.AddProgram('xmipp_ctf_group',
               src=['applications/programs/ctf_group/'],
               incs=[Dir('.').path],
               libs=_libsrecons, 
               deps=_depsrecons)
env.AddProgram('xmipp_ctf_phase_flip',
               src=['applications/programs/ctf_phase_flip/'],
               incs=[Dir('.').path],
               libs=_libsrecons, 
               deps=_depsrecons)
env.AddProgram('xmipp_ctf_show',
               src=['applications/programs/ctf_show/'],
               incs=[Dir('.').path],
               libs=_libsrecons, 
               deps=_depsrecons)
env.AddProgram('xmipp_ctf_sort_psds',
               src=['applications/programs/ctf_sort_psds/'],
               incs=[Dir('.').path],
               libs=_libsrecons, 
               deps=_depsrecons)
if not int(env['release']):
    env.AddProgram('xmipp_idr_xray_tomo',
                   src=['applications/programs/idr_xray_tomo/'],
                   incs=[Dir('.').path],
                   libs=_libsrecons, 
                   deps=_depsrecons)
env.AddProgram('xmipp_image_align',
               src=['applications/programs/image_align/'],
               incs=[Dir('.').path],
               libs=_libsrecons, 
               deps=_depsrecons)
env.AddProgram('xmipp_image_align_tilt_pairs',
               src=['applications/programs/image_align_tilt_pairs/'],
               incs=[Dir('.').path],
               libs=_libsrecons, 
               deps=_depsrecons)
env.AddProgram('xmipp_image_common_lines',
               src=['applications/programs/image_common_lines/'],
               incs=[Dir('.').path],
               libs=_libsrecons, 
               deps=_depsrecons)
env.AddProgram('xmipp_image_convert',
               src=['applications/programs/image_convert/'],
               incs=[Dir('.').path],
               libs=_libsclassif, 
               deps=_depsclassif)
env.AddProgram('xmipp_find_center',
               src=['applications/programs/image_find_center/'],
               incs=[Dir('.').path],
               libs=_libsclassif, 
               deps=_depsclassif)
env.AddProgram('xmipp_image_header',
               src=['applications/programs/image_header/'],
               incs=[Dir('.').path],
               libs=_libsclassif, 
               deps=_depsclassif)
env.AddProgram('xmipp_image_histogram',
               src=['applications/programs/image_histogram/'],
               incs=[Dir('.').path],
               libs=_libsclassif, 
               deps=_depsclassif)
env.AddProgram('xmipp_image_operate',
               src=['applications/programs/image_operate/'],
               incs=[Dir('.').path],
               libs=_libsclassif, 
               deps=_depsclassif)
env.AddProgram('xmipp_image_rotational_pca',
               src=['applications/programs/image_rotational_pca/'],
               incs=[Dir('.').path],
               libs=_libsrecons, 
               deps=_depsrecons)
env.AddProgram('xmipp_image_residuals',
               src=['applications/programs/image_residuals/'],
               incs=[Dir('.').path],
               libs=_libsrecons, 
               deps=_depsrecons)
env.AddProgram('xmipp_image_resize',
               src=['applications/programs/image_resize/'],
               incs=[Dir('.').path],
               libs=_libsrecons, 
               deps=_depsrecons)
env.AddProgram('xmipp_image_rotational_spectra',
               src=['applications/programs/image_rotational_spectra/'],
               incs=[Dir('.').path],
               libs=_libsrecons, 
               deps=_depsrecons)
env.AddProgram('xmipp_image_sort_by_statistics',
               src=['applications/programs/image_sort_by_statistics/'],
               incs=[Dir('.').path],
               libs=_libsrecons, 
               deps=_depsrecons)
env.AddProgram('xmipp_image_separate_objects',
               src=['applications/programs/image_separate_objects/'],
               incs=[Dir('.').path],
               libs=_libsclassif, 
               deps=_depsclassif)
env.AddProgram('xmipp_image_statistics',
               src=['applications/programs/image_statistics/'],
               incs=[Dir('.').path],
               libs=_libsclassif, 
               deps=_depsclassif)
env.AddProgram('xmipp_image_vectorize',
               src=['applications/programs/image_vectorize/'],
               incs=[Dir('.').path],
               libs=_libsclassif, 
               deps=_depsclassif)
env.AddProgram('xmipp_matrix_dimred',
               src=['applications/programs/matrix_dimred/'],
               incs=[Dir('.').path],
               libs=_libsdimred, 
               deps=_depsdimred)
env.AddProgram('xmipp_metadata_convert_to_spider',
               src=['applications/programs/metadata_convert_to_spider/'],
               incs=[Dir('.').path],
               libs=_libsinterf, 
               deps=_depsinterf)
env.AddProgram('xmipp_metadata_histogram',
               src=['applications/programs/metadata_histogram/'],
               incs=[Dir('.').path],
               libs=_libsnothing)
env.AddProgram('xmipp_metadata_import',
               src=['applications/programs/metadata_import/'],
               incs=[Dir('.').path],
               libs=_libsnothing)
env.AddProgram('xmipp_metadata_split',
               src=['applications/programs/metadata_split/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_metadata_split_3D',
               src=['applications/programs/metadata_split_3D/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_metadata_utilities',
               src=['applications/programs/metadata_utilities/'],
               incs=[Dir('.').path],
               libs=_libsnothing)
env.AddProgram('xmipp_metadata_xml',
               src=['applications/programs/metadata_xml/'],
               incs=[Dir('.').path],
               libs=_libsnothing)
env.AddProgram('xmipp_micrograph_scissor',
               src=['applications/programs/micrograph_scissor/'],
               incs=[Dir('.').path],
               libs=_libsnothing)
env.AddProgram('xmipp_micrograph_automatic_picking',
               src=['applications/programs/micrograph_automatic_picking/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_ml_align2d',
               src=['applications/programs/ml_align2d/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_mlf_align2d',
               src=['applications/programs/mlf_align2d/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_ml_refine3d',
               src=['applications/programs/ml_refine3d/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_mlf_refine3d',
               src=['applications/programs/mlf_refine3d/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_ml_tomo',
               src=['applications/programs/ml_tomo/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_mrc_create_metadata',
               src=['applications/programs/mrc_create_metadata/'],
               incs=[Dir('.').path],
               libs=_libsnothing)
env.AddProgram('xmipp_nma_alignment',
               src=['applications/programs/nma_alignment/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_nma_alignment_vol',
               src=['applications/programs/nma_alignment_vol/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_flexible_alignment',
               src=['applications/programs/flexible_alignment/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_pdb_analysis',
               src=['applications/programs/pdb_analysis/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_pdb_construct_dictionary',
               src=['applications/programs/pdb_construct_dictionary/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_pdb_nma_deform',
               src=['applications/programs/pdb_nma_deform/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_pdb_restore_with_dictionary',
               src=['applications/programs/pdb_restore_with_dictionary/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_phantom_create',
               src=['applications/programs/phantom_create/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_phantom_project',
               src=['applications/programs/phantom_project/'],
               incs=[Dir('.').path],
               libs=_libsinterf,
               deps=['lib/libXmippInterface.so', 'lib/libXmippRecons.so', 'lib/libXmippClassif.so'])
env.AddProgram('xmipp_phantom_simulate_microscope',
               src=['applications/programs/phantom_simulate_microscope/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_phantom_transform',
               src=['applications/programs/phantom_transform/'],
               incs=[Dir('.').path],
               libs=_libsinterf,
               deps=['lib/libXmippInterface.so', 'lib/libXmippRecons.so', 'lib/libXmippClassif.so'])
env.AddProgram('xmipp_reconstruct_art',
               src=['applications/programs/reconstruct_art/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_reconstruct_art_pseudo',
               src=['applications/programs/reconstruct_art_pseudo/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
if not int(env['release']):
    env.AddProgram('xmipp_reconstruct_art_xray',
                   src=['applications/programs/reconstruct_art_xray/'],
                   incs=[Dir('.').path],
                   libs=_libsrecons,
                   deps=_depsrecons)
env.AddProgram('xmipp_reconstruct_fourier',
               src=['applications/programs/reconstruct_fourier/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_reconstruct_significant',
               src=['applications/programs/reconstruct_significant/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_reconstruct_wbp',
               src=['applications/programs/reconstruct_wbp/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
#env.AddProgram('xmipp_reconstruct_fsc',
#               src=['applications/programs/reconstruct_fsc/'],
#               incs=[Dir('.').path],
#               libs=_libsnothing)
if not int(env['release']):
    env.AddProgram('xmipp_resolution_ibw',
                   src=['applications/programs/resolution_ibw/'],
                   incs=[Dir('.').path],
                   libs=_libsrecons,
                   deps=_depsrecons)
env.AddProgram('xmipp_resolution_ssnr',
               src=['applications/programs/resolution_ssnr/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_score_micrograph',
               src=['applications/programs/score_micrograph/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_transform_add_noise',
               src=['applications/programs/transform_add_noise/'],
               incs=[Dir('.').path],
               libs=_libsnothing)
env.AddProgram('xmipp_transform_adjust_volume_grey_levels',
               src=['applications/programs/transform_adjust_volume_grey_levels/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_transform_center_image',
               src=['applications/programs/transform_center_image/'],
               incs=[Dir('.').path],
               libs=_libsnothing)
env.AddProgram('xmipp_transform_dimred',
               src=['applications/programs/transform_dimred/'],
               incs=[Dir('.').path],
               libs=BASIC_LIBS+['XmippDimred', 'XmippClassif', 'XmippData', 'XmippExternal'], 
               deps=_depsdimred)
env.AddProgram('xmipp_transform_downsample',
               src=['applications/programs/transform_downsample/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_transform_filter',
               src=['applications/programs/transform_filter/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_transform_geometry',
               src=['applications/programs/transform_geometry/'],
               incs=[Dir('.').path],
               libs=_libsnothing)
env.AddProgram('xmipp_transform_mask',
               src=['applications/programs/transform_mask/'],
               incs=[Dir('.').path],
               libs=_libsnothing)
env.AddProgram('xmipp_transform_mirror',
               src=['applications/programs/transform_mirror/'],
               incs=[Dir('.').path],
               libs=_libsnothing)
env.AddProgram('xmipp_transform_morphology',
               src=['applications/programs/transform_morphology/'],
               incs=[Dir('.').path],
               libs=_libsnothing)
env.AddProgram('xmipp_transform_normalize',
               src=['applications/programs/transform_normalize/'],
               incs=[Dir('.').path],
               libs=_libsnothing)
env.AddProgram('xmipp_transform_randomize_phases',
               src=['applications/programs/transform_randomize_phases/'],
               incs=[Dir('.').path],
               libs=_libsnothing)
env.AddProgram('xmipp_transform_range_adjust',
               src=['applications/programs/transform_range_adjust/'],
               incs=[Dir('.').path],
               libs=_libsnothing)
env.AddProgram('xmipp_transform_symmetrize',
               src=['applications/programs/transform_symmetrize/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_transform_threshold',
               src=['applications/programs/transform_threshold/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_transform_window',
               src=['applications/programs/transform_window/'],
               incs=[Dir('.').path],
               libs=_libsnothing)
if not int(env['release']):
    env.AddProgram('xmipp_tomo_align_dual_tilt_series',
                   src=['applications/programs/tomo_align_dual_tilt_series/'],
                   incs=[Dir('.').path],
                   libs=_libsrecons,
                   deps=_depsrecons)
    env.AddProgram('xmipp_tomo_align_refinement',
                   src=['applications/programs/tomo_align_refinement/'],
                   incs=[Dir('.').path],
                   libs=_libsrecons,
                   deps=_depsrecons)
env.AddProgram('xmipp_tomo_align_tilt_series',
               src=['applications/programs/tomo_align_tilt_series/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons,
               cuda=int(env['cuda']))
env.AddProgram('xmipp_tomo_detect_missing_wedge',
               src=['applications/programs/tomo_detect_missing_wedge/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_tomo_project',
               src=['applications/programs/tomo_project/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_tomo_remove_fluctuations',
               src=['applications/programs/tomo_remove_fluctuations/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_tomo_extract_subvolume',
               src=['applications/programs/tomo_extract_subvolume/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_volume_align_prog',
               src=['applications/programs/volume_align_prog/'],
               incs=[Dir('.').path, Dir('#software/lib/python2.7/site-packages/numpy/core/include/').abspath],
               libs=_libsinterf,
               deps=_depsinterf)
env.AddProgram('xmipp_volume_center',
               src=['applications/programs/volume_center/'],
               incs=[Dir('.').path],
               libs=_libsnothing)
env.AddProgram('xmipp_volume_correct_bfactor',
               src=['applications/programs/volume_correct_bfactor/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_volume_enhance_contrast',
               src=['applications/programs/volume_enhance_contrast/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_volume_find_symmetry',
               src=['applications/programs/volume_find_symmetry/'],
               incs=[Dir('.').path],
               libs=_libsnothing)
env.AddProgram('xmipp_volume_from_pdb',
               src=['applications/programs/volume_from_pdb/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_volume_initial_simulated_annealing',
               src=['applications/programs/volume_initial_simulated_annealing/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_volume_validate_pca',
               src=['applications/programs/volume_validate_pca/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_volume_reslice',
               src=['applications/programs/volume_reslice/'],
               incs=[Dir('.').path],
               libs=_libsnothing)
env.AddProgram('xmipp_volume_segment',
               src=['applications/programs/volume_segment/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_volume_structure_factor',
               src=['applications/programs/volume_structure_factor/'],
               incs=[Dir('.').path],
               libs=_libsnothing)
env.AddProgram('xmipp_volume_to_pseudoatoms',
               src=['applications/programs/volume_to_pseudoatoms/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_volume_to_web',
               src=['applications/programs/volume_to_web/'],
               incs=[Dir('.').path],
               libs=_libsnothing)
env.AddProgram('xmipp_xray_import',
               src=['applications/programs/xray_import/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
env.AddProgram('xmipp_xray_psf_create',
               src=['applications/programs/xray_psf_create/'],
               incs=[Dir('.').path],
               libs=_libsnothing)
env.AddProgram('xmipp_xray_project',
               src=['applications/programs/xray_project/'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons)
if not int(env['release']):
    env.AddProgram('xmipp_xray_volume_correct',
                   src=['applications/programs/xray_volume_correct/'],
                   incs=[Dir('.').path],
                   libs=_libsrecons,
                   deps=_depsrecons)
# MPI program
_libsrecons += ['XmippParallel']
_depsrecons += ['lib/libXmippParallel.so']

env.AddProgram('xmipp_mpi_angular_class_average', 
               src=['applications/programs/mpi_angular_class_average'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons,
               mpi=True)
env.AddProgram('xmipp_mpi_angular_continuous_assign', 
               src=['applications/programs/mpi_angular_continuous_assign'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons,
               mpi=True)
env.AddProgram('xmipp_mpi_angular_discrete_assign', 
               src=['applications/programs/mpi_angular_discrete_assign'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons,
               mpi=True)
env.AddProgram('xmipp_mpi_angular_projection_matching', 
               src=['applications/programs/mpi_angular_projection_matching'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons,
               mpi=True)
env.AddProgram('xmipp_mpi_angular_project_library', 
               src=['applications/programs/mpi_angular_project_library'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons,
               mpi=True)
xmippMPIClassifyCL2D = env.AddProgram('xmipp_mpi_classify_CL2D', 
                                      src=['applications/programs/mpi_classify_CL2D'],
                                      incs=[Dir('.').path],
                                      libs=_libsrecons,
                                      deps=_depsrecons,
                                      mpi=True)
xmippClassifyCL2D = env.SymLink('bin/xmipp_classify_CL2D', 'bin/xmipp_mpi_classify_CL2D')
Depends(xmippClassifyCL2D, xmippMPIClassifyCL2D)
xmippMPIClassifyCLTomoProg = env.AddProgram('xmipp_mpi_classify_CLTomo_prog', 
                                            src=['applications/programs/mpi_classify_CLTomo_prog'],
                                            incs=[Dir('.').path, Dir('#software/lib/python2.7/site-packages/numpy/core/include/').abspath],
                                            libs=_libsrecons+_libsinterf,
                                            deps=_depsrecons+_depsinterf,
                                            mpi=True)
xmippClassifyCLTomo = env.SymLink('bin/xmipp_classify_CLTomo', 'bin/xmipp_mpi_classify_CLTomo')
Depends(xmippClassifyCLTomo, xmippMPIClassifyCLTomoProg)
env.AddProgram('xmipp_mpi_classify_CL2D_core_analysis', 
               src=['applications/programs/mpi_classify_CL2D_core_analysis'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons,
               mpi=True)
env.AddProgram('xmipp_mpi_ctf_correct_idr', 
               src=['applications/programs/mpi_ctf_correct_idr'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons,
               mpi=True)
env.AddProgram('xmipp_mpi_ctf_sort_psds', 
               src=['applications/programs/mpi_ctf_sort_psds'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons,
               mpi=True)
env.AddProgram('xmipp_mpi_image_operate', 
               src=['applications/programs/mpi_image_operate'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons,
               mpi=True)
env.AddProgram('xmipp_mpi_image_rotational_pca', 
               src=['applications/programs/mpi_image_rotational_pca'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons,
               mpi=True)
env.AddProgram('xmipp_mpi_performance_test', 
               src=['applications/programs/mpi_performance_test'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons,
               mpi=True)
env.AddProgram('xmipp_mpi_image_resize', 
               src=['applications/programs/mpi_image_resize'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons,
               mpi=True)
xmippMPIImageSort = env.AddProgram('xmipp_mpi_image_sort', 
                                   src=['applications/programs/mpi_image_sort'],
                                   incs=[Dir('.').path],
                                   libs=_libsrecons,
                                   deps=_depsrecons,
                                   mpi=True)
xmippImageSort = env.SymLink('bin/xmipp_image_sort', 'bin/xmipp_mpi_image_sort')
Depends(xmippImageSort, xmippMPIImageSort)
env.AddProgram('xmipp_mpi_ml_align2d', 
               src=['applications/programs/mpi_ml_align2d'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons,
               mpi=True)
env.AddProgram('xmipp_mpi_ml_tomo', 
               src=['applications/programs/mpi_ml_tomo'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons,
               mpi=True)
env.AddProgram('xmipp_mpi_mlf_align2d', 
               src=['applications/programs/mpi_mlf_align2d'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons,
               mpi=True)
env.AddProgram('xmipp_mpi_ml_refine3d', 
               src=['applications/programs/mpi_ml_refine3d'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons,
               mpi=True)
env.AddProgram('xmipp_mpi_mlf_refine3d', 
               src=['applications/programs/mpi_mlf_refine3d'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons,
               mpi=True)
env.AddProgram('xmipp_mpi_nma_alignment', 
               src=['applications/programs/mpi_nma_alignment'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons,
               mpi=True)
env.AddProgram('xmipp_mpi_xray_project', 
               src=['applications/programs/mpi_xray_project'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons,
               mpi=True)
env.AddProgram('xmipp_mpi_reconstruct_art', 
               src=['applications/programs/mpi_reconstruct_art'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons,
               mpi=True)
env.AddProgram('xmipp_mpi_reconstruct_fourier', 
               src=['applications/programs/mpi_reconstruct_fourier'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons,
               mpi=True)
env.AddProgram('xmipp_mpi_reconstruct_wbp', 
               src=['applications/programs/mpi_reconstruct_wbp'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons,
               mpi=True)
env.AddProgram('xmipp_mpi_reconstruct_significant', 
               src=['applications/programs/mpi_reconstruct_significant'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons,
               mpi=True)
env.AddProgram('xmipp_mpi_run', 
               src=['applications/programs/mpi_run'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons,
               mpi=True)
env.AddProgram('xmipp_mpi_tomo_extract_subvolume', 
               src=['applications/programs/mpi_tomo_extract_subvolume'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons,
               mpi=True)
env.AddProgram('xmipp_mpi_transform_filter', 
               src=['applications/programs/mpi_transform_filter'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons,
               mpi=True)
env.AddProgram('xmipp_mpi_transform_symmetrize', 
               src=['applications/programs/mpi_transform_symmetrize'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons,
               mpi=True)
env.AddProgram('xmipp_mpi_transform_geometry', 
               src=['applications/programs/mpi_transform_geometry'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons,
               mpi=True)
env.AddProgram('xmipp_mpi_transform_mask', 
               src=['applications/programs/mpi_transform_mask'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons,
               mpi=True)
env.AddProgram('xmipp_mpi_transform_normalize', 
               src=['applications/programs/mpi_transform_normalize'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons,
               mpi=True)
env.AddProgram('xmipp_mpi_transform_threshold', 
               src=['applications/programs/mpi_transform_threshold'],
               incs=[Dir('.').path],
               libs=_libsrecons,
               deps=_depsrecons,
               mpi=True)
if not int(env['release']):
    env.AddProgram('xmipp_mpi_write_test', 
                   src=['applications/programs/mpi_write_test'],
                   incs=[Dir('.').path],
                   libs=_libsrecons,
                   deps=_depsrecons,
                   mpi=True)

# Batches (apps)
xmippApropos = env.SymLink('bin/xmipp_apropos', 'applications/scripts/apropos/batch_apropos.py')
Depends(xmippApropos, packageDeps)
xmippCompile = env.SymLink('bin/xmipp_compile', 'applications/scripts/compile/batch_compile.py')
Depends(xmippCompile, packageDeps)
xmippExportEMX = env.SymLink('bin/xmipp_export_emx', 'applications/scripts/export_emx/batch_export_emx.py')
Depends(xmippExportEMX, packageDeps)
xmippImportBox = env.SymLink('bin/xmipp_import_box', 'applications/scripts/import_box/batch_import_box.py')
Depends(xmippImportBox, packageDeps)
xmippImportCTFparam = env.SymLink('bin/xmipp_import_ctfparam', 'applications/scripts/import_ctfparam/batch_import_ctfparam.py')
Depends(xmippImportCTFparam, packageDeps)
xmippImportCTFdat = env.SymLink('bin/xmipp_import_ctfdat', 'applications/scripts/import_ctfdat/batch_import_ctfdat.py')
Depends(xmippImportCTFdat, packageDeps)
xmippImportEMX = env.SymLink('bin/xmipp_import_emx', 'applications/scripts/import_emx/batch_import_emx.py')
Depends(xmippImportEMX, packageDeps)
xmippMetadataPlot = env.SymLink('bin/xmipp_metadata_plot', 'applications/scripts/metadata_plot/batch_metadata_plot.py')
Depends(xmippMetadataPlot, packageDeps)
xmippMetadataSelfileCreate = env.SymLink('bin/xmipp_metadata_selfile_create', 'applications/scripts/metadata_selfile_create/batch_metadata_selfile_create.py')
Depends(xmippMetadataSelfileCreate, packageDeps)
xmippProtocols = env.SymLink('bin/xmipp_protocols', 'protocols/batch_protocols.py')
Depends(xmippProtocols, packageDeps)
xmippBrowser = env.SymLink('bin/xmipp_browser', 'applications/scripts/browser/batch_browser.py')
Depends(xmippBrowser, packageDeps)
xmippMicrographParticlePicking = env.SymLink('bin/xmipp_micrograph_particle_picking', 'applications/scripts/micrograph_particle_picking/batch_micrograph_particle_picking.py')
Depends(xmippMicrographParticlePicking, packageDeps)
xmippChimeraClient = env.SymLink('bin/xmipp_chimera_client', 'applications/scripts/chimera_client/batch_chimera_client.py')
Depends(xmippChimeraClient, packageDeps)
xmippMicrographTiltpairPicking = env.SymLink('bin/xmipp_micrograph_tiltpair_picking', 'applications/scripts/micrograph_tiltpair_picking/batch_micrograph_tiltpair_picking.py')
Depends(xmippMicrographTiltpairPicking, packageDeps)
xmippMPIClassifyCLTomo = env.SymLink('bin/xmipp_mpi_classify_CLTomo', 'applications/scripts/mpi_classify_CLTomo/batch_mpi_classify_CLTomo.sh')
Depends(xmippMPIClassifyCLTomo, packageDeps)
xmippMPIStepsRunner = env.SymLink('bin/xmipp_mpi_steps_runner', 'protocols/batch_protocols.py')
Depends(xmippMPIStepsRunner, packageDeps)
xmippProjectionsExplorer = env.SymLink('bin/xmipp_projections_explorerj', 'applications/scripts/projections_explorerj/batch_projections_explorerj.py')
Depends(xmippProjectionsExplorer, packageDeps)
xmippShowj = env.SymLink('bin/xmipp_showj', 'applications/scripts/showj/batch_showj.py')
Depends(xmippShowj, packageDeps)
xmippTomoj = env.SymLink('bin/xmipp_tomoj', 'applications/scripts/tomoj/batch_tomoj.py')
Depends(xmippTomoj, packageDeps)
xmippVisualizePreprocessingMicrographj = env.SymLink('bin/xmipp_visualize_preprocessing_micrographj', 'applications/scripts/visualize_preprocessing_micrograph/batch_visualize_preprocessing_micrographj.py')
Depends(xmippVisualizePreprocessingMicrographj, packageDeps)
xmippVolumeAlign = env.SymLink('bin/xmipp_volume_align', 'applications/scripts/volume_align/batch_volume_align.sh')
Depends(xmippVolumeAlign, packageDeps)
xmippImagej = env.SymLink('bin/xmipp_imagej', 'external/runImageJ')
Depends(xmippImagej, packageDeps)
xmippPython = env.SymLink('bin/xmipp_python', File('#software/bin/python').abspath)
Depends(xmippPython, packageDeps)

def AddXmippTest(name, testprog, command):
    #testprog = AddXmippProgram(name, ['gtest'], 'tests')
    testname = 'xmipp_' + name
    xmlFileName = join('applications','tests','OUTPUT', testname) + ".xml"
    if  os.path.exists(xmlFileName):
       os.remove(xmlFileName)
    env['ENV']['LD_LIBRARY_PATH'] = ':'.join([Dir('lib').abspath, Dir('#software/lib').abspath])
    env['ENV']['XMIPP_HOME'] = Dir('.').abspath
    testcase = env.Alias('run_' + name , env.Command(xmlFileName, 'bin/xmipp_%s' % name, command))
    env.Depends(testcase, testprog)
    test = env.Alias('run_tests', testcase)
    AlwaysBuild(testcase)
    return testcase

_libstests = BASIC_LIBS+['XmippRecons', 'XmippClassif', 'XmippData', 'XmippExternal', 'XmippDimred', 'gtest']
_depstests = ['lib/libXmippRecons.so', 'lib/libXmippClassif.so', 'lib/libXmippDimred.so']
if int(env['gtest']):
    testCTF = env.AddProgram('xmipp_test_ctf', 
                             src=['applications/tests/test_ctf'],
                             incs=[Dir('.').path],
                             libs=_libstests,
                             deps=_depstests)
    AddXmippTest('test_ctf', testCTF, "$SOURCE --gtest_output=xml:$TARGET")
    testEuler = env.AddProgram('xmipp_test_euler', 
                               src=['applications/tests/test_euler'],
                               incs=[Dir('.').path],
                               libs=_libstests,
                               deps=_depstests)
    AddXmippTest('test_euler', testEuler, "$SOURCE --gtest_output=xml:$TARGET")
    testFftw = env.AddProgram('xmipp_test_fftw', 
                               src=['applications/tests/test_fftw'],
                               incs=[Dir('.').path],
                               libs=_libstests,
                               deps=_depstests)
    AddXmippTest('test_fftw', testFftw, "$SOURCE --gtest_output=xml:$TARGET")
    testFilters = env.AddProgram('xmipp_test_filters', 
                                 src=['applications/tests/test_filters'],
                                 incs=[Dir('.').path],
                                 libs=_libstests,
                                 deps=_depstests)
    AddXmippTest('test_filters', testFilters, "$SOURCE --gtest_output=xml:$TARGET")
    testFringeProcessing = env.AddProgram('xmipp_test_fringe_processing', 
                                          src=['applications/tests/test_fringe_processing'],
                                          incs=[Dir('.').path],
                                          libs=_libstests,
                                          deps=_depstests)
    AddXmippTest('test_fringe_processing', testFringeProcessing, "$SOURCE --gtest_output=xml:$TARGET")
    testFuncs = env.AddProgram('xmipp_test_funcs', 
                               src=['applications/tests/test_funcs'],
                               incs=[Dir('.').path],
                               libs=_libstests,
                               deps=_depstests)
    AddXmippTest('test_funcs', testFuncs, "$SOURCE --gtest_output=xml:$TARGET")
    testGeometry = env.AddProgram('xmipp_test_geometry', 
                                  src=['applications/tests/test_geometry'],
                                  incs=[Dir('.').path],
                                  libs=_libstests,
                                  deps=_depstests)
    AddXmippTest('test_geometry', testGeometry, "$SOURCE --gtest_output=xml:$TARGET")
    testImage = env.AddProgram('xmipp_test_image', 
                               src=['applications/tests/test_image'],
                               incs=[Dir('.').path],
                               libs=_libstests,
                               deps=_depstests)
    AddXmippTest('test_image', testImage, "$SOURCE --gtest_output=xml:$TARGET")
    testImageGeneric = env.AddProgram('xmipp_test_image_generic', 
                                      src=['applications/tests/test_image_generic'],
                                      incs=[Dir('.').path],
                                      libs=_libstests,
                                      deps=_depstests)
    AddXmippTest('test_image_generic', testImageGeneric, "$SOURCE --gtest_output=xml:$TARGET")
    testMatrix = env.AddProgram('xmipp_test_matrix', 
                               src=['applications/tests/test_matrix'],
                               incs=[Dir('.').path],
                               libs=_libstests,
                               deps=_depstests)
    AddXmippTest('test_matrix', testMatrix, "$SOURCE --gtest_output=xml:$TARGET")
    testMetadata = env.AddProgram('xmipp_test_metadata', 
                                  src=['applications/tests/test_metadata'],
                                  incs=[Dir('.').path],
                                  libs=_libstests,
                                  deps=_depstests)
    AddXmippTest('test_metadata', testMetadata, "$SOURCE --gtest_output=xml:$TARGET")
    testMultidim = env.AddProgram('xmipp_test_multidim', 
                                  src=['applications/tests/test_multidim'],
                                  incs=[Dir('.').path],
                                  libs=_libstests,
                                  deps=_depstests)
    AddXmippTest('test_multidim', testMultidim, "$SOURCE --gtest_output=xml:$TARGET")
    testPolar = env.AddProgram('xmipp_test_polar', 
                               src=['applications/tests/test_polar'],
                               incs=[Dir('.').path],
                               libs=_libstests,
                               deps=_depstests)
    AddXmippTest('test_polar', testPolar, "$SOURCE --gtest_output=xml:$TARGET")
    testPolynomials = env.AddProgram('xmipp_test_polynomials', 
                                     src=['applications/tests/test_polynomials'],
                                     incs=[Dir('.').path],
                                     libs=_libstests,
                                     deps=_depstests)
    AddXmippTest('test_polynomials', testPolynomials, "$SOURCE --gtest_output=xml:$TARGET")
    testSampling = env.AddProgram('xmipp_test_sampling', 
                                  src=['applications/tests/test_sampling'],
                                  incs=[Dir('.').path],
                                  libs=_libstests,
                                  deps=_depstests)
    AddXmippTest('test_sampling', testSampling, "$SOURCE --gtest_output=xml:$TARGET")
    testSymmetries = env.AddProgram('xmipp_test_symmetries', 
                                    src=['applications/tests/test_symmetries'],
                                    incs=[Dir('.').path],
                                    libs=_libstests,
                                    deps=_depstests)
    AddXmippTest('test_symmetries', testSymmetries, "$SOURCE --gtest_output=xml:$TARGET")
    testTransformation = env.AddProgram('xmipp_test_transformation', 
                                        src=['applications/tests/test_transformation'],
                                        incs=[Dir('.').path],
                                        libs=_libstests,
                                        deps=_depstests)
    AddXmippTest('test_transformation', testTransformation, "$SOURCE --gtest_output=xml:$TARGET")
    testDimred = env.AddProgram('xmipp_test_dimred', 
                                src=['applications/tests/test_dimred'],
                                incs=[Dir('.').path],
                                libs=_libstests,
                                deps=_depstests)
    AddXmippTest('test_dimred', testDimred, "$SOURCE --gtest_output=xml:$TARGET")
    testWavelets = env.AddProgram('xmipp_test_wavelets', 
                                  src=['applications/tests/test_wavelets'],
                                  incs=[Dir('.').path],
                                  libs=_libstests,
                                  deps=_depstests)
    AddXmippTest('test_wavelets', testWavelets, "$SOURCE --gtest_output=xml:$TARGET")
    testFilename = env.AddProgram('xmipp_test_filename', 
                                  src=['applications/tests/test_filename'],
                                  incs=[Dir('.').path],
                                  libs=_libstests,
                                  deps=_depstests)
    AddXmippTest('test_filename', testFilename, "$SOURCE --gtest_output=xml:$TARGET")

packageDeps = xmippParallel

Default(packageDeps)
Return('packageDeps')
