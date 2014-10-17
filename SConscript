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
     buildDir=['hdf5-1.8.10/src', 'hdf5-1.8.10/c++'],
     configDir=['hdf5-1.8.10', 'hdf5-1.8.10'],
     flags=[[''], ['--enable-cxx']],
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

# XMIPP SHARED LIBRARIES
# XmippExternal
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
              prefix='xmipp',
              deps=tiff + hdf5 + lastTarget)

xmippSqliteExt = env.AddPackageLibrary(
               'XmippSqliteExt',
               dirs=['external/sqliteExt'],
               patterns=['extension-functions.c'],
               libs=['m'],
               deps=tiff + hdf5 + lastTarget)

# XmippData
xmippData = env.AddPackageLibrary(
          'XmippData',
          dirs=['libraries/data'],
          patterns=['*.cpp'],
          incs=[Dir('.').path, Dir('libraries').path],
          libs=['XmippExternal'],
          deps=xmippExternal + xmippSqliteExt)

# XmippClassif
xmippClassif = env.AddPackageLibrary(
             'XmippClassif',
             dirs=['libraries/classification'],
             patterns=['*.cpp'],
             incs=[Dir('.').path, Dir('libraries').path],
             libs=['XmippExternal', 'XmippData'],
             deps=xmippExternal + xmippSqliteExt + xmippData)

# XmippDimred
xmippDimred = env.AddPackageLibrary(
            'XmippDimred',
            dirs=['libraries/dimred'],
            incs=[Dir('.').path, Dir('libraries').path],
            libs=['XmippExternal', 'XmippData'],
            deps=xmippExternal + xmippSqliteExt + xmippData)

# XmippRecons
xmippRecons = env.AddPackageLibrary(
            'XmippRecons',
            dirs=['libraries/reconstruction'],
            patterns=['*.cpp'],
            incs=[Dir('.').path, Dir('libraries').path, Dir('libraries/reconstruction').path, Dir('external').path],
            libs=['XmippExternal', 'XmippData', 'XmippClassif'],
            deps=xmippExternal + xmippSqliteExt + xmippData + xmippClassif)
lastTarget = xmippRecons

# XMIPP PROGRAMS
# 

Default(lastTarget)
Return('lastTarget')