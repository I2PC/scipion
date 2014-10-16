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

# Jpeg
# not needed anymore since Scipion already has it
#jpeg = env.AddLibrary(
#     'jpeg',
#     tar='jpegsrc.v8c.tgz',
#     buildDir='jpeg-8c',
#     targets=[File('#software/lib/libjpeg.so').abspath],
#     clean=[Dir('#software/tmp/jpeg-8c').abspath],
#     url='http://scipionwiki.cnb.csic.es/files/xmipp/software/external/jpegsrc.v8c.tgz')
#Depends(jpeg, lastTarget)

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

# XmippExternal
xmippExternal = env.AddPackageLibrary(
              'XmippExternal',
              tars=['external/bilib.tgz', 
                    'external/condor.tgz', 
                    'external/alglib-3.8.0.cpp.tgz'],
              untarTargets=['bilib/configs.h',
                            'condor/Matrix.h',
                            'alglib-3.8.0.cpp/gpl3.txt'],
              incs=['external/bilib' + s for s in ['', '/headers', '/types']],
              src=['bilib/sources/*.cc', 
                   'condor/*.cpp',
                   'alglib-3.8.0.cpp/src/*.cpp'],
              dirs=['external',
                    'external',
                    'external'],
              prefix='xmipp',
              deps=tiff + hdf5 + lastTarget)
lastTarget = xmippExternal

xmippSqliteExt = env.AddPackageLibrary(
               'XmippSqliteExt',
               lib='XmippSqlite')

# XmippData

# XmippClassif

# XmippDimred

# XmippRecons

Default(xmippExternal)
Return('lastTarget')