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

import os
from os.path import join, abspath, dirname
import sys
import platform

# for acceding SCIPION dict easily
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
PACKAGES_FOLDER =      DIR =                         4 
                    # folder name in temporal directory  
LIB_FOLDER =           TAR =                         5
                    # tarfile name in temporal directory
MAN_FOLDER =           DEPS =                        6 
                    # explicit dependencies              
TMP_FOLDER =           URL =                         7
                    # URL to download from               
INCLUDE_FOLDER =       FLAGS =                       8
                    # Other flags for the compiler
LOG_FOLDER =                                         9 #

# We start importing the environment
Import('env', 'SCIPION')

# Printing scipion Logo
env.ScipionLogo()

####################################
# BUILDING SCIPION MAIN DICTIONARY #
####################################

#################
# PREREQUISITES #
#################

USER_COMMANDS = False

# Python is needed. Any kind of python, to be able to execute this,
# but 2.7.6 version os python will be the one to execute Scipion. If
# the person doesn't have 2.7.6, SCons will compile it. Otherwise, a
# virtual environment will be used on top of system one, just if the
# user selects it. By default, python will be built
COMPILE_PYTHON = (env['PYVERSION'] != env['MANDATORY_PYVERSION']) and not USER_COMMANDS

# If we don't compile python, then we need a virtual environment on top of the system python
BUILD_VIRTUALENV = not COMPILE_PYTHON


############################
# FIRST LEVEL DEPENDENCIES #
############################

# Tcl/Tk
env.AddLibrary('tcl', 
               tar='tcl8.6.1-src.tar.gz', 
               dirname='tcl8.6.1', 
               src=os.path.join('tcl8.6.1', 'unix'))
# TODO: add flags here (and the same for all the other AddLibrary()s

# env.AddLibrary('tk',  
#                tar='tk8.6.1-src.tar.gz', 
#                dir='tk8.6.1',
#                src=os.path.join('tk8.6.1', 'unix'),
#                deps=['tcl'])

# sqlite
env.AddLibrary('sqlite',
               tar='sqlite-3.6.23.tgz',
               libs=['libsqlite3.so'])

# # zlib
# env.AddLibrary('zlib',
#                tar='zlib-1.2.8.tar.gz',
#                dir='zlib-1.2.8',
#                libs=['libz.so'])

# #############################
# # SECOND LEVEL DEPENDENCIES #
# #############################

# # python 2.7.8
# env.AddLibrary('python',
#                tar='Python-2.7.8.tgz',
#                url='http://scipionwiki.cnb.csic.es/files/scipion/software/python/Python-2.7.8.tgz',
#                libs=['libpython2.7.so'],
#                deps=['sqlite', 'tcl', 'tk'])

# ############################
# # THIRD LEVEL DEPENDENCIES #
# ############################

# # numpy
# env.AddLibrary('numpy',  
#                tar='numpy-1.8.1.tar.gz', 
#                dir='numpy-1.8.1', 
#                url='http://scipionwiki.cnb.csic.es/files/scipion/software/python/numpy-1.8.1.tar.gz',
#                deps=['python'])

# # matplotlib
# env.AddLibrary('matplotlib',  
#                tar='matplotlib-1.3.1.tar.gz', 
#                dir='matplotlib-1.3.1', 
#                url='http://scipionwiki.cnb.csic.es/files/scipion/software/python/matplotlib-1.3.1.tar.gz',
#                deps=['python'])

# # psutil
# env.AddLibrary('psutil', 
#                tar='psutil-2.1.1.tar.gz', 
#                dir='psutil-2.1.1', 
#                url='http://scipionwiki.cnb.csic.es/files/scipion/software/python/psutil-2.1.1.tar.gz',
#                deps=['python'])

# # mpi4py
# env.AddLibrary('mpi4py', 
#                tar='mpi4py-1.3.1.tar.gz', 
#                dir='mpi4py-1.3.1', 
#                url='http://scipionwiki.cnb.csic.es/files/scipion/software/python/mpi4py-1.3.1.tar.gz',
#                deps=['python'])

# # scipy
# env.AddLibrary('scipy', 
#                dft=False, 
#                tar='scipy-0.14.0.tar.gz', 
#                dir='scipy-0.14.0', 
#                url='http://scipionwiki.cnb.csic.es/files/scipion/software/python/scipy-0.14.0.tar.gz',
#                deps=['python'])

# # bibtex
# env.AddLibrary('bibtexparser', 
#                tar='bibtexparser-0.5.tgz',  
#                url='http://scipionwiki.cnb.csic.es/files/scipion/software/python/bibtexparser-0.5.tgz',
#                deps=['python'])

# # django
# env.AddLibrary('django', 
#                tar='Django-1.5.5.tgz', 
#                url='http://scipionwiki.cnb.csic.es/files/scipion/software/python/Django-1.5.5.tgz',
#                deps=['python'])

# # paramiko
# env.AddLibrary('paramiko', 
#                dft=False, 
#                tar='paramiko-1.14.0.tar.gz', 
#                dir='paramiko-1.14.0', 
#                url='http://scipionwiki.cnb.csic.es/files/scipion/software/python/paramiko-1.14.0.tar.gz',
#                deps=['python'])

# # PIL
# env.AddLibrary('pil',
#                tar='Imaging-1.1.7.tar.gz',
#                dir='Imaging-1.1.7',
#                url='http://scipionwiki.cnb.csic.es/files/scipion/software/python/Imaging-1.1.7_xmipp_setup.tgz',
#                deps=['python', 'zlib'])

# ######################
# # CONFIGURATION FILE #
# ######################
# # TODO: At this point, it is time to read the configuration file in order to alter (or not) the previously hard-coded libraries


# #############
# # DOWNLOADS #
# #############

# # for lib in SCIPION['LIBS']:
# #     env.Download('software/tmp/%s' % SCIPION['LIBS'][lib][TAR],
# #                  Value(SCIPION['LIBS'][lib][URL]))

# env.Download('software/tmp/%s' % SCIPION['LIBS']['sqlite'][TAR],
#              Value(SCIPION['LIBS']['sqlite'][URL]))

# #########
# # UNTAR #
# #########

# # def untar(pkg, target=None):
# #     return env.UntarLibrary(target if target else pkg,
# #                             'software/tmp/%s' % SCIPION['LIBS'][pkg][TAR])

# def untar(pkg):
#     env.UntarLibrary(target='software/tmp/%s/configure' % SCIPION['LIBS'][pkg][DIR],
#                      source='software/tmp/%s' % SCIPION['LIBS'][pkg][TAR])

# #sqliteUntar = untar('sqlite', 'software/tmp/sqlite-3.6.23/Makefile.in')
# sqliteUntar = untar('sqlite')
# # tclUntar = untar('tcl', 'software/tmp/tcl8.6.1/unix/Makefile.in')
# # tkUntar = untar('tk', 'software/tmp/tk8.6.1/unix/Makefile.in')
# # zlibUntar = untar('zlib', 'software/tmp/zlib-1.2.8/configure')
# # pythonUntar = untar('python', 'software/tmp/Python-2.7.8/Makefile.pre.in')
# # untar('numpy')
# # untar('matplotlib')
# # untar('psutil')
# # untar('mpi4py')
# # untar('scipy')
# # untar('bibtexparser')
# # untar('django')
# # untar('paramiko')
# # untar('pil')


# ##########################
# # EXECUTING COMPILATIONS #
# ##########################

# # EXTERNAL LIBRARIES
    
# t= env.AutoConfig(source=Dir('software/tmp/%s' % SCIPION['LIBS']['sqlite'][DIR]),
#                AutoConfigTarget='Makefile',
#                AutoConfigSource='configure',
#                AutoConfigParams=['CPPFLAGS=-w', 
#                                  'CFLAGS=-DSQLITE_ENABLE_UPDATE_DELETE_LIMIT=1', 
#                                  '--prefix=%s' % os.path.abspath(SCIPION['FOLDERS'][SOFTWARE_FOLDER])])

# # env.AutoConfig(source=Dir('software/tmp/%s' % SCIPION['LIBS']['sqlite'][DIR]),
# #                target='makefile-sqlite',
# #                AutoConfigTarget='Makefile',
# #                AutoConfigSource='configure',
# #                AutoConfigParams=['CPPFLAGS=-w', 
# #                                  'CFLAGS=-DSQLITE_ENABLE_UPDATE_DELETE_LIMIT=1', 
# #                                  '--prefix=%s' % os.path.abspath(SCIPION['FOLDERS'][SOFTWARE_FOLDER])])
# env['CROSS_BUILD'] = False
# env.Make(source=t,
#          target='software/lib/libsqlite3.so',
#          MakePath='software/tmp/%s' % SCIPION['LIBS']['sqlite'][DIR],
#          MakeEnv=os.environ,
#          MakeTargets='install',
#          MakeStdOut='aaaaa')

# # env.CompileLibrary('sqlite',
# #                    source=sqliteUntar,
# #                    flags=['CPPFLAGS=-w', 
# #                           'CFLAGS=-DSQLITE_ENABLE_UPDATE_DELETE_LIMIT=1', 
# #                           '--prefix=%s' % os.path.join(Dir(SCIPION['FOLDERS'][SOFTWARE_FOLDER]).abspath)], 
# #                    target='libsqlite3.so',
# #                    autoSource='configure') #SCIPION['LIBS']['sqlite'][DIR])

# # env.CompileLibrary('tcl', 
# #                    source=tclUntar,
# #                    flags=['--enable-threads', 
# #                           '--prefix=%s' % os.path.join(Dir(SCIPION['FOLDERS'][SOFTWARE_FOLDER]).abspath)],
# #                    target='libtcl.so',
# #                    autoSource=os.path.join('unix','Makefile.in'),
# #                    autoTarget=os.path.join('unix','Makefile'),
# #                    makePath='unix')

# # env.CompileLibrary('tk', 
# #                    source=tkUntar,
# #                    flags=['--enable-threads', 
# # #                          '--with-tcl="%s"' % os.path.join(Dir(os.path.join(SCIPION['FOLDERS'][SOFTWARE_FOLDER], 'lib64')).abspath),
# #                           '--prefix=%s' % os.path.join(Dir(SCIPION['FOLDERS'][SOFTWARE_FOLDER]).abspath)], 
# #                    target='libtk.so',
# #                    autoSource=os.path.join('unix','Makefile.in'),
# #                    autoTarget=os.path.join('unix','Makefile'),
# #                    makePath='unix')

# # env.CompileLibrary('zlib',
# #                    source=zlibUntar,
# #                    autoSource='configure',
# #                    autoTarget='configure.log',
# #                    flags=['--prefix=%s' % os.path.join(Dir(SCIPION['FOLDERS'][SOFTWARE_FOLDER]).abspath)],
# #                    target='libz.so')

# # # PYTHON

# # pythonMake = env.CompileLibrary('python',
# #                                 source=pythonUntar,
# #                                 flags=['--prefix=%s' % os.path.join(Dir(SCIPION['FOLDERS'][SOFTWARE_FOLDER]).abspath),
# #                                        'CFLAGS=-I/usr/include/ncurses'],
# #                                 target='libpython2.7.so',
# #                                 autoSource='Makefile.pre.in')

# # EMPackagesDeps = env.CompileWithSetupPy('python', deps=pythonMake)


# # # PYTHON MODULES

# # EMPackagesDeps += env.CompileWithSetupPy('numpy')
# # EMPackagesDeps += env.CompileWithSetupPy('matplotlib')
# # EMPackagesDeps += env.CompileWithSetupPy('psutil')
# # EMPackagesDeps += env.CompileWithSetupPy('mpi4py')
# # EMPackagesDeps += env.CompileWithSetupPy('scipy')
# # EMPackagesDeps += env.CompileWithSetupPy('bibtexparser')
# # EMPackagesDeps += env.CompileWithSetupPy('django')
# # EMPackagesDeps += env.CompileWithSetupPy('paramiko')
# # EMPackagesDeps += env.CompileWithSetupPy('pil')


# # # EM PACKAGES

# # # Xmipp3.1
# # env.AddPackage('xmipp', 
# #                dft=False,
# #                tar='Xmipp-3.1-src.tgz',
# #                dir='xmipp',
# #                url='http://xmipp.cnb.csic.es/Downloads/Xmipp-3.1-src.tar.gz')

# # # Bsoft
# # env.AddPackage('bsoft',
# #                dft=False,
# #                tar='bsoft1_8_8_Fedora_12.tgz',
# #                dir='bsoft')

# # # CtfFind
# # env.AddPackage('ctffind',
# #                dft=False,
# #                tar='ctffind_V3.5.tgz',
# #                dir='ctf')

# # # EMAN2
# # env.AddPackage('eman2',
# #                dft=False,
# #                tar='eman2.1beta3.linux64.tar.gz',
# #                dir='EMAN2')

# # # frealign
# # env.AddPackage('frealign',
# #                dft=False,
# #                tar='frealign_v9.07.tgz')

# # # relion
# # env.AddPackage('relion',
# #                dft=False,
# #                tar='relion-1.2.tgz')
# # # spider
# # env.AddPackage('spider',
# #                dft=False,
# #                tar='spider-web-21.13.tgz',
# #                dir='spider-web')

# # if not GetOption('clean'):
# #     env.Download('xmipp', type='EMPackage')
# #     env.Download('bsoft', type='EMPackage')
# #     env.Download('ctffind', type='EMPackage')
# #     env.Download('eman2', type='EMPackage')
# #     env.Download('frealign', type='EMPackage')
# #     env.Download('relion', type='EMPackage')
# #     env.Download('spider', type='EMPackage')

# # # Purge option
# # if GetOption('purge'):
# #     env.RemoveInstallation()

