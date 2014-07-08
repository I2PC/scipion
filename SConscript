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

# Python is needed. Any kind of python, to be able to execute this, but 2.7.6 version os python will be the one to execute Scipion. If the person doesn't have 2.7.6, SCons will compile it. Otherwise, a virtual environment will be used on top of system one, just if the user selects it. By default, python will be built
COMPILE_PYTHON = (env['PYVERSION'] != env['MANDATORY_PYVERSION']) and not USER_COMMANDS
# If python is not needed to be compiled, then a virtual environment is needed on top of the system python
BUILD_VIRTUALENV = not COMPILE_PYTHON

#Already compiled scons (using install.sh)
#Python at any version (if 2.7.7, this will be used, otherwise, a new one will be compiled)

############################
# FIRST LEVEL DEPENDENCIES #
############################

# Tcl/Tk
env.AddLibrary('tcl', 
               tar='tcl8.6.1-src.tar.gz', 
               dir='tcl8.6.1', 
               src=os.path.join('tcl8.6.1', 'unix'))
env.AddLibrary('tk',  
               tar='tk8.6.1-src.tar.gz', 
               dir='tk8.6.1',
               src=os.path.join('tk8.6.1', 'unix'),
               deps=['tcl'])

# sqlite
env.AddLibrary('sqlite',
               tar='sqlite-3.6.23.tgz',
               libs=['libsqlite3.so'])


#############################
# SECOND LEVEL DEPENDENCIES #
#############################

# python 2.7.7
env.AddLibrary('python',  
               tar='Python-2.7.7.tgz',
               url='http://scipionwiki.cnb.csic.es/files/scipion/software/python/Python-2.7.7.tgz',
               libs=['libpython2.7.so'],
               deps=['sqlite', 'tcl', 'tk'])

############################
# THIRD LEVEL DEPENDENCIES #
############################

# numpy
env.AddLibrary('numpy',  
               tar='numpy-1.8.1.tar.gz', 
               dir='numpy-1.8.1', 
               url='http://scipionwiki.cnb.csic.es/files/scipion/software/python/numpy-1.8.1.tar.gz',
               deps=['python'])

# matplotlib
env.AddLibrary('matplotlib',  
               tar='matplotlib-1.3.1.tar.gz', 
               dir='matplotlib-1.3.1', 
               url='http://scipionwiki.cnb.csic.es/files/scipion/software/python/matplotlib-1.3.1.tar.gz',
               deps=['python'])

# psutil
env.AddLibrary('psutil', 
               tar='psutil-2.1.1.tar.gz', 
               dir='psutil-2.1.1', 
               url='http://scipionwiki.cnb.csic.es/files/scipion/software/python/psutil-2.1.1.tar.gz',
               deps=['python'])

# mpi4py
env.AddLibrary('mpi4py', 
               tar='mpi4py-1.3.1.tar.gz', 
               dir='mpi4py-1.3.1', 
               url='http://scipionwiki.cnb.csic.es/files/scipion/software/python/mpi4py-1.3.1.tar.gz',
               deps=['python'])

# scipy
env.AddLibrary('scipy', 
               dft=False, 
               tar='scipy-0.14.0.tar.gz', 
               dir='scipy-0.14.0', 
               url='http://scipionwiki.cnb.csic.es/files/scipion/software/python/scipy-0.14.0.tar.gz',
               deps=['python'])

# bibtex
env.AddLibrary('bibtexparser', 
               tar='bibtexparser-0.5.tgz',  
               url='http://scipionwiki.cnb.csic.es/files/scipion/software/python/bibtexparser-0.5.tgz',
               deps=['python'])

# django
env.AddLibrary('django', 
               tar='Django-1.5.5.tgz', 
               url='http://scipionwiki.cnb.csic.es/files/scipion/software/python/Django-1.5.5.tgz',
               deps=['python'])

# paramiko
env.AddLibrary('paramiko', 
               dft=False, 
               tar='paramiko-1.14.0.tar.gz', 
               dir='paramiko-1.14.0', 
               url='http://scipionwiki.cnb.csic.es/files/scipion/software/python/paramiko-1.14.0.tar.gz',
               deps=['python'])

# PIL
env.AddLibrary('pil',
               tar='Imaging-1.1.7.tar.gz',
               dir='Imaging-1.1.7',
               url='http://scipionwiki.cnb.csic.es/files/scipion/software/python/Imaging-1.1.7.tar.gz',
               deps=['python'])

######################
# CONFIGURATION FILE #
######################
# TODO: At this point, it is time to read the configuration file in order to alter (or not) the previously hard-coded libraries


#############
# DOWNLOADS #
#############

if not GetOption('clean'):
    env.Download('sqlite')
    env.Download('tcl')
    env.Download('tk')
    env.Download('python')
    env.Download('numpy')
    env.Download('matplotlib')
    env.Download('psutil')
    env.Download('mpi4py')
    env.Download('scipy')
    env.Download('bibtexparser')
    env.Download('django')
    env.Download('paramiko')
    env.Download('pil')

#########
# UNTAR #
#########
sqliteUntar = env.UntarLibrary('sqlite')
tclUntar = env.UntarLibrary('tcl')
tkUntar = env.UntarLibrary('tk')
pythonUntar = env.UntarLibrary('python')
env.UntarLibrary('numpy')
env.UntarLibrary('matplotlib')
env.UntarLibrary('psutil')
env.UntarLibrary('mpi4py')
env.UntarLibrary('scipy')
env.UntarLibrary('bibtexparser')
env.UntarLibrary('django')
env.UntarLibrary('paramiko')
env.UntarLibrary('pil')

##########################
# EXECUTING COMPILATIONS #
##########################

# EXTERNAL LIBRARIES

env.CompileLibrary('sqlite',
                   source=sqliteUntar,
                   flags=['CPPFLAGS=-w', 
                          'CFLAGS=-DSQLITE_ENABLE_UPDATE_DELETE_LIMIT=1', 
                          '--prefix=%s' % os.path.join(Dir(SCIPION['FOLDERS'][SOFTWARE_FOLDER]).abspath)], 
                   target='libsqlite3.so')

env.CompileLibrary('tcl', 
                   source=tclUntar,
                   flags=['--enable-threads', 
                          '--prefix=%s' % os.path.join(Dir(SCIPION['FOLDERS'][SOFTWARE_FOLDER]).abspath)],
                   target='libtcl.so',
                   autoSource=os.path.join('unix','Makefile.in'),
                   autoTarget=os.path.join('unix','Makefile'),
                   makePath='unix')

env.CompileLibrary('tk', 
                   source=tkUntar,
                   flags=['--enable-threads', 
#                          '--with-tcl="%s"' % os.path.join(Dir(os.path.join(SCIPION['FOLDERS'][SOFTWARE_FOLDER], 'lib64')).abspath),
                          '--prefix=%s' % os.path.join(Dir(SCIPION['FOLDERS'][SOFTWARE_FOLDER]).abspath)], 
                   target='libtk.so',
                   autoSource=os.path.join('unix','Makefile.in'),
                   autoTarget=os.path.join('unix','Makefile'),
                   makePath='unix')


# PYTHON

pythonMake = env.CompileLibrary('python',
                                source=pythonUntar,
                                flags=['--prefix=%s' % os.path.join(Dir(SCIPION['FOLDERS'][SOFTWARE_FOLDER]).abspath),
                                       'CFLAGS=-I/usr/include/ncurses'],
                                target='libpython2.7.so',
                                autoSource='Makefile.pre.in')

EMPackagesDeps = env.CompileWithSetupPy('python', deps=pythonMake)


# PYTHON MODULES

EMPackagesDeps += env.CompileWithSetupPy('numpy')
EMPackagesDeps += env.CompileWithSetupPy('matplotlib')
EMPackagesDeps += env.CompileWithSetupPy('psutil')
EMPackagesDeps += env.CompileWithSetupPy('mpi4py')
EMPackagesDeps += env.CompileWithSetupPy('scipy')
EMPackagesDeps += env.CompileWithSetupPy('bibtexparser')
EMPackagesDeps += env.CompileWithSetupPy('django')
EMPackagesDeps += env.CompileWithSetupPy('paramiko')
EMPackagesDeps += env.CompileWithSetupPy('pil')


# EM PACKAGES

# Xmipp3.1
env.AddPackage('xmipp', 
               dft=False,
               tar='Xmipp-3.1-src.tgz',
               dir='xmipp',
               url='http://xmipp.cnb.csic.es/Downloads/Xmipp-3.1-src.tar.gz')

# Bsoft
env.AddPackage('bsoft',
               dft=False,
               tar='bsoft1_8_8_Fedora_12.tgz',
               dir='bsoft')

# CtfFind
env.AddPackage('ctffind',
               dft=False,
               tar='ctffind_V3.5.tgz',
               dir='ctf')

# EMAN2
env.AddPackage('eman2',
               dft=False,
               tar='eman2.1beta3.linux64.tar.gz',
               dir='EMAN2')

# frealign
env.AddPackage('frealign',
               dft=False,
               tar='frealign_v9.07.tgz')

# relion
env.AddPackage('relion',
               dft=False,
               tar='relion-1.2.tgz')
# spider
env.AddPackage('spider',
               dft=False,
               tar='spider-web-21.13.tgz',
               dir='spider-web')

if not GetOption('clean'):
    env.Download('xmipp', type='EMPackage')
    env.Download('bsoft', type='EMPackage')
    env.Download('ctffind', type='EMPackage')
    env.Download('eman2', type='EMPackage')
    env.Download('frealign', type='EMPackage')
    env.Download('relion', type='EMPackage')
    env.Download('spider', type='EMPackage')

# Purge option
if GetOption('purge'):
    env.RemoveInstallation()

