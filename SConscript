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
env.AddLibrary('tcl', dft=True, tar='tcl8.6.1-src.tar.gz', dir='tcl8.6.1', url='http://scipionwiki.cnb.csic.es/files/scipion/software/external/tcl8.6.1-src.tar.gz')
env.AddLibrary('tk', dft=True, tar='tk8.6.1-src.tar.gz', dir='tk8.6.1', url='http://scipionwiki.cnb.csic.es/files/scipion/software/external/tk8.6.1-src.tar.gz')

# sqlite
env.AddLibrary('sqlite', dft=True, tar='sqlite-3.6.23.tgz', dir='sqlite-3.6.23', url='http://scipionwiki.cnb.csic.es/files/scipion/software/external/sqlite-3.6.23.tgz')


#############################
# SECOND LEVEL DEPENDENCIES #
#############################

# python 2.7.7
env.AddLibrary('python', dft=True, tar='Python-2.7.7.tgz', dir='Python-2.7.7', url='http://scipionwiki.cnb.csic.es/files/scipion/software/python/Python-2.7.7.tgz')

############################
# THIRD LEVEL DEPENDENCIES #
############################

# numpy
env.AddLibrary('numpy', dft=True, tar='numpy-1.8.1.tar.gz', dir='numpy-1.8.1', url='http://scipionwiki.cnb.csic.es/files/scipion/software/python/numpy-1.8.1.tar.gz')

# matplotlib
env.AddLibrary('matplotlib', dft=True, tar='matplotlib-1.3.1.tar.gz', dir='matplotlib-1.3.1', url='http://scipionwiki.cnb.csic.es/files/scipion/software/python/matplotlib-1.3.1.tar.gz')

# psutil
env.AddLibrary('psutil', dft=True, tar='psutil-2.1.1.tar.gz', dir='psutil-2.1.1', url='http://scipionwiki.cnb.csic.es/files/scipion/software/python/psutil-2.1.1.tar.gz')

# mpi4py
env.AddLibrary('mpi4py', dft=True, tar='mpi4py-1.3.1.tar.gz', dir='mpi4py-1.3.1', url='http://scipionwiki.cnb.csic.es/files/scipion/software/python/mpi4py-1.3.1.tar.gz')

# scipy
env.AddLibrary('scipy', dft=False, tar='scipy-0.14.0.tar.gz', dir='scipy-0.14.0', url='http://scipionwiki.cnb.csic.es/files/scipion/software/python/scipy-0.14.0.tar.gz')

# bibtex
env.AddLibrary('bibtexparser', dft=True, tar='bibtexparser-0.5.tgz', dir='bibtexparser-0.5', url='http://scipionwiki.cnb.csic.es/files/scipion/software/python/bibtexparser-0.5.tgz')

# django
env.AddLibrary('django', dft=True, tar='Django-1.5.5.tgz', dir='Django-1.5.5', url='http://scipionwiki.cnb.csic.es/files/scipion/software/python/Django-1.5.5.tgz')

# paramiko
env.AddLibrary('paramiko', dft=False, tar='paramiko-1.14.0.tar.gz', dir='paramiko-1.14.0', url='http://scipionwiki.cnb.csic.es/files/scipion/software/python/paramiko-1.14.0.tar.gz')

######################
# CONFIGURATION FILE #
######################
# TODO: At this point, it is time to read the configuration file in order to alter (or not) the previously hard-coded libraries


#############
# DOWNLOADS #
#############

sqliteDownload = env.DownloadLibrary('sqlite')
tclDownload = env.DownloadLibrary('tcl')
tkDownload = env.DownloadLibrary('tk')
pythonDownload = env.DownloadLibrary('python')
numpyDownload = env.DownloadLibrary('numpy')
matplotlibDownload = env.DownloadLibrary('matplotlib')
psutilDownload = env.DownloadLibrary('psutil')
mpi4pyDownload = env.DownloadLibrary('mpi4py')
scipyDownload = env.DownloadLibrary('scipy')
bibtexparserDownload = env.DownloadLibrary('bibtexparser')
djangoDownload = env.DownloadLibrary('django')
paramikoDownload = env.DownloadLibrary('paramiko')

#########
# UNTAR #
#########
sqliteUntar = env.UntarLibrary('sqlite')
tclUntar = env.UntarLibrary('tcl')
tkUntar = env.UntarLibrary('tk')
pythonUntar = env.UntarLibrary('python')
numpyUntar = env.UntarLibrary('numpy')
matplotlibUntar = env.UntarLibrary('matplotlib')
psutilUntar = env.UntarLibrary('psutil')
mpi4pyUntar = env.UntarLibrary('mpi4py')
scipyUntar = env.UntarLibrary('scipy')
bibtexparserUntar = env.UntarLibrary('bibtexparser')
djangoUntar = env.UntarLibrary('django')
paramikoUntar = env.UntarLibrary('paramiko')

# Explicitly write depends
#Depends(sqliteUntar, sqliteDownload)
#Depends(tclUntar, tclDownload)
#Depends(tkUntar, tkDownload)
#Depends(pythonUntar, pythonDownload)
#Depends(numpyUntar, numpyDownload)
#Depends(matplotlibUntar, matplotlibDownload)
#Depends(psutilUntar, psutilDownload)
#Depends(mpi4pyUntar, mpi4pyDownload)
#Depends(scipyUntar, scipyDownload)
#Depends(bibtexparserUntar, bibtexparserDownload)
#Depends(djangoUntar, djangoDownload)
#Depends(paramikoUntar, paramikoDownload)

##########################
# EXECUTING COMPILATIONS #
##########################

print "\t target %s \n \tdepends on source %s" % (os.path.join(SCIPION['FOLDERS'][TMP_FOLDER], SCIPION['LIBS']['sqlite'][DIR], '.libs', 'libsqlite3.so'), Glob(os.path.join(SCIPION['FOLDERS'][TMP_FOLDER], SCIPION['LIBS']['sqlite'][DIR]), '*.c')[0])

sqliteCompile = env.CompileLibrary('sqlite',
                                   flags=['CPPFLAGS=-w', 
                                          'CFLAGS=-DSQLITE_ENABLE_UPDATE_DELETE_LIMIT=1', 
                                          '--prefix=%s' % os.path.join(Dir(SCIPION['FOLDERS'][SOFTWARE_FOLDER]).abspath)], 
                                   source=Glob(os.path.join(SCIPION['FOLDERS'][TMP_FOLDER], SCIPION['LIBS']['sqlite'][DIR]), '*.c'),
                                   target=File('libsqlite3.so'),
                                   autoSource='Makefile.in',
                                   autoTarget='Makefile')

tclCompile = env.CompileLibrary('tcl', 
                                flags=['--enable-threads', 
                                       '--prefix=%s' % os.path.join(Dir(SCIPION['FOLDERS'][SOFTWARE_FOLDER]).abspath)],
                                source=Glob(os.path.join(SCIPION['FOLDERS'][TMP_FOLDER], SCIPION['LIBS']['tcl'][DIR]), 'unix', '*.c'),
                                target=File('libtcl.so'),
                                autoSource=os.path.join('unix','Makefile.in'),
                                autoTarget=os.path.join('unix','Makefile'),
                                makePath='unix')
#tkCompile = env.CompileLibrary('tk', 
#                               flags=['--enable-threads', 
#                                      '--with-tcl="%s"' % os.path.join(Dir(os.path.join(SCIPION['FOLDERS'][TMP_FOLDER], SCIPION['LIBS']['tcl'][DIR],'unix')).abspath),
#                                      '--prefix=%s' % os.path.join(Dir(SCIPION['FOLDERS'][SOFTWARE_FOLDER]).abspath)], 
#                               libs=['.libs'],
#                               deps=tkUntar+tclUntar,
#                               autoSource=os.path.join('unix','Makefile.in'),
#                               autoTarget=os.path.join('unix','Makefile'),
#                               makePath='unix')


