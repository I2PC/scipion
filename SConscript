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

# We start importing the environment
Import('env')

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
#Python at any version (if 2.7.6, this will be used, otherwise, a new one will be compiled)

############################
# FIRST LEVEL DEPENDENCIES #
############################

# Tcl/Tk
env.AddLibrary('tcl',    dft=True, tar='tcl8.6.1-src.tar.gz', dir='tcl8.6.1',      url='http://scipionwiki.cnb.csic.es/files/scipion/software/external/tcl8.6.1-src.tar.gz')
env.AddLibrary('tk',     dft=True, tar='tk8.6.1-src.tar.gz',  dir='tk8.6.1',       url='http://scipionwiki.cnb.csic.es/files/scipion/software/external/tk8.6.1-src.tar.gz')

# sqlite
env.AddLibrary('sqlite', dft=True, tar='sqlite-3.6.23.tgz',   dir='sqlite-3.6.23', url='http://scipionwiki.cnb.csic.es/files/scipion/software/external/sqlite-3.6.23.tgz')


#############################
# SECOND LEVEL DEPENDENCIES #
#############################

# python 2.7.6
env.AddLibrary('python', dft=True, tar='Python-2.7.7.tgz',    dir='Python-2.7.7',  url='http://scipionwiki.cnb.csic.es/files/scipion/software/python/Python-2.7.7.tgz')

############################
# THIRD LEVEL DEPENDENCIES #
############################

# numpy

# matplotlib

# psutils

# mpi4py

# scipy

# bibtex

# django

# paramiko


#############
# DOWNLOADS #
#############

sqliteDownload = env.DownloadLibrary('sqlite')
#tclDownload = env.DownloadLibrary('tcl')
#tkDownload = env.DownloadLibrary('tk')
#pythonDownload = env.DownloadLibrary('python')


#########
# UNTAR #
#########
if not os.path.exists(os.path.join('software', 'tmp', 'sqlite-2.6.23.tgz')):
    print "Doesn't exist"
sqliteUntar = env.UntarLibrary('sqlite', tar=sqliteDownload)
#tclUntar = env.UntarLibrary('tcl', tar=tclDownload)
#tkUntar = env.UntarLibrary('tk', tar=tkDownload)
#pythonUntar = env.UntarLibrary('python', tar=pythonDownload)

# Explicitly write depends
Depends(sqliteUntar, sqliteDownload)
#Depends(tclUntar, tclDownload)
#Depends(tkUntar, tkDownload)
#Depends(pythonUntar, pythonDownload)

##########################
# EXECUTING COMPILATIONS #
##########################

env.Alias('sqlitePackage', sqliteDownload)
env.Alias('sqliteDir', sqliteUntar)
