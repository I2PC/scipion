# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
This module contains the PATH related utilities 
inside the utils module
"""

import os
import shutil
from os.path import exists, join, splitext, isdir, expanduser, basename
import pyworkflow as pw


def findFile(filename, *paths):
    """ Search if the file is present in
    some path in the *paths provided.
    Return None if not found.
    """
    if filename:
        for p in paths:
            fn = join(p, filename)
            if exists(fn):
                return fn
    return None

def findResource(filename):
    """ This function will search for a give
    resource filename in the paths specified
    in pyworkflow.RESOURCES path list.
    """
    return findFile(filename, *pw.RESOURCES)

def replaceExt(filename, newExt):
    """ Replace the current path extension(from last .)
    with a new one. The new one should not contains the ."""
    return splitext(filename)[0] + '.' + newExt


def replaceBaseExt(filename, newExt):
    """ Replace the current basename extension(from last .)
    with a new one. The new one should not contains the .
    """
    return replaceExt(basename(filename), newExt)

def removeBaseExt(filename):
    """Take the basename of the filename and remove extension"""
    return removeExt(basename(filename))

def removeExt(filename):
    """ Remove extension from basename """
    return splitext(filename)[0]

def joinExt(*extensions):
    """Join several path parts with a ."""
    return '.'.join(extensions)

def cleanPath(*paths):
    """Remove a list of paths, either folders or files"""
    for p in paths:
        if isdir(p):
            shutil.rmtree(p)
        else:
            os.remove(p)
            
def makePath(*paths):
    """Create a list of paths if they don't exists"""
    for p in paths:
        if not exists(p):
            os.makedirs(p)
            
def existsPath(*paths):
    """Check if the list of paths exists.
    Will return the list of missing files,
    if the list is empty means that all path exists
    """
    return [p for p in paths if not exists(p)]

def getHomePath(user=''):
    """Return the home path of a give user."""
    return expanduser("~" + user)

def createFolderForFile(filePath):
    """
    Create folder for file in local host if it does not exist.
    filePath -- File path which parent directory we must create (/file path/...).
    """
    filePathParentDirectory = os.path.dirname(filePath)
    makePath(filePathParentDirectory)
