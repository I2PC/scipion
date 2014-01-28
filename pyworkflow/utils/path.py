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
from os.path import exists, join, splitext, isdir, isfile, expanduser, expandvars, basename, dirname, split, relpath
import pyworkflow as pw
from glob import glob


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
    """ Join several path parts with a ."""
    return '.'.join(extensions)

def cleanPath(*paths):
    """ Remove a list of paths, either folders or files"""
    for p in paths:
        if exists(p):
            if isdir(p):
                shutil.rmtree(p)
            else:
                os.remove(p)
                
def cleanPattern(pattern):
    """ Remove all files that match the pattern. """
    files = glob(pattern)
    cleanPath(*files)
            
def makePath(*paths):
    """ Create a list of paths if they don't exists.
    Recursively create all folder needed in a path.
    If a path passed is a file, only the directory will be created.
    """
    for p in paths:
        if not exists(p) and len(p):
            os.makedirs(p)

def makeFilePath(*files):
    """ Make the path to ensure that files can be written. """
    makePath(*[dirname(f) for f in files])    
            
def missingPaths(*paths):
    """ Check if the list of paths exists.
    Will return the list of missing files,
    if the list is empty means that all path exists
    """
    return [p for p in paths if not exists(p)]

def getHomePath(user=''):
    """Return the home path of a give user."""
    return expanduser("~" + user)

def expandPattern(pattern, vars=True, user=True):
    """ Expand environment vars and user from a given pattern. """
    if vars:
        pattern = expandvars(pattern)
    if user:
        pattern = expanduser(pattern)
    return pattern

def getFiles(folderPath):
    """
    Gets all files of given folder and it subfolders.
    folderPath -- Folder path to get files.
    returns -- Set with all folder files.
    """
    filePaths = set()
    for path, dirs, files in os.walk(folderPath):
        for f in files:
            filePaths.add(join(path, f))
    return filePaths

def copyTree(source, dest):
    """
    Wrapper arount the shutil.copytree, but allowing
    that the dest folder also exists.
    """
    if not exists(dest):
        shutil.copytree(source, dest, symlinks=True)
    else:
        for f in os.listdir(source):
            fnPath = os.path.join(source, f)
            if isfile(fnPath):
                shutil.copy(fnPath, dest)
            elif os.path.isdir(fnPath):
                copyTree(fnPath, join(dest, f))
                
def copyFile(source, dest):
    """ Shortcut to shutil.copy. """
    shutil.copy(source, dest)
    
def moveFile(source, dest):
    """ Move file from source to dest. """
    copyFile(source, dest)
    cleanPath(source)
    
def createLink(source, dest):
    """ Creates a link to a given file path. """
    if exists(dest):
        os.remove(dest)
    destDir = split(dest)[0]
    os.symlink(relpath(source,destDir),dest)
    
def getLastFile(pattern):
    """ Return the last file matching the pattern. """
    files = glob(pattern)
    if len(files):
        files.sort()
        return files[-1]
    return None


