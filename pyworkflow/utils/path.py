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
from glob import glob

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
    """ Join several path parts with a ."""
    return '.'.join(extensions)

def getFile(filePath):
    """ Given a path, remove the folders and return the file """
    return split(filePath)[1]

def getExt(filePath):
    """ Return the extesion given a file. """
    return splitext(filePath)[1]

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
    Wrapper around the shutil.copytree, but allowing
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

def moveTree(src, dest):
    copyTree(src, dest)
    cleanPath(src)
                
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


# Console (and XMIPP) escaped colors, and the related tags that we create
# with Text.tag_config(). This dict is used in OutputText:addLine()
colorName = {'30': 'gray',
             '31': 'red',
             '32': 'green',
             '33': 'yellow',
             '34': 'blue',
             '35': 'magenta',
             '36': 'cyan',
             '37': 'white'}

def renderTextFile(fname, add, offset=0, lineNo=0, numberLines=True,
                   maxSize=100, headSize=20, tailSize=None):
    """
    Call callback function add() on each fragment of text from file fname,
    delimited by lines and/or color codes.
      add: callback function add(txt, tag='normal')
      offset: byte offset - we start reading the file from there
      lineNo: lines will be numbered from this value on
      numberLines: whether to prepend the line numbers
    """
    textfile = open(fname)
    size = (os.stat(fname).st_size - offset) / 1024  # in kB

    for line in iterBigFile(textfile, offset, size,
                            maxSize, headSize, tailSize):
        if line is not None:
            lineNo += 1
            renderLine(line, add, lineNo, numberLines)
        else:
            add("""\n
    ==> Too much data to read (%d kB) -- %d kB omitted
    ==> Click on """ % (size, size - headSize - (tailSize or headSize)))
            add(fname, 'link:%s' % fname)
            add(' to open it with the default viewer\n\n')
            if numberLines:
                add('    ==> Line numbers below are not '
                         'in sync with the input data\n\n')

    offset = textfile.tell()  # save last position in file
    textfile.close()

    return offset, lineNo


def renderLine(line, add, lineNo=1, numberLines=True):
    """
    Find all the fragments of formatted text in line and call
    add(fragment, tag) for each of them.
    """
    # Prepend line number
    if numberLines and lineNo:
        add('%05d:' % lineNo, 'cyan')
        add('   ')

    # iter 1\riter 2\riter 3  -->  iter 3
    if '\r' in line:
        line = line[line.rfind('\r')+1:]  # overwriting!

    # Find all console escape codes and use the appropriate tag instead.
    pos = 0  # current position in the line we are parsing
    while True:
        # line looks like: 'blah blah \x1b[XX;YYmremark\x1b[0m blah blah'
        # where YY is the color code (31 is red, for example).
        start = line.find('\x1b[', pos)
        if start < 0:
            add(line[pos:])
            break
        add(line[pos:start])
        if line[start+2:start+4] == '0m':
            start += 4
        colorCode = line[start+5:start+7]
        end = line.find('\x1b[0m', start + 7)
        if end < 0:  # what, no end for color formatting?
            # maybe we should also warn the user...
            add(line[start+8:])
            break
        add(line[start+8:end], colorName.get(colorCode, 'normal'))
        pos = end + 4


def iterBigFile(textfile, offset=0, size=None,
                maxSize=100, headSize=20, tailSize=None):
    """
    Yield lines from file textfile. If the size to read is bigger
    than maxSize then yield the first lines until headSize bytes, then
    yield None, then yield the last lines from tailSize bytes to the end.
    """
    if size is None:
        # Size in kB of the part of the file that we will read
        textfile.seek(0, 2)
        sizeKb = (textfile.tell() - offset) / 1024
    else:
        sizeKb = size

    headSizeB = headSize * 1024
    tailSizeB = (tailSize or headSize) * 1024

    textfile.seek(offset)
    if sizeKb > maxSize:  # max size that we want to read (in kB)
        for line in textfile.read(headSizeB).split('\n'):
            yield line + '\n'
        yield None  # Special result to mark omitting lines
        textfile.seek(-tailSizeB, 2)  # ready to show the last bytes

    # Add the remaining lines (from our last offset)
    for line in textfile:
        yield line
