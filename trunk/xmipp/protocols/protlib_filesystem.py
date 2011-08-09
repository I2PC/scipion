'''
/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
 '''

#---------------------------------------------------------------------------
# Filesystem utilities
#---------------------------------------------------------------------------
import os
from protlib_utils import printLog
from shutil import copyfile
#from xmipp import *

# The following are Wrappers to be used from Protocols
# provinding filesystem utitities

def createDir(log, path):
    """ Create directory, does not add workingdir"""
    from distutils.dir_util import mkpath
    from distutils.errors import DistutilsFileError
    try:
        if not os.path.exists(path):
            mkpath(path, 0777, True)
            printLog("Created dir " + path, log)
    except DistutilsFileError, e:
        printLog(redStr("Couldn't create dir: '%(path)s': %(e)s" % locals()), 
                 log, err=True, isError=True)

def changeDir(log, path):
    """ Change to Directory """
    try:
        os.chdir(path)
        printLog("Changed to dir " + path,log)
    except os.error, (errno, errstr):
        printLog("Could not change to directory '%s': Error (%d): %s" % (path, errno, errstr),log,err=True,isError=True)

def deleteDir(log, path):
    from distutils.dir_util import remove_tree
    if os.path.exists(path):
        remove_tree(path, True)
        printLog("Deleted directory " + path, log)
           
def deleteFile(log, filename, verbose):
    if os.path.exists(filename):
        os.remove(filename)
        if verbose:
            printLog('Deleted file %s' % filename,log)
    else:
        if verbose:
            printLog('Do not need to delete %s; already gone' % filename,log)
            
def copyFile(log, source, dest):
    try:
        copyfile(source, dest)
        printLog("Copied '%s' to '%s'" % (source, dest))
    except Exception, e:
        printLog("Could not copy '%s' to '%s'. Error: %s" % (source, dest, str(e)), log, err=True, isError=True)
    

def deleteFiles(log, filelist, verbose):
    for file in filelist:
        deleteFile(log, file, verbose)
        
#--------------------------- Xmipp specific tools ---------------------------------
def getXmippPath(subfolder=''):
    '''Return the path the the Xmipp installation folder
    if a subfoder is provided, will be concatenated to the path'''
    protdir = os.popen('which xmipp_protocols', 'r').read()
    xmippdir = os.path.dirname(os.path.dirname(protdir))
    return os.path.join(xmippdir, subfolder)

def includeProtocolsDir():
    protDir = getXmippPath('protocols')
    import sys
    sys.path.append(protDir)

def getProtocolTemplate(prot):
    protDir = getXmippPath('protocols')
    srcProtName = '%s.py' % prot.key
    #srcProtDir = getXmippPath('protocols')
    srcProtAbsPath = os.path.join(protDir, srcProtName)
    return srcProtAbsPath
