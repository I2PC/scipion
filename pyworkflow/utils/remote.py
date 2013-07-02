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
This modules contains some wrappers over the ssh paramiko library
to provide handling functions on remote hosts through ssh.

Main classes are:

    RemotePath: this class use an ssh/sftp connection to mimic the 
        functions in utils.path but remotely and provides basic file transfer
        
        
"""

import os, shutil, errno
import socket, paramiko, hashlib

from pyworkflow.utils.path import *
from pyworkflow.utils.log import *

LOCAL_USER_AND_HOST = ''
SSH_PORT = 22
PAIRS_SEPARATOR = ':'

log = getGeneralLogger('pyworkflow.utils.remote')



def testHostConfig(host):
    """ Test the connection to a remote host give its configuration.(see HostConfig class)
    Params:
        host: configuration of the remote host, should contains hostName, userName and password.
    Returns: True if the host could be reached.
    """    
    try:
        ssh = sshConnectFromHost(host)
        rpath = RemotePath(ssh)
        rpath.listdir('.')
        rpath.close()
        return True
    except Exception, ex:
        return False

def sshConnect(hostName, userName, password, port=SSH_PORT, **args):
    """ Common way to create a ssh connection.
    Params:
        hostName: Remote host name.
        userName: User name.
        password: Password.
        port: port to establish connection (usually 22)
    Returns: ssh connection handler.
    """
    ssh = paramiko.SSHClient()
    ssh.load_system_host_keys()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(hostName, port, userName, password, **args)
    return ssh

def sshConnectFromHost(host):
    """ Establish a connection receiving a HostConfig class. """
    return sshConnect(host.getHostName(), host.getUserName(), host.getPassword())
    
class RemotePath(object):
    """ This class will server as a Wrapper to the
    paramiko sftp protocol througth a ssh connection.
    This class will implement useful methods for remote
    path handling such as: creating files, deleting folders...
    """
    
    @classmethod
    def fromCredentials(cls, hostName, userName, password, port=SSH_PORT, **args):
        ssh = sshConnect(hostName, userName, password, port, **args)
        rpath = RemotePath(ssh)        
        return rpath
    
    def __init__(self, ssh):
        self.ssh = ssh
        self.sftp = ssh.open_sftp()
        # Shortcut some of the sftp methods
        self.listdir = self.sftp.listdir
        
    def exists(self, path):
        """ Check if a remote path exists(like os.path.exists remotely). """
        try:
            self.sftp.stat(path)
        except IOError, e:
            if e.errno == errno.ENOENT:
                return False
        else:
            return True
        
    def isdir(self, path):
        """ Check if a remote path is a directory (like os.path.isdir remotely). """
        try:
            from stat import S_ISDIR
            return S_ISDIR(self.sftp.stat(path).st_mode)
        except IOError:
            #Path does not exist, so by definition not a directory
            return False
        
    def makedirs(self, remoteFolder):
        """ Like os.makedirs remotely. """
        if len(remoteFolder) and not self.exists(remoteFolder):
            log.info('RemotePath.makedirs, path: %s' % remoteFolder)
            parent = dirname(remoteFolder)
            # if have parent and it doen't exist, create it recursively
            if len(parent) and not self.exists(parent): 
                self.makedirs(parent)
            self.sftp.mkdir(remoteFolder)
        
    def getFile(self, remoteFile, localFile):
        """ Wrapper around sftp.get that ensures
        path exists for localFile.
        """
        makeFilePath(localFile)
        self.sftp.get(remoteFile, localFile)
        
    def putFile(self, localFile, remoteFile):
        """ Wrapper around sftp.put that ensures
        the remote path exists for put the file.
        """
        self.makeFilePath(remoteFile)
        log.info('RemotePath.putFile, local=%s, remote=%s' % (localFile, remoteFile))
        self.sftp.put(localFile, remoteFile)
        
    def makeFilePath(self, *remoteFiles):
        """ Create the remote folder path for remoteFiles. """
        self.makePath(*[dirname(r) for r in remoteFiles])
        
    def makePath(self, *remoteFolders):
        """ Make all path in remoteFolders list. """
        for p in remoteFolders:
            if len(p):
                self.makedirs(p)
        
    def close(self):
        """ Close both ssh and sftp connections. """
        self.sftp.close()
        self.ssh.close()
        
    def copyTree(self, localFolder, remoteFolder):
        """
        Same as shutil.copytree, but allowing
        that the dest folder also exists.
        dest is a remote folder.
        """
        if not self.exists(remoteFolder):
            self.makePath(remoteFolder)
            
        for path in os.listdir(localFolder):
            localPath = join(localFolder, path)
            remotePath = join(remoteFolder, path)
            # Create folders
            if isdir(localPath):
                self.copyTree(localPath, remotePath)
            elif isfile(localPath):
                self.putFile(localPath, remotePath)
                
    def cleanPath(self, *remotePaths):
        """
        Same as pyworkflow.utils.path.cleanPath but for remote folders or files.
        """
        for path in remotePaths:
            if self.exists(path):
                # Create folders
                if self.isdir(path):
                    files = [join(path, f) for f in self.sftp.listdir(path)]
                    self.cleanPath(*files)
                    log.info('RemotePath.cleanPath, rmdir path: %s' % path)
                    self.sftp.rmdir(path)
                else:# isfile(path):
                    log.info('RemotePath.cleanPath, remove path: %s' % path)
                    self.sftp.remove(path)
    
if __name__ == '__main__':
    pass
