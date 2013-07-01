'''
Created on Apr 8, 2013

@author: antonio
'''
import os
from os.path import join
import shutil
import socket
import paramiko
import hashlib

from pyworkflow.utils.path import *
from pyworkflow.utils.log import *

LOCAL_USER_AND_HOST = ''
SSH_PORT = 22
PAIRS_SEPARATOR = ':'

log = getGeneralLogger('pyworkflow.utils.file_transfer')

class FileTransfer():
    
    ssh = None   
    sftp = None 
    
    def __init__(self):
        # Default ssh session options.
        self.ssh = paramiko.SSHClient()
        self.ssh.load_system_host_keys()
        self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        
        
        
    def transferFiles(self,
                      filePaths, 
                      hostsPasswords, 
                      gatewayHosts = None, 
                      numberTrials = 1,                        
                      forceOperation = False,
                      operationId = 1):
        """
        filePaths -- Files dictionary with this format: "userName@hostName:absolute_file_path": ["userName1@hostName1:absolute_file_path1", "userName2@hostName2:absolute_file_path2"]
        Key is the source file path and value the target file paths.
        gatewayHosts -- Gateway hosts dictionary with this format: "{userName1@hostName1:userName2@hostName2":["userName3@hostName3","userName4@hostName4"]}
        hostsPasswords -- Passwords needed to connect to involved hosts with this format: "userName@hostName":"hostPassword"
        numberTrials -- Number of trials in error cases.
        forceOperation -- Flag to indicate if, when an error happens and number of trials is exceeded, the operation must continue with the rest of files.        
        operationId -- Operation identifier.
        """
        classifiedFiles = self.__classifyFilePaths(filePaths)
        
        for userAndHostPairs, cursorFilePaths in classifiedFiles.iteritems():
            gatewayHostsCredentials = None
            if (gatewayHosts is not None and userAndHostPairs in gatewayHosts):
                gatewayHostsCredentials = gatewayHosts[userAndHostPairs]
            self.__sendFilesBetweenPairs(userAndHostPairs, 
                                         cursorFilePaths,
                                         hostsPasswords, 
                                         gatewayHostsCredentials, 
                                         numberTrials,                        
                                         forceOperation,
                                         operationId) 
    
    
    def copyFiles(self,
                  filePaths,                        
                  numberTrials = 1,                        
                  forceOperation = False,
                  operationId = 1):
        """
        filePaths -- Files dictionary with this format: "source_file_path": "target_file_path"
        forceOperation -- Flag to indicate if, when an error happens and number of trials is exceeded, the operation must continue with the rest of files.        
        operationId -- Operation identifier.
        """
        for sourceFilePath, targetFilePath in filePaths.iteritems(): 
            self.__copyLocalFile(sourceFilePath, targetFilePath)
    
    def transferFilesTo(self,
                        filePaths,
                        hostName,
                        userName,
                        hostPassword,
                        gatewayHosts = None, 
                        numberTrials = 1,                        
                        forceOperation = False,
                        operationId = 1):
        """
        filePaths -- Files dictionary with this format: "source_file_path": "target_file_path"
        hostName -- Remote host to transfer files.
        username -- User name for remote host.
        hostsPassword -- Passwords needed to connect to involved host.
        gatewayHosts -- Gateway hosts List with this format: ["userName1@hostName1","userName2@hostName2"]
        numberTrials -- Number of trials in error cases.
        forceOperation -- Flag to indicate if, when an error happens and number of trials is exceeded, the operation must continue with the rest of files.        
        operationId -- Operation identifier.
        """
        # Create ssh session to remote host
        log.info("Connecting to: " + userName + "@" + hostName)
        self.ssh.connect(hostName, SSH_PORT, userName, hostPassword)
        self.sftp = self.ssh.open_sftp()
        for sourceFilePath, targetFilePath in filePaths.iteritems():
            self.__sendLocalFile(sourceFilePath, targetFilePath, gatewayHosts, self.sftp)
        self.ssh.close()
        self.sftp.close()  
            
    def transferFilesFrom(self,
                        filePaths,
                        hostName,
                        userName,
                        hostPassword,
                        gatewayHosts = None, 
                        numberTrials = 1,                        
                        forceOperation = False,
                        operationId = 1):
        """
        filePaths -- Files dictionary with this format: "source_file_path": "target_file_path"
        hostName -- Remote host to transfer files.
        username -- User name for remote host.
        hostsPassword -- Passwords needed to connect to involved host.
        gatewayHosts -- Gateway hosts List with this format: ["userName1@hostName1","userName2@hostName2"]
        numberTrials -- Number of trials in error cases.
        forceOperation -- Flag to indicate if, when an error happens and number of trials is exceeded, the operation must continue with the rest of files.        
        operationId -- Operation identifier.
        """
        # Create ssh session to remote host
        log.info("Connecting to: " + userName + "@" + hostName)
        self.ssh.connect(hostName, SSH_PORT, userName, hostPassword)
        self.sftp = self.ssh.open_sftp()
        for sourceFilePath, targetFilePath in filePaths.iteritems():
            self.__getRemoteFile(sourceFilePath, targetFilePath, gatewayHosts, self.sftp)
        self.ssh.close()
        self.sftp.close()  
                
        
    def deleteFiles(self, 
                    filePaths,                    
                    hostsPasswords, 
                    gatewayHosts=None, 
                    numberTrials=1, 
                    forceOperation=False, 
                    operationId=1):
        """
        Delete a list of file paths.
        filepaths -- List of file paths to remove with this format: "userName@hostName:absolute_file_path"
        gatewayHosts -- Gateway hosts dictionary with this format: "userName1@hostName1:userName2@hostName2":"userName@hostName"
        hostsPasswords -- Passwords needed to connect to involved hosts with this format: "userName@hostName":"hostPassword"
        numberTrials -- Number of trials in error cases.
        forceOperation -- Flag to indicate if, when an error happens and number of trials is exceeded, the operation must continue with the rest of files.
        operationId -- Operation identifier.
        """
        
        # As we are going to create a session for each target host we must get different target hosts.
        userAndHosts = self.__getDiferentHostsFromFilePaths(filePaths)
        for userAndHost in userAndHosts:
            if (not self.__isLocalCredential(userAndHost)):
                resultFilePaths = self.__getFilePaths(filePaths, userAndHost)
                # Recover host credentials to remove files
                userName = self.__getUserAndHost(userAndHost)[0]
                hostName = self.__getUserAndHost(userAndHost)[1]
                hostPassword = hostsPasswords[userAndHost]
                # Create ssh session to remote host
                log.info("Connecting to: " + userName + "@" + hostName)
                self.ssh.connect(hostName, SSH_PORT, userName, hostPassword)
                self.sftp = self.ssh.open_sftp()
                for resultFilePath in resultFilePaths:
                    filePath = self.__getLocationAndFilePath(resultFilePath)[1]
                    log.info("Deleting file " + filePath)
                    self.sftp.remove(filePath)
                self.ssh.close()
                self.sftp.close()  
            else:
                pass
            
    def deleteDirectories(self, 
                          directoryPaths,                    
                          hostsPasswords, 
                          gatewayHosts=None, 
                          numberTrials=1, 
                          forceOperation=False, 
                          operationId=1):
        """
        Delete a list of directory paths.
        directoryPaths -- List of file paths to remove with this format: "userName@hostName:absolute_directory_path"
        gatewayHosts -- Gateway hosts dictionary with this format: "userName1@hostName1:userName2@hostName2":"userName@hostName"
        hostsPasswords -- Passwords needed to connect to involved hosts with this format: "userName@hostName":"hostPassword"
        numberTrials -- Number of trials in error cases.
        forceOperation -- Flag to indicate if, when an error happens and number of trials is exceeded, the operation must continue with the rest of files.
        operationId -- Operation identifier.
        """
        
        # As we are going to create a session for each target host we must get different target hosts.
        userAndHosts = self.__getDiferentHostsFromFilePaths(directoryPaths)
        for userAndHost in userAndHosts:
            if (not self.__isLocalCredential(userAndHost)):
                resultDirectoryPaths = self.__getFilePaths(directoryPaths, userAndHost)
                # Recover host credentials to remove files
                userName = self.__getUserAndHost(userAndHost)[0]
                hostName = self.__getUserAndHost(userAndHost)[1]
                hostPassword = hostsPasswords[userAndHost]
                # Create ssh session to remote host
                log.info("Connecting to: " + userName + "@" + hostName)
                self.ssh.connect(hostName, SSH_PORT, userName, hostPassword)
                self.sftp = self.ssh.open_sftp()
                for resultDirectoryPath in resultDirectoryPaths:
                    directoryName = self.__getLocationAndFilePath(resultDirectoryPath)[1]
                    log.info("Deleting directory " + directoryName)
                    self.sftp.rmdir()
                self.ssh.close()
                self.sftp.close()  
            else:
                pass
        
    def checkFiles(self, 
                   filePaths,                   
                   hostsPasswords,
                   gatewayHosts=None, 
                   numberTrials=1, 
                   forceOperation = False,  
                   operationId=1):
        """
        Check if file paths exists.
        filepaths -- List of file paths to check with this format: "userName@hostName:absolute_file_path"
        gatewayHosts -- Gateway hosts dictionary with this format: "userName1@hostName1:userName2@hostName2":"userName@hostName"
        hostsPasswords -- Passwords needed to connect to involved hosts with this format: "userName@hostName":"hostPassword"
        numberTrials -- Number of trials in error cases.
        forceOperation -- Flag to indicate if, when an error happens and number of trials is exceeded, the operation must continue with the rest of files.
        operationId -- Operation identifier.
        returns -- List of not located file paths.
        """
        returnFilePaths = []
        # As we are going to create a session for each target host we must get different target hosts.
        userAndHosts = self.__getDiferentHostsFromFilePaths(filePaths)
        for userAndHost in userAndHosts:
            resultFilePaths = self.__getFilePaths(filePaths, userAndHost)
            if (self.__isLocalCredential(userAndHost)):
                for resultFilePath in resultFilePaths:
                    filePath = self.__getLocationAndFilePath(resultFilePath)[1]
                    log.info("Checking: " + filePath)
                    if (len (missingPaths(filePath)) != 0):
                        returnFilePaths.append(filePath)
                        log.info("Check fail!!")
            else:
                
                # Recover host credentials
                userName = self.__getUserAndHost(userAndHost)[0]
                hostName = self.__getUserAndHost(userAndHost)[1]
                hostPassword = hostsPasswords[userAndHost]
                # Create ssh session to remote host
                log.info("Connecting to: " + userName + "@" + hostName)
                self.ssh.connect(hostName, SSH_PORT, userName, hostPassword)
                self.sftp = self.ssh.open_sftp()
                for resultFilePath in resultFilePaths:
                    filePath = self.__getLocationAndFilePath(resultFilePath)[1]
                    log.info("Checking: " + filePath)
                    try:
                        self.sftp.lstat(filePath)                    
                    except IOError:
                        returnFilePaths.append(resultFilePath)
                        log.info("Check fail!!")
                self.ssh.close()
                self.sftp.close()  
        return returnFilePaths
    
    def checkOneHostFiles(self, 
                   filePaths,
                   hostName,
                   userName,                   
                   hostsPassword,
                   gatewayHosts=None, 
                   numberTrials=1, 
                   forceOperation = False,  
                   operationId=1):
        """
        Check if file paths exists.
        filepaths -- List of file paths to check with this format: ["absolute_file_path1","absolute_file_path2"]
        gatewayHosts -- Gateway hosts dictionary with this format: "userName1@hostName1:userName2@hostName2":"userName@hostName"
        hostsPasswords -- Passwords needed to connect to involved hosts with this format: "userName@hostName":"hostPassword"
        numberTrials -- Number of trials in error cases.
        forceOperation -- Flag to indicate if, when an error happens and number of trials is exceeded, the operation must continue with the rest of files.
        operationId -- Operation identifier.
        returns -- List of not located file paths.
        """
        returnFilePaths = []
        print("**************************************** CHECKING******************************************")
        log.info("Connecting to: " + userName + "@" + hostName)
        self.ssh.connect(hostName, SSH_PORT, userName, hostsPassword)
        self.sftp = self.ssh.open_sftp()
        
        isLocalHost = self.__isLocalHost(hostName)
        
        for fileName in filePaths:
            log.info("Checking: " + fileName)
            if (isLocalHost):            
                if (len (missingPaths(fileName)) != 0):
                    returnFilePaths.append(fileName)
                    log.info("Check fail!!")
            else:
                try:
                    self.sftp.lstat(fileName)                    
                except IOError:
                    returnFilePaths.append(fileName)
                    log.info("Check fail!!")
        self.ssh.close()
        self.sftp.close()
        return returnFilePaths
    
    def __classifyFilePaths(self, filePaths):
        """
        Classify file paths to send depending on the source and target user/host involved credentials.
        filePaths -- Files dictionary with this format: "userName@hostName:absolute_file_path": [userName1@hostName1:absolute_file_path1, userName2@hostName2:absolute_file_path2]
        Key is the source file path and value the target file paths.
        returns -- Dictionary with file paths ordered by the two user/host credentials that are involved. The dictionary will have this structure:
        {"userName1@hostName1:userName2@hostName2" : {"filePath1":"filePath2"}
        """
        result = {}
        
        for sourcePath, targetPaths in filePaths.iteritems():
            for targetPath in targetPaths:
                sourceParts = self.__getLocationAndFilePath(sourcePath)
                targetParts = self.__getLocationAndFilePath(targetPath) 
                sourceUserAndHost = sourceParts[0]
                targetUserAndHost = targetParts[0]
                sourceFilePath = sourceParts[1]
                targetFilePath = targetParts[1]
                resultKey = sourceUserAndHost + PAIRS_SEPARATOR + targetUserAndHost
                auxDict = {}
                if (resultKey in result):     
                    auxDict = result[resultKey]
                    if sourceFilePath in auxDict:
                        # This will not happen in Scipion because one source path is not going
                        # to be sent to different file paths in the same machine.
                        msg = 'File ' + sourcePath + ' can not be sent to ' + auxDict[sourcePath] + " and to " + targetPath + ' because they are in the same machine.'
                        raise Exception(msg)
                auxDict[sourceFilePath] = targetFilePath
                result[resultKey] = auxDict
        
        return result
            
    def __sendFilesBetweenPairs(self, 
                                userAndHostPairs, 
                                filePaths,
                                hostsPasswords, 
                                gatewayHosts = None, 
                                numberTrials = 1,                        
                                forceOperation = False,
                                operationId = 1):  
        """
        filePaths -- Dictionary with this structure: {"filePath1":"filePath2"}
        gatewayHosts -- Gateway hosts List with this format: ["userName1@hostName1","userName2@hostName2"] 
        TODO: Get the target directories to check it existence and create them only once.
        """     
        pairParts = self.__getUserAndHosts(userAndHostPairs)
        sourceCredentials = pairParts[0]
        targetCredentials = pairParts[1]
        # We see what type of sending operation is.        
        if self.__isLocalCredential(sourceCredentials):
            if self.__isLocalCredential(targetCredentials):
                self.copyFiles(filePaths, numberTrials, forceOperation, operationId)
            else:
                targetUserAndHost = self.__getUserAndHost(targetCredentials)
                self.transferFilesTo(filePaths, targetUserAndHost[1], targetUserAndHost[0], hostsPasswords[targetCredentials], gatewayHosts, numberTrials, forceOperation, operationId)
        else:
            sourceUserAndHost = self.__getUserAndHost(sourceCredentials)
            if self.__isLocalCredential(targetCredentials):
                self.transferFilesFrom(filePaths, sourceUserAndHost[1], sourceUserAndHost[0], hostsPasswords[sourceCredentials], gatewayHosts, numberTrials, forceOperation, operationId)
            else:
                pass            
            
    def __getUserAndHosts(self, userAndHostPairs):
        """
        Separate user and host pair in their individuals.
        userAndHostPairs -- User and host pair: "userName1@hostName1:userName2@hostName2"
        returns -- Spplited pairs: ["userName1@hostName1", "userName2@hostName2"]
        """
        return userAndHostPairs.split(PAIRS_SEPARATOR)
    
    def __getUserAndHost(self, userAndHost):
        """
        Function to get the user and the host name from 'userName@hostName' string
        Returns: Tuple with ('userName', 'hostName')
        """
        return userAndHost.split("@")
    
    def __getLocationAndFilePath(self, locationAndFile):
        """
        Function to get the user and the userName@hostName and file path from 'userName@hostname:filePath' string
        Returns: Tuple with ('userName@hostName', 'filePath')
        """
        if (":" in locationAndFile):
            auxLocationAndFile = locationAndFile.split(":", 1)
            if (self.__isLocalCredential(auxLocationAndFile[0])):
                auxLocationAndFile[0] = LOCAL_USER_AND_HOST # Ease classification and other operations
        else:
            auxLocationAndFile = []
            auxLocationAndFile.append(LOCAL_USER_AND_HOST)
            auxLocationAndFile.append(locationAndFile)
        
        return auxLocationAndFile
    
    def __isLocalHost(self, hostName):
        """
        Checks if one host name is local machine name.
        No name or empty name means local host.
        """
        if (hostName is None or
            hostName == LOCAL_USER_AND_HOST):
            return True
        elif (socket.gethostname() == hostName):
            return True
        else:
            return False

    def __isLocalCredential(self, userAndHost):
        """
        Checks if one userName@hostName credential is about local machine.
        """
        if (userAndHost is None or
            userAndHost == LOCAL_USER_AND_HOST):
            return True
        else:
            hostName = self.__getUserAndHost(userAndHost)[1]
            return self.__isLocalHost(hostName);
        
    def __getDiferentHostsFromFilePaths (self, filePaths):
        """
        Gets different userName@hostName credentials of files list.
        files -- File list with this format: userName1@hostName1:absolute_file_path1
        returns -- List with different files userName@hostName credentials.
        """
        diferentHosts = []
        for filePath in filePaths:
            userAndHost = self.__getLocationAndFilePath(filePath)[0]
            if (self.__isLocalCredential(userAndHost)):
                userAndHost = LOCAL_USER_AND_HOST; # To simplify code for local credentials.
            if (userAndHost not in diferentHosts):
                diferentHosts.append(userAndHost)
        return diferentHosts;
    
    def __getFilePaths(self, filePaths, userAndHost):
        """
        Gets a file path from a list of target file paths for the given user and host.
        filePaths -- List of files with this format: "useName@hostName:/filePath"
        userAndHost -- String with this format. "userName@hostName" 
        returns -- Target file path list ("useName@hostName:/filePath") for the given user and host.
        """  
        resultFilePaths = []             
        for filePath in filePaths:
            locationAndFilePath = self.__getLocationAndFilePath(filePath)
            if (locationAndFilePath[0] == userAndHost):
                resultFilePaths.append(filePath);
        return resultFilePaths
    
    def __copyLocalFile(self, sourceFilePath, targetFilePath):
        """
        Send local file to remote machine
        sourceFilePath -- Source file path (/file path/...).
        targetFilePath -- Target file path (/file path/...).
        """
        log.info("Copying " + sourceFilePath + " to " + targetFilePath)
        # Check if file already existsFilePath and it is up to date
        existsFilePath = False
        if (exists(targetFilePath)):
            if ( self.__getLocalSHA1(sourceFilePath) ==  self.__getLocalSHA1(targetFilePath)):
                existsFilePath = True
                log.info(targetFilePath + " already existed")
        if (not existsFilePath):                    
            makeFilePath(targetFilePath)
            shutil.copy2(sourceFilePath, targetFilePath)
        
    def __sendLocalFile(self, sourceFilePath, targetFilePath, gatewayHosts, sftp):
        """
        Send local file to remote machine
        sourceFilePath -- Source file path (/file path/...).
        targetFilePath -- Target file path (/file path/...).
        sftp -- sftp connection.
        """
        log.info("Sending " + sourceFilePath + " to " + targetFilePath)
        # Check if file already existsFilePath and it is up to date
        existsFilePath = False
        if (self.__existsRemotePath(targetFilePath, sftp)):
            if ( self.__getLocalSHA1(sourceFilePath) ==  self.__getRemoteSHA1(targetFilePath, self.ssh)):
                existsFilePath = True
                log.info(targetFilePath + " already existed")
        if (not existsFilePath):        
            self.__createRemoteFolderForFile(targetFilePath, sftp)
            try:
                sftp.put(sourceFilePath, targetFilePath)
            except IOError as err:
                log.error("Fail sending local file " + sourceFilePath + " to remote file " + targetFilePath + " - " + str(err))
                raise
        
        
    def __getRemoteFile(self, sourceFilePath, targetFilePath, gatewayHosts, sftp):
        """
        Send local file to remote machine
        sourceFilePath -- Source file path (/file path/...).
        targetFilePath -- Target file path (/file path/...).
        sftp -- sftp connection.
        """
        log.info("Getting " + sourceFilePath + " to " + targetFilePath)
        # Check if file already existsFilePath and it is up to date
        existsFilePath = False
        if (exists(targetFilePath)):
            if ( self.__getRemoteSHA1(sourceFilePath, self.ssh) ==  self.__getLocalSHA1(targetFilePath)):
                existsFilePath = True
                log.info(targetFilePath + " already existed")
        if (not existsFilePath):        
            makeFilePath(targetFilePath)
            try:
                sftp.get(sourceFilePath, targetFilePath)
            except IOError as err:
                log.error("Fail getting remote file " + sourceFilePath + " to local file " + targetFilePath + " - " + str(err))
                raise
    
    def __createRemoteFolderForFile(self, filePath, sftp):
        """
        Create folder for file in remote host defined by sfpt.
        filePath -- File path which parent directory we must create (/file path/...).
        sftp -- Remote sftp session.
        """
        filePathParentDirectory = os.path.dirname(filePath)
        self.__mkdirP(filePathParentDirectory, sftp)

    def __mkdirP(self, remoteDirectory, sftp):
        """
        Create remote folder structure creating all non-existent folders.
        remoteDirectory -- Remote directory to create.
        sftp -- Remote sftp session.
        """
        if not self.__existsRemotePath(remoteDirectory, sftp):
            self.__mkdirP(os.path.dirname(remoteDirectory), sftp)
            sftp.mkdir(remoteDirectory) 
            
    def __getLocalSHA1(self, filePath):
        return hashlib.sha1(file(filePath, 'r').read()).hexdigest()

    def __getRemoteSHA1(self, filePath, ssh):
        stdin, stdout, stderr = ssh.exec_command("sha1sum '" + filePath + "'")
        return stdout.readlines()[0].split()[0] 
            
    def __existsRemotePath(self, path, sftp):
        try:
            sftp.lstat(path)
            return True            
        except IOError:
            return False
    
################################################################

#                AUXILIARY FUNCTIONS                           #

################################################################        

def isRemoteDir(sftp, path):
    """ Check if one remote directory exists
    Params:
        sftp: Sftp session.
        path: Path to check.
    Returns: True if the given path is a directory, false if it is not a directory.
    """
    try:
        sftp.chdir(path)
        return True
    except (IOError, paramiko.SFTPError):
        return False

def getRemoteFolderFiles(hostName, userName, password, folderPath, recursive=True):
    """ Recover all files in the given folder.
    Params:
        hostName: Remote host name.
        userName: User name.
        password: Password.
        folderPath: Folder to get files.
        recursive: if True go recursively inside other subfolders.
    Returns: List of files.
    """    
    # Default ssh session options.
    ssh = paramiko.SSHClient()
    ssh.load_system_host_keys()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(hostName, 22, userName, password)
    sftp = ssh.open_sftp()
    remoteFiles = getRemoteFiles(sftp, folderPath, recursive)
    sftp.close()
    ssh.close()
    return remoteFiles
    
def getRemoteFiles(sftp, folderPath, recursive=True):
    """ Recover all files in the given folder.
    Params:
        sftp: Sftp session.
        folderPath: Folder to get files.
        recursive: if True go recursively inside other subfolders.
    Returns: List of files.
    """    
    filePathList = sftp.listdir(folderPath)
    resultFilePathList = []
    for filePath in filePathList:
        if isRemoteDir(sftp, join(folderPath, filePath)) and recursive:
            resultFilePathList += getRemoteFiles(sftp, join(folderPath, filePath))
        else:
            resultFilePathList.append(join(folderPath, filePath))
    return resultFilePathList

def removeRemoteFolder (hostName, userName, password, folderPath):
    """ Removes a remote folder and all it content.
    Params:
        hostName: Remote host name.
        userName: User name.
        password: Password.
        folderPath: Folder to delete.
    Returns: 
        Tuple with standard input, standard output and error output.
    """
    ssh = paramiko.SSHClient()
    ssh.load_system_host_keys()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(hostName, 22, userName, password)    
    return removeTree(ssh, folderPath)


def removeTree(ssh, folderPath):
    """ Removes a remote folder and all it content.
    Params:
        ssh: Ssh session.
        folderPath: Folder to delete.
    Returns: 
        Tuple with standard input, standard output and error output.
    """
    stdin, stdout, stderr = ssh.exec_command("rm -dfr " + folderPath)
    ssh.close()
    return stdin, stdout, stderr


def getFilePathList(filePaths):
    """
    Get target file paths list from target file path dictionary.
    filePaths -- Files dictionary with this format: "userName@hostName:absolute_file_path": ["userName1@hostName1:absolute_file_path1", "userName2@hostName2:absolute_file_path2"]
    returns -- List of target files with this format: "useName@hostName:/filePath".
    """
    resultFilePathList = []
    for filePathList in filePaths.values():
        for filePath in filePathList:
            resultFilePathList.append(filePath)
    return resultFilePathList


if __name__ == '__main__':
    pass
