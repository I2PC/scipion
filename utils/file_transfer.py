'''
Created on Apr 8, 2013

@author: antonio
'''
import socket
import paramiko
import os
import shutil
import hashlib
from pyworkflow.utils.path import *
from pyworkflow.utils.log import *

LOCAL_USER_AND_HOST = ''
SSH_PORT = 22;
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
                      operationId = 1, 
                      numberTrials = 1,                        
                      forceOperation = False):
        """
        filePaths -- Files dictionary with this format: "userName@hostName:absolute_file_path": ["userName1@hostName1:absolute_file_path1", "userName2@hostName2:absolute_file_path2"]
        Key is the source file path and value the target file paths.
        """
        classifiedFiles = self.__classifyFilePaths(filePaths)
        #print ('Classified files to transfer: ' + str(classifiedFiles))
        
        for userAndHostPairs, cursorFilePaths in classifiedFiles.iteritems():
            self.__sendFilesBetweenPairs(userAndHostPairs, 
                                         cursorFilePaths,
                                         hostsPasswords, 
                                         gatewayHosts,
                                         operationId, 
                                         numberTrials,                        
                                         forceOperation)
        
        
    def close(self):
        if (self.sftp is not None):
            self.sftp.close()
        if (self.ssh is not None):
            self.ssh.close()

    
                
        
    def deleteFiles(self, 
                    filePaths,                    
                    hostsPasswords, 
                    gatewayHosts=None, 
                    operationId=1, 
                    numberTrials=1, 
                    forceOperation=False):
        """
        Delete a list of file paths.
        filepaths -- List of file paths to remove with this format: "userName@hostName:absolute_file_path"
        gatewayHosts -- Gateway hosts dictionary with this format: "userName1@hostName1:userName2@hostName2":"userName@hostName"
        hostsPasswords -- Passwords needed to connect to involved hosts with this format: "userName@hostName":"hostPassword"
        operationId -- Operation identifier.
        numberTrials -- Number of trials in error cases.
        forceOperation -- Flag to indicate if, when an error happens and number of trials is exceeded, the operation must continue with the rest of files.
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
            else:
                pass
            
    def deleteDirectories(self, 
                          directoryPaths,                    
                          hostsPasswords, 
                          gatewayHosts=None, 
                          operationId=1, 
                          numberTrials=1, 
                          forceOperation=False):
        """
        Delete a list of directory paths.
        directoryPaths -- List of file paths to remove with this format: "userName@hostName:absolute_directory_path"
        gatewayHosts -- Gateway hosts dictionary with this format: "userName1@hostName1:userName2@hostName2":"userName@hostName"
        hostsPasswords -- Passwords needed to connect to involved hosts with this format: "userName@hostName":"hostPassword"
        operationId -- Operation identifier.
        numberTrials -- Number of trials in error cases.
        forceOperation -- Flag to indicate if, when an error happens and number of trials is exceeded, the operation must continue with the rest of files.
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
            else:
                pass
        
    def checkFiles(self, 
                   filePaths,                   
                   hostsPasswords,
                   gatewayHosts=None,  
                   operationId=1, 
                   numberTrials=1, 
                   forceOperation = False):
        """
        Check if file paths exists.
        filepaths -- List of file paths to remove with this format: "userName@hostName:absolute_file_path"
        gatewayHosts -- Gateway hosts dictionary with this format: "userName1@hostName1:userName2@hostName2":"userName@hostName"
        hostsPasswords -- Passwords needed to connect to involved hosts with this format: "userName@hostName":"hostPassword"
        operationId -- Operation identifier.
        numberTrials -- Number of trials in error cases.
        forceOperation -- Flag to indicate if, when an error happens and number of trials is exceeded, the operation must continue with the rest of files.
        returns -- List of not located file paths.
        """
        returnFilePaths = []
        # As we are going to create a session for each target host we must get different target hosts.
        userAndHosts = self.__getDiferentHostsFromFilePaths(filePaths)
        for userAndHost in userAndHosts:
            if (not self.__isLocalCredential(userAndHost)):
                resultFilePaths = self.__getFilePaths(filePaths, userAndHost)
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
                        log.info("Check fail!!")
                        returnFilePaths.append(resultFilePath)
            else:
                pass
        return returnFilePaths
    
    def __classifyFilePaths(self, filePaths):
        """
        Classify file paths to send depending on the source and target user/host involved credentials.
        filePaths -- Files dictionary with this format: "userName@hostName:absolute_file_path": [userName1@hostName1:absolute_file_path1, userName2@hostName2:absolute_file_path2]
        Key is the source file path and value the target file paths.
        returns -- Dictionary with file paths ordered by the two user/host credentials that are involved. The dictionary will have this structure:
        {"userName1@hostName1:userName2@hostName2" : {"userName@hostName:filePath":["userName3@hostName3:filePath3", "userName4@hostName4:filePath4"]}
        """
        result = {}
        
        for sourcePath, targetPaths in filePaths.iteritems():
            for targetPath in targetPaths:
                sourceParts = self.__getLocationAndFilePath(sourcePath)
                targetParts = self.__getLocationAndFilePath(targetPath) 
                sourceUserAndHost = sourceParts[0]
                targetUserAndHost = targetParts[0]
                tempKey = sourceUserAndHost + PAIRS_SEPARATOR + targetUserAndHost
                resultKey = self.__getKey(result, tempKey)
                auxDict = {}
                auxList = []
                if (resultKey is not None):     
                    auxDict = result[resultKey]
                    if sourcePath in auxDict:
                        # This will not happen in Scipion because one source path is not going
                        # to be sent to different file paths in the same machine.
                        auxList = auxDict[sourcePath]
                else:
                    resultKey = tempKey
                auxList.append(targetPath)
                auxDict[sourcePath] = auxList
                result[resultKey] = auxDict
        
        return result
                   
    def __equalUserAndHost(self, 
                         userAndHost,
                         filePath):
        """
        Check if file path has the given credentials.
        userAndHost -- String with "userName@hostName" format.
        """
        fileUserAndHost = self.__getLocationAndFilePath(filePath)
        if (fileUserAndHost == userAndHost):
            return True
        else:
            return False
        
    def __getKey(self,
               dictionary,
               userAndHostPair):
        """
        Get key for user and host credential in dictionary.
        dictionary -- Dictionary where make checking operation.
        userAndHostPair -- User and host pair credential to check with "userName1@hostName1:userName2@hostName2" format.
        returns -- Key for userAndHost in dictionary if "userName1@hostName1:userName2@hostName2" or "userName2@hostName2:userName1@hostName1" are
        in given dictionary. None in other cases.
        """
        if userAndHostPair in dictionary:
            return userAndHostPair
        else:
            sppliter = self.__getUserAndHosts(userAndHostPair)
            newUserAndHost = (sppliter[1] + ":" + sppliter[0])
            if newUserAndHost in dictionary:
                return newUserAndHost
            else:
                return None
        
    def __equalUserAndHostPairs(self, userAndHostPair1, userAndHostPair2):
        """
        Check if two credentials pairs are equals.
        userAndHostPair1: Host credential pair with this format: "userName1@hostName1:userName2@hostName2"
        userAndHostPair2: Host credential pair with this format: "userName1@hostName1:userName2@hostName2"
        Two credentials pairs are equals if they are identical or if they are in reverse order: 
                                "userName1@hostName1:userName2@hostName2" == "userName1@hostName1:userName2@hostName2"
                                "userName1@hostName1:userName2@hostName2" == "userName2@hostName2:userName1@hostName1"
        """
        if userAndHostPair1 == userAndHostPair2:
            return True
        else:
            pairParts1 = self.__getUserAndHosts(userAndHostPair1)
            pairParts2 = self.__getUserAndHosts(userAndHostPair2)
            return (pairParts1[0] == pairParts2[1] and
                    pairParts1[1] == pairParts2[0])
            
    def __sendFilesBetweenPairs(self, 
                                userAndHostPairs, 
                                filePaths,
                                hostsPasswords, 
                                gatewayHosts = None,
                                operationId = 1, 
                                numberTrials = 1,                        
                                forceOperation = False):  
        """
        TODO: Get the target directories to check it existence and create them only once.
        """     
        pairParts = self.__getUserAndHosts(userAndHostPairs)
        userAndHost1 = pairParts[0]
        userAndHost2 = pairParts[1]
        # We see what type of sending operation is.
        localAndLocal = False
        localAndRemote = False
        remoteAndRemote = False
        remoteUserAndHost = None        
        if self.__isLocalCredential(userAndHost1):
            if self.__isLocalCredential(userAndHost2):
                localAndLocal = True
            else:
                localAndRemote = True
                remoteUserAndHost = userAndHost2
        else:
            if self.__isLocalCredential(userAndHost2):
                localAndRemote = True
                remoteUserAndHost = userAndHost1
            else:
                remoteAndRemote = True
        
        # If the operation involves local and remote machine we create ssh session to remote host
        if (localAndRemote):
            remoteUserAndHostParts = self.__getUserAndHost(remoteUserAndHost)
            log.info("Connecting to: " + remoteUserAndHostParts[0] + "@" + remoteUserAndHostParts[1])
            self.ssh.connect(remoteUserAndHostParts[1], SSH_PORT, remoteUserAndHostParts[0], hostsPasswords[remoteUserAndHost])
            self.sftp = self.ssh.open_sftp()
                
                
        for sourcePath, targetPaths in filePaths.iteritems():
            sourceParts = self.__getLocationAndFilePath(sourcePath)
            sourceUserAndHost = sourceParts[0]
            sourceFilePath = sourceParts[1]
            for targetPath in targetPaths:
                targetParts = self.__getLocationAndFilePath(targetPath)
                targetUserAndHost = targetParts[0]
                targetFilePath = targetParts[1]
                if localAndLocal:
                    self.__copyLocalFile(sourceFilePath, targetFilePath)
                elif localAndRemote:                     
                    # Check if source or target is remote
                    if self.__isLocalCredential(sourceUserAndHost):
                        # The source path is local
                        self.__sendLocalFile(sourceFilePath, targetFilePath, gatewayHosts, self.sftp)                        
                    else:
                        self.__getRemoteFile(sourceFilePath, targetFilePath, self.sftp)
                elif remoteAndRemote:
                    pass
                else:
                    raise Exception("There was a problem with the sending type")          
            
    def __getUserAndHosts(self, userAndHostPairs):
        """
        Separate user and host pair in their individuals.
        userAndHostPairs -- User and host pair: "userName1@hostName1:userName2@hostName2"
        returns -- Spplited pairs: ["userName1@hostName1", "userName2@hostName2"]
        """
        return userAndHostPairs.split(':')
    
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
            createFolderForFile(targetFilePath)
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
        
        
    def __getRemoteFile(self, sourceFilePath, targetFilePath, sftp):
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
            createFolderForFile(targetFilePath)
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
        stdin, stdout, stderr = ssh.exec_command("sha1sum " + filePath)
        return stdout.readlines()[0].split()[0]
    
    def __isRemoteDir(self, sftp, path):
        try:
            sftp.chdir(path)
            return True
        except IOError:
            return False 
    
    def __existsRemotePath(self, path, sftp):
        try:
            sftp.lstat(path)
            return True            
        except IOError:
            return False
    
################################################################

#                AUXILIARY FUNCTIONS                           #

################################################################        


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
