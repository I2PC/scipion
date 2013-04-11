'''
Created on Apr 8, 2013

@author: antonio
'''
import socket
import paramiko
import os
import shutil

LOCAL_USER_AND_HOST = ""
SSH_PORT = "22";


class FileTransferProcess():
    
    ssh = None   
    sftp = None 
    
    def __init__(self):
        # Default ssh session options.
        ssh = paramiko.SSHClient()
        ssh.load_system_host_keys()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    def transferFiles(self,
                      sourceFiles, 
                      targetFiles,
                      gatewayHosts = None, 
                      hostsPaswords, 
                      operationId, 
                      numberTrials, 
                      forceOperation = False):
        # Order source files by it host.
        filesByHosts = getSourceFilesByHosts(sourceFiles)
        # First we send local files.
        localFiles = filesByHosts[LOCAL_USER_AND_HOST]    
        if (localFiles is not None
            and len(localFiles) > 0):
            self.trasferLocalFiles(localFiles, getTargetFilesById(localFiles.keys(), targetFiles), gatewayHosts, 
                              hostsPaswords, operationId, numberTrials, forceOperation)
            
        self.sftp.close()
        self.ssh.close()
        
        
    def __copyLocalFiles(self, 
                         sourceFiles, 
                         targetFiles, 
                         operationId, 
                         numberTrials, 
                         forceOperation=False):
        
        """
        sourceFiles -- Local source files dictionary with this format: "id":"filePath"
        targetFiles -- Target source files dictionary with this format: "id":"filePath"
        operationId -- Operation identifier.
        numberTrials -- Number of trials in error cases.
        forceOperation -- Flag to indicate if, when an error happens and number of trials is exceeded, the operation must continue with the rest of files.
        """
        for sourceFileId, sourceFile in sourceFiles.iterItems():
            targetFile = targetFiles[sourceFileId]
            copyLocalFile(sourceFile, targetFile)
            
    
    def __trasferLocalFilesToRemoteHost(self, 
                                        sourceFiles, 
                                        targetFiles,
                                        gatewayHosts, 
                                        hostsPaswords, 
                                        operationId, 
                                        numberTrials, 
                                        forceOperation=False):
        """
        sourceFiles -- Local source files dictionary with this format: "id":"filePath"
        targetFiles -- Target files dictionary with this format: "id": [userName1@hostName1:absolute_file_path1, userName2@hostName2:absolute_file_path2]
        gatewayHosts -- Gateway hosts dictionary with this format: "userName1@hostName1:userName2@hostName2":"userName@hostName"
        hostsPasswords -- Passwords needed to connect to involved hosts with this format: "userName@hostName":"hostPassword"
        operationId -- Operation identifier.
        numberTrials -- Number of trials in error cases.
        forceOperation -- Flag to indicate if, when an error happens and number of trials is exceeded, the operation must continue with the rest of files.
        """
        # As we are going to create a session for each target host we must get different target hosts.
        targetUserAndHosts = getDiferentTargetHosts(targetFiles)
        for targetUserAndHost in targetUserAndHosts:
            # Recover target credentials to send files
            targetUserName = getUserAndHost(targetUserAndHost)[0]
            targetHostName = getUserAndHost(targetUserAndHost)[1]
            targetPassword = hostsPaswords[targetUserAndHost]
            # Create ssh session to remote host
            self.ssh.connect(targetHostName, SSH_PORT, targetUserName, targetPassword)
            self.sftp = self.ssh.open_sftp()
            for sourceFileId, sourceFile in sourceFiles.iteritems():
                targetFilePaths = targetFiles[sourceFileId]
                #One file is not going to be copy to  more than one target file path in the same host.
                resultTargetFilePaths = getFilePaths(targetFilePaths, targetUserAndHost)
                if (resultTargetFilePaths is not None and len(resultTargetFilePaths) != 0):
                    targetFilePath = resultTargetFilePaths[0]
                    # Send local file
                    targetFile = getLocationAndFilePath(targetFilePath)[1]
                    sendLocalFile(sourceFile, targetFile, gatewayHosts, self.sftp)
                    
    def __trasferRemoteFilesToLocalHost(self, 
                          sourceFiles, 
                          targetFiles,
                          gatewayHosts, 
                          hostsPaswords, 
                          operationId, 
                          numberTrials, 
                          forceOperation=False):
        """
        sourceFiles -- Remote source files dictionary with this format: "id":"sourceFilePath"
        targetFiles -- Target files dictionary with this format: "id":"targetFilePath"
        gatewayHosts -- Gateway hosts dictionary with this format: "userName1@hostName1:userName2@hostName2":"userName@hostName"
        hostsPasswords -- Passwords needed to connect to involved hosts with this format: "userName@hostName":"hostPassword"
        operationId -- Operation identifier.
        numberTrials -- Number of trials in error cases.
        forceOperation -- Flag to indicate if, when an error happens and number of trials is exceeded, the operation must continue with the rest of files.
        """
        # As we are going to create a session for each target host we must get different target hosts.
        sourceUserAndHosts = getDiferentHostsFromFilePaths(sourceFiles)
        for sourceUserAndHost in sourceUserAndHosts:
            sourceUserName = getUserAndHost(sourceUserAndHost)[0]
            sourceHostName = getUserAndHost(sourceUserAndHost)[1]
            sourcePassword = hostsPaswords[sourceUserAndHost]
            # Create ssh session to remote host
            self.ssh.connect(sourceHostName, SSH_PORT, sourceUserName, sourcePassword)
            self.sftp = self.ssh.open_sftp()
            for sourceFileId, sourceFilePath in sourceFiles.iteritems():
                targetFilePath = targetFiles[sourceFileId]
                sourceFile = getLocationAndFilePath(sourceFilePath)[1]
                getRemoteFile(sourceFile, targetFilePath, self.sftp)
                
        
    def deleteFiles(self, 
                    filePaths,
                    gatewayHosts = None, 
                    hostsPaswords, 
                    operationId, 
                    numberTrials, 
                    forceOperation = False):
        """
        filepaths -- List of file paths to remove with this format: "userName@hostName:absolute_file_path"
        gatewayHosts -- Gateway hosts dictionary with this format: "userName1@hostName1:userName2@hostName2":"userName@hostName"
        hostsPasswords -- Passwords needed to connect to involved hosts with this format: "userName@hostName":"hostPassword"
        operationId -- Operation identifier.
        numberTrials -- Number of trials in error cases.
        forceOperation -- Flag to indicate if, when an error happens and number of trials is exceeded, the operation must continue with the rest of files.
        """
        
        # As we are going to create a session for each target host we must get different target hosts.
        userAndHosts = getDiferentHostsFromFilePaths(filePaths)
        for userAndHost in userAndHosts:
            resultFilePaths = getFilePaths(filePaths, userAndHost)
            # Recover host credentials to remove files
            userName = getUserAndHost(userAndHost)[0]
            hostName = getUserAndHost(userAndHost)[1]
            hostPassword = hostsPaswords[userAndHost]
            # Create ssh session to remote host
            self.ssh.connect(hostName, SSH_PORT, userName, hostPassword)
            self.sftp = self.ssh.open_sftp()
            for resultFilePath in resultFilePaths:
                self.sftp.remove(getLocationAndFilePath(resultFilePath)[1])
        
    def checkFiles(self, 
                   filePaths,
                   gatewayHosts = None, 
                   hostsPaswords, 
                   operationId, 
                   numberTrials, 
                   forceOperation = False):
        """
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
        userAndHosts = getDiferentHostsFromFilePaths(filePaths)
        for userAndHost in userAndHosts:
            resultFilePaths = getFilePaths(filePaths, userAndHost)
            # Recover host credentials to remove files
            userName = getUserAndHost(userAndHost)[0]
            hostName = getUserAndHost(userAndHost)[1]
            hostPassword = hostsPaswords[userAndHost]
            # Create ssh session to remote host
            self.ssh.connect(hostName, SSH_PORT, userName, hostPassword)
            self.sftp = self.ssh.open_sftp()
            for resultFilePath in resultFilePaths:
                try:
                    self.sftp.lstat(getLocationAndFilePath(resultFilePath)[1])                    
                except IOError:
                    returnFilePaths.append(resultFilePath)
        return returnFilePaths
    
    def __trasferRemoteFiles(self,
                           sourceFiles, 
                           targetFiles,
                           gatewayHosts, 
                           hostsPaswords, 
                           operationId, 
                           numberTrials, 
                           forceOperation=False):
        pass


################################################################

#                AUXILIARY FUNCTIONS                           #

################################################################

def getUserAndHost(userAndHost):
    """
    Function to get the user and the host name from 'userName@hostName' string
    Returns: Tuple with ('userName', 'hostName')
    """
    return userAndHost.split("@")

def getLocationAndFilePath(localtionAndFile):
    """
    Function to get the user and the userName@hostName and file path from 'userName@hostname:filePath' string
    Returns: Tuple with ('userName@hostName', 'filePath')
    """
    return localtionAndFile.split(":", 1)

def isLocalHost(hostName):
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

def isLocalCredential(userAndHost):
    """
    Checks if one userName@hostName credential is about local machine.
    """
    if (userAndHost is None or
        userAndHost == LOCAL_USER_AND_HOST):
        return True
    else:
        hostName = getUserAndHost(userAndHost)[1]
        return isLocalHost(hostName);
        

def getSourceFilesByHosts(files):
    """    
    This function order files by it userName@hostName credentials insert them in a new dictionary 
    where data is stored with this structure: "userName@hostName" : {"1": "filepath1", "2" : "filePath2"}
    files -- Dictionary with "id":"useName@hostName:filePath" structure.
    returns -- Dictionary with files ordered by it host credentials.
    """
    filesByHosts = {}
    registeredFiles = {}
    for fileId, fileCredentials in files:
        fileCredentialsTuple = getLocationAndFilePath(fileCredentials)
        userAndHost = fileCredentialsTuple[0]
        filePath = fileCredentialsTuple[1]
        if (isLocalCredential(userAndHost)):
            userAndHost = LOCAL_USER_AND_HOST; # To simplify code for local credentials.
        if (userAndHost in filesByHosts):
            registeredFiles = filesByHosts[userAndHost]
            registeredFiles[fileId] = filePath
        else:
            registeredFiles = {fileId : filePath} 
        filesByHosts[userAndHost] = registeredFiles;
        
    return filesByHosts; 

def getTargetFilesById(sourceFilesIds, targetFiles):
    """
    This functon gets all target files from a target files dictionary that match with a source file id list.
    sourceFilesIds -- List of file ids.
    targetFiles -- Target files dictionary.
    returns -- Dictionary with the target files that match with the given source file id list. 
    """
    resultTargetFiles = {}
    for fileId in sourceFilesIds :
        resultTargetFiles[fileId, targetFiles[fileId]]
    return resultTargetFiles;
    
def getDiferentTargetHosts (targetFiles):
    """
    This function gets different userName@hostName credentials of target files dictionary.
    targetFiles -- Target files dictionary with this format: "id": [userName1@hostName1:absolute_file_path1, userName2@hostName2:absolute_file_path2]
    returns -- List with different target userName@hostName credentials.
    """
    diferentTargetHosts = []
    for targetFileList in targetFiles.values():
        for targetFile in targetFileList:
            userAndHost = getLocationAndFilePath(targetFile)
            if (isLocalCredential(userAndHost)):
                userAndHost = LOCAL_USER_AND_HOST; # To simplify code for local credentials.
            if (userAndHost not in diferentTargetHosts):
                diferentTargetHosts.append(userAndHost)
    return diferentTargetHosts;

def getDiferentHostsFromFilePaths (filePaths):
    """
    Gets different userName@hostName credentials of files list.
    files -- File list with this format: userName1@hostName1:absolute_file_path1
    returns -- List with different files userName@hostName credentials.
    """
    diferentHosts = []
    for filePath in filePaths:
            userAndHost = getLocationAndFilePath(filePath)
            if (isLocalCredential(userAndHost)):
                userAndHost = LOCAL_USER_AND_HOST; # To simplify code for local credentials.
            if (userAndHost not in diferentHosts):
                diferentHosts.append(userAndHost)
    return diferentHosts;
    
 
def getFilePaths(filePaths, userAndHost):
    """
    Gets a file path from a list of target file paths for the given user and host.
    filePaths -- List of files with this format: "useName@hostName:/filePath"
    userAndHost -- String with this format. "userName@hostName" 
    returns -- Target file path list ("useName@hostName:/filePath") for the given user and host.
    """  
    resultFilePaths = []             
    for filePath in filePaths:
        locationAndFilePath = getLocationAndFilePath(filePath)
        if (locationAndFilePath[0] == userAndHost):
            resultFilePaths.append(filePath);
    return resultFilePaths

def copyLocalFile(sourceFile, targetFile):
    """
    Send local file to remote machine
    sourceFile -- Source file path (/file path/...).
    targetFile -- Target file path (/file path/...).
    """
    createLocalFolderForFile(targetFile)
    shutil.copy2(sourceFile, targetFile)

def sendLocalFile(sourceFile, targetFile, sftp):
    """
    Send local file to remote machine
    sourceFile -- Source file path (/file path/...).
    targetFile -- Target file path (/file path/...).
    ssh -- ssh session.
    """
    createRemoteFolderForFile(targetFile, sftp)
    sftp.put(sourceFile, targetFile)
    
def getRemoteFile(sourceFile, targetFile, sftp):
    """
    Send local file to remote machine
    sourceFile -- Source file path (/file path/...).
    targetFile -- Target file path (/file path/...).
    ssh -- ssh session.
    """
    createLocalFolderForFile(targetFile)
    sftp.get(sourceFile, targetFile)
                

def createRemoteFolderForFile(filePath, sftp):
    """
    Create folder for file in remote host defined by sfpt.
    filePath -- File path which parent directory we must create (/file path/...).
    sftp -- Remote sftp session.
    """
    filePathParentDirectory = os.path.dirname(filePath)
    # We check if this file path exist
    try:
        sftp.lstat(filePathParentDirectory)
    except IOError:
        # Directory does not exist
        sftp.mkdir(filePathParentDirectory)
        
def createLocalFolderForFile(filePath):
    """
    Create folder for file in local host.
    filePath -- File path which parent directory we must create (/file path/...).
    """
    filePathParentDirectory = os.path.dirname(filePath)
    # We check if this file path exist
    if not os.path.exists(filePathParentDirectory):
        os.makedirs(filePathParentDirectory)

def getTargetFilePathList(targetFilePaths):
    """
    Get target file paths list from target file path dictionary.
    targetFiles -- Target files dictionary with this format: "id": [userName1@hostName1:absolute_file_path1, userName2@hostName2:absolute_file_path2]
    returns -- List of target files with this format: "useName@hostName:/filePath".
    """
    filePaths = []
    for targetFilePath in targetFilePaths:
        for filePath in  targetFilePath:
            filePaths.append(filePath)
    return filePaths

def getLocalToLocalFiles(localSourceFilePaths, targetFiles):
    for localSourceFileId in localSourceFilePaths.keys():
        targetFilePaths = targetFiles[localSourceFileId]
        for targetFilePath in targetFilePaths:
            getLocationAndFilePath(targetFilePath)
            pass

if __name__ == '__main__':
    pass