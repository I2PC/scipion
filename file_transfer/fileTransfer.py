'''
Created on Apr 8, 2013

@author: antonio
'''
import socket
import paramiko
import os
import shutil

LOCAL_USER_AND_HOST = ''
SSH_PORT = '22';
PAIRS_SEPARATOR = ':'


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
                      hostsPaswords, 
                      gatewayHosts = None,
                      operationId = 1, 
                      numberTrials = 1,                        
                      forceOperation = False):
        """
        filePaths -- Files dictionary with this format: "userName@hostName:absolute_file_path": ["userName1@hostName1:absolute_file_path1", "userName2@hostName2:absolute_file_path2"]
        Key is the source file path and value the target file paths.
        """
        classifiedFiles = self.__classifyFilePaths(filePaths)
        
        for userAndHostPairs, cursorFilePaths in classifiedFiles.iteritems():
            self.__sendFilesBetweenPairs(userAndHostPairs, 
                                         cursorFilePaths,
                                         hostsPaswords, 
                                         gatewayHosts,
                                         operationId, 
                                         numberTrials,                        
                                         forceOperation)
        
        self.sftp.close()
        self.ssh.close()
    

    
                
        
    def deleteFiles(self, 
                    filePaths,                    
                    hostsPaswords, 
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
                hostPassword = hostsPaswords[userAndHost]
                # Create ssh session to remote host
                self.ssh.connect(hostName, SSH_PORT, userName, hostPassword)
                self.sftp = self.ssh.open_sftp()
                for resultFilePath in resultFilePaths:
                    self.sftp.remove(self.__getLocationAndFilePath(resultFilePath)[1])
            else:
                pass
            
    def deleteDirectories(self, 
                          directoryPaths,                    
                          hostsPaswords, 
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
                hostPassword = hostsPaswords[userAndHost]
                # Create ssh session to remote host
                self.ssh.connect(hostName, SSH_PORT, userName, hostPassword)
                self.sftp = self.ssh.open_sftp()
                for resultDirectoryPath in resultDirectoryPaths:
                    self.sftp.rmdir(self.__getLocationAndFilePath(resultDirectoryPath)[1])
            else:
                pass
        
    def checkFiles(self, 
                   filePaths,                   
                   hostsPaswords,
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
                hostPassword = hostsPaswords[userAndHost]
                # Create ssh session to remote host
                self.ssh.connect(hostName, SSH_PORT, userName, hostPassword)
                self.sftp = self.ssh.open_sftp()
                for resultFilePath in resultFilePaths:
                    try:
                        self.sftp.lstat(self.__getLocationAndFilePath(resultFilePath)[1])                    
                    except IOError:
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
                                hostsPaswords, 
                                gatewayHosts = None,
                                operationId = 1, 
                                numberTrials = 1,                        
                                forceOperation = False):       
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
            self.ssh.connect(remoteUserAndHostParts[1], SSH_PORT, remoteUserAndHostParts[0], hostsPaswords[remoteUserAndHost])
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
    
    def __copyLocalFile(self, sourceFile, targetFile):
        """
        Send local file to remote machine
        sourceFile -- Source file path (/file path/...).
        targetFile -- Target file path (/file path/...).
        """
        self.__createLocalFolderForFile(targetFile)
        shutil.copy2(sourceFile, targetFile)
        
    def __sendLocalFile(self, sourceFile, targetFile, gatewayHosts, sftp):
        """
        Send local file to remote machine
        sourceFile -- Source file path (/file path/...).
        targetFile -- Target file path (/file path/...).
        sftp -- sftp connection.
        """
        self.__createRemoteFolderForFile(targetFile, sftp)
        sftp.put(sourceFile, targetFile)
        
        
    def __getRemoteFile(self, sourceFile, targetFile, sftp):
        """
        Send local file to remote machine
        sourceFile -- Source file path (/file path/...).
        targetFile -- Target file path (/file path/...).
        sftp -- sftp connection.
        """
        self.__createLocalFolderForFile(targetFile)
        sftp.get(sourceFile, targetFile)
    
    def __createRemoteFolderForFile(self, filePath, sftp):
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
            
    def __createLocalFolderForFile(self, filePath):
        """
        Create folder for file in local host.
        filePath -- File path which parent directory we must create (/file path/...).
        """
        filePathParentDirectory = os.path.dirname(filePath)
        # We check if this file path exist
        if not os.path.exists(filePathParentDirectory):
            os.makedirs(filePathParentDirectory)
    
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