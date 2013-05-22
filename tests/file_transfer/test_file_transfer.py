"""
Created on Apr 9, 2013

@author: antonio
"""
import unittest
#from fileTransfer.fileTransfer import FileTransfer
from utils.file_transfer import *
import xml.etree.ElementTree as ET
import time

class TestFileTransfer(unittest.TestCase):

    def setUp(self):
        self.fileTransfer = FileTransfer()
    
    def testLocalFilesToOneRemoteHostTranfer(self):
        print "Local files to one remote host"
        passTest = transferFiles('localToRemote.xml', self.fileTransfer, "localToOneRemoteHost")
        deleteFiles('localToRemote.xml', self.fileTransfer, "localToOneRemoteHost")
        self.assertTrue(passTest)
    
    def testOneRemoteHostFilesToLocalTranfer(self):
        print "Remote files from one remote host to local machine"
        passTest = transferFiles('remoteToLocal.xml', self.fileTransfer, "oneRemoteHostToLocal")
        deleteFiles('remoteToLocal.xml', self.fileTransfer, "oneRemoteHostToLocal")
        self.assertTrue(passTest)
     
    def testLocalFilesToSeveralRemoteHostTranfer(self):
        print "Local files to several remote hosts"
        passTest = transferFiles('localToRemote.xml', self.fileTransfer, "localToSeveralRemoteHosts")
        deleteFiles('localToRemote.xml', self.fileTransfer, "localToSeveralRemoteHosts")
        self.assertTrue(passTest)
     
    def testSeveralRemoteHostFilesToLocalTranfer(self):
        print "Remote files from several remote hosts to local machine"
        passTest = transferFiles('remoteToLocal.xml', self.fileTransfer, "severalRemoteHostToLocal")
        deleteFiles('remoteToLocal.xml', self.fileTransfer, "severalRemoteHostToLocal")
        self.assertTrue(passTest)
    
    def testAllTranfer(self):
        print "All types"
        passTest = transferFiles('allTransfers.xml', self.fileTransfer, "all")
        deleteFiles('allTransfers.xml', self.fileTransfer, "all")
        self.assertTrue(passTest) 
     
    def testLocalToLocal(self):
        print "Local to local"
        passTest = transferFiles('localToLocal.xml', self.fileTransfer, "localToLocal")
        deleteFiles('localToLocal.xml', self.fileTransfer, "localToLocal")
        self.assertTrue(passTest) 
      
    def testLocalFilesTo(self):
        print "Local files to one remote host with method fileTransferTo"
        passTest = transferFilesOneHost('localToRemote.xml', self.fileTransfer, "fileTransferTo", True)
        deleteFiles('localToRemote.xml', self.fileTransfer, "localToOneRemoteHost")
        self.assertTrue(passTest)
     
    def tearDown(self):
        pass

def transferFiles(filesMetadata, fileTransfer, transferOperationName):
    """
    Transfer files taking the data from an XML file
    filesMetaData -- Xml file with all the transfer information.
    fileTransfer -- File transfer instance.
    transferOperationName -- Transfer operation name into de metadata.
    returns -- True if all files was copied correctly, false otherwise
    """
    tree = ET.parse(filesMetadata)
    root = tree.getroot()
    transferNode = root.find('fileTransfer[@name="' + transferOperationName + '"]')
    operationId = transferNode.attrib['operationId']
    trials = int(transferNode.attrib['trials'])
    forceOperation = bool(transferNode.attrib['forceOperation'])
    hostPasswords = getHostPasswords(transferNode)
    filePaths = getFilePaths(transferNode)
    start = time.time()
    fileTransfer.transferFiles(filePaths, hostPasswords, gatewayHosts=None, numberTrials=trials, forceOperation=forceOperation, operationId=operationId)
    print "Transfer time: " + str(time.time()- start)
    filePathList = getFilePathList(filePaths)
    passTest = len(fileTransfer.checkFiles(filePathList, hostPasswords, gatewayHosts=None, numberTrials=trials, forceOperation=forceOperation, operationId=operationId)) == 0
    return passTest 

def transferFilesOneHost(filesMetadata, fileTransfer, transferOperationName, isTo):
    """
    Transfer files taking the data from an XML file
    filesMetaData -- Xml file with all the transfer information.
    fileTransfer -- File transfer instance.
    transferOperationName -- Transfer operation name into de metadata.
    returns -- True if all files was copied correctly, false otherwise
    """
    #filesMetadata = 'localToRemote.xml'
    #transferOperationName = "fileTransferTo"
    tree = ET.parse(filesMetadata)
    root = tree.getroot()
    transferNode = root.find('fileTransfer[@name="' + transferOperationName + '"]')
    operationId = transferNode.attrib['operationId']
    trials = int(transferNode.attrib['trials'])
    forceOperation = bool(transferNode.attrib['forceOperation'])
    hostPasswords = getHostPasswords(transferNode)    
    #userAndHost = fileTransfer.__getUserAndHost(hostPasswords.keys()[0])
    userAndHost = getUserAndHost(hostPasswords.keys()[0])
    hostPassword = hostPasswords.values()[0]    
    tempFilePaths = getFilePaths(transferNode)    
    filePaths = {}
    for sourcePath, targetPaths in tempFilePaths.iteritems():
        filePaths[sourcePath] = targetPaths[0]    
    start = time.time()
    if (isTo):
        fileTransfer.transferFilesTo(filePaths, userAndHost[1], userAndHost[0], hostPassword, gatewayHosts = None, numberTrials = trials, forceOperation = forceOperation, operationId = operationId)
    else:
        fileTransfer.transferFilesFrom(filePaths, userAndHost[1], userAndHost[0], hostPassword, gatewayHosts = None, numberTrials = trials, forceOperation = forceOperation, operationId = operationId)
    print "Transfer time: " + str(time.time()- start)
    filePathList = getFilePathList(filePaths)
    passTest = len(fileTransfer.checkFiles(filePathList, hostPasswords, gatewayHosts=None, numberTrials=trials, forceOperation=forceOperation, operationId=operationId)) == 0
    return passTest 

def deleteFiles(filesMetadata, fileTransfer, transferOperationName):
    """
    Delete files taking the data from an XML file
    filesMetaData -- Xml file with all the transfer information.
    fileTransfer -- File transfer instance.
    transferOperationName -- Tranfer operation name into de metadata.
    """
    tree = ET.parse(filesMetadata)
    root = tree.getroot()
    transferNode = root.find('fileTransfer[@name="' + transferOperationName + '"]')
    operationId = transferNode.attrib['operationId']
    trials = int(transferNode.attrib['trials'])
    forceOperation = bool(transferNode.attrib['forceOperation'])
    hostPasswords = getHostPasswords(transferNode)
    filePaths = getFilePaths(transferNode)
    filePathList = getFilePathList(filePaths)
    fileTransfer.deleteFiles(filePathList, hostPasswords, gatewayHosts=None, numberTrials=trials, forceOperation=forceOperation, operationId=operationId)        
    
def getHostPasswords(transferNode):
    """
    Get host passwords from xml metadata node.
    transferNode -- Node from xml metadata file.
    """
    hostPasswords = {}
    for hostPassword in transferNode.iter('hostPassword'):
        userAndHost = hostPassword.find('userAndHost').text
        password = hostPassword.find('password').text
        hostPasswords[userAndHost] = password
    return hostPasswords;

def getFilePaths(transferNode):
    """
    Get file paths from xml metadata node.
    transferNode -- Node from xml metadata file.
    """
    filePaths = {}
    for transfer in transferNode.iter('transfer'):
        sourceFilePath = transfer.find('sourceFilePath').text
        targetFilePathList = []
        for target in transfer.iter('targetFilePath'):
            targetFilePath = target.text
            targetFilePathList.append(targetFilePath)
        filePaths[sourceFilePath] = targetFilePathList
    return filePaths

def getUserAndHost(userAndHost):
        """
        Function to get the user and the host name from 'userName@hostName' string
        Returns: Tuple with ('userName', 'hostName')
        """
        return userAndHost.split("@")

if __name__ == "__main__":
    #import sys;sys.argv = ["", "Test.testName"]
    unittest.main()