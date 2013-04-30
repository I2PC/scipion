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
        
    def tearDown(self):
        self.fileTransfer.close()

def transferFiles(filesMetadata, fileTransfer, transferOperationName):
    """
    Transfer files taking the data from an XML file
    filesMetaData -- Xml file with all the transfer information.
    fileTransfer -- File transfer instance.
    transferOperationName -- Tranfer operation name into de metadata.
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
    fileTransfer.transferFiles(filePaths, hostPasswords, gatewayHosts=None, operationId=operationId, numberTrials=trials, forceOperation=forceOperation)
    print "Transfer time: " + str(time.time()- start)
    filePathList = getFilePathList(filePaths)
    passTest = len(fileTransfer.checkFiles(filePathList, hostPasswords, gatewayHosts=None, operationId=operationId, numberTrials=trials, forceOperation=forceOperation)) == 0
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
    fileTransfer.deleteFiles(filePathList, hostPasswords, gatewayHosts=None, operationId=operationId, numberTrials=trials, forceOperation=forceOperation)        
    
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

if __name__ == "__main__":
    #import sys;sys.argv = ["", "Test.testName"]
    unittest.main()