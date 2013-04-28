"""
Created on Apr 9, 2013

@author: antonio
"""
import unittest
#from fileTransfer.fileTransfer import FileTransfer
from file_transfer.fileTransfer import *
import xml.etree.ElementTree as ET
import time

class TestFileTransfer(unittest.TestCase):

    def setUp(self):
        self.fileTransferProces = FileTransfer()
    
    def testLocalFilesToOneRemoteHostTranfer(self):
        passTest = transferFiles('localToRemote.xml', self.fileTransferProces, "localToOneRemoteHost")
        self.assertTrue(passTest)
    
    def testOneRemoteHostFilesToLocalTranfer(self):
        passTest = transferFiles('remoteToLocal.xml', self.fileTransferProces, "oneRemoteHostToLocal")
        self.assertTrue(passTest)
    def testAllTranfer(self):
        passTest = transferFiles('allTransfers.xml', self.fileTransferProces, "all")
        self.assertTrue(passTest) 
       
    def testLocalFilesToSeveralRemoteHostTranfer(self):
        passTest = transferFiles('localToRemote.xml', self.fileTransferProces, "localToSeveralRemoteHosts")
        self.assertTrue(passTest)
    
    def tearDown(self):
        pass



def transferFiles(filesMetadata, fileTransferProces, transferOperationName):
    tree = ET.parse(filesMetadata)
    root = tree.getroot()
    transferNode = root.find('fileTransfer[@name="' + transferOperationName + '"]')
    operationId = transferNode.attrib['operationId']
    trials = int(transferNode.attrib['trials'])
    forceOperation = bool(transferNode.attrib['forceOperation'])
    hostPasswords = getHostPasswords(transferNode)
    filePaths = getFilePaths(transferNode)
    start = time.time()
    fileTransferProces.transferFiles(filePaths, hostPasswords, gatewayHosts=None, operationId=operationId, numberTrials=trials, forceOperation=forceOperation);
    print "Transfer time: " + str(time.time()- start)
    filePathList = getFilePathList(filePaths)
    passTest = len(fileTransferProces.checkFiles(filePathList, hostPasswords, gatewayHosts=None, operationId=operationId, numberTrials=trials, forceOperation=forceOperation)) == 0
    return passTest        
    
def getHostPasswords(transferNode):
    hostPasswords = {}
    for hostPassword in transferNode.iter('hostPassword'):
        userAndHost = hostPassword.find('userAndHost').text
        password = hostPassword.find('password').text
        hostPasswords[userAndHost] = password
    return hostPasswords;

def getFilePaths(transferNode):
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