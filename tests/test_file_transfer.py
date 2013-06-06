"""
Created on Apr 9, 2013

@author: antonio
"""
import unittest
import argparse
import sys
from utils.path import getFolderFiles, cleanPath
from os.path import join, basename
from pyworkflow.utils.file_transfer import *
from pyworkflow.utils.utils import executeRemote, getLocalUserName, getLocalHostName

remoteHostName1 = None
remoteUserName1 = None
remotePassword1 = None
remoteSourceFolder1 = None
remoteTargetFolder1 = None

remoteHostName2 = None
remoteUserName2 = None
remotePassword2 = None
remoteSourceFolder2 = None
remoteTargetFolder2 = None

localSourceFolder = None
localTargetFolder = None

class TestFileTransfer(unittest.TestCase):  
    
    def setUp(self):
        self.fileTransfer = FileTransfer()
        self.gatewayHosts = None 
        self.numberTrials = 1                      
        self.forceOperation = False
        self.operationId = 1
        self.hostsRefs = self.getHostsRefs()
        self.hostPasswords = self.getHostsPasswords()
      
    def testLocalToSeveralRemoteHosts(self):
        tempFolder = "localToSeveralRemote"
        filePaths = {}
        sourceFilesPathList = getFolderFiles(localSourceFolder)
        for sourceFilePath in sourceFilesPathList:
            sourceFileName = basename(sourceFilePath)
            targetFilePath1 = join(remoteTargetFolder1, tempFolder, sourceFileName)
            targetFilePath2 = join(remoteTargetFolder2, tempFolder, sourceFileName)
            targetFilePathList = []
            targetFilePathList.append(self.hostsRefs[0] + ":" + targetFilePath1)
            targetFilePathList.append(self.hostsRefs[1] + ":" + targetFilePath2)
            filePaths[sourceFilePath] = targetFilePathList
        self.fileTransfer.transferFiles(filePaths, self.hostPasswords, gatewayHosts=self.gatewayHosts, numberTrials=self.numberTrials, forceOperation=self.forceOperation, operationId=self.operationId)
        checkPathList = getFilePathList(filePaths)
        passTest = len(self.fileTransfer.checkFiles(checkPathList, self.hostPasswords, gatewayHosts=self.gatewayHosts, numberTrials=self.numberTrials, forceOperation=self.forceOperation, operationId=self.operationId)) == 0
#         self.fileTransfer.deleteFiles(checkPathList, self.hostPasswords, gatewayHosts=self.gatewayHosts, numberTrials=self.numberTrials, forceOperation=self.forceOperation, operationId=self.operationId)
        removeRemoteFolder (remoteHostName1, remoteUserName1, remotePassword1, join(remoteTargetFolder1, tempFolder))
        removeRemoteFolder (remoteHostName2, remoteUserName2, remotePassword2, join(remoteTargetFolder2, tempFolder))
        self.assertTrue(passTest)
         
    def testSeveralRemoteHostsToLocal(self):
        tempFolder = "severalRemoteToLocal"
        filePaths = {}
        remoteSourceFilePathList1 = getRemoteFolderFiles(remoteHostName1, remoteUserName1, remotePassword1, remoteSourceFolder1)
        remoteSourceFilePathList2 = getRemoteFolderFiles(remoteHostName2, remoteUserName2, remotePassword2, remoteSourceFolder2)
        for remoteSourceFilePath in remoteSourceFilePathList1:
            targetFilePathList = []
            targetFilePathList.append(join(localTargetFolder, tempFolder, "test1", basename(remoteSourceFilePath)))
            filePaths[self.hostsRefs[0] + ":" + remoteSourceFilePath] = targetFilePathList
        for remoteSourceFilePath in remoteSourceFilePathList2:
            targetFilePathList = []
            targetFilePathList.append(join(localTargetFolder, tempFolder, "test2", basename(remoteSourceFilePath)))
            filePaths[self.hostsRefs[1] + ":" + remoteSourceFilePath] = targetFilePathList
        self.fileTransfer.transferFiles(filePaths, self.hostPasswords, gatewayHosts=self.gatewayHosts, numberTrials=self.numberTrials, forceOperation=self.forceOperation, operationId=self.operationId)
        checkPathList = getFilePathList(filePaths)
        passTest = len(self.fileTransfer.checkFiles(checkPathList, self.hostPasswords, gatewayHosts=self.gatewayHosts, numberTrials=self.numberTrials, forceOperation=self.forceOperation, operationId=self.operationId)) == 0
#         self.fileTransfer.deleteFiles(checkPathList, self.hostPasswords, gatewayHosts=self.gatewayHosts, numberTrials=self.numberTrials, forceOperation=self.forceOperation, operationId=self.operationId)
        cleanPath(join(localTargetFolder, tempFolder))
        self.assertTrue(passTest) 
           
    def testLocalToLocal(self):
        tempFolder = "localToLocal"
        filePaths = {}
        sourceFilesPathList = getFolderFiles(localSourceFolder)
        for sourceFilePath in sourceFilesPathList:
            sourceFileName = basename(sourceFilePath)
            targetFilePath = join(localTargetFolder, tempFolder, sourceFileName)
            targetFilePathList = []
            targetFilePathList.append(targetFilePath)
            filePaths[sourceFilePath] = targetFilePathList
        self.fileTransfer.transferFiles(filePaths, self.hostPasswords, gatewayHosts=self.gatewayHosts, numberTrials=self.numberTrials, forceOperation=self.forceOperation, operationId=self.operationId)
        checkPathList = getFilePathList(filePaths)
        passTest = len(self.fileTransfer.checkFiles(checkPathList, self.hostPasswords, gatewayHosts=self.gatewayHosts, numberTrials=self.numberTrials, forceOperation=self.forceOperation, operationId=self.operationId)) == 0
#         self.fileTransfer.deleteFiles(checkPathList, self.hostPasswords, gatewayHosts=self.gatewayHosts, numberTrials=self.numberTrials, forceOperation=self.forceOperation, operationId=self.operationId)
        cleanPath(join(localTargetFolder, tempFolder))
        self.assertTrue(passTest)
          
    def testLocalToOneRemoteHost(self):
        tempFolder = "localToOneRemote"
        filePaths = {}
        checkPathList = []
        sourceFilesPathList = getFolderFiles(localSourceFolder)
        for sourceFilePath in sourceFilesPathList:
            sourceFileName = basename(sourceFilePath)
            targetFilePath = join(remoteTargetFolder1, tempFolder, sourceFileName)
            filePaths[sourceFilePath] = targetFilePath
            checkPathList.append(self.hostsRefs[0] + ":" + targetFilePath)
        self.fileTransfer.transferFilesTo(filePaths, remoteHostName1, remoteUserName1, remotePassword1, gatewayHosts=self.gatewayHosts, numberTrials=self.numberTrials, forceOperation=self.forceOperation, operationId=self.operationId)
        passTest = len(self.fileTransfer.checkFiles(checkPathList, self.hostPasswords, gatewayHosts=self.gatewayHosts, numberTrials=self.numberTrials, forceOperation=self.forceOperation, operationId=self.operationId)) == 0
        removeRemoteFolder (remoteHostName1, remoteUserName1, remotePassword1, join(remoteTargetFolder1, tempFolder))
        self.assertTrue(passTest)
     
    def testOneRemoteHostToLocal(self):
        tempFolder = "oneRemoteHostToLocal"
        filePaths = {}
        checkPathList = []
        remoteSourceFilePathList2 = getRemoteFolderFiles(remoteHostName2, remoteUserName2, remotePassword2, remoteSourceFolder2)
        for remoteSourceFilePath in remoteSourceFilePathList2:
            remoteSourceFileName = basename(remoteSourceFilePath)
            targetFilePath = join(localTargetFolder, tempFolder, "test2", remoteSourceFileName)
            filePaths[remoteSourceFilePath] = targetFilePath
            checkPathList.append(targetFilePath)
          
        self.fileTransfer.transferFilesFrom(filePaths, remoteHostName2, remoteUserName2, remotePassword2, gatewayHosts=self.gatewayHosts, numberTrials=self.numberTrials, forceOperation=self.forceOperation, operationId=self.operationId)
        passTest = len(self.fileTransfer.checkFiles(checkPathList, self.hostPasswords, gatewayHosts=self.gatewayHosts, numberTrials=self.numberTrials, forceOperation=self.forceOperation, operationId=self.operationId)) == 0
        cleanPath(join(localTargetFolder, tempFolder))
        self.assertTrue(passTest)

    def testAll(self):
        tempFolder = "allTransfer"
        filePaths = {}
        
        sourceFilesPathList = getFolderFiles(localSourceFolder)
        remoteSourceFilePathList1 = getRemoteFolderFiles(remoteHostName1, remoteUserName1, remotePassword1, remoteSourceFolder1)
        remoteSourceFilePathList2 = getRemoteFolderFiles(remoteHostName2, remoteUserName2, remotePassword2, remoteSourceFolder2)
        
        localUserName = getLocalUserName()
        localHostname = getLocalHostName()
        
        for sourceFilePath in sourceFilesPathList:
            sourceFileName = basename(sourceFilePath)
            targetFilePath1 = join(remoteTargetFolder1, tempFolder, sourceFileName)
            targetFilePath2 = join(remoteTargetFolder2, tempFolder, sourceFileName)
            targetFilePath3 = join (localTargetFolder, tempFolder, sourceFileName)
            targetFilePathList = []
            targetFilePathList.append(self.hostsRefs[0] + ":" + targetFilePath1)
            targetFilePathList.append(self.hostsRefs[1] + ":" + targetFilePath2)
            targetFilePathList.append(localUserName + "@" + localHostname + ":" + targetFilePath3)
            filePaths[sourceFilePath] = targetFilePathList
        
        for remoteSourceFilePath in remoteSourceFilePathList1:
            targetFilePathList = []
            targetFilePathList.append(join(localTargetFolder, tempFolder, "test1", basename(remoteSourceFilePath)))
            filePaths[self.hostsRefs[0] + ":" + remoteSourceFilePath] = targetFilePathList
        for remoteSourceFilePath in remoteSourceFilePathList2:
            targetFilePathList = []
            targetFilePathList.append(join(localTargetFolder, tempFolder, "test2", basename(remoteSourceFilePath)))
            filePaths[self.hostsRefs[1] + ":" + remoteSourceFilePath] = targetFilePathList
        
        self.fileTransfer.transferFiles(filePaths, self.hostPasswords, gatewayHosts=self.gatewayHosts, numberTrials=self.numberTrials, forceOperation=self.forceOperation, operationId=self.operationId)
        checkPathList = getFilePathList(filePaths)
        passTest = len(self.fileTransfer.checkFiles(checkPathList, self.hostPasswords, gatewayHosts=self.gatewayHosts, numberTrials=self.numberTrials, forceOperation=self.forceOperation, operationId=self.operationId)) == 0
#         self.fileTransfer.deleteFiles(checkPathList, self.hostPasswords, gatewayHosts=self.gatewayHosts, numberTrials=self.numberTrials, forceOperation=self.forceOperation, operationId=self.operationId)
        cleanPath(join(localTargetFolder, tempFolder))
        removeRemoteFolder (remoteHostName1, remoteUserName1, remotePassword1, join(remoteTargetFolder1, tempFolder))
        removeRemoteFolder (remoteHostName2, remoteUserName2, remotePassword2, join(remoteTargetFolder2, tempFolder))
        self.assertTrue(passTest) 
         
    def getHostsPasswords(self):
        hostPasswords = {}
        hostPasswords[self.hostsRefs[0]] = remotePassword1
        hostPasswords[self.hostsRefs[1]] = remotePassword2
        return hostPasswords
    
    def getHostsRefs(self):
        hostRefs = []
        hostRefs.append(remoteUserName1 + "@" + remoteHostName1)
        hostRefs.append(remoteUserName2 + "@" + remoteHostName2)
        return hostRefs

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('remoteHostName1')
    parser.add_argument('remoteUserName1')
    parser.add_argument('remotePassword1')
    parser.add_argument('remoteSourceFolder1')
    parser.add_argument('remoteTargetFolder1')
    parser.add_argument('remoteHostName2')
    parser.add_argument('remoteUserName2')
    parser.add_argument('remotePassword2')
    parser.add_argument('remoteSourceFolder2')
    parser.add_argument('remoteTargetFolder2')
    parser.add_argument('localSourceFolder')
    parser.add_argument('localTargetFolder')
    parser.add_argument('unittest_args', nargs='*')
    args = parser.parse_args()
    remoteHostName1 = args.remoteHostName1
    remoteUserName1 = args.remoteUserName1
    remotePassword1 = args.remotePassword1
    remoteSourceFolder1 = args.remoteSourceFolder1
    remoteTargetFolder1 = args.remoteTargetFolder1
    remoteHostName2 = args.remoteHostName2
    remoteUserName2 = args.remoteUserName2
    remotePassword2 = args.remotePassword2
    remoteSourceFolder2 = args.remoteSourceFolder2
    remoteTargetFolder2 = args.remoteTargetFolder2   
    localSourceFolder = args.localSourceFolder
    localTargetFolder = args.localTargetFolder 
    sys.argv[1:] = args.unittest_args
    unittest.main()