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

import sys
import os
import time
import unittest

from pyworkflow.utils.utils import getLocalUserName, getLocalHostName
from pyworkflow.utils.remote import join, RemotePath
import pyworkflow.utils.file_transfer as ft


try:
    remoteHostName1 = os.environ['remoteHostName1']
    remoteUserName1 = os.environ['remoteUserName1']
    remotePassword1 = os.environ['remotePassword1']
    remoteSourceFolder1 = os.environ['remoteSourceFolder1']
    remoteTargetFolder1 = os.environ['remoteTargetFolder1']

    remoteHostName2 = os.environ['remoteHostName2']
    remoteUserName2 = os.environ['remoteUserName2']
    remotePassword2 = os.environ['remotePassword2']
    remoteSourceFolder2 = os.environ['remoteSourceFolder2']
    remoteTargetFolder2 = os.environ['remoteTargetFolder2']

    localSourceFolder = os.environ['localSourceFolder']
    localTargetFolder = os.environ['localTargetFolder']
except KeyError as e:
    print e

class TestFileTransfer(unittest.TestCase):

    def setUp(self):
        self.fileTransfer = ft.FileTransfer()
        self.gatewayHosts = None
        self.numberTrials = 1
        self.forceOperation = False
        self.operationId = 1
        self.hostsRefs = self.getHostsRefs()
        self.hostPasswords = self.getHostsPasswords()

    def testLocalToSeveralRemoteHosts(self):
        tempFolder = "localToSeveralRemote"
        filePaths = {}
        sourceFilesPathList = ft.getFiles(localSourceFolder)
        for sourceFilePath in sourceFilesPathList:
            sourceFileName = ft.basename(sourceFilePath)
            targetFilePath1 = join(remoteTargetFolder1, tempFolder, sourceFileName)
            targetFilePath2 = join(remoteTargetFolder2, tempFolder, sourceFileName)
            targetFilePathList = []
            targetFilePathList.append(self.hostsRefs[0] + ":" + targetFilePath1)
            targetFilePathList.append(self.hostsRefs[1] + ":" + targetFilePath2)
            filePaths[sourceFilePath] = targetFilePathList
        self.fileTransfer.transferFiles(filePaths, self.hostPasswords, gatewayHosts=self.gatewayHosts, numberTrials=self.numberTrials, forceOperation=self.forceOperation, operationId=self.operationId)
        checkPathList = ft.getFilePathList(filePaths)
        passTest = len(self.fileTransfer.checkFiles(checkPathList, self.hostPasswords, gatewayHosts=self.gatewayHosts, numberTrials=self.numberTrials, forceOperation=self.forceOperation, operationId=self.operationId)) == 0
#         self.fileTransfer.deleteFiles(checkPathList, self.hostPasswords, gatewayHosts=self.gatewayHosts, numberTrials=self.numberTrials, forceOperation=self.forceOperation, operationId=self.operationId)
        ft.removeRemoteFolder(remoteHostName1, remoteUserName1, remotePassword1, join(remoteTargetFolder1, tempFolder))
        ft.removeRemoteFolder(remoteHostName2, remoteUserName2, remotePassword2, join(remoteTargetFolder2, tempFolder))
        self.assertTrue(passTest)

    def testSeveralRemoteHostsToLocal(self):
        tempFolder = "severalRemoteToLocal"
        filePaths = {}
        remoteSourceFilePathList1 = ft.getRemoteFolderFiles(remoteHostName1, remoteUserName1, remotePassword1, remoteSourceFolder1)
        remoteSourceFilePathList2 = ft.getRemoteFolderFiles(remoteHostName2, remoteUserName2, remotePassword2, remoteSourceFolder2)
        for remoteSourceFilePath in remoteSourceFilePathList1:
            targetFilePathList = []
            targetFilePathList.append(join(localTargetFolder, tempFolder, "test1", ft.basename(remoteSourceFilePath)))
            filePaths[self.hostsRefs[0] + ":" + remoteSourceFilePath] = targetFilePathList
        for remoteSourceFilePath in remoteSourceFilePathList2:
            targetFilePathList = []
            targetFilePathList.append(join(localTargetFolder, tempFolder, "test2", ft.basename(remoteSourceFilePath)))
            filePaths[self.hostsRefs[1] + ":" + remoteSourceFilePath] = targetFilePathList
        self.fileTransfer.transferFiles(filePaths, self.hostPasswords, gatewayHosts=self.gatewayHosts, numberTrials=self.numberTrials, forceOperation=self.forceOperation, operationId=self.operationId)
        checkPathList = ft.getFilePathList(filePaths)
        passTest = len(self.fileTransfer.checkFiles(checkPathList, self.hostPasswords, gatewayHosts=self.gatewayHosts, numberTrials=self.numberTrials, forceOperation=self.forceOperation, operationId=self.operationId)) == 0
#         self.fileTransfer.deleteFiles(checkPathList, self.hostPasswords, gatewayHosts=self.gatewayHosts, numberTrials=self.numberTrials, forceOperation=self.forceOperation, operationId=self.operationId)
        ft.cleanPath(join(localTargetFolder, tempFolder))
        self.assertTrue(passTest)

    def testLocalToLocal(self):
        tempFolder = "localToLocal"
        filePaths = {}
        sourceFilesPathList = ft.getFiles(localSourceFolder)
        for sourceFilePath in sourceFilesPathList:
            sourceFileName = ft.basename(sourceFilePath)
            targetFilePath = join(localTargetFolder, tempFolder, sourceFileName)
            targetFilePathList = []
            targetFilePathList.append(targetFilePath)
            filePaths[sourceFilePath] = targetFilePathList
        self.fileTransfer.transferFiles(filePaths, self.hostPasswords, gatewayHosts=self.gatewayHosts, numberTrials=self.numberTrials, forceOperation=self.forceOperation, operationId=self.operationId)
        checkPathList = ft.getFilePathList(filePaths)
        passTest = len(self.fileTransfer.checkFiles(checkPathList, self.hostPasswords, gatewayHosts=self.gatewayHosts, numberTrials=self.numberTrials, forceOperation=self.forceOperation, operationId=self.operationId)) == 0
#         self.fileTransfer.deleteFiles(checkPathList, self.hostPasswords, gatewayHosts=self.gatewayHosts, numberTrials=self.numberTrials, forceOperation=self.forceOperation, operationId=self.operationId)
        ft.cleanPath(join(localTargetFolder, tempFolder))
        self.assertTrue(passTest)

    def testLocalToOneRemoteHost(self):
        tempFolder = "localToOneRemote"
        filePaths = {}
        checkPathList = []
        sourceFilesPathList = ft.getFiles(localSourceFolder)
        for sourceFilePath in sourceFilesPathList:
            sourceFileName = ft.basename(sourceFilePath)
            targetFilePath = join(remoteTargetFolder1, tempFolder, sourceFileName)
            filePaths[sourceFilePath] = targetFilePath
            checkPathList.append(self.hostsRefs[0] + ":" + targetFilePath)
        self.fileTransfer.transferFilesTo(filePaths, remoteHostName1, remoteUserName1, remotePassword1, gatewayHosts=self.gatewayHosts, numberTrials=self.numberTrials, forceOperation=self.forceOperation, operationId=self.operationId)
        passTest = len(self.fileTransfer.checkFiles(checkPathList, self.hostPasswords, gatewayHosts=self.gatewayHosts, numberTrials=self.numberTrials, forceOperation=self.forceOperation, operationId=self.operationId)) == 0
        ft.removeRemoteFolder(remoteHostName1, remoteUserName1, remotePassword1, join(remoteTargetFolder1, tempFolder))
        self.assertTrue(passTest)

    def testOneRemoteHostToLocal(self):
        tempFolder = "oneRemoteHostToLocal"
        filePaths = {}
        checkPathList = []
        remoteSourceFilePathList2 = ft.getRemoteFolderFiles(remoteHostName2, remoteUserName2, remotePassword2, remoteSourceFolder2)
        for remoteSourceFilePath in remoteSourceFilePathList2:
            remoteSourceFileName = ft.basename(remoteSourceFilePath)
            targetFilePath = join(localTargetFolder, tempFolder, "test2", remoteSourceFileName)
            filePaths[remoteSourceFilePath] = targetFilePath
            checkPathList.append(targetFilePath)

        self.fileTransfer.transferFilesFrom(filePaths, remoteHostName2, remoteUserName2, remotePassword2, gatewayHosts=self.gatewayHosts, numberTrials=self.numberTrials, forceOperation=self.forceOperation, operationId=self.operationId)
        passTest = len(self.fileTransfer.checkFiles(checkPathList, self.hostPasswords, gatewayHosts=self.gatewayHosts, numberTrials=self.numberTrials, forceOperation=self.forceOperation, operationId=self.operationId)) == 0
        ft.cleanPath(join(localTargetFolder, tempFolder))
        self.assertTrue(passTest)

    def testAll(self):
        tempFolder = "allTransfer"
        filePaths = {}

        sourceFilesPathList = ft.getFiles(localSourceFolder)
        remoteSourceFilePathList1 = ft.getRemoteFolderFiles(remoteHostName1, remoteUserName1, remotePassword1, remoteSourceFolder1)
        remoteSourceFilePathList2 = ft.getRemoteFolderFiles(remoteHostName2, remoteUserName2, remotePassword2, remoteSourceFolder2)

        localUserName = getLocalUserName()
        localHostname = getLocalHostName()

        for sourceFilePath in sourceFilesPathList:
            sourceFileName = ft.basename(sourceFilePath)
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
            targetFilePathList.append(join(localTargetFolder, tempFolder, "test1", ft.basename(remoteSourceFilePath)))
            filePaths[self.hostsRefs[0] + ":" + remoteSourceFilePath] = targetFilePathList
        for remoteSourceFilePath in remoteSourceFilePathList2:
            targetFilePathList = []
            targetFilePathList.append(join(localTargetFolder, tempFolder, "test2", ft.basename(remoteSourceFilePath)))
            filePaths[self.hostsRefs[1] + ":" + remoteSourceFilePath] = targetFilePathList

        self.fileTransfer.transferFiles(filePaths, self.hostPasswords, gatewayHosts=self.gatewayHosts, numberTrials=self.numberTrials, forceOperation=self.forceOperation, operationId=self.operationId)
        checkPathList = ft.getFilePathList(filePaths)
        passTest = len(self.fileTransfer.checkFiles(checkPathList, self.hostPasswords, gatewayHosts=self.gatewayHosts, numberTrials=self.numberTrials, forceOperation=self.forceOperation, operationId=self.operationId)) == 0
#         self.fileTransfer.deleteFiles(checkPathList, self.hostPasswords, gatewayHosts=self.gatewayHosts, numberTrials=self.numberTrials, forceOperation=self.forceOperation, operationId=self.operationId)
        ft.cleanPath(join(localTargetFolder, tempFolder))
        ft.removeRemoteFolder (remoteHostName1, remoteUserName1, remotePassword1, join(remoteTargetFolder1, tempFolder))
        ft.removeRemoteFolder (remoteHostName2, remoteUserName2, remotePassword2, join(remoteTargetFolder2, tempFolder))
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




class TestRedirect(unittest.TestCase):


    def testRedirect(self):
        fOut = open('/home/antonio/Documents/Scipion/testFileTransfer/stdout.txt', 'a')
#         fErr = open('/home/antonio/Documents/Scipion/stderr.txt', 'a')
        sys.stdout = fOut
#         sys.stderr = fErr
        numberList = list(range(20))
        for num in numberList:
            print "WACKA " + str(num)
            #self.copyFile('/home/antonio/Documents/Scipion/testFileTransfer/stdout.txt', '/home/antonio/Documents/Scipion/testFileTransfer/stdoutCopy.txt')
            fOut.flush()
            time.sleep(3)
        fOut.close()

    def copyFile(self, file1, file2):
        f = open(file1)
        g = open(file2, 'a')
        linea = f.readline()
        while linea != "":
            g.write(linea)
            linea = f.readline()
        g.close()
        f.close()




class TestRemotePath(unittest.TestCase):

    def setUp(self):
        self.rpath = RemotePath.fromCredentials('crunchy', 'josem', 'kkk')
        self.localFolder = "/home/josem/Scipion"
        self.remoteFolder = "/gpfs/fs1/home/bioinfo/josem/Scipion"

    def test_copyTree(self):
        self.rpath.copyTree(self.localFolder, self.remoteFolder)

    def test_cleanPath(self):
        self.rpath.cleanPath(self.remoteFolder)
