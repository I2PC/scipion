'''
Created on Apr 9, 2013

@author: antonio
'''
import unittest
import fileTransferProcess.fileTransferProcess.FileTransferProcess as FileTransferProcess

class TestFileTransferProcess(unittest.TestCase):

    sourceFilesCase1 = None
    targetFilesCase1 = None
    gatewayHostsCase1 = None
    hostsPaswordCase1s = None
    operationIdCase1 = None 
    numberTrialsCase1 = None
    fileTransferProcess = None
    forceOperationCase1 = None

    def setUp(self):
        self.fileTransferProces = FileTransferProcess()
        self.sourceFilesCase1 = {'1':'/home/antonio/Documents/Scipion/ScipionTestFiles/1212_00029.xmp', '2':"/home/antonio/Documents/Scipion/ScipionTestFiles/1212_00036.xmp"}
        self.targetFilesCase1 = {'1':['apoza@glassfishdev.cnb.csic.es:/home/apoza/testFileTransfer/testFile1.xmp'], '2':['apoza@glassfishdev.cnb.csic.es:/home/apoza/testFileTransfer/testFile2.xmp']}
        self.hostsPaswordsCase1 = {'apoza@glassfishdev.cnb.csic.es':'BF6fYiFYiYD'}
        self.operationIdCase1 = 1 
        self.numberTrialsCase1 = 1
        self.forceOperationCase1 = False
    
    def tearDownModule(self):
        filePaths = self.fileTransferProces.getTargetFilePathList(self.targetFilesCase1)
        self.fileTransferProces.deleteFiles(filePaths, self.gatewayHostsCase1, self.hostsPaswordsCase1, 1, 1, False)
        
    def testSimpleLocalToRemoteFileTranfer(self):
        self.fileTransferProces.transferFiles(self.sourceFilesCase1, self.targetFilesCase1, self.gatewayHostsCase1, self.hostsPaswordsCase1, self.operationIdCase1, self.numberTrialsCase1);
        filePaths = self.fileTransferProces.getTargetFilePathList(self.targetFilesCase1)
        passTest = len(self.fileTransferProces.checkFiles(filePaths, self.gatewayHostsCase1, self.hostsPaswordCase1s, self.operationIdCase1, self.numberTrialsCase1, self.forceOperationCase1)) == 0
        self.assertTrue(passTest)
    


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()