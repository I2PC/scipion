'''
Created on May 8, 2013

@author: antonio
'''
import sys
from pyworkflow.manager import Manager
from pyworkflow.project import *
from pyworkflow.em import *
from em.packages.xmipp3.data import XmippSetOfMicrographs
import os
from os.path import split
from pyworkflow.utils.file_transfer import FileTransfer
from utils.file_transfer import FileTransfer

projName = sys.argv[1]
hostName = sys.argv[2]
folder = sys.argv[3]
userName = sys.argv[4]
password = sys.argv[5]

manager = Manager()
proj = manager.createProject(projName)  # Now it will be loaded if exists

projectFolder = split(proj.path)[1]
print ('Project folder: ' + projectFolder) 
targetProjectFolder = join(folder, projectFolder)

filePaths = {}

result = proj.mapper.selectByClass('EmanProtBoxing')
if len(result):    
    for emanProtBoxing in result:
        print ("EmanProtBoxing filePaths: ")
        for filePath in emanProtBoxing.getFiles():
            sourceFilePath = os.path.join(proj.path, filePath)
            targetFilePathList = []
            targetFilePath = userName + '@' + hostName + ':' + os.path.join(targetProjectFolder, filePath)
            targetFilePathList.append(targetFilePath)
            filePaths[sourceFilePath] = targetFilePathList
else:
    print "Not EmanProtBoxing found"

print "File paths: " + str(filePaths)
hostsPasswords = {}
hostsPasswords[userName + '@' + hostName] = password

fileTransfer = FileTransfer()
fileTransfer.transferFiles(filePaths, hostsPasswords, gatewayHosts = None, numberTrials = 1, forceOperation = False, operationId = 1)
