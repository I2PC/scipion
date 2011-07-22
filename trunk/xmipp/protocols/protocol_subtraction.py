#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Xmipp protocol for subtraction
#
# Example use:
# ./xmipp_partial_projection_subtraction.py
#
# Authors: 
#

from xmipp import *
from protlib_base import *
from protlib_utils import *
#from protlib_utils import getListFromVector
#import os,sys
from protocol_partial_projection_subtraction_header import *

def checkErrors():    
    '''This function will be used to validate the protocols
    should be implemented in all derived protocols'''
    errors = []        
    # Check if there is workingdir 
    # move this to gui
    if not os.path.exists(ProtocolName):
        errors.append("Refered protocol named %s does not exist"%ProtocolName)
        
    return errors 
        
class ProtPartialProjectionSubtraction(XmippProtocol):

    def __init__(self, scriptname,project=None):
        
	XmippProtocol.__init__(self, protDict.projsubs.key, scriptname, project)
        print "0"
	self.Import = 'from protocol_partial_projection_subtraction import *'
	print "1"
	#import config
        #super(ProtPartialProjectionSubtraction,self).__init__(protDict.projsubs.key, scriptname, project)
        print "2"	
	self.myName='partial_projection_subtraction'
        self.subtractionDir ='Subtraction'
        self.referenceDir   ='Refs'
        self.referenceStack ='ref'
        self.subImgsDir  = 'SubImgs'
        self.subtractedStack ='subtracted'
        self.volsDir ='Vols'
        self.tempFileName=''
        self.current_angles='current_angles.doc'
        self.Import = 'from protlib_partial_projection_subtraction import *'
        self.scaledImages = 'scaled'
        self.runName = RunName
        
        # Check these params
        self.DoDeleteWorkingDir = False
        self.DoParallel = True
        
    def preRun(self):

        self.pmprotWorkingDir = WorkingDir
        self.Iteration_Working_Directory = os.path.join(self.pmprotWorkingDir,'Iter_'+ str(iterationNo))
        self.subtractionDir = os.path.join(self.WorkingDir,RunName)
        self.volsDir = os.path.join(self.WorkingDir,self.volsDir)
        self.referenceDir = os.path.join(self.WorkingDir,self.referenceDir)
        self.subImgsDir = os.path.join(self.WorkingDir,self.subImgsDir)
        
        self.scaledImages = os.path.join(self.WorkingDir,self.scaledImages)
        
        self.doScaleImages = doScaleImages
        
        self.dimX = dimX
        self.dimY = dimY
        
        self.dRradiusMax = dRradiusMax
        self.dRradiusMin = dRradiusMin
        
        if(MaxChangeInAngles > 100):
            self.MaxChangeInAngles=-1
        else:
            self.MaxChangeInAngles = MaxChangeInAngles
            
        #new vwersion need zfill
        tmpFilename = 'Iter_'+ str(iterationNo) + '_' + self.current_angles #zfill()
        self.filename_currentAngles = os.path.join(self.Iteration_Working_Directory,tmpFilename)
        
        if(self.doScaleImages):
            
            md = MetaData(self.filename_currentAngles)
            img = Image()
            img.readApplyGeo(md,       1, False, DATA, ALL_IMAGES,False)
            x=y=z=n=0
            (x,y,z,n) = img.getDimensions()
            factorX = x / self.dimX
            if (self.dimY<0):
                factorY = y / self.dimX
            else:
                factorY = y / self.dimY

            self.dRradiusMax = round(self.dRradiusMax/factorX)
            self.dRradiusMin = round(self.dRradiusMin/factorY)
            
        print "dRradiusMax: ", self.dRradiusMax
        print "dRradiusMin: ", self.dRradiusMin
        
        #if not ctf info available use original sel FIXME
        tmpFileName = os.path.join(self.pmprotWorkingDir,'CtfGroups/ctf_group??????.sel')
        self.defGroups=['']
        import glob
        self.defGroups += glob.glob(tmpFileName)
        #if not CTF groups available
        if  len(self.defGroups)<2:
            self.defGroups.append(self.filename_currentAngles) 
            
        self.defocusGroupNo = len(self.defGroups)
        
        self.DocFileExp=['']
        if DocFileExp!='':
            tmpFileName = DocFileExp
        else:
            tmpFileName = FileName(self.filename_currentAngles)
        tmpFileName = tmpFileName.withoutExtension()
        tmpFileName = tmpFileName + '_ctfgroups.doc'
        if os.path.exists(tmpFileName):
            os.remove(tmpFileName)
        for iterN in range(1, self.defocusGroupNo ):
            tmpDocFileExp = 'ctfgroup_' + str(iterN).zfill(6) + '@' + tmpFileName
            self.DocFileExp.append(tmpDocFileExp)
            
        self.reconstructedVolume = [""]
        for iterN in range(1, self.defocusGroupNo ):
            tmpReconstruct = os.path.join(self.volsDir, 'rec_ctfg' + str(iterN).zfill(6) + '.vol')
            self.reconstructedVolume.append(tmpReconstruct)
            
        self.maskReconstructedVolume = [""]
        for iterN in range(1, self.defocusGroupNo ):
            tmpReconstruct = os.path.join(self.volsDir, 'rec_ctfg' + str(iterN).zfill(6) + '_mask.vol')
            self.maskReconstructedVolume.append(tmpReconstruct)
    
        tmp = self.referenceStack  #'ref'
        self.referenceStackDoc = [""]
        for iterN in range(1, self.defocusGroupNo ):
            tmpReference = os.path.join(self.referenceDir, tmp + '_' + str(iterN).zfill(6) + '.doc')
            self.referenceStackDoc.append(tmpReference)
    
        tmp = self.referenceStack
        self.referenceStack = [""]
        for iterN in range(1, self.defocusGroupNo ):
            tmpReference = os.path.join(self.referenceDir, tmp + '_' + str(iterN).zfill(6) + '.stk')
            self.referenceStack.append(tmpReference)
    
        tmp = self.subtractedStack
        self.subtractedStack = [""]
        for iterN in range(1, self.defocusGroupNo ):
            tmpSubtract = os.path.join(self.subImgsDir, tmp + '_' + str(iterN).zfill(6) + '.stk')
            self.subtractedStack.append(tmpSubtract)
        #FIXME
        #DO we need to import this paths or are already in pythonpath
        #import shutil,time
        scriptdir = os.path.split(os.path.dirname(os.popen('which xmipp_protocols', 'r').read()))[0] + '/protocols'
        sys.path.append(scriptdir)

        #from pysqlite2 import dbapi2 as sqlite
        self.Db.setPrintWrapperParameters(PrintWrapperParameters)
        self.Db.setPrintWrapperCommand(PrintWrapperCommand)
        self.Db.setVerify(Verify,ViewVerifyedFiles)
        
        
    def otherActionsToBePerformedBeforeLoop(self):
        _dataBase = self.Db
        #Create directories
        _dataBase.insertAction('createDir', None, path = self.volsDir)
        _dataBase.insertAction('createDir', None, path = self.referenceDir)
        _dataBase.insertAction('createDir', None, path = self.subImgsDir)
        
        #Create auxiliary metadata with image names , angles and CTF
        
        if(doScaleImages):
            _VerifyFiles = [self.scaledImages+".stk"]
            _VerifyFiles.append(self.scaledImages+".xmd")
            id = _dataBase.insertAction('scaleImages', _VerifyFiles, None, None, None
                                       , dimX = self.dimX
                                       , dimY = self.dimY
                                       , DoParallel = self.DoParallel
                                       , filename_currentAngles = self.filename_currentAngles
                                       , MpiJobSize = MpiJobSize
                                       , NumberOfMpiProcesses = NumberOfMpiProcesses
                                       , NumberOfThreads = NumberOfThreads
                                       , scaledImages = self.scaledImages
                                       , SystemFlavour = SystemFlavour
                                       )            
    def actionsToBePerformedInsideLoop(self):
        
        _dataBase = self.Db
        for iterN in range(1, self.defocusGroupNo):
            #Create auxiliary metadata with image names , angles and CTF
            if(self.doScaleImages):
                inputSelfile = self.scaledImages +'.xmd'
            else:
                inputSelfile = self.filename_currentAngles
                
            print 'self.defGroups[' + str(iterN) + ']: ', self.defGroups[iterN]
            print 'self.defocusGroupNo: ' + str(self.defocusGroupNo)
            
            if(self.defocusGroupNo > 2 and self.doScaleImages):
                _VerifyFiles = []
                auxFilename = FileName(self.DocFileExp[iterN])
                _VerifyFiles.append(auxFilename.removeBlockName())
                id = self.Db.insertAction('joinImageCTFscale', _VerifyFiles
                                        , CTFgroupName = self.defGroups[iterN]
                                        , DocFileExp = self.DocFileExp[iterN]
                                        , inputSelfile = inputSelfile
                                        )
            elif (self.defocusGroupNo > 2 and not self.doScaleImages):
                _VerifyFiles = []
                auxFilename = FileName(self.DocFileExp[iterN])
                _VerifyFiles.append(auxFilename.removeBlockName())
                id = self.Db.insertAction('joinImageCTF', _VerifyFiles
                                        , CTFgroupName = self.defGroups[iterN]
                                        , DocFileExp = self.DocFileExp[iterN]
                                        , inputSelfile = inputSelfile
                                        )                
            else:
                self.DocFileExp[iterN] = self.defGroups[iterN]
            #reconstruct each CTF group
            _VerifyFiles = []
            _VerifyFiles.append(self.reconstructedVolume[iterN])
            id = self.Db.insertAction('reconstructVolume', _VerifyFiles
                                        , DocFileExp = self.DocFileExp[iterN]
                                        , DoParallel = self.DoParallel
                                        , MpiJobSize = MpiJobSize
                                        , NumberOfMpiProcesses = NumberOfMpiProcesses
                                        , NumberOfThreads = NumberOfThreads
                                        , reconstructedVolume = self.reconstructedVolume[iterN]
                                        , SymmetryGroup = SymmetryGroup
                                        , SystemFlavour = SystemFlavour)
            
            #mask volume before projection
            _VerifyFiles = []
            auxFilename = FileName(self.DocFileExp[iterN])
            _VerifyFiles.append(self.maskReconstructedVolume[iterN])
            id = self.Db.insertAction('maskVolume', _VerifyFiles
                                        , dRradiusMax = self.dRradiusMax
                                        , dRradiusMin = self.dRradiusMin
                                        , maskReconstructedVolume = self.maskReconstructedVolume[iterN]
                                        , reconstructedVolume = self.reconstructedVolume[iterN])
    
            #project reconstructe4d volumes
            _VerifyFiles = []
            _VerifyFiles.append(self.referenceStack[iterN])
            tmp = self.referenceStack[iterN]
            _VerifyFiles.append(tmp.replace('.stk','.doc'))
            _VerifyFiles.append(tmp.replace('.stk','_sampling.xmd'))
            
            id = self.Db.insertAction('createProjections', _VerifyFiles
                                        , AngSamplingRateDeg = AngSamplingRateDeg
                                        , DocFileExp = self.DocFileExp[iterN]
                                        , DoParallel = self.DoParallel
                                        , maskReconstructedVolume = self.maskReconstructedVolume[iterN]
                                        , MaxChangeInAngles = self.MaxChangeInAngles
                                        , MpiJobSize = MpiJobSize
                                        , NumberOfMpiProcesses = NumberOfMpiProcesses
                                        , NumberOfThreads = NumberOfThreads
                                        , referenceStack = self.referenceStack[iterN]
                                        , SymmetryGroup = SymmetryGroup
                                        , SystemFlavour = SystemFlavour)
                          
                     
                     
    #project reconstructe4d volumes
            _Parameters = {
                          'DocFileExp':self.DocFileExp[iterN]
                          ,'referenceStackDoc':self.referenceStackDoc[iterN]
                          ,'subtractedStack':self.subtractedStack[iterN]
                          }
            command = "subtractionScript"
            _VerifyFiles = []
            _VerifyFiles.append(self.subtractedStack[iterN])
            _VerifyFiles.append(self.subtractedStack[iterN]+'ref')
            _VerifyFiles.append(self.subtractedStack[iterN]+'exp')
            
            id = self.Db.insertAction('subtractionScript', _VerifyFiles
                                        , DocFileExp = self.DocFileExp[iterN] 
                                        , referenceStackDoc = self.referenceStackDoc[iterN]
                                        , subtractedStack = self.subtractedStack[iterN])
                          
    def defineActions(self):
        self.preRun()
        self.otherActionsToBePerformedBeforeLoop()
        self.actionsToBePerformedInsideLoop()
        
    def validate(self):
        return checkErrors()
    
def ImportProtocol():
    print "Import Protocol"
    scriptdir = os.getcwd()
    sys.path.append(scriptdir)

    from protlib_utils import unique_filename
    from protlib_utils import loadModule
    pardir=os.path.abspath(os.getcwd())
    tempFileName=unique_filename(ProtocolName)
    fn = FileName(tempFileName)
    fn=fn.getBaseName()
    shutil.copy(ProtocolName,fn+'.py')
    
    #FIXME use load module??
    exec  "import " + fn
    #loadModule(fn)
    #import arg
    global ProjectDir
    ProjectDir = eval(fn +'.ProjectDir')
    global WorkingDir
    WorkingDir = eval(fn +'.WorkingDir')
    global LogDir
    LogDir = eval(fn +'.LogDir')
    global SymmetryGroup
    SymmetryGroup = eval(fn +'.SymmetryGroup')
    global AngSamplingRateDeg
    if(len(AngSamplingRateDeg) < 1):
        AngSamplingRateDeg=getComponentFromVector(eval(fn +'.AngSamplingRateDeg'),iterationNo - 1)
    global MaxChangeInAngles
    if(len(MaxChangeInAngles) < 1):
        MaxChangeInAngles=getComponentFromVector(eval(fn +'.MaxChangeInAngles'),iterationNo - 1)
    global refDirNameAngSamplingRateDeg
    #refDirName=  eval(fn +'.LibraryDir')
    
           
if __name__ == '__main__':
    ImportProtocol()
    protocolMain(ProtPartialProjectionSubtraction)
    
