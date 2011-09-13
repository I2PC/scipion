#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Xmipp protocol for partial projection subtraction
#
# Example use:
# ./protocol_subtraction_header.py
#
# Authors: Roberto Marabini,
#          Alejandro Echeverria Rey.
#
#from xmipp import *
#from protlib_base import *
#from protlib_utils import *
        
import os, shutil, sys
from xmipp import FileName, MetaData, Image
from xmipp import DATA, ALL_IMAGES
#from xmipp import MetaData, FILENAMENUMBERLENGTH, AGGR_COUNT, MDL_CTFMODEL,MDL_COUNT
from protlib_base import XmippProtocol, protocolMain
#from protlib_utils import getListFromVector, getBoolListFromVector, getComponentFromVector
from protlib_utils import getComponentFromVector
from protlib_sql import XmippProjectDb
from config_protocols import protDict

class ProtPartialProjectionSubtraction(XmippProtocol):

    def __init__(self, scriptname,project=None):
        
        XmippProtocol.__init__(self, protDict.projsubs.name, scriptname, project)
        self.Import = 'from protocol_subtraction_before_loop import *;\
                       from protocol_subtraction_in_loop import *;'        
        
        self.myName='partial_projection_subtraction'
        self.subtractionDir ='Subtraction'
        self.referenceDir   ='Refs'
        self.referenceStack ='ref'
        self.subImgsDir  = 'SubImgs'
        self.subtractedStack ='subtracted'
        self.volsDir ='Vols'
        self.tempFileName=''
        self.current_angles='current_angles.doc'
        self.scaledImages = 'scaled'
        
        # Check these params
        self.DoDeleteWorkingDir = False
        
    def ImportProtocol(self):

        from protlib_utils import unique_filename
        from protlib_utils import loadModule
        pardir=os.path.abspath(os.getcwd())
        tempFileName=unique_filename(self.ProtocolName)
        fn = FileName(tempFileName)
        fn=fn.getBaseName()
        shutil.copy(self.ProtocolName,fn+'.py')
        
        #FIXME use load module??
        exec  "import " + fn
    
        self.pmprotWorkingDir = eval(fn +'.WorkingDir')
        if (self.SymmetryGroup == ''):
            self.SymmetryGroup    = eval(fn +'.SymmetryGroup')
        AngSamplingRateDeg=getComponentFromVector(eval(fn +'.AngSamplingRateDeg'),self.iterationNo - 1)
        if (len(self.AngSamplingRateDeg) <1):
            self.AngSamplingRateDeg    = AngSamplingRateDeg
        MaxChangeInAngles=float(getComponentFromVector(eval(fn +'.MaxChangeInAngles'),self.iterationNo - 1))
        if (len(self.MaxChangeInAngles) <1):
            self.MaxChangeInAngles    = MaxChangeInAngles
        
    def validate(self):    
        '''This function will be used to validate the protocols
        should be implemented in all derived protocols'''
        errors = []        
        # Check if there is workingdir 
        # move this to gui
        if not os.path.exists(self.ProtocolName):
            errors.append("Refered protocol named %s does not exist"%ProtocolName)
            
        return errors 

    def summary(self):
        super(ProtPartialProjectionSubtraction, self).summary()
        summary = ['Summary'] 
        #summary = ['Performed %d iterations with angular sampling rate %s' 
        #           % (self.NumberOfIterations, self.AngSamplingRateDeg)]
        #summary += ['Final Resolution is %s'%'not yet implemented']
        #summary += ['Number of CTFgroups and References is %d %d respectively'
        #                %(self.NumberOfCtfGroups,self.numberOfReferences)]

        return summary
    
    def preRun(self):

        self.ImportProtocol()
        
        self.Iteration_Working_Directory = os.path.join(self.pmprotWorkingDir,'Iter_'+ str(self.iterationNo))
        self.subtractionDir = self.workingDirPath(self.RunName)
        self.volsDir = self.workingDirPath(self.volsDir)
        self.referenceDir = self.workingDirPath(self.referenceDir)
        self.subImgsDir = self.workingDirPath(self.subImgsDir)
        self.scaledImages = self.workingDirPath(self.scaledImages)
                
        if(self.MaxChangeInAngles > 100):
            self.MaxChangeInAngles=-1
            
        #new vwersion need zfill
        tmpFilename = 'Iter_'+ str(self.iterationNo) + '_' + self.current_angles #zfill()
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
            
        #if not ctf info available use original sel FIXME
        tmpFileName = os.path.join(self.pmprotWorkingDir,'CtfGroups/ctf_group??????.sel')
        self.defGroups=['']
        import glob
        self.defGroups += glob.glob(tmpFileName)
        #if not CTF groups available
        if  len(self.defGroups)<2:
            self.defGroups.append(self.filename_currentAngles) 
            
        self.defocusGroupNo = len(self.defGroups)
        
        if self.DocFileExp!='':
            tmpFileName = self.DocFileExp
        else:
            tmpFileName = FileName(self.filename_currentAngles)
        self.DocFileExp=['']
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
        self.Db.setPrintWrapperParameters(self.PrintWrapperParameters)
        self.Db.setPrintWrapperCommand(self.PrintWrapperCommand)
        self.Db.setVerify(self.Verify,self.ViewVerifyedFiles)
        
        
    def otherActionsToBePerformedBeforeLoop(self):
        _log = self.Log
        _dataBase = self.Db
        #Create directories
        _dataBase.insertAction('createDir', path = self.volsDir)
        _dataBase.insertAction('createDir', path = self.referenceDir)
        _dataBase.insertAction('createDir', path = self.subImgsDir)
        
        if(self.doScaleImages):
            _VerifyFiles = [self.scaledImages+".stk"]
            _VerifyFiles.append(self.scaledImages+".xmd")
            id = _dataBase.insertAction('scaleImages', verifyfiles = _VerifyFiles
                                       , dimX = self.dimX
                                       , dimY = self.dimY
                                       , filename_currentAngles = self.filename_currentAngles
                                       , MpiJobSize = self.MpiJobSize
                                       , NumberOfMpi = self.NumberOfMpi
                                       , NumberOfThreads = self.NumberOfThreads
                                       , scaledImages = self.scaledImages
                                       , SystemFlavour = self.SystemFlavour
                                       )            
    def actionsToBePerformedInsideLoop(self):
        _log = self.Log
        _dataBase = self.Db
        for iterN in range(1, self.defocusGroupNo):
            #Create auxiliary metadata with image names , angles and CTF
            if(self.doScaleImages):
                inputSelfile = self.scaledImages +'.xmd'
            else:
                inputSelfile = self.filename_currentAngles
                            
            if(self.defocusGroupNo > 2 and self.doScaleImages):
                _VerifyFiles = []
                auxFilename = FileName(self.DocFileExp[iterN])
                _VerifyFiles.append(auxFilename.removeBlockName())
                id = self.Db.insertAction('joinImageCTFscale', verifyfiles = _VerifyFiles
                                        , CTFgroupName = self.defGroups[iterN]
                                        , DocFileExp = self.DocFileExp[iterN]
                                        , inputSelfile = inputSelfile
                                        )
            elif (self.defocusGroupNo > 2 and not self.doScaleImages):
                _VerifyFiles = []
                auxFilename = FileName(self.DocFileExp[iterN])
                _VerifyFiles.append(auxFilename.removeBlockName())
                id = self.Db.insertAction('joinImageCTF', verifyfiles = _VerifyFiles
                                        , CTFgroupName = self.defGroups[iterN]
                                        , DocFileExp = self.DocFileExp[iterN]
                                        , inputSelfile = inputSelfile
                                        )                
            else:
                self.DocFileExp[iterN] = self.defGroups[iterN]
            #reconstruct each CTF group
            _VerifyFiles = []
            _VerifyFiles.append(self.reconstructedVolume[iterN])
            id = self.Db.insertAction('reconstructVolume', verifyfiles = _VerifyFiles
                                        , DocFileExp = self.DocFileExp[iterN]
                                        , MpiJobSize = self.MpiJobSize
                                        , NumberOfMpi = self.NumberOfMpi
                                        , NumberOfThreads = self.NumberOfThreads
                                        , reconstructedVolume = self.reconstructedVolume[iterN]
                                        , SymmetryGroup = self.SymmetryGroup
                                        , SystemFlavour = self.SystemFlavour)
            
            #mask volume before projection
            _VerifyFiles = []
            auxFilename = FileName(self.DocFileExp[iterN])
            _VerifyFiles.append(self.maskReconstructedVolume[iterN])
            id = self.Db.insertAction('maskVolume', verifyfiles = _VerifyFiles
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
            
            id = self.Db.insertAction('createProjections', verifyfiles = _VerifyFiles
                                        , AngSamplingRateDeg = self.AngSamplingRateDeg
                                        , DocFileExp = self.DocFileExp[iterN]
                                        , maskReconstructedVolume = self.maskReconstructedVolume[iterN]
                                        , MaxChangeInAngles = self.MaxChangeInAngles
                                        , MpiJobSize = self.MpiJobSize
                                        , NumberOfMpi = self.NumberOfMpi
                                        , NumberOfThreads = self.NumberOfThreads
                                        , referenceStack = self.referenceStack[iterN]
                                        , SymmetryGroup = self.SymmetryGroup
                                        , SystemFlavour = self.SystemFlavour)
                          
                     
                     
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
            
            id = self.Db.insertAction('subtractionScript', verifyfiles = _VerifyFiles
                                        , DocFileExp = self.DocFileExp[iterN] 
                                        , referenceStackDoc = self.referenceStackDoc[iterN]
                                        , subtractedStack = self.subtractedStack[iterN])
                          
    def defineActions(self):
        self.preRun()
        self.otherActionsToBePerformedBeforeLoop()
        self.actionsToBePerformedInsideLoop()
    
