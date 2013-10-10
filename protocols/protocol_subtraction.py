#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Xmipp protocol for partial projection subtraction
#
# Example use:
# ./protocol_subtraction_header.py
#
# Authors: Roberto Marabini,
#          Alejandro Echeverria Rey.
#
        
import os, shutil, sys
from os.path import join, exists
from xmipp import FileName, MetaData, Image
from xmipp import DATA, ALL_IMAGES
from xmipp import MDL_COUNT

from protlib_base import XmippProtocol, protocolMain
#from protlib_utils import getListFromVector, getBoolListFromVector, getComponentFromVector
from protlib_utils import getComponentFromVector, runShowJ
from protlib_sql import XmippProjectDb
from config_protocols import protDict
from protlib_filesystem import copyFile
from protlib_gui_ext import showError

class ProtPartialProjectionSubtraction(XmippProtocol):
    def __init__(self, scriptname,project=None):
        
        super(ProtPartialProjectionSubtraction, self).__init__(protDict.subtraction.name, scriptname, project)

        self.Import = 'from protocol_subtraction_before_loop import *;\
                       from protocol_subtraction_in_loop import *;'        
        
        self.setPreviousRun(self.ImportRun)
        
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
        self.resultsImagesName = 'results_images'
        
        self.localCurrentAngles='sub_current_angles.doc'
        
        # Check these params
        self.DoDeleteWorkingDir = False
        
        self.CtfGroupDirectoryName = "CtfGroups"
        self.CtfGroupRootName =  "ctf"
        self.localStackCTFs = "ctf_ctf.stk"
        self.localDocCTFs = "ctf_images.sel"
        

    def ImportProtocol(self):
        
        self.pmprotWorkingDir = self.PrevRun.WorkingDir
#        if (self.SymmetryGroup == ''):
#            self.SymmetryGroup = self.PrevRun.SymmetryGroup
            
#        if (len(self.AngSamplingRateDeg) <1):
#            self.AngSamplingRateDeg    = getComponentFromVector(self.PrevRun.AngSamplingRateDeg,self.iterationNo - 1)
            
#        if (len(self.MaxChangeInAngles) <1):
#            self.MaxChangeInAngles    = float(getComponentFromVector(self.PrevRun.MaxChangeInAngles,self.iterationNo - 1))
            
        file_name_tmp = join(self.CtfGroupDirectoryName, self.CtfGroupRootName) +'Info.xmd'
        file_name = join(self.PrevRun.WorkingDir, file_name_tmp)
                         
        if exists(file_name):
            auxMD = MetaData("numberGroups@"+file_name)
            self.defocusGroupNo = auxMD.getValue(MDL_COUNT,auxMD.firstObject())
        else:
            self.defocusGroupNo = 1

        
    def validate(self):    
        errors = []
        #1 call base class, checks if project exists
        super(ProtPartialProjectionSubtraction, self).validate()
        #create auxiliary file names and check if the files exists
        #ReferenceFileName
        prot = self.project.getProtocolFromRunName(self.ImportRun)
        remoteWorkingDir  = join(prot.WorkingDir,"Iter_%03d"%(self.iterationNo))
        referenceName     = "reconstruction_filtered_Ref3D_%03d.vol"%(self.reference3DNo)
        self.ReferenceFileNames =  join(remoteWorkingDir,referenceName)
        if not exists(self.ReferenceFileNames):
            errors.append('Reference FileName = %s does not exists'%self.ReferenceFileName)
        #DocFileExp: Iter_(X-1)/Iter_(X-1)_current_angles.doc
        #remember to filter the metadata for ref3D 
        docFileName     = "current_angles.doc"
        self.DocFileExp =  join(remoteWorkingDir,docFileName)
        if not exists(self.DocFileExp):
            errors.append('DocFileExp FileName = %s does not exists'%self.DocFileExp)
        return errors 

    def summary(self):
        super(ProtPartialProjectionSubtraction, self).summary()
        summary = [] 
        
        file_name_tmp = join(self.CtfGroupDirectoryName, self.CtfGroupRootName) +'Info.xmd'
        file_name = ""+join(self.PrevRun.WorkingDir, file_name_tmp)
                         
        if exists(file_name):
            auxMD = MetaData("numberGroups@"+file_name)
            self.NumberOfCtfGroups = auxMD.getValue(MDL_COUNT,auxMD.firstObject())
        else:
            self.NumberOfCtfGroups = 1

        summary += ['Number of CTFgroups is %d ' %(self.NumberOfCtfGroups)]

        return summary
    
    
    def visualize(self):
    
        plots = [k for k in ['DisplayReference', 'DisplayExperimental', 'DisplaySubtracted'] if self.ParamsDict[k]]
        
        if len(plots):
            self.launchProjSubPlots(plots)
            
            
    def visualizeVar(self, varName):

        self.launchProjSubPlots([varName])
        
        
    def launchProjSubPlots(self, selectedPlots):
        ''' Launch some plots for a Projection Matching protocol run '''
        import numpy as np
        _log = self.Log

        xplotter=None
        self._plot_count = 0
        
        def doPlot(plotName):
            return plotName in selectedPlots
        
        file_name_tmp = join(self.CtfGroupDirectoryName, self.CtfGroupRootName) +'Info.xmd'
        file_name = join(self.PrevRun.WorkingDir, file_name_tmp)
                         
        if exists(file_name):
            auxMD = MetaData("numberGroups@"+file_name)
            self.NumberOfCtfGroups = auxMD.getValue(MDL_COUNT,auxMD.firstObject())
        else:
            self.NumberOfCtfGroups = 1
            
        
        if doPlot('DisplayReference'):
            for indexCtf in range(1, self.NumberOfCtfGroups+1): 
                file_name = self.getFilename('SubtractedStack', ctf=indexCtf)+'ref'
                    
                if exists(file_name):
                    try:
                        runShowJ(file_name)
                    except Exception, e:
                        showError("Error launching java app", str(e), self.master)

        if doPlot('DisplayExperimental'):
            for indexCtf in range(1, self.NumberOfCtfGroups+1): 
                file_name = self.getFilename('SubtractedStack', ctf=indexCtf)+'exp'
                    
                if exists(file_name):
                    try:
                        runShowJ(file_name)
                    except Exception, e:
                        showError("Error launching java app", str(e), self.master)

        if doPlot('DisplaySubtracted'):
            for indexCtf in range(1, self.NumberOfCtfGroups+1): 
                file_name = self.getFilename('SubtractedStack', ctf=indexCtf)
                print "DisplaySubtracted file_name: ", file_name
                if exists(file_name):
                    try:
                        runShowJ(file_name)
                    except Exception, e:
                        showError("Error launching java app", str(e), self.master)

        if xplotter:
            xplotter.show()


    def createFilenameTemplates(self):  
        
        extraParams = {'ReferenceVolumeName': 'reference_volume.vol',
        'ProjSubDir': "ProjSub",
         'LocalCurrentAngles': 'sub_current_angles.doc',
         'LocalCurrentAnglesCtfGroups': 'sub_current_angles_ctfgroups.doc',
        'CtfGroupDirectory': "CtfGroups",
        'CtfGroupRootName': "ctf",
        'CtfGroupRecRootName':'rec_ctfg',
        'VolsDir': 'Vols',
        'ReferencesDir': 'Refs',
        'SubtractionsDir': 'SubImgs',
        'ScaledXmd': 'scaled.xmd',
        'ScaledStk': 'scaled.stk',
        'CtfGroupSubsetFileName': "ctf_images.sel"
        }
        
        self.ParamsDict.update(extraParams)
        for k, v in extraParams.iteritems():
            setattr(self, k, v)
#                              
        Ref3D = 'refGroup%(ref)03d'
        Ctf = 'ctfGroup%(ctf)06d'
        Vol = 'Rec_' + Ctf
        StackCTFs = join('%(CtfGroupDirectory)s', 'ctf_ctf.stk')
        DocCTFs = join('%(CtfGroupDirectory)s', 'ctf_images.sel')

        return {
                'ReconstructedFileNamesIters': self.workingDirPath(join('%(VolsDir)s', Vol + '.vol')),
                'ReconstructedMaskedFileNamesIters': self.workingDirPath(join('%(VolsDir)s', Vol + '_masked.vol')),
                'DocCTFsBlocks' : Ctf + '@' + self.workingDirPath(DocCTFs),
                'StackCTFsBlocks' : Ctf + '@' + self.workingDirPath(StackCTFs),
                'ScaledXmdBlocks': self.workingDirPath('%(ScaledXmd)s'),
                'ScaledXmdBlocks': self.workingDirPath('%(ScaledXmd)s'),
                'SubCurrentAngles': Ctf + '@' + self.workingDirPath('%(LocalCurrentAngles)s'),
                'SubCurrentAnglesAllExpImgs': 'all_exp_images@' + self.workingDirPath('%(LocalCurrentAngles)s'),
                'SubCurrentAnglesCftGroups': Ctf + '@' + self.workingDirPath('%(LocalCurrentAnglesCtfGroups)s'),
                'SubCurrentAnglesCftGroupsAllExpImgs': 'all_exp_images@' + self.workingDirPath('%(LocalCurrentAnglesCtfGroups)s'),
                'ReferenceStack'   : self.workingDirPath(join('%(ReferencesDir)s','ref_'+ Ctf + '.stk')),
                'ReferenceStackDoc': self.workingDirPath(join('%(ReferencesDir)s','ref_'+ Ctf + '.doc')),
                'ReferenceStackSampling': self.workingDirPath(join('%(ReferencesDir)s','ref_'+ Ctf + '_sampling.xmd')),
                'SubtractedStack'  : self.workingDirPath(join('%(SubtractionsDir)s','sub_'+ Ctf + '.stk')),
                'SubtractedDoc'    : self.workingDirPath(join('%(SubtractionsDir)s','sub_'+ Ctf + '.doc'))
                }
     
    def preRun(self):

        self.ImportProtocol()
        
        self.Iteration_Working_Directory = os.path.join(self.pmprotWorkingDir,'Iter_00'+ str(self.iterationNo))
        self.subtractionDir = self.workingDirPath(self.RunName)
        self.volsDir        = self.workingDirPath(self.volsDir)
        self.referenceDir   = self.workingDirPath(self.referenceDir)
        self.subImgsDir     = self.workingDirPath(self.subImgsDir)
        self.scaledImages   = self.workingDirPath(self.scaledImages)
        self.resultsImagesName = self.workingDirPath(self.resultsImagesName)
        
        self.localFilenameCurrentAngles = self.workingDirPath(self.localCurrentAngles)
        
        self.CtfGroupDirectory  = self.workingDirPath(self.CtfGroupDirectoryName)
        tmpCTFname              = join(self.CtfGroupDirectoryName,self.localStackCTFs)
        self.projmatchStackCTFs = join(self.pmprotWorkingDir,tmpCTFname)
        self.localStackCTFs     = self.workingDirPath(tmpCTFname)
        
        tmpCTFname            = join(self.CtfGroupDirectoryName,self.localDocCTFs)
        self.projmatchDocCTFs = join(self.pmprotWorkingDir,tmpCTFname)
        self.localDocCTFs     = self.workingDirPath(tmpCTFname)
                
        if(self.MaxChangeInAngles > 100):
            self.MaxChangeInAngles=-1
            
        tmpFilename = self.current_angles
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
            
        
    def otherActionsToBePerformedBeforeLoop(self):
        _log = self.Log
        _dataBase = self.Db
        #Create directories
        _dataBase.insertStep('createDir', path = self.volsDir)
        _dataBase.insertStep('createDir', path = self.referenceDir)
        _dataBase.insertStep('createDir', path = self.subImgsDir)
        _dataBase.insertStep('createDir', path = self.CtfGroupDirectory)
#        copyFile(_log, self.filename_currentAngles, self.localFilenameCurrentAngles)
        _dataBase.insertStep('copyFile',source=self.filename_currentAngles, dest=self.localFilenameCurrentAngles)
        
        if(self.defocusGroupNo > 1):
            _dataBase.insertStep('copyFile',source=self.projmatchStackCTFs, dest=self.localStackCTFs)
            _dataBase.insertStep('copyFile',source=self.projmatchDocCTFs, dest=self.localDocCTFs)

        if(self.doScaleImages):
            _VerifyFiles = [self.scaledImages+".stk"]
            _VerifyFiles.append(self.scaledImages+".xmd")
            id = _dataBase.insertStep('scaleImages', verifyfiles = _VerifyFiles
                                       , dimX = self.dimX
                                       , dimY = self.dimY
                                       , filename_currentAngles = self.getFilename('SubCurrentAnglesAllExpImgs')
                                       , MpiJobSize = self.MpiJobSize
                                       , NumberOfMpi = self.NumberOfMpi
                                       , NumberOfThreads = self.NumberOfThreads
                                       , scaledImages = self.scaledImages
                                       )            
        _dataBase.connection.commit()


    def actionsToBePerformedInsideLoop(self):
        _log = self.Log
        _dataBase = self.Db
        for indexCTFInsideLoop in range(1, self.defocusGroupNo+1):
            
            if(self.doScaleImages):
                inputSelfile = self.scaledImages +'.xmd'
            else:
                inputSelfile = self.getFilename('SubCurrentAnglesAllExpImgs')
                            
            if(self.DocFileExp == ''):
                DocFileExpAux = self.getFilename('SubCurrentAnglesCftGroups', ctf=indexCTFInsideLoop)
            else:
                DocFileExpAux = self.DocFileExp
                    
            if(self.defocusGroupNo > 1):
                
                
                if(self.doScaleImages):
                    _VerifyFiles = []
                    auxFilename = FileName(self.getFilename('SubCurrentAnglesCftGroups', ctf=indexCTFInsideLoop))
                    _VerifyFiles.append(auxFilename.removeBlockName())
                    id = self.Db.insertStep('joinImageCTFscale', verifyfiles = _VerifyFiles
                                            , CTFgroupName = self.getFilename('DocCTFsBlocks', ctf=indexCTFInsideLoop)
                                            , DocFileExp = DocFileExpAux
                                            , inputSelfile = inputSelfile
                                            )
                else: #No scale
                    _VerifyFiles = []
                    auxFilename = FileName(self.getFilename('SubCurrentAnglesCftGroups', ctf=indexCTFInsideLoop))
                    _VerifyFiles.append(auxFilename.removeBlockName())
                    id = self.Db.insertStep('joinImageCTF', verifyfiles = _VerifyFiles
                                            , CTFgroupName =  self.getFilename('DocCTFsBlocks', ctf=indexCTFInsideLoop)
                                            , DocFileExp = DocFileExpAux
                                            , inputSelfile = inputSelfile
                                            )                
            else: # No Ctf correction in ProjMatch 
                if(self.doScaleImages):
                    _VerifyFiles = []
                    id = self.Db.insertStep('joinImageCTFscale', verifyfiles = _VerifyFiles
                                            , CTFgroupName = self.localFilenameCurrentAngles
                                            , DocFileExp = DocFileExpAux
                                            , inputSelfile = inputSelfile
                                            )
                else:
                    _VerifyFiles = []
                    id = self.Db.insertStep('joinImageCTF', verifyfiles = _VerifyFiles
                                            , CTFgroupName = self.localFilenameCurrentAngles
                                            , DocFileExp = DocFileExpAux
                                            , inputSelfile = inputSelfile
                                            )

                
                
            #reconstruct each CTF group
            _VerifyFiles = []
            _VerifyFiles.append(self.getFilename('ReconstructedFileNamesIters', ctf=indexCTFInsideLoop))
            id = self.Db.insertStep('reconstructVolume', verifyfiles = _VerifyFiles
                                        , DocFileExp = DocFileExpAux
                                        , MpiJobSize = self.MpiJobSize
                                        , NumberOfMpi = self.NumberOfMpi
                                        , NumberOfThreads = self.NumberOfThreads
                                        , reconstructedVolume = self.getFilename('ReconstructedFileNamesIters', ctf=indexCTFInsideLoop)
                                        , SymmetryGroup = self.SymmetryGroup)
            
            
            #mask volume before projection
            _VerifyFiles = []
            auxFilename = FileName(self.getFilename('SubCurrentAnglesCftGroups', ctf=indexCTFInsideLoop))
            _VerifyFiles.append(self.getFilename('ReconstructedMaskedFileNamesIters', ctf=indexCTFInsideLoop))
            
            id = self.Db.insertStep('maskVolume', verifyfiles = _VerifyFiles
                                        , dRradiusMax = self.dRradiusMax
                                        , dRradiusMin = self.dRradiusMin
                                        , maskReconstructedVolume = self.getFilename('ReconstructedMaskedFileNamesIters', ctf=indexCTFInsideLoop)
                                        , reconstructedVolume = self.getFilename('ReconstructedFileNamesIters', ctf=indexCTFInsideLoop)
                                        , NumberOfMpi = self.NumberOfMpi
                                        , NumberOfThreads = self.NumberOfThreads)
            
            
            #project reconstructe4d volumes
            _VerifyFiles = []
            _VerifyFiles.append(self.getFilename('ReferenceStack', ctf=indexCTFInsideLoop))
            _VerifyFiles.append(self.getFilename('ReferenceStackDoc', ctf=indexCTFInsideLoop))
            _VerifyFiles.append(self.getFilename('ReferenceStackSampling', ctf=indexCTFInsideLoop))
            
            id = self.Db.insertStep('createProjections', verifyfiles = _VerifyFiles
                                        , AngSamplingRateDeg = self.AngSamplingRateDeg
                                        , DocFileExp = DocFileExpAux
                                        , maskReconstructedVolume = self.getFilename('ReconstructedMaskedFileNamesIters', ctf=indexCTFInsideLoop)
                                        , MaxChangeInAngles = self.MaxChangeInAngles
                                        , MpiJobSize = self.MpiJobSize
                                        , NumberOfMpi = self.NumberOfMpi
                                        , NumberOfThreads = self.NumberOfThreads
                                        , referenceStack = self.getFilename('ReferenceStack', ctf=indexCTFInsideLoop)
                                        , SymmetryGroup = self.SymmetryGroup)
                     
                     
            #project reconstructe4d volumes
            _VerifyFiles = []
            _VerifyFiles.append(self.getFilename('SubtractedStack', ctf=indexCTFInsideLoop))
            _VerifyFiles.append(self.getFilename('SubtractedStack', ctf=indexCTFInsideLoop)+'ref')
            _VerifyFiles.append(self.getFilename('SubtractedStack', ctf=indexCTFInsideLoop)+'exp')
            
            id = self.Db.insertStep('subtractionScript', verifyfiles = _VerifyFiles
                                        , DocFileExp = DocFileExpAux
                                        , referenceStackDoc = self.getFilename('ReferenceStackDoc', ctf=indexCTFInsideLoop)
                                        , subtractedStack = self.getFilename('SubtractedStack', ctf=indexCTFInsideLoop)
                                        , resultsImagesName = self.resultsImagesName)
                          
            self.Db.connection.commit()

    def defineSteps(self):
        self.preRun()
        self.otherActionsToBePerformedBeforeLoop()
        self.actionsToBePerformedInsideLoop()
    
