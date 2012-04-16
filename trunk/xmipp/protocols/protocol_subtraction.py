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
from protlib_utils import getComponentFromVector
from protlib_sql import XmippProjectDb
from config_protocols import protDict
from protlib_filesystem import copyFile

class ProtPartialProjectionSubtraction(XmippProtocol):

    def __init__(self, scriptname,project=None):
        
        super(ProtPartialProjectionSubtraction, self).__init__(protDict.subtraction.name, scriptname, project)
        
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
        self.resultsImagesName = 'results_images'
        
        self.localCurrentAngles='sub_current_angles.doc'
        
        # Check these params
        self.DoDeleteWorkingDir = False
        
        self.CtfGroupDirectoryName = "CtfGroups"
        self.CtfGroupRootName =  "ctf"
        self.localStackCTFs = "ctf_ctf.stk"
        self.localDocCTFs = "ctf_images.sel"

        

    def ImportProtocol(self):
        
        print "ImportProtocol"
        
        importProt = self.getProtocolFromRunName(self.ProtocolName) 
#        self.importDir = importProt.WorkingDir
#        self.TiltPairs = importProt.TiltPairs
#        self.importMicroscope = importProt.getFilename('microscope')
#        self.importMicrographs = importProt.getFilename('micrographs')
        
        self.pmprotWorkingDir = importProt.WorkingDir
        if (self.SymmetryGroup == ''):
            self.SymmetryGroup = importProt.SymmetryGroup
            
        if (len(self.AngSamplingRateDeg) <1):
            self.AngSamplingRateDeg    = getComponentFromVector(importProt.AngSamplingRateDeg,self.iterationNo - 1)
            
        if (len(self.MaxChangeInAngles) <1):
            self.MaxChangeInAngles    = float(getComponentFromVector(importProt.MaxChangeInAngles,self.iterationNo - 1))
            

#        self.defocusGroupNo = importProt.NumberOfCtfGroups
        file_name_tmp = join(self.CtfGroupDirectoryName, self.CtfGroupRootName) +'Info.xmd'
        file_name = join(self.pmprotWorkingDir, file_name_tmp)
        print "file_name: ", file_name
                         
        if exists(file_name):
            auxMD = MetaData("numberGroups@"+file_name)
            self.defocusGroupNo = auxMD.getValue(MDL_COUNT,auxMD.firstObject())
        else:
            self.defocusGroupNo = 1


#        from protlib_filesystem import uniqueFilename
#        from protlib_utils import loadModule
#        pardir=os.path.abspath(os.getcwd())
#        tempFileName=uniqueFilename(self.ProtocolName)
#        fn = FileName(tempFileName)
#        fn=fn.getBaseName()
#        shutil.copy(self.ProtocolName,fn+'.py')
        
        #FIXME use load module??
        #exec  "import " + fn
        
#        self.pmprotWorkingDir = eval(fn +'.WorkingDir')
#        if (self.SymmetryGroup == ''):
#            self.SymmetryGroup    = eval(fn +'.SymmetryGroup')
#        AngSamplingRateDeg=getComponentFromVector(eval(fn +'.AngSamplingRateDeg'),self.iterationNo - 1)
#        if (len(self.AngSamplingRateDeg) <1):
#            self.AngSamplingRateDeg    = AngSamplingRateDeg
#        MaxChangeInAngles=float(getComponentFromVector(eval(fn +'.MaxChangeInAngles'),self.iterationNo - 1))
#        if (len(self.MaxChangeInAngles) <1):
#            self.MaxChangeInAngles    = MaxChangeInAngles
        
    def validate(self):    
        errors = []
        #1 call base class, checks if project exists
        super(ProtPartialProjectionSubtraction, self).validate()
        
#        # Check if there is workingdir 
#        # move this to gui
#        if not os.path.exists(self.ProtocolName):
#            errors.append("Refered protocol named %s does not exist" % self.ProtocolName)
            
        return errors 

    def summary(self):
        super(ProtPartialProjectionSubtraction, self).summary()
        summary = ['Summary'] 
        #summary = ['Performed %d iterations with angular sampling rate %s' 
        #           % (self.NumberOfIterations, self.AngSamplingRateDeg)]
        #summary += ['Final Resolution is %s'%'not yet implemented']
        return summary
    
    
    def visualize(self):
    
        plots = [k for k in ['DisplayReference', 'DisplayExperimental', 'DisplaySubtracted'] if self.ParamsDict[k]]
        
        if len(plots):
            self.launchProjSubPlots(plots)
            
            
    def visualizeVar(self, varName):
#        if varName == 'DoShowReferences':
#            self.visualizeReferences()
#        else:
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
        file_name = join(self.pmprotWorkingDir, file_name_tmp)
        print "file_name: ", file_name
                         
        if exists(file_name):
            auxMD = MetaData("numberGroups@"+file_name)
            self.defocusGroupNo = auxMD.getValue(MDL_COUNT,auxMD.firstObject())
        else:
            self.defocusGroupNo = 1
            
        
        if doPlot('DisplayReference'):
            for indexCtf in range(1, self.defocusGroupNo+1): 
                file_name = self.getFilename('SubtractedStack', ctf=indexCtf)+'ref'
                    
                if exists(file_name):
                    try:
                        runShowJ(file_name)
                    except Exception, e:
                        from protlib_gui_ext import showError
                        showError("Error launching java app", str(e))

        if doPlot('DisplayExperimental'):
            for indexCtf in range(1, self.defocusGroupNo+1): 
                file_name = self.getFilename('SubtractedStack', ctf=indexCtf)+'exp'
                    
                if exists(file_name):
                    try:
                        runShowJ(file_name)
                    except Exception, e:
                        from protlib_gui_ext import showError
                        showError("Error launching java app", str(e))

        if doPlot('DisplaySubtracted'):
            for indexCtf in range(1, self.defocusGroupNo+1): 
                file_name = self.getFilename('SubtractedStack', ctf=indexCtf)
                    
                if exists(file_name):
                    try:
                        runShowJ(file_name)
                    except Exception, e:
                        from protlib_gui_ext import showError
                        showError("Error launching java app", str(e))

        if xplotter:
            xplotter.show()


    def createFilenameTemplates(self):  
        
        #Some class variables
#        LibraryDir = "ReferenceLibrary"
        extraParams = {'ReferenceVolumeName': 'reference_volume.vol',
#        'LibraryDir': LibraryDir,
#        'ProjectLibraryRootName': join(LibraryDir, "gallery"),
        'ProjSubDir': "ProjSub",
#        'ProjMatchName': self.Name,
#        'ClassAverageName': 'class_average',
#        #ProjMatchRootName = ProjMatchDir + "/" + ProjMatchName
#        'ForReconstructionSel': "reconstruction.sel",
#        'ForReconstructionDoc': "reconstruction.doc",
#        'MultiAlign2dSel': "multi_align2d.sel",
#        'BlockWithAllExpImages' : 'all_exp_images',
#        'DocFileWithOriginalAngles': 'original_angles.doc',
#        'Docfile_with_current_angles': 'current_angles',
         'LocalCurrentAngles': 'sub_current_angles.doc',
         'LocalCurrentAnglesCtfGroups': 'sub_current_angles_ctfgroups.doc',
#        'FilteredReconstruction': "filtered_reconstruction",
#        'ReconstructedVolume': "reconstruction",
#        'MaskReferenceVolume': "masked_reference",
#        'OutputFsc': "resolution.fsc",
        'CtfGroupDirectory': "CtfGroups",
        'CtfGroupRootName': "ctf",
        #
        'CtfGroupRecRootName':'rec_ctfg',
        'VolsDir': 'Vols',
        'ReferencesDir': 'Refs',
        'SubtractionsDir': 'SubImgs',
        'ScaledXmd': 'scaled.xmd',
        'ScaledStk': 'scaled.stk',
        'CtfGroupSubsetFileName': "ctf_images.sel"
        }
        
        self.ParamsDict.update(extraParams)
        # Set as protocol variables
        for k, v in extraParams.iteritems():
            setattr(self, k, v)
#                              
        #Iter = 'Iter_%(iter)03d'
        Ref3D = 'refGroup%(ref)03d'
        Ctf = 'ctfGroup%(ctf)06d'
        Vol = 'Rec_' + Ctf
#        IterDir = self.workingDirPath(Iter)
#        
#        #ProjMatchDirs = join(IterDir, '%(ProjMatchDir)s.doc')
        #ProjSubDirs = join(IterDir, '%(ProjSubDir)s')
        StackCTFs = join('%(CtfGroupDirectory)s', 'ctf_ctf.stk')
        DocCTFs = join('%(CtfGroupDirectory)s', 'ctf_images.sel')
#        _OutClassesXmd = join(ProjMatchDirs, '%(ProjMatchName)s_' + Ref3D + '.xmd')
#        _OutClassesXmdS1 = join(ProjMatchDirs, '%(ProjMatchName)s_split_1_' + Ref3D + '.xmd')
#        _OutClassesXmdS2 = join(ProjMatchDirs, '%(ProjMatchName)s_split_2_' + Ref3D + '.xmd')
#        CtfGroupBase = join(self.workingDirPath(), self.CtfGroupDirectory, '%(CtfGroupRootName)s')
#        ProjLibRootNames = join(IterDir, '%(ProjectLibraryRootName)s_' + Ref3D)
        return {
                # Global filenames templates
#                'IterDir': IterDir,
#                'ProjSubDirs': ProjSubDirs,
#                'DocfileInputAnglesIters': join(IterDir, '%(Docfile_with_current_angles)s.doc'),
#                'LibraryDirs': join(IterDir, '%(LibraryDir)s'),
#                'ProjectLibraryRootNames': ProjLibRootNames,
#                'ProjMatchRootNames': join(ProjMatchDirs, '%(ProjMatchName)s_' + Ref3D + '.doc'),
#                'ProjMatchRootNamesWithoutRef': join(ProjMatchDirs, '%(ProjMatchName)s.doc'),
#                'OutClasses': join(ProjMatchDirs, '%(ProjMatchName)s'),
#                'OutClassesXmd': _OutClassesXmd,
#                'OutClassesStk': join(ProjMatchDirs, '%(ProjMatchName)s_' + Ref3D + '.stk'),
#                'OutClassesDiscarded': join(ProjMatchDirs, '%(ProjMatchName)s_discarded.xmd'),
#                'ReconstructionXmd': Ref3D + '@' +_OutClassesXmd,
#                'ReconstructionXmdSplit1': Ref3D + '@' +_OutClassesXmdS1,
#                'ReconstructionXmdSplit2': Ref3D + '@' +_OutClassesXmdS2,
#                'MaskedFileNamesIters': join(IterDir, '%(MaskReferenceVolume)s_' + Ref3D + '.vol'),
#                'ReconstructedFileNamesIters': join(IterDir, '%(ReconstructedVolume)s_' + Ref3D + '.vol'),
#            tmpReconstruct = os.path.join(self.volsDir, 'rec_ctfg' + str(iterN).zfill(6) + '.vol')

                'ReconstructedFileNamesIters': self.workingDirPath(join('%(VolsDir)s', Vol + '.vol')),
                'ReconstructedMaskedFileNamesIters': self.workingDirPath(join('%(VolsDir)s', Vol + '_masked.vol')),

#                'ReconstructedFileNamesItersSplit1': join(IterDir, '%(ReconstructedVolume)s_split_1_' + Ref3D + '.vol'),
#                'ReconstructedFileNamesItersSplit2': join(IterDir, '%(ReconstructedVolume)s_split_2_' + Ref3D + '.vol'),
#                'ReconstructedFilteredFileNamesIters': join(IterDir, '%(ReconstructedVolume)s_filtered_' + Ref3D + '.vol'),
#                'ResolutionXmdFile': join(IterDir, '%(ReconstructedVolume)s_' + Ref3D + '_frc.xmd'),
#                'ResolutionXmd': 'resolution@' + join(IterDir, '%(ReconstructedVolume)s_' + Ref3D + '_frc.xmd'),
#                'ResolutionXmdMax': 'resolution_max@' + join(IterDir, '%(ReconstructedVolume)s_' + Ref3D + '_frc.xmd'),
#                'MaskedFileNamesIters': join(IterDir, '%(MaskReferenceVolume)s_' + Ref3D + '.vol'),
#                # Particular templates for executeCtfGroups  
#                'ImageCTFpairs': CtfGroupBase + '_images.sel',
#                'CTFGroupSummary': CtfGroupBase + 'Info.xmd',
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
                'SubtractedStack'  : self.workingDirPath(join('%(SubtractionsDir)s','sub_'+ Ctf + '.stk')),
                'SubtractedDoc'    : self.workingDirPath(join('%(SubtractionsDir)s','sub_'+ Ctf + '.doc'))
#                'StackWienerFilters': CtfGroupBase + '_wien.stk',
#                'SplitAtDefocus': CtfGroupBase + '_split.doc',
#                # Particular templates for angular_project_library 
#                'ProjectLibraryStk': ProjLibRootNames + '.stk',
#                'ProjectLibraryDoc': ProjLibRootNames + '.doc',
#                'ProjectLibrarySampling': ProjLibRootNames + '_sampling.xmd',
#                'ProjectLibraryGroupSampling': ProjLibRootNames + '_group%(group)06d_sampling.xmd',
                }
     
#        file_name = join(self.CtfGroupDirectory, self.CtfGroupRootName) +'Info.xmd'
#        if exists(file_name):
#            auxMD = MetaData("numberGroups@"+file_name)
#            self.NumberOfCtfGroups = auxMD.getValue(MDL_COUNT,auxMD.firstObject())
#        else:        
#            #if not CTF groups available
#            self.NumberOfCtfGroups = 1
#            self.defGroups.append(self.filename_currentAngles)
#        
#        
#        #if not CTF groups available
#        if  len(self.defGroups)<2:
#            self.defGroups.append(self.filename_currentAngles) 
#            
#        self.defocusGroupNo = len(self.defGroups)
#        print 'self.defocusGroupNo: ', self.defocusGroupNo
#        
#        if self.DocFileExp!='':
#            tmpFileName = self.DocFileExp
#        else:
#            tmpFileName = FileName(self.filename_currentAngles)
#        self.DocFileExp=['']
#        tmpFileName = tmpFileName.withoutExtension()
#        tmpFileName = tmpFileName + '_ctfgroups.doc'
#        if os.path.exists(tmpFileName):
#            os.remove(tmpFileName)
#        for iterN in range(1, self.defocusGroupNo ):
#            tmpDocFileExp = 'ctfgroup_' + str(iterN).zfill(6) + '@' + tmpFileName
#            self.DocFileExp.append(tmpDocFileExp)
#            
#        self.reconstructedVolume = [""]
#        for iterN in range(1, self.defocusGroupNo ):
#            tmpReconstruct = os.path.join(self.volsDir, 'rec_ctfg' + str(iterN).zfill(6) + '.vol')
#            self.reconstructedVolume.append(tmpReconstruct)
#            
#        self.maskReconstructedVolume = [""]
#        for iterN in range(1, self.defocusGroupNo ):
#            tmpReconstruct = os.path.join(self.volsDir, 'rec_ctfg' + str(iterN).zfill(6) + '_mask.vol')
#            self.maskReconstructedVolume.append(tmpReconstruct)
#    
#        tmp = self.referenceStack  #'ref'
#        self.referenceStackDoc = [""]
#        for iterN in range(1, self.defocusGroupNo ):
#            tmpReference = os.path.join(self.referenceDir, tmp + '_' + str(iterN).zfill(6) + '.doc')
#            self.referenceStackDoc.append(tmpReference)
#    
#        tmp = self.referenceStack
#        self.referenceStack = [""]
#        for iterN in range(1, self.defocusGroupNo ):
#            tmpReference = os.path.join(self.referenceDir, tmp + '_' + str(iterN).zfill(6) + '.stk')
#            self.referenceStack.append(tmpReference)
#    
#        tmp = self.subtractedStack
#        self.subtractedStack = [""]
#        for iterN in range(1, self.defocusGroupNo ):
#            tmpSubtract = os.path.join(self.subImgsDir, tmp + '_' + str(iterN).zfill(6) + '.stk')
#            self.subtractedStack.append(tmpSubtract)
     
     
     
    def preRun(self):

        self.ImportProtocol()
        
        self.Iteration_Working_Directory = os.path.join(self.pmprotWorkingDir,'Iter_00'+ str(self.iterationNo))
        self.subtractionDir = self.workingDirPath(self.RunName)
        self.volsDir = self.workingDirPath(self.volsDir)
        self.referenceDir = self.workingDirPath(self.referenceDir)
        self.subImgsDir = self.workingDirPath(self.subImgsDir)
        self.scaledImages = self.workingDirPath(self.scaledImages)
        self.resultsImagesName = self.workingDirPath(self.resultsImagesName)
        
        self.localFilenameCurrentAngles = self.workingDirPath(self.localCurrentAngles)
        
        self.CtfGroupDirectory = self.workingDirPath(self.CtfGroupDirectoryName)
        tmpCTFname = join(self.CtfGroupDirectoryName,self.localStackCTFs)
        self.projmatchStackCTFs = join(self.pmprotWorkingDir,tmpCTFname)
        self.localStackCTFs = self.workingDirPath(tmpCTFname)
        
        tmpCTFname = join(self.CtfGroupDirectoryName,self.localDocCTFs)
        self.projmatchDocCTFs = join(self.pmprotWorkingDir,tmpCTFname)
        self.localDocCTFs = self.workingDirPath(tmpCTFname)
                
        if(self.MaxChangeInAngles > 100):
            self.MaxChangeInAngles=-1
            
        #new vwersion need zfill
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
            
        #if not ctf info available use original sel FIXME
#        tmpFileName = os.path.join(self.pmprotWorkingDir,'CtfGroups/ctf_group??????.sel')
#        self.defGroups=['']
#        import glob
#        self.defGroups += glob.glob(tmpFileName)
#        print 'self.defGroups: ', self.defGroups
        
        
#        file_name = join(self.CtfGroupDirectory, self.CtfGroupRootName) +'Info.xmd'
#        if exists(file_name):
#            auxMD = MetaData("numberGroups@"+file_name)
#            self.NumberOfCtfGroups = auxMD.getValue(MDL_COUNT,auxMD.firstObject())
#        else:        
#            #if not CTF groups available
#            self.NumberOfCtfGroups = 1
#            self.defGroups.append(self.filename_currentAngles)
#        
#        
#        #if not CTF groups available
#        if  len(self.defGroups)<2:
#            self.defGroups.append(self.filename_currentAngles) 
#            
#        self.defocusGroupNo = len(self.defGroups)
#        print 'self.defocusGroupNo: ', self.defocusGroupNo
#        
#        if self.DocFileExp!='':
#            tmpFileName = self.DocFileExp
#        else:
#            tmpFileName = FileName(self.filename_currentAngles)
#        self.DocFileExp=['']
#        tmpFileName = tmpFileName.withoutExtension()
#        tmpFileName = tmpFileName + '_ctfgroups.doc'
#        if os.path.exists(tmpFileName):
#            os.remove(tmpFileName)
#        for iterN in range(1, self.defocusGroupNo ):
#            tmpDocFileExp = 'ctfgroup_' + str(iterN).zfill(6) + '@' + tmpFileName
#            self.DocFileExp.append(tmpDocFileExp)
#            
#        self.reconstructedVolume = [""]
#        for iterN in range(1, self.defocusGroupNo ):
#            tmpReconstruct = os.path.join(self.volsDir, 'rec_ctfg' + str(iterN).zfill(6) + '.vol')
#            self.reconstructedVolume.append(tmpReconstruct)
#            
#        self.maskReconstructedVolume = [""]
#        for iterN in range(1, self.defocusGroupNo ):
#            tmpReconstruct = os.path.join(self.volsDir, 'rec_ctfg' + str(iterN).zfill(6) + '_mask.vol')
#            self.maskReconstructedVolume.append(tmpReconstruct)
#    
#        tmp = self.referenceStack  #'ref'
#        self.referenceStackDoc = [""]
#        for iterN in range(1, self.defocusGroupNo ):
#            tmpReference = os.path.join(self.referenceDir, tmp + '_' + str(iterN).zfill(6) + '.doc')
#            self.referenceStackDoc.append(tmpReference)
#    
#        tmp = self.referenceStack
#        self.referenceStack = [""]
#        for iterN in range(1, self.defocusGroupNo ):
#            tmpReference = os.path.join(self.referenceDir, tmp + '_' + str(iterN).zfill(6) + '.stk')
#            self.referenceStack.append(tmpReference)
#    
#        tmp = self.subtractedStack
#        self.subtractedStack = [""]
#        for iterN in range(1, self.defocusGroupNo ):
#            tmpSubtract = os.path.join(self.subImgsDir, tmp + '_' + str(iterN).zfill(6) + '.stk')
#            self.subtractedStack.append(tmpSubtract)
        #FIXME
        #DO we need to import this paths or are already in pythonpath
        #import shutil,time
#        scriptdir = os.path.split(os.path.dirname(os.popen('which xmipp_protocols', 'r').read()))[0] + '/protocols'
#        sys.path.append(scriptdir)

        # Configure dabase
        ##self.Db.setVerify(self.Verify,self.ViewVerifyedFiles)
        ##self.Db.setParentDefault(XmippProjectDb.lastStep)
        
        
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
        print "self.defocusGroupNo: ", self.defocusGroupNo
        for iterN in range(1, self.defocusGroupNo+1):
            print "iterN: ", iterN
            #Create auxiliary metadata with image names , angles and CTF
            if(self.doScaleImages):
                inputSelfile = self.scaledImages +'.xmd'
            else:
                inputSelfile = self.filename_currentAngles
                            
            if(self.defocusGroupNo > 1 and self.doScaleImages):
                _VerifyFiles = []
                auxFilename = FileName(self.getFilename('SubCurrentAnglesCftGroups', ctf=iterN))
                _VerifyFiles.append(auxFilename.removeBlockName())
                id = self.Db.insertStep('joinImageCTFscale', verifyfiles = _VerifyFiles
                                        , CTFgroupName = self.getFilename('DocCTFsBlocks', ctf=iterN)
                                        , DocFileExp = self.getFilename('SubCurrentAnglesCftGroups', ctf=iterN)
                                        , inputSelfile = inputSelfile
                                        )
            elif (self.defocusGroupNo > 1 and not self.doScaleImages):
                _VerifyFiles = []
                auxFilename = FileName(self.getFilename('SubCurrentAnglesCftGroups', ctf=iterN))
                _VerifyFiles.append(auxFilename.removeBlockName())
                id = self.Db.insertStep('joinImageCTF', verifyfiles = _VerifyFiles
                                        , CTFgroupName =  self.getFilename('DocCTFsBlocks', ctf=iterN)
                                        , DocFileExp = self.getFilename('SubCurrentAnglesCftGroups', ctf=iterN)
                                        , inputSelfile = inputSelfile
                                        )                
#            else: # No Ctf correction in ProjMatch
#                self.DocFileExp[iterN] =  self.getFilename('StackCTFs', ctf=iterN)

            #reconstruct each CTF group
            _VerifyFiles = []
            _VerifyFiles.append(self.getFilename('ReconstructedFileNamesIters', ctf=iterN))
            id = self.Db.insertStep('reconstructVolume', verifyfiles = _VerifyFiles
                                        , DocFileExp = self.getFilename('SubCurrentAnglesCftGroups', ctf=iterN)
                                        , MpiJobSize = self.MpiJobSize
                                        , NumberOfMpi = self.NumberOfMpi
                                        , NumberOfThreads = self.NumberOfThreads
                                        , reconstructedVolume = self.getFilename('ReconstructedFileNamesIters', ctf=iterN)
                                        , SymmetryGroup = self.SymmetryGroup)
            
            #mask volume before projection
            _VerifyFiles = []
            auxFilename = FileName(self.getFilename('SubCurrentAnglesCftGroups', ctf=iterN))
#            _VerifyFiles.append(self.maskReconstructedVolume[iterN])
            id = self.Db.insertStep('maskVolume', verifyfiles = _VerifyFiles
                                        , dRradiusMax = self.dRradiusMax
                                        , dRradiusMin = self.dRradiusMin
                                        , maskReconstructedVolume = self.getFilename('ReconstructedMaskedFileNamesIters', ctf=iterN)
                                        , reconstructedVolume = self.getFilename('ReconstructedFileNamesIters', ctf=iterN)
                                        , NumberOfMpi = self.NumberOfMpi
                                        , NumberOfThreads = self.NumberOfThreads)
    
            #project reconstructe4d volumes
            _VerifyFiles = []
#            _VerifyFiles.append(self.getFilename('ReferenceStack', ctf=iterN))
#            tmp = self.referenceStack[iterN]
#            _VerifyFiles.append(tmp.replace('.stk','.doc'))
#            _VerifyFiles.append(tmp.replace('.stk','_sampling.xmd'))
            
            id = self.Db.insertStep('createProjections', verifyfiles = _VerifyFiles
                                        , AngSamplingRateDeg = self.AngSamplingRateDeg
                                        , DocFileExp = self.getFilename('SubCurrentAnglesCftGroups', ctf=iterN)
                                        , maskReconstructedVolume = self.getFilename('ReconstructedMaskedFileNamesIters', ctf=iterN)
                                        , MaxChangeInAngles = self.MaxChangeInAngles
                                        , MpiJobSize = self.MpiJobSize
                                        , NumberOfMpi = self.NumberOfMpi
                                        , NumberOfThreads = self.NumberOfThreads
                                        , referenceStack = self.getFilename('ReferenceStack', ctf=iterN)
                                        , SymmetryGroup = self.SymmetryGroup)
                          
                     
                     
    #project reconstructe4d volumes
            _Parameters = {
                          'DocFileExp'        : self.getFilename('SubCurrentAnglesCftGroups', ctf=iterN)
                          ,'referenceStackDoc': self.getFilename('ReferenceStackDoc', ctf=iterN)
                          ,'subtractedStack'  : self.getFilename('SubtractedStack', ctf=iterN)
                          }
            command = "subtractionScript"
            _VerifyFiles = []
            _VerifyFiles.append(self.getFilename('SubtractedStack', ctf=iterN))
            _VerifyFiles.append(self.getFilename('SubtractedStack', ctf=iterN)+'ref')
            _VerifyFiles.append(self.getFilename('SubtractedStack', ctf=iterN)+'exp')
            
            id = self.Db.insertStep('subtractionScript', verifyfiles = _VerifyFiles
                                        , DocFileExp = self.getFilename('SubCurrentAnglesCftGroups', ctf=iterN) 
                                        , referenceStackDoc = self.getFilename('ReferenceStackDoc', ctf=iterN)
                                        , subtractedStack = self.getFilename('SubtractedStack', ctf=iterN)
                                        , resultsImagesName = self.resultsImagesName)
                          
            self.Db.connection.commit()

#    def defineActions(self):
    def defineSteps(self):
        self.preRun()
        self.otherActionsToBePerformedBeforeLoop()
        self.actionsToBePerformedInsideLoop()
    
