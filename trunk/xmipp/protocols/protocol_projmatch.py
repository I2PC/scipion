#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Xmipp protocol for projection matching
#
# Example use:
# ./xmipp_protocol_projmatch.py
#
# Authors: Roberto Marabini,
#          Sjors Scheres,    March 2008
#        Rewritten by Roberto Marabini
#

from protlib_base import *
from xmipp import *
from protlib_utils import getListFromVector
import os


class ProtProjMatch(XmippProtocol):

#    def __init__(self, scriptname, workingdir, projectdir=None, logdir='Logs', restartStep=1, isIter=True):
    def __init__(self, scriptname, project=None):
        super(ProtProjMatch,self).__init__(protDict.projmatch.key,scriptname, RunName, project,NumberOfMpiProcesses)
        #Some class variables
        self.ReferenceVolumeName = 'reference_volume.vol'
        self.LibraryDir = "ReferenceLibrary"
        self.ProjectLibraryRootName = self.LibraryDir + "/gallery"
        self.ProjMatchDir = "ProjMatchClasses"
        self.ProjMatchName = 'proj_match'
        self.ClassAverageName = 'class_average'
        #ProjMatchRootName = ProjMatchDir + "/" + ProjMatchName
        self.ForReconstructionSel = "reconstruction.sel"
        self.ForReconstructionDoc = "reconstruction.doc"
        self.MultiAlign2dSel = "multi_align2d.sel"
        self.DocFileWithOriginalAngles = 'original_angles.doc'
        self.docfile_with_current_angles = 'current_angles'
        self.FilteredReconstruction = "filtered_reconstruction"
        
        self.ReconstructedVolume = "reconstruction"#
        self.maskReferenceVolume = "masked_reference"#
        
        self.OutputFsc = "resolution.fsc"
        self.CtfGroupDirectory = "CtfGroups"
        self.CtfGroupRootName = "ctf"
        self.CtfGroupSubsetFileName = self.CtfGroupRootName + "_images.sel"
        
        self.reconstructedFileNamesIters = []# names for reconstructed volumes
        #maskedFileNamesIter = []# names masked volumes used as reference
        self.numberOfReferences = 1#number of references
        self.createAuxTable = False
        self.NumberOfCtfGroups = 1
        self.Import = 'from protocol_projmatch_before_loop import *;\
                       from protocol_projmatch_in_loop import *;'
        self.runName=RunName
        self.DoDeleteWorkingDir = DoDeleteWorkingDir
        #self.WorkingDir = os.path.join(self.Name,_runName)

        
    def validate(self):
        errors = []
        #1 call base class, checks if project exists
        super(ProtProjMatch, self).validate()
        
        #2 Check reference and projection size match
        _ReferenceFileNames = getListFromVector(ReferenceFileNames)
        _Parameters = {
              'ReferenceFileNames':_ReferenceFileNames
            , 'SelFileName':SelFileName
            }
        from protocol_projmatch_before_loop import checkVolumeProjSize
        _retval, _error_message = checkVolumeProjSize(None,**_Parameters)
        if(not _retval):
            errors.append(_error_message)
    
    
        # 3 Never allow DoAlign2D and DoCtfCorrection together
        if (int(getComponentFromVector(DoAlign2D,1))==1 and DoCtfCorrection):
            errors.append("You cannot realign classes AND perform CTF-correction. Switch either of them off!")
    
        #4N outter radius is compulsory
        _OuterRadius = getComponentFromVector(OuterRadius,1)
        _InnerRadius = getComponentFromVector(InnerRadius,1)
        if _OuterRadius <= _InnerRadius:
            errors.append("OuterRadius must be larger than InnerRadius")
        #not clear if this a need to return anything since is part of the class
        #FIXME
        return errors 
    
    def summary(self):
        super(ProtProjMatch, self).summary()
        summary = ['Performed %d iterations with angular sampling rate %s' 
                   % (NumberOfIterations, AngSamplingRateDeg)]
        summary += ['Final Resolution is %s'%'not yet implemented']
        summary += ['Number of CTFgroups and References is %d %d respectively'
                        %(self.NumberOfCtfGroups,self.numberOfReferences)]

        return summary
    
    def preRun(self):
        print "in PRERUN"
        #Convert directories/files  to absolute path from projdir
        self.CtfGroupDirectory = os.path.join(self.WorkingDir, self.CtfGroupDirectory)
        self.CtfGroupSubsetFileName = os.path.join(self.CtfGroupDirectory, self.CtfGroupSubsetFileName)
    #vector for iterations??????
    #    global ProjMatchDir
    #    ProjMatchDir = WorkingDir +'/' + ProjMatchDir
        self.DocFileWithOriginalAngles = os.path.join(self.WorkingDir, self.DocFileWithOriginalAngles)
    
        
        # Convert vectors to list
        self.ReferenceFileNames = getListFromVector(ReferenceFileNames)
        self.numberOfReferences = len(self.ReferenceFileNames)
        #directory with ProjMatchClasses
        self.ProjMatchDirs=[" "]
        self.LibraryDirs=[" "]
        self.DocFileInputAngles=[self.DocFileWithOriginalAngles]
        #ProjMatchRootName=[" "]
        
        for iterN in range(NumberOfIterations):
            fnBaseIter = "%s/Iter_%02d/" % (self.WorkingDir, iterN + 1)
            self.ProjMatchDirs.append(fnBaseIter + self.ProjMatchDir)
            self.LibraryDirs.append( fnBaseIter + self.LibraryDir)
            self.DocFileInputAngles.append("%s%s.doc" % (fnBaseIter, self.docfile_with_current_angles))
        
        auxList = (self.numberOfReferences + 1) * [None]
        self.ProjectLibraryRootNames=[[None]]
        for iterN in range(NumberOfIterations):
            fnBaseIter = "%s/Iter_%02d/" % (self.WorkingDir, iterN + 1)
            for refN in range(self.numberOfReferences):                
                auxList[refN + 1]= "%s%s_ref_%02d.stk" % (fnBaseIter, self.ProjectLibraryRootName, refN)
            self.ProjectLibraryRootNames.append(list(auxList))
                    
        self.ProjMatchRootNames=[[None]]
        for iterN in range(NumberOfIterations):
            for refN in range(self.numberOfReferences):
                auxList[refN + 1]="%s/%s_ref_%02d.doc" % (self.ProjMatchDirs[iterN + 1], self.ProjMatchName, refN + 1)
            self.ProjMatchRootNames.append(list(auxList))
    
    
        #name of masked volumes
        #add dummy name so indexes start a 1
        self.maskedFileNamesIters = [[None]]
        for iterN in range(NumberOfIterations):
            fnBaseIter = "%s/Iter_%02d/" % (self.WorkingDir, iterN + 1)
            for refN in range(self.numberOfReferences):
                auxList[refN + 1] = "%s%s_ref_%02d.vol" % (fnBaseIter, self.maskReferenceVolume, refN + 1)
            self.maskedFileNamesIters.append(list(auxList))
    
        ####################################################################
        #add initial reference, useful for many routines
        #NOTE THAT INDEXES START AT 1
        self.reconstructedFileNamesIters.append([None] + self.ReferenceFileNames)
        for iterN in range(NumberOfIterations):
            fnBaseIter = "%s/Iter_%02d/" % (self.WorkingDir, iterN + 1)
            for refN in range(self.numberOfReferences):
                auxList[refN + 1] = "%s%s_ref_%02d.vol" % (fnBaseIter, self.ReconstructedVolume, refN + 1)
            self.reconstructedFileNamesIters.append(list(auxList))
    
        self.docfile_with_current_anglesList=[None]
        for iterN in range(NumberOfIterations):
            fnBaseIter = "%s/Iter_%02d/%s.doc" % (self.WorkingDir, iterN + 1, self.docfile_with_current_angles)
            self.docfile_with_current_anglesList.append(fnBaseIter)
    
        #parameter for projection matching
        self.Align2DIterNr          = [-1]+getListFromVector(Align2DIterNr,NumberOfIterations)
        self.Align2dMaxChangeOffset = [-1]+getListFromVector(Align2dMaxChangeOffset,NumberOfIterations)
        self.Align2dMaxChangeRot    = [-1]+getListFromVector(Align2dMaxChangeRot,NumberOfIterations)
        self.AngSamplingRateDeg     = [-1]+getListFromVector(AngSamplingRateDeg,NumberOfIterations)
        self.DiscardPercentage      = [-1]+getListFromVector(DiscardPercentage,NumberOfIterations)
        self.DoAlign2D              = [False]+getBoolListFromVector(DoAlign2D,NumberOfIterations)
        self.DoComputeResolution    = [False]+getBoolListFromVector(DoComputeResolution,NumberOfIterations)
        self.DoSplitReferenceImages = [False]+getBoolListFromVector(DoSplitReferenceImages,NumberOfIterations)
        self.InnerRadius            = [False]+getListFromVector(InnerRadius,NumberOfIterations)
        self.MaxChangeInAngles      = [-1]+getListFromVector(MaxChangeInAngles,NumberOfIterations)
        self.MaxChangeOffset        = [-1]+getListFromVector(MaxChangeOffset,NumberOfIterations)
        self.MinimumCrossCorrelation= [-1]+getListFromVector(MinimumCrossCorrelation,NumberOfIterations)
        self.OnlyWinner             = [False]+getBoolListFromVector(OnlyWinner,NumberOfIterations)
        self.OuterRadius            = [False]+getListFromVector(OuterRadius,NumberOfIterations)
        self.PerturbProjectionDirections = [False]+getBoolListFromVector(PerturbProjectionDirections,NumberOfIterations)
        self.ReferenceIsCtfCorrected     = [-1]+getListFromVector(str(ReferenceIsCtfCorrected) + " True", NumberOfIterations)
        self.ScaleNumberOfSteps          = [-1]+getListFromVector(ScaleNumberOfSteps,NumberOfIterations)
        self.ScaleStep              = [-1]+getListFromVector(ScaleStep,NumberOfIterations)
        self.Search5DShift          = [-1]+getListFromVector(Search5DShift,NumberOfIterations)
        self.Search5DStep           = [-1]+getListFromVector(Search5DStep,NumberOfIterations)
        self.SymmetryGroup          = [-1]+getListFromVector(SymmetryGroup,NumberOfIterations)
         
        # Configure dabase
        self.Db.setPrintWrapperParameters(PrintWrapperParameters)
        self.Db.setPrintWrapperCommand(PrintWrapperCommand)
        self.Db.setVerify(Verify,ViewVerifyedFiles)
        self.Db.setParentDefault(XmippProtocolDbStruct.lastStep)

    def otherActionsToBePerformedBeforeLoop(self):
        print "in otherActionsToBePerformedBeforeLoop"

        _dataBase = self.Db
        if DoCtfCorrection:
            auxMD1 = MetaData(CTFDatName)
            auxMD2 = MetaData()
            auxMD2.aggregate(auxMD1, AGGR_COUNT,MDL_CTFMODEL,MDL_CTFMODEL,MDL_COUNT)
            self.NumberOfCtfGroups = auxMD2.size()
            
        #create dir for iteration 1 (This need to be 0 or 1? ROB FIXME
        _dataBase.insertAction('createDir', path = self.getIterDirName(0))
    
        #Check references and projections size match
        #Is already done in preconditions but I like to
        #run protocols from command line bypassing the gui
        _dataBase.insertAction('checkVolumeProjSize', 
                                                         ReferenceFileNames = self.ReferenceFileNames
                                                       , SelFileName = SelFileName)
    
        #Check Option compatibility
        _dataBase.insertAction('checkOptionsCompatibility', DoAlign2D       = self.DoAlign2D[1]
                                                          , DoCtfCorrection = DoCtfCorrection)
    
        #7 make CTF groups
        suffixes = ['_images.sel']
        if DoCtfCorrection:
            suffixes += ['Info.xmd', '_ctf.stk', '_wien.stk', '_split.doc']        
        fnBase = os.path.join(self.CtfGroupDirectory, self.CtfGroupRootName)
        _VerifyFiles = [fnBase + s for s in suffixes]
        _dataBase.insertAction('executeCtfGroups', _VerifyFiles, None, None,  CTFDatName = CTFDatName
                                                                            , CtfGroupDirectory = self.CtfGroupDirectory
                                                                            , CtfGroupMaxDiff = CtfGroupMaxDiff
                                                                            , CtfGroupMaxResol = CtfGroupMaxResol
                                                                            , CtfGroupRootName = self.CtfGroupRootName
                                                                            , DataArePhaseFlipped = DataArePhaseFlipped
                                                                            , DoAutoCtfGroup = DoAutoCtfGroup
                                                                            , DoCtfCorrection = DoCtfCorrection
                                                                            , PaddingFactor = PaddingFactor
                                                                            , SelFileName = SelFileName
                                                                            , SplitDefocusDocFile = SplitDefocusDocFile
                                                                            , WienerConstant = WienerConstant)
        #Create Initial angular file. Either fill it with zeros or copy input
        _VerifyFiles = [self.DocFileWithOriginalAngles]
        _dataBase.insertAction('initAngularReferenceFile', _VerifyFiles, None, None
                                                                , DocFileName =DocFileName
                                                                , DocFileWithOriginalAngles =self.DocFileWithOriginalAngles
                                                                , SelFileName =SelFileName)
    
        #Save all parameters in dict for future runs (this is done by database)
        #so far no parameter is being saved, but dummy=0
        self.Db.setIteration(XmippProtocolDbStruct.doAlways)
        _dataBase.insertAction('self.saveParameters', SystemFlavour =SystemFlavour)
        _dataBase.insertAction('self.loadParameters', None, None)
        self.Db.setIteration(1)
        #no entries will be save untill this commit
        print "commit databse"
        _dataBase.connection.commit()
    
    def getIterDirName(self, iterN):
        return os.path.join(self.projectDir, self.WorkingDir, 'Iter%02d' % iterN)
    
    def actionsToBePerformedInsideLoop(self):
        _log = self.Log
        _dataBase = self.Db
        
        for iterN in range(1, NumberOfIterations + 1):
            _dataBase.setIteration(iterN)
            # create IterationDir
            _dataBase.insertAction('createDir', path = self.getIterDirName(iterN))
    
            #Create directory with classes
            _dataBase.insertAction('createDir', path = self.ProjMatchDirs[iterN])
        
            #Create directory with image libraries
            id = _dataBase.insertAction('createDir', path = self.LibraryDirs[iterN])

            for refN in range(1, self.numberOfReferences + 1):
                # Mask reference volume
                _VerifyFiles = [self.maskedFileNamesIters[iterN][refN]]
                _dataBase.insertAction('executeMask', _VerifyFiles, id, None
                                    , DoMask =             DoMask
                                    , DoSphericalMask =    DoSphericalMask
                                    , maskedFileName =     self.maskedFileNamesIters[iterN][refN]
                                    , maskRadius =         MaskRadius
                                    , reconstructedFileName = self.reconstructedFileNamesIters[iterN - 1][refN]
                                    , userSuppliedMask =   MaskFileName)

                # angular_project_library
                #file with projections
                auxFn=self.ProjectLibraryRootNames[iterN][refN]
                _VerifyFiles=[auxFn]
                auxFn=auxFn[:-4]#remove extension
                #file with projection angles angles 
                _VerifyFiles.append(auxFn + ".doc")
                #file with sampling point neighbourhood 
                _VerifyFiles.append(auxFn + "_sampling.xmd")
                #file with sampling point neighbourhood for each ctf group, this is reduntant but useful
                for i in range (1,self.NumberOfCtfGroups+1):
                    _VerifyFiles.append(auxFn + "_group" + str(i).zfill(FILENAMENUMBERLENGTH) +"_sampling.xmd")
                            
                _dataBase.insertAction('angular_project_library', _VerifyFiles, None, None,
                                      AngSamplingRateDeg = self.AngSamplingRateDeg[iterN]
                                    , CtfGroupSubsetFileName = self.CtfGroupSubsetFileName
                                    , DoCtfCorrection = DoCtfCorrection
                                    , DocFileInputAngles = self.DocFileInputAngles[iterN-1]
                                    , DoParallel = self.DoParallel
                                    , DoRestricSearchbyTiltAngle = DoRestricSearchbyTiltAngle
                                    , MaxChangeInAngles = self.MaxChangeInAngles[iterN]
                                    , maskedFileNamesIter = self.maskedFileNamesIters[iterN][refN]
                                    , MpiJobSize = MpiJobSize
                                    , NumberOfMpiProcesses = NumberOfMpiProcesses
                                    , NumberOfThreads = NumberOfThreads
                                    , OnlyWinner = self.OnlyWinner[iterN]
                                    , PerturbProjectionDirections = self.PerturbProjectionDirections[iterN]
                                    , ProjectLibraryRootName = self.ProjectLibraryRootNames[iterN][refN]
                                    , SystemFlavour = SystemFlavour
                                    , SymmetryGroup = self.SymmetryGroup[iterN]
                                    , SymmetryGroupNeighbourhood = SymmetryGroupNeighbourhood
                                    , Tilt0  = Tilt0
                                    , TiltF = TiltF)
                # projectionMatching    
                #File with list of images and references
                _VerifyFiles=[self.ProjMatchRootNames[iterN][refN]]
                for i in range (1,self.NumberOfCtfGroups+1):
                    _VerifyFiles.append(auxFn + "_group" + str(i).zfill(6) +"_sampling.xmd")
                    
                _dataBase.insertAction('projection_matching', _VerifyFiles, None, None,
                                      AvailableMemory =AvailableMemory
                                    , CtfGroupRootName = self.CtfGroupRootName
                                    , CtfGroupDirectory = self.CtfGroupDirectory
                                    , DoComputeResolution =self.DoComputeResolution[iterN]
                                    , DoCtfCorrection = DoCtfCorrection
                                    , DoScale =DoScale
                                    , DoParallel = self.DoParallel
                                    , InnerRadius =self.InnerRadius[iterN]
                                    , MaxChangeOffset =self.MaxChangeOffset[iterN]
                                    , MpiJobSize =MpiJobSize
                                    , NumberOfCtfGroups =self.NumberOfCtfGroups
                                    , NumberOfMpiProcesses =NumberOfMpiProcesses
                                    , NumberOfThreads =NumberOfThreads
                                    , OuterRadius =self.OuterRadius[iterN]
                                    , PaddingFactor =PaddingFactor
                                    , ProjectLibraryRootName =self.ProjectLibraryRootNames[iterN][refN]
                                    , ProjMatchRootName =self.ProjMatchRootNames[iterN][refN]
                                    , ReferenceIsCtfCorrected =self.ReferenceIsCtfCorrected[iterN]
                                    , ScaleStep =self.ScaleStep[iterN]
                                    , ScaleNumberOfSteps =self.ScaleNumberOfSteps[iterN]
                                    , Search5DShift=self.Search5DShift[iterN]
                                    , Search5DStep=self.Search5DStep[iterN]
                                    , SystemFlavour=SystemFlavour)
                
    
            
            #assign the images to the different references based on the crosscorrelation coheficient
            #if only one reference it just copy the docfile generated in the previous step
            _VerifyFiles = [self.DocFileInputAngles[iterN]]
            _dataBase.insertAction('assign_images_to_references', [self.DocFileInputAngles[iterN]],
                                    None, None 
                                     , DocFileInputAngles = self.DocFileInputAngles[iterN]#Output file with angles
                                     , NumberOfCtfGroups  = self.NumberOfCtfGroups
                                     , ProjMatchRootName  = self.ProjMatchRootNames[iterN]#LIST
                                     , NumberOfReferences = self.numberOfReferences
                         )
    
            #align images, not possible for ctf groups
            for refN in range(1, self.numberOfReferences + 1):
                _VerifyFiles = []
                id = _dataBase.insertAction('angular_class_average', _VerifyFiles, None, None
                         , Align2DIterNr = self.Align2DIterNr[iterN]#
                         , Align2dMaxChangeRot = self.Align2dMaxChangeRot[iterN]#
                         , Align2dMaxChangeOffset = self.Align2dMaxChangeOffset[iterN]#
                         , CtfGroupDirectory = self.CtfGroupDirectory#
                         , CtfGroupRootName  = self.CtfGroupRootName#
                         , DiscardPercentage = self.DiscardPercentage[iterN]#
                         , DoAlign2D         = self.DoAlign2D[iterN]#
                         , DoCtfCorrection = DoCtfCorrection#
                         , DocFileInputAngles = self.DocFileInputAngles[iterN]#
                         , DoParallel = self.DoParallel#
                         , DoSplitReferenceImages =self.DoSplitReferenceImages[iterN]#
                         , InnerRadius =self.InnerRadius[iterN]#
                         , MaxChangeOffset =self.MaxChangeOffset[iterN]#
                         , MinimumCrossCorrelation =self.MinimumCrossCorrelation[iterN]#
                         , NumberOfReferences =self.numberOfReferences#
                         , NumberOfCtfGroups = self.NumberOfCtfGroups#
                         , NumberOfMpiProcesses =NumberOfMpiProcesses#
                         , NumberOfThreads =NumberOfThreads#
                         , PaddingFactor =PaddingFactor#
                         , ProjectLibraryRootName =self.ProjectLibraryRootNames[iterN][refN]#
                         , ProjMatchRootName =self.ProjMatchRootNames[iterN][refN]#
                         , refN =refN
                         , SystemFlavour =SystemFlavour#
                         )
                
                ##############REMOVE SHUTIL.COPY
                # Mask reference volume
                
                id = _dataBase.insertAction('executeMask', [self.maskedFileNamesIters[iterN][refN]],
                                            None,None
                                                , DoMask =             DoMask
                                                , DoSphericalMask =    DoSphericalMask
                                                , maskedFileName =     self.maskedFileNamesIters[iterN][refN]
                                                , maskRadius =         MaskRadius
                                                , reconstructedFileName = self.reconstructedFileNamesIters[iterN - 1][refN]
                                                , userSuppliedMask =   MaskFileName)
    
    #reconstruct
    #resolution
    #
                
    ######################
    ######################
    ########################            
                #REMOVE
                #Create directory
            _dataBase.insertAction('createDir', path =self.getIterDirName(iterN + 1))
            command = "shutil.copy('%s','%s');dummy" % (self.ReferenceFileNames[0], self.reconstructedFileNamesIters[iterN][refN])
            id = _dataBase.insertAction(command, self.reconstructedFileNamesIters[iterN][refN])
        command = "print 'ALL DONE';dummy"
        self.Db.setIteration(XmippProtocolDbStruct.doAlways)
        _dataBase.setIteration(XmippProtocolDbStruct.doAlways)
        self.Db.setIteration(iterN)

        id = _dataBase.insertAction(command)
    
        _dataBase.connection.commit()

    def defineActions(self):
        self.preRun()
        self.otherActionsToBePerformedBeforeLoop()
        self.actionsToBePerformedInsideLoop()
from protlib_utils import *

if __name__ == '__main__':
    protocolMain(ProtProjMatch)