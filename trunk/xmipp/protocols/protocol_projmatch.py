#!/usr/bin/env xmipp_python
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


import os
from os.path import join
from xmipp import MetaData, FILENAMENUMBERLENGTH, AGGR_COUNT, MDL_CTFMODEL, MDL_COUNT
from protlib_base import XmippProtocol, protocolMain
from protlib_utils import getListFromVector, getBoolListFromVector, getComponentFromVector
from protlib_sql import XmippProjectDb, SqliteDb
from config_protocols import protDict

class ProtProjMatch(XmippProtocol):
#    def __init__(self, scriptname, workingdir, projectdir=None, logdir='Logs', restartStep=1, isIter=True):
    def __init__(self, scriptname, project=None):
        super(ProtProjMatch, self).__init__(protDict.projmatch.name, scriptname, project)
        self.CtfGroupDirectory = self.workingDirPath(self.CtfGroupDirectory)
        self.CtfGroupSubsetFileName = join(self.CtfGroupDirectory, self.CtfGroupSubsetFileName)
        self.DocFileWithOriginalAngles = self.workingDirPath(self.DocFileWithOriginalAngles)
        
        #maskedFileNamesIter = []# names masked volumes used as reference
        self.numberOfReferences = 1#number of references
#        self.createAuxTable = False
        self.NumberOfCtfGroups = 1
        self.Import = 'from protocol_projmatch_before_loop import *;\
                       from protocol_projmatch_in_loop import *;' 
                #Convert directories/files  to absolute path from projdir
    def validate(self):
        errors = []
        #1 call base class, checks if project exists
        super(ProtProjMatch, self).validate()
        
        #2 Check reference and projection size match
        _ReferenceFileNames = getListFromVector(self.ReferenceFileNames)
        _Parameters = {
              'ReferenceFileNames':_ReferenceFileNames
            , 'SelFileName':self.SelFileName
            }
        from protocol_projmatch_before_loop import checkVolumeProjSize
        # a!!!
        #_retval, _error_message = checkVolumeProjSize(None,**_Parameters)
        #if(not _retval):
        #    errors.append(_error_message)
    
    
        # 3 Never allow DoAlign2D and DoCtfCorrection together
        if (int(getComponentFromVector(self.DoAlign2D, 1)) == 1 and self.DoCtfCorrection):
            errors.append("You cannot realign classes AND perform CTF-correction. Switch either of them off!")
    
        #4N outter radius is compulsory
        _OuterRadius = getComponentFromVector(self.OuterRadius, 1)
        _InnerRadius = getComponentFromVector(self.InnerRadius, 1)
        if _OuterRadius <= _InnerRadius:
            errors.append("OuterRadius must be larger than InnerRadius")

        return errors 
    
    def summary(self):
        from protlib_sql import XmippProtocolDb

        super(ProtProjMatch, self).summary()
        
        auxMD1 = MetaData(self.CTFDatName)
        auxMD2 = MetaData()
        auxMD2.aggregate(auxMD1, AGGR_COUNT, MDL_CTFMODEL, MDL_CTFMODEL, MDL_COUNT)
        self.NumberOfCtfGroups = auxMD2.size()

        self.ReferenceFileNames = getListFromVector(self.ReferenceFileNames)
        self.numberOfReferences = len(self.ReferenceFileNames)
        
        
        
        #runId = self.Db.getRunId(protDict.projmatch.name, protDict.projmatch.RunName)
        #print "runId", runId
        #THIS NEED TO BE FIXED. dO NOT OPEN THE DATA BASE JUST TO GET A NUMBERs
        _dataBase  = XmippProtocolDb(self, self.scriptName, False)
        
        iteration = _dataBase.getRunIter()
        summary = ['Performed <%d/%d> iterations with angular sampling rate <%s>' 
                   % (iteration, self.NumberOfIterations, self.AngSamplingRateDeg)]
        summary += ['Final Resolution is <%s>' % 'not yet implemented']
        summary += ['Number of CTFgroups and References is <%d> and <%d> respectively'
                        % (self.NumberOfCtfGroups, self.numberOfReferences)]

        return summary
    
    def createFilenameTemplates(self):  
        #Some class variables
        LibraryDir = "ReferenceLibrary"
        extraParams = {'ReferenceVolumeName': 'reference_volume.vol',
        'LibraryDir': LibraryDir,
        'ProjectLibraryRootName': join(LibraryDir, "gallery"),
        'ProjMatchDir': "ProjMatchClasses",
        'ProjMatchName': self.Name,
        'ClassAverageName': 'class_average',
        #ProjMatchRootName = ProjMatchDir + "/" + ProjMatchName
        'ForReconstructionSel': "reconstruction.sel",
        'ForReconstructionDoc': "reconstruction.doc",
        'MultiAlign2dSel': "multi_align2d.sel",
        'BlockWithAllExpImages' : 'all_exp_images',
        'DocFileWithOriginalAngles': 'original_angles.doc',
        'Docfile_with_current_angles': 'current_angles',
        'FilteredReconstruction': "filtered_reconstruction",
        'ReconstructedVolume': "reconstruction",
        'MaskReferenceVolume': "masked_reference",
        'OutputFsc': "resolution.fsc",
        'CtfGroupDirectory': "CtfGroups",
        'CtfGroupRootName': "ctf",
        'CtfGroupSubsetFileName': "ctf_images.sel"
        }
        self.ParamsDict.update(extraParams);
        # Set as protocol variables
        for k, v in extraParams.iteritems():
            setattr(self, k, v)
                              
        Iter = 'Iter_%(iter)03d'
        Ref3D = 'Ref3D_%(ref)03d'
        Ctf = 'CtfGroup_%(ctf)06d'
        IterDir = self.workingDirPath(Iter)
        
        #ProjMatchDirs = join(IterDir, '%(ProjMatchDir)s.doc')
        ProjMatchDirs = join(IterDir, '%(ProjMatchDir)s')
        _OutClassesXmd = join(ProjMatchDirs, '%(ProjMatchName)s_' + Ref3D + '.xmd')
        _OutClassesXmdS1 = join(ProjMatchDirs, '%(ProjMatchName)s_split_1_' + Ref3D + '.xmd')
        _OutClassesXmdS2 = join(ProjMatchDirs, '%(ProjMatchName)s_split_2_' + Ref3D + '.xmd')
        CtfGroupBase = join(self.workingDirPath(), self.CtfGroupDirectory, '%(CtfGroupRootName)s')
        ProjLibRootNames = join(IterDir, '%(ProjectLibraryRootName)s_' + Ref3D)
        return {
                # Global filenames templates
                'IterDir': IterDir,
                'ProjMatchDirs': ProjMatchDirs,
                'DocfileInputAnglesIters': join(IterDir, '%(Docfile_with_current_angles)s.doc'),
                #'LibraryDirs': join(IterDir, '%(LibraryDir)s.doc'),
                'LibraryDirs': join(IterDir, '%(LibraryDir)s'),
                'ProjectLibraryRootNames': ProjLibRootNames,
                'ProjMatchRootNames': join(ProjMatchDirs, '%(ProjMatchName)s_' + Ref3D + '.doc'),
                'ProjMatchRootNamesWithoutRef': join(ProjMatchDirs, '%(ProjMatchName)s.doc'),
                'OutClasses': join(ProjMatchDirs, '%(ProjMatchName)s'),
                'OutClassesXmd': _OutClassesXmd,
                'OutClassesStk': join(ProjMatchDirs, '%(ProjMatchName)s_' + Ref3D + '.stk'),
                'OutClassesDiscarded': join(ProjMatchDirs, '%(ProjMatchName)s_discarded.xmd'),
                'ReconstructionXmd': Ref3D + '@' +_OutClassesXmd,
                'ReconstructionXmdSplit1': Ref3D + '@' +_OutClassesXmdS1,
                'ReconstructionXmdSplit2': Ref3D + '@' +_OutClassesXmdS2,
                'MaskedFileNamesIters': join(IterDir, '%(MaskReferenceVolume)s_' + Ref3D + '.vol'),
                'ReconstructedFileNamesIters': join(IterDir, '%(ReconstructedVolume)s_' + Ref3D + '.vol'),
                'ReconstructedFileNamesItersSplit1': join(IterDir, '%(ReconstructedVolume)s_split_1_' + Ref3D + '.vol'),
                'ReconstructedFileNamesItersSplit2': join(IterDir, '%(ReconstructedVolume)s_split_2_' + Ref3D + '.vol'),
                'ReconstructedFilteredFileNamesIters': join(IterDir, '%(ReconstructedVolume)s_filtered_' + Ref3D + '.vol'),
                'ResolutionXmdFile': join(IterDir, '%(ReconstructedVolume)s_' + Ref3D + '.frc'),
                'ResolutionXmd': 'resolution@' + join(IterDir, '%(ReconstructedVolume)s_' + Ref3D + '.frc'),
                'ResolutionXmdMax': 'resolution_max@' + join(IterDir, '%(ReconstructedVolume)s_' + Ref3D + '.frc'),
                'MaskedFileNamesIters': join(IterDir, '%(MaskReferenceVolume)s_' + Ref3D + '.vol'),
                # Particular templates for executeCtfGroups  
                'ImageCTFpairs': CtfGroupBase + '_images.sel',
                'CTFGroupSummary': CtfGroupBase + 'Info.xmd',
                'StackCTFs': CtfGroupBase + '_ctf.stk',
                'StackWienerFilters': CtfGroupBase + '_wien.stk',
                'SplitAtDefocus': CtfGroupBase + '_split.doc',
                # Particular templates for angular_project_library 
                'ProjectLibraryStk': ProjLibRootNames + '.stk',
                'ProjectLibraryDoc': ProjLibRootNames + '.doc',
                'ProjectLibrarySampling': ProjLibRootNames + '_sampling.xmd',
                'ProjectLibraryGroupSampling': ProjLibRootNames + '_group%(group)06d_sampling.xmd',
                }
         
    
    def preRun(self):
        print "in PRERUN"
    #vector for iterations??????
    #    global ProjMatchDir
    #    ProjMatchDir = WorkingDir +'/' + ProjMatchDir
    
        
        # Convert vectors to list
        self.ReferenceFileNames = getListFromVector(self.ReferenceFileNames)
        self.numberOfReferences = len(self.ReferenceFileNames)
        #directory with ProjMatchClasses
#        self.ProjMatchDirs = [" "]
#        self.LibraryDirs = [" "]
#        self.DocFileInputAngles = [self.DocFileWithOriginalAngles]
#        #ProjMatchRootName=[" "]
#        
#        for iterN in range(self.NumberOfIterations):
#            fnBaseIter = "%s/Iter_%02d/" % (self.WorkingDir, iterN + 1)
#            self.ProjMatchDirs.append(fnBaseIter + self.ProjMatchDir)
#            self.LibraryDirs.append(fnBaseIter + self.LibraryDir)
#            self.DocFileInputAngles.append("%s%s.doc" % (fnBaseIter, self.docfile_with_current_angles))
#        
#        auxList = (self.numberOfReferences + 1) * [None]
#        self.ProjectLibraryRootNames = [[None]]
#        for iterN in range(self.NumberOfIterations):
#            fnBaseIter = "%s/Iter_%02d/" % (self.WorkingDir, iterN + 1)
#            for refN in range(self.numberOfReferences):                
#                auxList[refN + 1] = "%s%s_ref_%02d.stk" % (fnBaseIter, self.ProjectLibraryRootName, refN)
#            self.ProjectLibraryRootNames.append(list(auxList))
#                    
#        self.ProjMatchRootNames = [[None]]
#        for iterN in range(self.NumberOfIterations):
#            for refN in range(self.numberOfReferences):
#                auxList[refN + 1] = "%s/%s_ref_%02d.doc" % (self.ProjMatchDirs[iterN + 1], self.ProjMatchName, refN + 1)
#            self.ProjMatchRootNames.append(list(auxList))
#    
#        self.ProjMatchRootNamesWithoutRef = [[None]]
#        for iterN in range(self.NumberOfIterations):
#            self.ProjMatchRootNamesWithoutRef.append(list("%s/%s.doc" % (self.ProjMatchDirs[iterN + 1], self.ProjMatchName)))
#    
#        #name of masked volumes
#        #add dummy name so indexes start a 1
#        self.maskedFileNamesIters = [[None]]
#        for iterN in range(self.NumberOfIterations):
#            fnBaseIter = "%s/Iter_%02d/" % (self.WorkingDir, iterN + 1)
#            for refN in range(self.numberOfReferences):
#                auxList[refN + 1] = "%s%s_ref_%02d.vol" % (fnBaseIter, self.maskReferenceVolume, refN + 1)
#            self.maskedFileNamesIters.append(list(auxList))
#    
#        ####################################################################
#        #add initial reference, useful for many routines
#        #NOTE THAT INDEXES START AT 1

        # Construct special filename list with zero special case
        self.DocFileInputAngles = [self.DocFileWithOriginalAngles] + [self.getFilename('DocfileInputAnglesIters', iter=i) for i in range(1, self.NumberOfIterations + 1)]
        print 'self.DocFileInputAngles: ', self.DocFileInputAngles
        self.reconstructedFileNamesIters = [[None] + self.ReferenceFileNames]
        for iterN in range(1, self.NumberOfIterations + 1):
            self.reconstructedFileNamesIters.append([None] + [self.getFilename('ReconstructedFileNamesIters', iter=iterN, ref=r) for r in range(1, self.numberOfReferences + 1)])

        self.reconstructedFilteredFileNamesIters = [[None] + self.ReferenceFileNames]
        for iterN in range(1, self.NumberOfIterations + 1):
            self.reconstructedFilteredFileNamesIters.append([None] + [self.getFilename('ReconstructedFilteredFileNamesIters', iter=iterN, ref=r) for r in range(1, self.numberOfReferences + 1)])

#            fnBaseIter = "%s/Iter_%02d/" % (self.WorkingDir, iterN + 1)
#            for refN in range(self.numberOfReferences):
#                auxList[refN + 1] = "%s%s_ref_%02d.vol" % (fnBaseIter, self.ReconstructedVolume, refN + 1)
    
#        self.docfile_with_current_anglesList = [None]
#        for iterN in range(self.NumberOfIterations):
#            fnBaseIter = "%s/Iter_%02d/%s.doc" % (self.WorkingDir, iterN + 1, self.docfile_with_current_angles)
#            self.docfile_with_current_anglesList.append(fnBaseIter)
        _tmp = self.FourierMaxFrequencyOfInterest
        self.FourierMaxFrequencyOfInterest = list(-1 for  k in range(0, self.NumberOfIterations + 1))
        self.FourierMaxFrequencyOfInterest[1] = _tmp
        #parameter for projection matching
        self.Align2DIterNr = [-1] + getListFromVector(self.Align2DIterNr, self.NumberOfIterations)
        self.Align2dMaxChangeOffset = [-1] + getListFromVector(self.Align2dMaxChangeOffset, self.NumberOfIterations)
        self.Align2dMaxChangeRot = [-1] + getListFromVector(self.Align2dMaxChangeRot, self.NumberOfIterations)
        self.AngSamplingRateDeg = [-1] + getListFromVector(self.AngSamplingRateDeg, self.NumberOfIterations)
        self.ConstantToAddToFiltration = [-1] + getListFromVector(self.ConstantToAddToFiltration, self.NumberOfIterations)
        self.ConstantToAddToMaxReconstructionFrequency = [-1] + getListFromVector(self.ConstantToAddToMaxReconstructionFrequency, self.NumberOfIterations)
        self.DiscardPercentage = [-1] + getListFromVector(self.DiscardPercentage, self.NumberOfIterations)
        self.DiscardPercentagePerClass = [-1] + getListFromVector(self.DiscardPercentagePerClass, self.NumberOfIterations)
        self.DoAlign2D = [False] + getBoolListFromVector(self.DoAlign2D, self.NumberOfIterations)
        self.DoComputeResolution = [False] + getBoolListFromVector(self.DoComputeResolution, self.NumberOfIterations)
        self.DoSplitReferenceImages = [False] + getBoolListFromVector(self.DoSplitReferenceImages, self.NumberOfIterations)
        self.InnerRadius = [False] + getListFromVector(self.InnerRadius, self.NumberOfIterations)
        self.MaxChangeInAngles = [-1] + getListFromVector(self.MaxChangeInAngles, self.NumberOfIterations)
        self.MaxChangeOffset = [-1] + getListFromVector(self.MaxChangeOffset, self.NumberOfIterations)
        self.MinimumCrossCorrelation = [-1] + getListFromVector(self.MinimumCrossCorrelation, self.NumberOfIterations)
        self.OnlyWinner = [False] + getBoolListFromVector(self.OnlyWinner, self.NumberOfIterations)
        self.OuterRadius = [False] + getListFromVector(self.OuterRadius, self.NumberOfIterations)
        self.PerturbProjectionDirections = [False] + getBoolListFromVector(self.PerturbProjectionDirections, self.NumberOfIterations)
        self.ReferenceIsCtfCorrected = [-1] + getListFromVector(str(self.ReferenceIsCtfCorrected) + " True", self.NumberOfIterations)
        self.ScaleNumberOfSteps = [-1] + getListFromVector(self.ScaleNumberOfSteps, self.NumberOfIterations)
        self.ScaleStep = [-1] + getListFromVector(self.ScaleStep, self.NumberOfIterations)
        self.Search5DShift = [-1] + getListFromVector(self.Search5DShift, self.NumberOfIterations)
        self.Search5DStep = [-1] + getListFromVector(self.Search5DStep, self.NumberOfIterations)
        self.SymmetryGroup = [-1] + getListFromVector(self.SymmetryGroup, self.NumberOfIterations)
         
        # Configure dabase
        ###############self.Db.setVerify(self.Verify,self.ViewVerifyedFiles)
        ###############self.Db.setParentDefault(XmippProjectDb.lastStep)
        

    def otherActionsToBePerformedBeforeLoop(self):
        print "in otherActionsToBePerformedBeforeLoop"
        _VerifyFiles = []

        _dataBase = self.Db
        if self.DoCtfCorrection:
            auxMD1 = MetaData(self.CTFDatName)
            auxMD2 = MetaData()
            auxMD2.aggregate(auxMD1, AGGR_COUNT, MDL_CTFMODEL, MDL_CTFMODEL, MDL_COUNT)
            self.NumberOfCtfGroups = auxMD2.size()
            print "self.NumberOfCtfGroups: ", self.NumberOfCtfGroups
            
        #create dir for iteration 1 (This need to be 0 or 1? ROB FIXME
        #!a _dataBase.insertStep('createDir', path = self.getIterDirName(0))
    
        #Check references and projections size match
        #Is already done in preconditions but I like to
        #run protocols from command line bypassing the gui
        _dataBase.insertStep('checkVolumeProjSize',
                                                         ReferenceFileNames=self.ReferenceFileNames
                                                       , SelFileName=self.SelFileName)
    
        #Check Option compatibility
        _dataBase.insertStep('checkOptionsCompatibility', DoAlign2D=self.DoAlign2D[1]
                                                          , DoCtfCorrection=self.DoCtfCorrection)
    
        #7 make CTF groups
        verifyFiles = [self.getFilename('ImageCTFpairs')]
        if self.DoCtfCorrection:
            verifyFiles += [self.getFilename(k) \
                            for k in ['CTFGroupSummary','StackCTFs','StackWienerFilters','SplitAtDefocus']]

        _dataBase.insertStep('executeCtfGroups', verifyfiles=verifyFiles
                                               , CTFDatName=self.CTFDatName
                                               , CtfGroupDirectory=self.CtfGroupDirectory
                                               , CtfGroupMaxDiff=self.CtfGroupMaxDiff
                                               , CtfGroupMaxResol=self.CtfGroupMaxResol
                                               , CtfGroupRootName=self.CtfGroupRootName
                                               , DataArePhaseFlipped=self.DataArePhaseFlipped
                                               , DoAutoCtfGroup=self.DoAutoCtfGroup
                                               , DoCtfCorrection=self.DoCtfCorrection
                                               , PaddingFactor=self.PaddingFactor
                                               , SelFileName=self.SelFileName
                                               , SplitDefocusDocFile=self.SplitDefocusDocFile
                                               , WienerConstant=self.WienerConstant)
        #Create Initial angular file. Either fill it with zeros or copy input
        _dataBase.insertStep('initAngularReferenceFile', verifyfiles=[self.DocFileWithOriginalAngles]
                                                                , BlockWithAllExpImages = self.BlockWithAllExpImages
                                                                , CtfGroupDirectory=self.CtfGroupDirectory
                                                                , CtfGroupRootName=self.CtfGroupRootName
                                                                , DocFileName=self.DocFileName
                                                                , DocFileWithOriginalAngles=self.DocFileWithOriginalAngles
                                                                , SelFileName=self.SelFileName)
    
        #Save all parameters in dict for future runs (this is done by database)
        #so far no parameter is being saved, but dummy=0
        #self.Db.setIteration(XmippProjectDb.doAlways)
        #_dataBase.insertStep('self.saveParameters', SystemFlavour = self.SystemFlavour)
        #_dataBase.insertStep('self.loadParameters', None, None)
        #self.Db.setIteration(1)
        #no entries will be save untill this commit
        #print "commit databse"
        _dataBase.connection.commit()
    
    def getIterDirName(self, iterN):
        return join(self.projectDir, self.WorkingDir, 'Iter_%02d' % iterN)
    
    def actionsToBePerformedInsideLoop(self):
        _log = self.Log
        _dataBase = self.Db

        for iterN in range(1, self.NumberOfIterations + 1):
            _dataBase.setIteration(iterN)
            # create IterationDir
            _dataBase.insertStep('createDir', path=self.getFilename('IterDir', iter=iterN))
    
            #Create directory with classes
            _dataBase.insertStep('createDir', path=self.getFilename('ProjMatchDirs', iter=iterN))
        
            #Create directory with image libraries
            id = _dataBase.insertStep('createDir', path=self.getFilename('LibraryDirs', iter=iterN))

            ProjMatchRootNameList = ['']
            for refN in range(1, self.numberOfReferences + 1):
                # Mask reference volume
                maskedFileName = self.getFilename('MaskedFileNamesIters', iter=iterN, ref=refN)
                _dataBase.insertStep('executeMask'
                                    , verifyfiles=[maskedFileName]
                                    , parent_step_id=id
                                    , DoMask=self.DoMask
                                    , DoSphericalMask=self.DoSphericalMask
                                    , maskedFileName=maskedFileName
                                    , maskRadius=self.MaskRadius
                                    , ReconstructedFilteredVolume = self.reconstructedFilteredFileNamesIters[iterN-1][refN]
                                    , userSuppliedMask=self.MaskFileName)

                # angular_project_library
                #file with projections
                auxFn = self.getFilename('ProjectLibraryRootNames', iter=iterN, ref=refN)
                #file with sampling point neighbourhood for each ctf group, this is reduntant but useful
                #Files: projections, projection_angles, sampling_points and neighbourhood
                _VerifyFiles = [self.getFilename('ProjectLibrary' + e, iter=iterN, ref=refN)
                                     for e in ['Stk', 'Doc', 'Sampling']]
                _VerifyFiles = _VerifyFiles + [self.getFilename('ProjectLibraryGroupSampling', iter=iterN, ref=refN, group=g) \
                                     for g in range (1, self.NumberOfCtfGroups + 1)]
                projLibFn =  self.getFilename('ProjectLibraryStk', iter=iterN, ref=refN)  
                         
                _dataBase.insertStep('angular_project_library', verifyfiles=_VerifyFiles
                                    , AngSamplingRateDeg=self.AngSamplingRateDeg[iterN]
                                    , BlockWithAllExpImages = self.BlockWithAllExpImages
                                    , CtfGroupSubsetFileName=self.CtfGroupSubsetFileName
                                    , DoCtfCorrection=self.DoCtfCorrection
                                    , DocFileInputAngles=self.DocFileInputAngles[iterN - 1]
                                    , DoParallel=self.DoParallel
                                    , DoRestricSearchbyTiltAngle=self.DoRestricSearchbyTiltAngle
                                    , MaxChangeInAngles=self.MaxChangeInAngles[iterN]
                                    , maskedFileNamesIter=maskedFileName
                                    , NumberOfMpi=self.NumberOfMpi
                                    , NumberOfThreads=self.NumberOfThreads
                                    , MpiJobSize=self.MpiJobSize
                                    , OnlyWinner=self.OnlyWinner[iterN]
                                    , PerturbProjectionDirections=self.PerturbProjectionDirections[iterN]
                                    , ProjectLibraryRootName=projLibFn
                                    , SymmetryGroup=self.SymmetryGroup[iterN]
                                    , SymmetryGroupNeighbourhood=self.SymmetryGroupNeighbourhood
                                    , Tilt0=self.Tilt0
                                    , TiltF=self.TiltF)
                # projectionMatching    
                #File with list of images and references
                ProjMatchRootName = self.getFilename('ProjMatchRootNames', iter=iterN, ref=refN)
                ProjMatchRootNameList.append(ProjMatchRootName)
#                for i in range (1, self.NumberOfCtfGroups + 1):
#                    _VerifyFiles.append(auxFn + "_group" + str(i).zfill(6) + "_sampling.xmd")
                    
                _dataBase.insertStep('projection_matching', verifyfiles=[ProjMatchRootName],
                                      AvailableMemory=self.AvailableMemory
                                    , CtfGroupRootName=self.CtfGroupRootName
                                    , CtfGroupDirectory=self.CtfGroupDirectory
                                    , DocFileInputAngles=self.DocFileInputAngles[iterN - 1]
                                    , DoComputeResolution=self.DoComputeResolution[iterN]
                                    , DoCtfCorrection=self.DoCtfCorrection
                                    , DoScale=self.DoScale
                                    , DoParallel=self.DoParallel
                                    , InnerRadius=self.InnerRadius[iterN]
                                    , MaxChangeOffset=self.MaxChangeOffset[iterN]
                                    , MpiJobSize=self.MpiJobSize
                                    , NumberOfCtfGroups=self.NumberOfCtfGroups
                                    , NumberOfMpi=self.NumberOfMpi
                                    , NumberOfThreads=self.NumberOfThreads
                                    , OuterRadius=self.OuterRadius[iterN]
                                    , PaddingFactor=self.PaddingFactor
                                    , ProjectLibraryRootName=projLibFn
                                    , ProjMatchRootName=ProjMatchRootName
                                    , ReferenceIsCtfCorrected=self.ReferenceIsCtfCorrected[iterN]
                                    , ScaleStep=self.ScaleStep[iterN]
                                    , ScaleNumberOfSteps=self.ScaleNumberOfSteps[iterN]
                                    , Search5DShift=self.Search5DShift[iterN]
                                    , Search5DStep=self.Search5DStep[iterN]
                                    )

            
            #assign the images to the different references based on the crosscorrelation coheficient
            #if only one reference it just copy the docfile generated in the previous step
            _dataBase.insertStep('assign_images_to_references', verifyfiles=[self.DocFileInputAngles[iterN]]
                                     , BlockWithAllExpImages = self.BlockWithAllExpImages
                                     , DocFileInputAngles=self.DocFileInputAngles[iterN]#Output file with angles
                                     , NumberOfCtfGroups=self.NumberOfCtfGroups
                                     , ProjMatchRootName=ProjMatchRootNameList#LIST
                                     , NumberOfReferences=self.numberOfReferences
                         )
    
            #align images, not possible for ctf groups
            _VerifyFiles = [self.getFilename('OutClassesDiscarded', iter=iterN)]
            _VerifyFiles = _VerifyFiles + [self.getFilename('OutClassesXmd', iter=iterN, ref=g) \
                                     for g in range (1, self.numberOfReferences + 1)]
            _VerifyFiles = _VerifyFiles + [self.getFilename('OutClassesStk', iter=iterN, ref=g) \
                                     for g in range (1, self.numberOfReferences + 1)]

            _dataBase.insertStep('angular_class_average', verifyfiles=_VerifyFiles
                             , Align2DIterNr=self.Align2DIterNr[iterN]#
                             , Align2dMaxChangeOffset=self.Align2dMaxChangeOffset[iterN]#
                             , Align2dMaxChangeRot=self.Align2dMaxChangeRot[iterN]#
                             , CtfGroupDirectory=self.CtfGroupDirectory#
                             , CtfGroupRootName=self.CtfGroupRootName#
                             , DiscardImages=self.DiscardImages#
                             , DiscardPercentage=self.DiscardPercentage[iterN]#
                             , DiscardPercentagePerClass=self.DiscardPercentagePerClass[iterN]#
                             , DoAlign2D=self.DoAlign2D[iterN]#
                             , DoComputeResolution=self.DoComputeResolution[iterN]
                             , DoCtfCorrection=self.DoCtfCorrection#
                             , DocFileInputAngles=self.DocFileInputAngles[iterN]#
                             , DoParallel=self.DoParallel
                             , DoSaveImagesAssignedToClasses=self.DoSaveImagesAssignedToClasses#
                             , DoSplitReferenceImages=self.DoSplitReferenceImages[iterN]#
                             , InnerRadius=self.InnerRadius[iterN]#
                             , MaxChangeOffset=self.MaxChangeOffset[iterN]#
                             , MinimumCrossCorrelation=self.MinimumCrossCorrelation[iterN]#
                             , MpiJobSize=self.MpiJobSize
                             , NumberOfMpi=self.NumberOfMpi
                             , NumberOfThreads=self.NumberOfThreads
                             , OutClasses=self.getFilename('OutClasses', iter=iterN)#
                             , PaddingFactor=self.PaddingFactor#
                             , ProjectLibraryRootName=self.getFilename('ProjectLibraryStk', iter=iterN, ref=refN)#
                             )
            

            for refN in range(1, self.numberOfReferences + 1):
    
                #self._ReconstructionMethod=arg.getComponentFromVector(_ReconstructionMethod, iterN-1)
                #self._ARTLambda=arg.getComponentFromVector(_ARTLambda, iterN-1)

                #if (DoReconstruction):
                _VerifyFiles = [self.getFilename('ReconstructedFileNamesIters', iter=iterN, ref=refN)]
                id = _dataBase.insertStep('reconstruction', verifyfiles=_VerifyFiles
                                              , ARTReconstructionExtraCommand = self.ARTReconstructionExtraCommand
                                              , WBPReconstructionExtraCommand = self.WBPReconstructionExtraCommand
                                              , FourierReconstructionExtraCommand = self.FourierReconstructionExtraCommand
                                              , Iteration_number  = iterN
                                              , DoParallel=self.DoParallel#
                                              , maskedFileNamesIter=maskedFileName
                                              , MpiJobSize=self.MpiJobSize
                                              , NumberOfMpi=self.NumberOfMpi#
                                              , NumberOfThreads=self.NumberOfThreads#
                                              , ReconstructionMethod = self.ReconstructionMethod
                                              , FourierMaxFrequencyOfInterest = self.FourierMaxFrequencyOfInterest[iterN]
                                              , ARTLambda = self.ARTLambda
                                              , SymmetryGroup = self.SymmetryGroup[iterN]
                                              , ReconstructionXmd = self.getFilename('ReconstructionXmd', iter=iterN, ref=refN)
                                              , ReconstructedVolume = self.getFilename('ReconstructedFileNamesIters', iter=iterN, ref=refN)
                                              , PaddingFactor = self.PaddingFactor
                                              , ResolSam = self.ResolSam
                                              , ResolutionXmdPrevIterMax = self.getFilename('ResolutionXmdMax', iter=iterN-1, ref=refN)
                                              , ConstantToAddToFiltration = self.ConstantToAddToMaxReconstructionFrequency[iterN]
                                              )
                    
                if(self.DoSplitReferenceImages[iterN]):
                    
                    _VerifyFiles = [self.getFilename('ReconstructedFileNamesItersSplit1', iter=iterN, ref=refN)]
                    id = _dataBase.insertStep('reconstruction', verifyfiles=_VerifyFiles
                                              , ARTReconstructionExtraCommand = self.ARTReconstructionExtraCommand
                                              , WBPReconstructionExtraCommand = self.WBPReconstructionExtraCommand
                                              , FourierReconstructionExtraCommand = self.FourierReconstructionExtraCommand
                                              , Iteration_number  = iterN
                                              , DoParallel=self.DoParallel#
                                              , maskedFileNamesIter=maskedFileName
                                              , MpiJobSize=self.MpiJobSize
                                              , NumberOfMpi=self.NumberOfMpi#
                                              , NumberOfThreads=self.NumberOfThreads#
                                              , ReconstructionMethod = self.ReconstructionMethod
                                              , FourierMaxFrequencyOfInterest = self.FourierMaxFrequencyOfInterest[iterN]
                                              , ARTLambda = self.ARTLambda
                                              , SymmetryGroup = self.SymmetryGroup[iterN]
                                              , ReconstructionXmd = self.getFilename('ReconstructionXmdSplit1', iter=iterN, ref=refN)
                                              , ReconstructedVolume = self.getFilename('ReconstructedFileNamesItersSplit1', iter=iterN, ref=refN)
                                              , PaddingFactor = self.PaddingFactor
                                              , ResolSam = self.ResolSam
                                              , ResolutionXmdPrevIterMax = self.getFilename('ResolutionXmdMax', iter=iterN-1, ref=refN)
                                              , ConstantToAddToFiltration = self.ConstantToAddToMaxReconstructionFrequency[iterN]
                                              )

                    _VerifyFiles = [self.getFilename('ReconstructedFileNamesItersSplit2', iter=iterN, ref=refN)]
                    id = _dataBase.insertStep('reconstruction', verifyfiles=_VerifyFiles
                                              , ARTReconstructionExtraCommand = self.ARTReconstructionExtraCommand
                                              , WBPReconstructionExtraCommand = self.WBPReconstructionExtraCommand
                                              , FourierReconstructionExtraCommand = self.FourierReconstructionExtraCommand
                                              , Iteration_number  = iterN
                                              , DoParallel=self.DoParallel#
                                              , MpiJobSize=self.MpiJobSize
                                              , maskedFileNamesIter=maskedFileName
                                              , NumberOfMpi=self.NumberOfMpi#
                                              , NumberOfThreads=self.NumberOfThreads#
                                              , ReconstructionMethod = self.ReconstructionMethod
                                              , FourierMaxFrequencyOfInterest = self.FourierMaxFrequencyOfInterest[iterN]
                                              , ARTLambda = self.ARTLambda
                                              , SymmetryGroup = self.SymmetryGroup[iterN]
                                              , ReconstructionXmd = self.getFilename('ReconstructionXmdSplit2', iter=iterN, ref=refN)
                                              , ReconstructedVolume = self.getFilename('ReconstructedFileNamesItersSplit2', iter=iterN, ref=refN)
                                              , PaddingFactor = self.PaddingFactor
                                              , ResolSam = self.ResolSam
                                              , ResolutionXmdPrevIterMax = self.getFilename('ResolutionXmdMax', iter=iterN-1, ref=refN)
                                              , ConstantToAddToFiltration = self.ConstantToAddToMaxReconstructionFrequency[iterN]
                                              )
                    
                    _VerifyFiles = [self.getFilename('ResolutionXmdFile', iter=iterN, ref=refN)]
                    id = _dataBase.insertStep('compute_resolution', verifyfiles=_VerifyFiles
                                                 , FourierMaxFrequencyOfInterest = self.FourierMaxFrequencyOfInterest[iterN]
                                                 , ReconstructionMethod = self.ReconstructionMethod
                                                 , ResolutionXmdCurrIter = self.getFilename('ResolutionXmd', iter=iterN, ref=refN)
                                                 , ResolutionXmdCurrIterMax = self.getFilename('ResolutionXmdMax', iter=iterN, ref=refN)
                                                 , ResolutionXmdPrevIterMax = self.getFilename('ResolutionXmdMax', iter=iterN-1, ref=refN)
                                                 , OuterRadius = self.OuterRadius[iterN]
                                                 , ReconstructedVolumeSplit1 = self.getFilename('ReconstructedFileNamesItersSplit1', iter=iterN, ref=refN)
                                                 , ReconstructedVolumeSplit2 = self.getFilename('ReconstructedFileNamesItersSplit2', iter=iterN, ref=refN)
                                                 , ResolSam = self.ResolSam
                                                  )

                    id = _dataBase.insertStep('filter_volume', verifyfiles=_VerifyFiles
                                              , FourierMaxFrequencyOfInterest = self.FourierMaxFrequencyOfInterest[iterN]
                                              , ReconstructedVolume = self.getFilename('ReconstructedFileNamesIters', iter=iterN, ref=refN)
                                              , ReconstructedFilteredVolume = self.reconstructedFilteredFileNamesIters[iterN][refN]
                                              , DoComputeResolution = self.DoComputeResolution
                                              , OuterRadius = self.OuterRadius
                                              , DoLowPassFilter = self.DoLowPassFilter
                                              , UseFscForFilter = self.UseFscForFilter
                                              , ConstantToAddToFiltration = self.ConstantToAddToFiltration[iterN]
                                              , ResolutionXmdPrevIterMax = self.getFilename('ResolutionXmdMax', iter=iterN-1, ref=refN)
                                              , ResolSam = self.ResolSam
                                              )
        _dataBase.connection.commit()

    def defineSteps(self):
        self.preRun()
        self.otherActionsToBePerformedBeforeLoop()
        self.actionsToBePerformedInsideLoop()

