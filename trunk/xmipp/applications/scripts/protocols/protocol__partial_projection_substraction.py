#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Xmipp protocol for subtraction
#
# Example use:
# ./xmipp_partial_projection_subtraction.py
#
# Authors: 
#
# {begin_of_header}
#-----------------------------------------------------------------------------
# {section} Global parameters
#-----------------------------------------------------------------------------
#Comment
Comment='Describe your project here...'

#Comment
""" Subtraction subdirectory. 'run_001' will create 'Subtraction/run_001' directory
"""
RunName='run_001'

# {file} Protocol Name
ProtocolName='ProjMatch/xmipp_2.4_subtraction_crunchy_shifts/xmipp_protocol_projmatch_shifts_backup.py'

# {expert}{file} CTFDat file with CTF data:
""" The input selfile may be a subset of the images in the CTFDat file, but all 
    images in the input selfile must be present in the CTFDat file. This field is 
    obligatory if CTF correction is to be performed. 
    Note that this file should be positioned in the project directory, and that the
    image names and ctf parameter filenames should be in absolute paths.
    Usually optained from the protocol file
"""

#Show results for iteration
""" Use data coming from iteration
"""
iterationNo=5

# Resume at iteration
""" Set to 1 to start a new run, set to -1 to continue the process (where you left it),
    set to a positive number N to restart at the begining of iteration N
    Note1: Do NOT delete working directory if this option is not set to 1
    Note2: Set this option to -1 if you want to perform extra iterations after
           successfully finish an execution
"""
ContinueAtIteration =1

# {expert} Resume at Iter (vs Step)
"""This option control how to resume a previously performed run.
    Set to TRUE to restart at the beginning of iteration N
    or FALSE to continue at step N. (N is set in the next parameter).
    NOTE:If you do not know what are you doing make it equal to False
"""
IsIter =False

#{expert}{file} Select the input volume
"""{expert DELETE THIS ENTRY}Name of the reference volume by default Iter_X_reconstruction.vol

"""
szInputVolumeName=''
#{expert}{file} Select docfile used to compute references
"""Name of the doc file used to compute reference library, by dfault
   ../Iter_(X-1)/Iter_(X-1)_current_angles.doc
"""
DocFileExp=''

# {expert} Mask reference volume?
doMask =True

#Crown Mask radius (inner)
dRradiusMin=39

#{expert}  iCrown mask radius center (outter)
dRradiusMax=64

# {expert} Backup Experimental Proj. angles
doBackupProjectionAngles =True

# {expert}Angular sampling rate
"""Angular distance (in degrees) between neighboring projection  points
   usually obtained from protocol file
"""   
AngSamplingRateDeg=''

# {expert}  Angular search range 
"""Maximum change in rot & tilt  (in +/- degrees)
   usually obtained from protocol file
"""   
MaxChangeInAngles =''

# {expert} Symmetry group
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry
    for a description of the symmetry groups format
    If no symmetry is present, give c1
    usually obtained from protocol file
"""
SymmetryGroup=''

CTFDatName=''
# {expert} Correct by CTF:
""" Set to True if you want to correct by CTF
"""
doCTFCorrection=False

# Scale images?
""" Set to True if you want to scale images (Using padding/windowing in Fourier space)
"""
doScaleImages=True

# New X dimension
""" New X dimension.
"""
dimX = 32

# {expert} New Y dimensions
""" New Y dimension. -1 means Y = X
"""
dimY = -1


#------------------------------------------------------------------------------------------------
# {section} Parallelization issues
#------------------------------------------------------------------------------------------------
# Number of (shared-memory) threads?
""" This option provides shared-memory parallelization on multi-core machines. 
    It does not require any additional software, other than xmipp
"""
NumberOfThreads = 2

# Number of MPI processes to use
NumberOfMpiProcesses = 2

#MPI job size 
"""Minimum size of jobs in mpi processes. Set to 1 for large images (e.g. 500x500) and to 10 for small images (e.g. 100x100)
"""
MpiJobSize ='1'

#Submmit to queue
"""Submmit to queue"""
SubmmitToQueue=False
#------------------------------------------------------------------------------------------------
# {section}{expert}{condition}(SubmmitToQueue=True) Queue 
#------------------------------------------------------------------------------------------------

# Queue name
"""Name of the queue to submit the job"""
QueueName="default"
# Queue hours
"""This establish a maximum number of hours the job will
be running, after that time it will be killed by the 
queue system"""
QueueHours=72

# minumum size of jobs in mpi processe. Set to 1 for large images (e.g. 500x500) and to 10 for small images (e.g. 100x100)
#MpiJobSize ='3'

#------------------------------------------------------------------------------------------------
# {hidden} Analysis of results
""" This script serves only for GUI-assisted visualization of the results
"""
AnalysisScript='visualize_partial_projection_subtraction.py'
#-----------------------------------------------------------------------------
# {section} Debug
#-----------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------

#Verify
"""Check that some output files are created. 
"""
Verify=True

# {expert} print wrapper name
PrintWrapperCommand=True

# {expert} print wrapper parameters
PrintWrapperParameters=True

# {expert} show file verification
ViewVerifyedFiles=True 

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
SystemFlavour = "TORQUE-OPENMPI"


import os,sys
def checkErrors():    
    '''This function will be used to validate the protocols
    should be implemented in all derived protocols'''
    errors = []        
    # Check if there is workingdir 
    # move this to gui
    if not os.path.exists(ProtocolName):
        errors.append("Refered protocol named %s does not exist"%ProtocolName)
        
    return errors 

from protlib_base import *
from xmipp import *

class ProtPartialProjectionSubtraction(XmippProtocol):

    def __init__(self, scriptname,project=None):
        from config import *
        super(ProtPartialProjectionSubtraction,self).__init__(ProtocolNames.subtraction, scriptname, RunName, project,NumberOfMpiProcesses)
        #super(ProtPartialProjectionSubtraction,self).__init__(scriptname, workingdir, ProjectDir, logdir, restartStep, isIter)
    #def __init__(self, scriptname, workingdir, projectdir=None, logdir='Logs', restartStep=1, isIter=True):
    #    super(ProtProjMatch,self).__init__(ProtocolNames.projmatch,scriptname, RunName, project)
        self.myName='partial_projection_subtraction'
        #self.run_file1='./readDocfileAndPairExperimentalAndReferenceImages_v3.sh'
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
        
    def preRun(self):

        print "in PRERUN"
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
        
        print "dRradiusMax: ", self.dRradiusMax
        print "dRradiusMin: ", self.dRradiusMin
            
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
        
        tmpFileName = os.path.join(self.pmprotWorkingDir,'CtfGroups/ctf_group??????.sel')
        self.defGroups=['']
        import glob
        self.defGroups = glob.glob(tmpFileName)
        
        print "self.defGroups: ", self.defGroups 
        
        self.defocusGroupNo = len(self.defGroups)
        
        self.DocFileExp=['']
        if DocFileExp!='':
            tmpFileName = DocFileExp
        else:
            tmpFileName = FileName(self.filename_currentAngles)
        tmpFileName = tmpFileName.withoutExtension()
        for iterN in range(1, self.defocusGroupNo ):
            tmpDocFileExp = 'ctfgroup_' + str(iterN).zfill(6) + '@' + tmpFileName + '_ctfgroups.doc'
            self.DocFileExp.append(tmpDocFileExp
                              )
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
        
    def actionsToBePerformedInsideLoop(self):
        
        print "actionsToBePerformedInsideLoop"
        print "self.defocusGroupNo: ", self.defocusGroupNo
        _dataBase = self.Db
        for iterN in range(1, self.defocusGroupNo):
            print "inside for"
            #Create auxiliary metadata with image names , angles and CTF
            if(self.doScaleImages):
                inputSelfile = self.scaledImages
            else:
                inputSelfile = self.filename_currentAngles

            _VerifyFiles = []
            auxFilename = FileName(self.DocFileExp[iterN])
            _VerifyFiles.append(auxFilename.removeBlockName())
            id = _dataBase.insertAction('joinImageCTF', _VerifyFiles, None, None
                                        , CTFgroupName = self.defGroups[iterN]
                                        , DocFileExp = self.DocFileExp[iterN]
                                        , inputSelfile = inputSelfile)
            
            if(self.doScaleImages):
                _VerifyFiles = []
                id = _dataBase.insertAction('updateImageLabel', _VerifyFiles, None, None
                                            , inputSelfile = self.DocFileExp[iterN])
                
            #reconstruct each CTF group
            _VerifyFiles = []
            _VerifyFiles.append(self.reconstructedVolume[iterN])
            id = _dataBase.insertAction('reconstructVolume', _VerifyFiles, None, None
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
            id = _dataBase.insertAction('maskVolume', _VerifyFiles, None, None
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
            
            id = _dataBase.insertAction('createProjections', _VerifyFiles, None, None
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
            
            id = _dataBase.insertAction('subtractionScript', _VerifyFiles, None, None
                                        , DocFileExp = self.DocFileExp[iterN] 
                                        , referenceStackDoc = self.referenceStackDoc[iterN]
                                        , subtractedStack = self.subtractedStack[iterN])
                          
            
        
        
    def otherActionsToBePerformedBeforeLoop(self):
        _dataBase = self.Db
        #Create directories
        _dataBase.insertAction('createDir', path = self.volsDir)
        _dataBase.insertAction('createDir', path = self.referenceDir)
        _dataBase.insertAction('createDir', path = self.subImgsDir)
        
        #Create auxiliary metadata with image names , angles and CTF
        
        if(doScaleImages):
            _VerifyFiles = [self.scaledImages+".stk"]
            _VerifyFiles.append(self.scaledImages+".xmd")
            id = _dataBase.insertAction('scaleImages', _VerifyFiles, None, None
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
            self.scaledImages = self.scaledImages + ".xmd"

    def defineActions(self):
        self.preRun()
        self.otherActionsToBePerformedBeforeLoop()
        self.actionsToBePerformedInsideLoop()
        
    def validate(self):
        return checkErrors()
    
def ImportProtocol():
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
    import arg
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
        AngSamplingRateDeg=arg.getComponentFromVector(eval(fn +'.AngSamplingRateDeg'),iterationNo)
    global MaxChangeInAngles
    if(len(MaxChangeInAngles) < 1):
        MaxChangeInAngles=arg.getComponentFromVector(eval(fn +'.MaxChangeInAngles'),iterationNo)
    global refDirNameAngSamplingRateDeg
    #refDirName=  eval(fn +'.LibraryDir')
        
if __name__ == '__main__':
    ImportProtocol()
    protocolMain(ProtPartialProjectionSubtraction)
    #import sys
    #ImportProtocol()
    #script  = sys.argv[0]
    #p = ProtPartialProjectionSubtraction(script, WorkingDir, ProjectDir, LogDir, ContinueAtIteration, IsIter)
    #p.run()
