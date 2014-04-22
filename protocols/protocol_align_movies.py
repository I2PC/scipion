#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Author: Carlos Oscar Sorzano, October 2011

#from config_protocols import protDict
from protlib_base import *
from protlib_utils import which, runJob, runShowJ,printLog
from protlib_filesystem import deleteFile, createLink2, exists, replaceFilenameExt, createDir, copyFile
import xmipp 
from protlib_gui_ext import showWarning

## The dictionary with specific filename templates 
## is defined here to allow use of it outside the protocol
_templateDict = {
        # This templates are relative to a micrographDir
        'movieAverage': join('%(movieDir)s','%(baseName)s_aligned.spi'),
        'moviesFlows': join('%(movieDir)s','%(baseName)s_flow.spi'),
        }

def _getFilename(key, **args):
    return _templateDict[key] % args

class ProtAlignMovies(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.align_movies.name, scriptname, project)
        self.Import = "from protocol_align_movies import *"
        self.setPreviousRun(self.ImportRun) 
        self.inputFilename('microscope', 'acquisition','micrographs')        
        self.MicrographsMD   = self.Input['micrographs']
        self.MicroscopeMD     = self.Input['microscope']
        self.AcquisitionMD   = self.Input['acquisition']
        self.micrographs = self.getFilename('micrographs',workingDir=self.workingDirPath())
        self.micrographsClass = 'classes@'+self.getFilename('micrographs',workingDir=self.workingDirPath())
        self.microscope  = self.getFilename('microscope',workingDir=self.workingDirPath())
        self.acquisition = self.getFilename('acquisition',workingDir=self.workingDirPath())

    def defineSteps(self):
        extraDir=self.workingDirPath('extra')
        parent_id = self.insertStep('createDir',verifyfiles=[extraDir],path=extraDir)
        # Create verifyFiles for the MPI and output directories
        MD = xmipp.MetaData(self.MicrographsMD)
        #if removed in import do not process them
        MD.removeDisabled()
    
        # Now the estimation actions
        for objId in MD:
            inputMovie = MD.getValue(xmipp.MDL_MICROGRAPH_MOVIE, objId)
            movieNameList = os.path.basename(inputMovie)
            movieName = os.path.splitext(movieNameList)[0]
            #micrographDir = os.path.join(extraDir,shortname)                    
            if self.DoGPU:
                func = self.insertStep
            else:
                func = self.insertParallelStep   
            func('alignSingleMovie1'
                                    , verifyfiles=[_getFilename('movieAverage',movieDir=extraDir,baseName=movieName)]
                                    , WorkingDir=self.WorkingDir
                                    , inputMovie=inputMovie 
                                    , movieAverage = _getFilename('movieAverage',movieDir=extraDir,baseName=movieName)                                   
                                    , WinSize=self.WinSize
                                    , DoGPU = self.DoGPU
                                    , GPUCore = self.GPUCore
                                    , parent_step_id=XmippProjectDb.FIRST_STEP
                                    )
        
        # Gather results after external actions
        self.insertStep('gatherResults',verifyfiles=[self.micrographs,self.microscope,self.acquisition],
                           WorkingDir=self.WorkingDir,
                           extraDir=extraDir,
                           currProtMicrographClass=self.micrographsClass,
                           currProtMicrograph=self.micrographs,
                           currProtMicroscope=self.microscope,
                           currProtAcquisitionInfo=self.acquisition,
                           prevProtMicrograph=self.MicrographsMD,
                           prevProtMicroscope=self.MicroscopeMD,
                           prevProtAcquisitionInfo=self.AcquisitionMD
                           )
    
    def createFilenameTemplates(self):
        return _templateDict
    
    def validate(self):
        #TODO check if GPU is available
        errors = []
        return errors

    def summary(self):
        message = []
        from protlib_xmipp import getMdSize
        size = getMdSize(self.Input['micrographs'])
        message.append("Movies aligned <%d>" % (size))
        message.append("Input directory: [%s]" % self.PrevRun.WorkingDir)
        return message
    
    def papers(self):
        papers=[]
        return papers
    
    def visualize(self):
        summaryFile = self.getFilename('micrographs')
        if exists(summaryFile):
            runShowJ(summaryFile, extraParams = "--mode metadata")
        else:
            showWarning('Warning', 'There are not results yet',self.master)
    

def alignSingleMovie1(log,WorkingDir
                     , inputMovie   
                     , movieAverage                                 
                     , WinSize
                     , DoGPU
                     , GPUCore
                     ):

        # Align estimation with Xmipp
        args = '-i %s -o %s --winSize %d'%(inputMovie,movieAverage,WinSize)
        if DoGPU:
            progName = 'xmipp_optical_alignment_gpu'
            args += ' --gpu %d'%(GPUCore)
        else:
            progName = 'xmipp_optical_alignment_cpu'

        runJob(log,progName, args)
                
def gatherResults(log, 
                  WorkingDir,
                  extraDir,
                  currProtMicrographClass,
                  currProtMicrograph,
                  currProtMicroscope,
                  currProtAcquisitionInfo,
                  prevProtMicrograph,
                  prevProtMicroscope,
                  prevProtAcquisitionInfo
                  ):
    
    import glob
    fnList=glob.glob(os.path.join(extraDir,'*_aligned.spi'))
    copyFile(log,prevProtMicrograph,currProtMicrograph) 
    copyFile(log,prevProtMicroscope,currProtMicroscope) 
    copyFile(log,prevProtAcquisitionInfo,currProtAcquisitionInfo) 
    md=xmipp.MetaData(currProtMicrograph)
    for id  in md:
        inputMovie = md.getValue(xmipp.MDL_MICROGRAPH, id)
        movieNameList = os.path.basename(inputMovie)
        movieName = os.path.splitext(movieNameList)[0]
        file = _getFilename('movieAverage',movieDir=extraDir,baseName=movieName) 
        md.setValue(xmipp.MDL_MICROGRAPH,file,id)
    md.sort(xmipp.MDL_MICROGRAPH)
    md.write(currProtMicrographClass,xmipp.MD_APPEND)
