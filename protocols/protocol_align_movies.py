#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Author: Carlos Oscar Sorzano, October 2011

#from config_protocols import protDict
from protlib_base import *
from protlib_utils import which, runJob, runShowJ,printLog
from protlib_filesystem import deleteFile, createLink2, exists, replaceFilenameExt, createDir
import xmipp 
from protlib_gui_ext import showWarning

## The dictionary with specific filename templates 
## is defined here to allow use of it outside the protocol
_templateDict = {
        # This templates are relative to a micrographDir
        'movieAverage': join('%(movieDir)s','%(baseName)s_aligned.spi'),
        'moviesFlows': join('%(movieDir)s','%(baseName)s_flow.spi'),
        'micrographs': join('%(workingDir)s','micrographs.xmd')
        }

def _getFilename(key, **args):
    return _templateDict[key] % args

class ProtAlignMovies(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.align_movies.name, scriptname, project)
        self.Import = "from protocol_align_movies import *"
        self.setPreviousRun(self.ImportRun) 
        self.inputFilename('microscope', 'acquisition','micrographs')        
        self.MicrographsMD = self.Input['micrographs']
        self.micrographs= _getFilename('micrographs',workingDir=self.workingDirPath())

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
            
            self.insertParallelStep('alignSingleMovie'
                                    , verifyfiles=[_getFilename('movieAverage',movieDir=extraDir,baseName=movieName)]
                                    , WorkingDir=self.WorkingDir
                                    , inputMovie=inputMovie 
                                    , movieAverage = [_getFilename('movieAverage',movieDir=extraDir,baseName=movieName)]                                   
                                    , WinSize=self.WinSize
                                    , parent_step_id=XmippProjectDb.FIRST_STEP
                                    , DoGPU = self.DoGPU
                                    )
        
        # Gather results after external actions
        self.insertStep('gatherResults',verifyfiles=[self.micrographs],
                           WorkingDir=self.WorkingDir,
                           extraDir=extraDir,
                           summaryFile=self.micrographs)
    
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
    
def alignSingleMovie(log,WorkingDir
                     , inputMovie   
                     , movieAverage                                 
                     , WinSize
                     , doGPU
                     ):

        # Align estimation with Xmipp
        args= ("%s %s %d")%(inputMovie,movieAverage,WinSize)
        
        print args
        if doGPU:
            args += ' '
        else:
            args += ' '
            
        runJob(log,'xmipp_optical_alignment', args)
                
def gatherResults(log, WorkingDir,extraDir,summaryFile):
    
    import glob
    fnList=glob.glob(os.path.join(WorkingDir,extraDir,'*_aligned.spi'))
    md=xmipp.MetaData()
    for file in fnList:
        id=md.addObject()
        md.setValue(xmipp.MDL_MICROGRAPH,file,id)
#TODO: ADD CTF
    md.sort(xmipp.MDL_MICROGRAPH)
    md.write(summaryFile)
