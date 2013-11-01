#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# General script for Xmipp-based pre-processing of micrographs
# Author: Joaquin Oton, Oct 2013

from glob import glob
from protlib_base import *
import xmipp
from protlib_filesystem import  join, createDir, createLink, removeBasenameExt
from protlib_utils import runJob
from protlib_xmipp import redStr, RowMetaData
import math


class ProtXrayImport(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.xray_import.name, scriptname, project)
        self.Import = "from protocol_xray_import import *"
        self.TomogramsMd = self.getFilename('tomograms')
           
    def defineSteps(self):
        tomograms = self.getTomograms()
        
        for t in tomograms:
            self.insertTomogramStep(t)
            
        self.insertStep('createResultMd',
                        WorkingDir=self.WorkingDir, 
                        TomogramList=tomograms, resultMd=self.TomogramsMd, verifyfiles=[self.TomogramsMd])


    def validate(self):
        errors = []
        
        tomogramList = self.getTomograms()
        if len(tomogramList) == 0:
            if self.Synchrotron == "Alba":
                errors.append("There are no tomograms to process in " + self.DirTomograms)   

        return errors

    def summary(self):
        message = []
        
        if self.Synchrotron == "Alba":
            if self.ImportFrom == "single file":
                message.append("Import of [%s] tomogram" % self.Tomogram)   
            else:
                message.append("Import of <%d> tomograms from [%s]" % (len(self.getTomograms()), self.DirTomograms))   
        else:
            message.append("Import of Bessy tomogram from [%s] and initial index %s" % (self.DirBessyData, self.TIni))

    
           
        
        return message
    
    def getTomograms(self):
        ''' Return a list with micrographs in WorkingDir'''
        if self.Synchrotron == "Alba":
            if self.ImportFrom == "single file":
                return [self.Tomogram]
            else:
                return glob(join(self.DirTomograms, '*.hdf*'))
        else:
            return [self.DirBessyData]    
    
    def visualize(self):
        resultMd = self.workingDirPath('tomograms.xmd')
        if os.path.exists(resultMd):
            from protlib_utils import runShowJ
            runShowJ(resultMd)
        pass
        
    def insertTomogramStep(self, tomogram):
        self.ParamsDict['tomogram'] = tomogram
        
        if self.Synchrotron == "Alba":
            params = '--mistral "%(tomogram)s" '
        else:
            params = '--bessy "%(tomogram)s" %(TIni)s %(TEnd)s %(FIni)s %(FEnd)s '
        print params

        if self.DoLog:
            params += '--log '
        if self.DoCorrect:
            params += '--correct '
        if self.DoCrop:
            params += '--crop %(Crop)d '
        if self.DoBadPixelsMask:
            params += '--bad_pixels_filter "%(BadPixelsMask)s" '

        # At this moment hdf5 library does not support parallel access
        if self.Synchrotron != "Alba":
            params += '--thr "%(NumberOfThreads)s" '
    
        _, tomoDir, tomoRoot = getTomoDirs(self.WorkingDir, tomogram)
        params += '--oroot "%s" ' % tomoRoot
        
        self.insertStep("importTomogram", 
                        TomogramDir=tomoDir,
                        TomogramRoot=tomoRoot,
                        params=params % self.ParamsDict, verifyfiles=[tomoRoot+'.mrc'])
       

def getTomoDirs(WorkingDir, tomogram):
    """ Give the tomogram filename, return the root and prefix
    to store the result files. """
    tomoBaseName = removeBasenameExt(tomogram)
    tomoDir = join(WorkingDir, tomoBaseName)
    tomoRoot = join(tomoDir, tomoBaseName)
    
    return tomoBaseName, tomoDir, tomoRoot
      

def importTomogram(log, TomogramDir, TomogramRoot, params):
    createDir(log, TomogramDir)
    runJob(log, "xmipp_xray_import", params, NumberOfMpi=1)
#     createLink(log, TomogramRoot + '.mrc', TomogramRoot + '.st')
    
def createResultMd(log, WorkingDir, TomogramList, resultMd):
    md = xmipp.MetaData()
    mdOut = xmipp.MetaData()
    
    for tomogram in TomogramList:
        tomoBaseName, _, tomoRoot = getTomoDirs(WorkingDir, tomogram)
        mdOut.setValue(xmipp.MDL_TOMOGRAMMD, tomoRoot + '.xmd', mdOut.addObject())
        
#         blockList = xmipp.getBlocksInMetaDataFile(tomoRoot+'.xmd')
#         for block in blockList:
#             md.read("%(block)s@%(tomoRoot)s.xmd" % locals())
#             md.write('%(block)s_%(tomoBaseName)s@%(resultMd)s' % locals(), xmipp.MD_APPEND)
    
    mdOut.write('tomograms@%s' % (resultMd), xmipp.MD_APPEND)
        
