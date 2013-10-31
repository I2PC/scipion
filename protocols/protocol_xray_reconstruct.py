#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# General script for Xmipp-based reconstruction of X-ray tomograms
# Author: Joaquin Oton, Oct 2013

from glob import glob
from protlib_base import *
import xmipp
from protlib_filesystem import  join, createDir, createLink, removeBasenameExt, removeFilenameExt
from protlib_utils import runJob
from protlib_xmipp import redStr, RowMetaData
import math


class ProtXrayRecon(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.xray_reconstruct.name, scriptname, project)
        self.Import = "from protocol_xray_reconstruct import *"
        self.setPreviousRun(self.ImportRun)
        self.inputFilename('tomograms')
        self.TomogramsMd = self.getFilename('tomograms')
 
    def getTomogramMd(self):
        """ Read and return the Metadata with input tomograms"""
        TomogramsMd = self.Input['tomograms']
        return xmipp.MetaData('tomograms@%s' % (TomogramsMd))
              
    def defineSteps(self):
        md = self.getTomogramMd()
        mdOut = xmipp.MetaData()        
        
        thickness = self.thickness
        recVolList = []
        
        for objId in md:
            fnRootIn = removeFilenameExt(md.getValue(xmipp.MDL_TOMOGRAMMD, objId))
            fnIn = fnRootIn + '.mrc'
            fnBaseName = basename(fnRootIn)
            tomoDir = self.workingDirPath(fnBaseName)
            fnRootOut = join(tomoDir, fnBaseName)
            fnTmpOut = self.tmpPath(fnBaseName + '.mrc')
            fnOut = fnRootOut + '.mrc'
            
            self.insertStep('createDir', path=tomoDir)      
            
            if self.DoInvertContrast:
                fnTmpIn = self.tmpPath(fnBaseName + '_inverted.mrc')
                params = '-i %(fnIn)s --mult -1 -o %(fnTmpIn)s'
                self.insertRunJobStep('xmipp_image_operate', params % locals(), verifyFiles=[fnTmpIn])
                fnIn = fnTmpIn
            
            if self.recProgram == 'imod':
                prog = 'tilt'
                params = '-LOG 0.0 -MODE 2 -PERPENDICULAR -FULLIMAGE 1324,1284 -SCALE 0.0,%(thickness)s '\
                         '-THICKNESS %(thickness)s -RADIAL 0.03,0.05  -ActionIfGPUFails 1,2 -UseGPU 0 '\
                         '-TILTFILE %(fnRootIn)s.tlt %(fnIn)s %(fnTmpOut)s'
            else:
                prog = 'tomo3d'
                params = '-i %(fnIn)s -a %(fnRootIn)s.tlt -o %(fnTmpOut)s -z %(thickness)s'
            
            self.insertRunJobStep(prog, params % locals(), verifyFiles=[fnTmpOut])
            
            params = '-i %(fnTmpOut)s -o %(fnOut)s --face top'
            self.insertRunJobStep('xmipp_volume_reslice', params % locals(), verifyFiles=[fnOut])
            
            recVolList.append(fnOut)
        
        self.insertStep('createResultMd', volList=recVolList, resultMd=self.TomogramsMd, verifyfiles=[self.TomogramsMd])     

        # Removing temporary files
        self.insertDeleteTmpDir()

    def validate(self):
        errors = []
        md = self.getTomogramMd()
        if md.size() < 1:
            errors.append("There are no tomograms to reconstruct in " + self.Input['tomograms'])
        
        return errors

    def summary(self):
        md = self.getTomogramMd()
        size = md.size()
        message = "Reconstruction with <%s> of " % (self.recProgram)
        if size > 1:
            message += "<%d> tomograms" % (size)
        elif size == 1:
            message += "<%s> tomogram" % removeBasenameExt(md.getValue(xmipp.MDL_TOMOGRAMMD, md.firstObject()))
        
        message += ' of thickness <%s>' % self.thickness
        return [message]
    
    def visualize(self):
        resultMd = self.workingDirPath('tomograms.xmd')
        if os.path.exists(resultMd):
            from protlib_utils import runShowJ
            runShowJ(resultMd, extraParams='--mode gallery')
        pass
        
def createResultMd(log, volList, resultMd):
    mdOut = xmipp.MetaData()
     
    for vol in volList:
        mdOut.setValue(xmipp.MDL_TOMOGRAM_VOLUME, vol, mdOut.addObject())
     
    mdOut.write('tomograms@%s' % (resultMd), xmipp.MD_APPEND)
