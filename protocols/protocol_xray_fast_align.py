#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# General script for Xmipp-based pre-processing of X-ray Tomograms
# Author: Joaquin Oton, Oct 2013

from glob import glob
from protlib_base import *
import xmipp
from protlib_filesystem import  join, basename, createDir, createLink, removeFilenamePrefix, removeBasenameExt, removeFilenameExt, copyFile
from protlib_utils import runJob, printLog
from protlib_xmipp import redStr, RowMetaData
import math


class ProtXrayFastAlign(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.xray_fast_align.name, scriptname, project)
        self.Import = "from protocol_xray_fast_align import *"
        self.setPreviousRun(self.ImportRun)
        self.inputFilename('tomograms')
        self.TomogramsMd = self.getFilename('tomograms')
        self.tomoAlignedList = []
  
    def getTomogramMd(self):
        """ Read and return the Metadata with input tomograms"""
        TomogramsMd = self.Input['tomograms']
        return xmipp.MetaData('tomograms@%s' % (TomogramsMd))
            
    def defineSteps(self):
        md = self.getTomogramMd()
        mdOut = xmipp.MetaData()
#         for tomoName in blockList:
        for objId in md:
#             md.read('%s@%s' % (tomoName, TomogramsMd))
#             objId = md.firstObject()
            fnRootIn = removeFilenameExt(md.getValue(xmipp.MDL_TOMOGRAMMD, objId))
            fnBaseName = basename(fnRootIn)
            tomoDir = self.workingDirPath(fnBaseName)
            fnRootOut = join(tomoDir, fnBaseName)
            fnRootTmp = self.tmpPath(fnBaseName)
            fnOut = fnRootOut + '.mrc'
            
            self.insertStep('createDir', verifyfiles=[tomoDir], path=tomoDir)
            
            params = '-input %(fnRootIn)s.mrc -output %(fnRootTmp)s.prexf  -tiltfile %(fnRootIn)s.tlt '\
                     '-rotation 0.0 -sigma1 0.03 -radius2 0.25 -sigma2 0.05'
            self.insertRunJobStep('tiltxcorr', params=params % locals(), verifyFiles=[fnRootTmp + '.prexf'])

            
            params = '-input %(fnRootTmp)s.prexf -nfit 0 -goutput %(fnRootTmp)s.prexg' 
            self.insertRunJobStep('xftoxg', params=params % locals(), verifyFiles=[fnRootTmp + '.prexg'])
            
            params = '-input %(fnRootIn)s.mrc -output %(fnRootTmp)s.preali -mode 0 -xform %(fnRootTmp)s.prexg -float 2' 
            self.insertRunJobStep('newstack', params=params % locals(), verifyFiles=[fnRootTmp + '.preali'])
            
            params = '-input %(fnRootTmp)s.preali -output %(fnRootTmp)s.fid -tiltfile %(fnRootIn)s.tlt '\
            '-prexf %(fnRootTmp)s.prexg -rotation 0.0 -sigma1 0.03 -radius2 0.25 -sigma2 0.05 -border 49,49 '\
            '-size 300,300 -overlap 0.33,0.33'
            self.insertRunJobStep('tiltxcorr', params=params % locals(), verifyFiles=[fnRootTmp + '.fid'])
   
   
            tiltAlignParams = {}
            tiltAlignParams['-tiltfile'] = '%(fnRootIn)s.tlt' 
            tiltAlignParams['-ModelFile'] = '%(fnRootTmp)s.fid'
            tiltAlignParams['-ImageFile'] = '%(fnRootTmp)s.preali'
            tiltAlignParams['-OutputTransformFile'] = '%(fnRootTmp)s.tltxf' 
            tiltAlignParams['-OutputLocalFile'] = '%(fnRootTmp)s_local.xf'
            tiltAlignParams['-OutputTiltFile'] = '%(fnRootOut)s.tlt' 
            tiltAlignParams['-RotationAngle'] = 0.0 
            tiltAlignParams['-AngleOffset'] = 0.0 
            tiltAlignParams['-RotOption'] = -1
            tiltAlignParams['-RotDefaultGrouping'] = 5 
            tiltAlignParams['-TiltOption'] = 0 
            tiltAlignParams['-MagReferenceView'] = 1 
            tiltAlignParams['-MagOption'] = 0 
            tiltAlignParams['-MagDefaultGrouping'] = 4 
            tiltAlignParams['-XStretchOption'] = 0 
            tiltAlignParams['-XStretchDefaultGrouping'] = 7 
            tiltAlignParams['-SkewOption'] = 0 
            tiltAlignParams['-SkewDefaultGrouping'] = 11 
            tiltAlignParams['-ResidualReportCriterion'] = 3.0  
            tiltAlignParams['-SurfacesToAnalyze'] = 1 
            tiltAlignParams['-MetroFactor'] = 0.25 
            tiltAlignParams['-MaximumCycles'] = 1000 
            tiltAlignParams['-AxisZShift'] = 0.0 
            tiltAlignParams['-LocalAlignments'] = 0 
            tiltAlignParams['-MinFidsTotalAndEachSurface'] = '8,3' 
            tiltAlignParams['-LocalOutputOptions'] = '0,0,0' 
            tiltAlignParams['-LocalRotOption'] = 3 
            tiltAlignParams['-LocalRotDefaultGrouping'] = 6 
            tiltAlignParams['-LocalTiltOption'] = 5 
            tiltAlignParams['-LocalTiltDefaultGrouping'] = 6 
            tiltAlignParams['-LocalMagReferenceView'] = 1 
            tiltAlignParams['-LocalMagOption'] = 3 
            tiltAlignParams['-LocalMagDefaultGrouping'] = 7 
            tiltAlignParams['-LocalXStretchOption'] = 0 
            tiltAlignParams['-LocalXStretchDefaultGrouping'] = 7 
            tiltAlignParams['-LocalSkewOption'] = 0
            tiltAlignParams['-LocalSkewDefaultGrouping'] = 11 
            tiltAlignParams['-BeamTiltOption'] = 0 
   
            params = ''
            for k, v in tiltAlignParams.iteritems():
                params += " %s %s " % (k, str(v))
   
            self.insertRunJobStep('tiltalign', params=params % locals(), verifyFiles=[fnRootTmp + '.tltxf', fnRootTmp + '_local.xf', fnRootOut + '.tlt'])
   
            params = '-in1 %(fnRootTmp)s.prexg -in2 %(fnRootTmp)s.tltxf -output %(fnRootTmp)s_fid.xf' 
            self.insertRunJobStep('xfproduct', params=params % locals(), verifyFiles=[fnRootTmp + '_fid.xf'])
            
            self.insertStep('copyFile', source=fnRootTmp + '_fid.xf', dest=fnRootTmp + '.xf', verifyfiles=[fnRootTmp + '.xf'])

            params = '-input %(fnRootIn)s.mrc -output %(fnOut)s -offset 0,0 -xform %(fnRootTmp)s_fid.xf -origin -taper 1,0' 
            self.insertRunJobStep('newstack', params=params % locals(), verifyFiles=[fnOut])
            
            # Create metadata 
            self.insertStep('createtMd', tomoOut=fnOut, tomoRoot=fnRootOut,fnTiltIn=fnRootIn+'.tlt',fnTiltOut=fnRootOut+'.tlt', verifyfiles=[fnRootOut + '.xmd'])
            # Add aligned tomo to list
            self.tomoAlignedList.append(fnRootOut)
            
        self.insertStep('createResultMd', TomogramList=self.tomoAlignedList, resultMd=self.TomogramsMd, verifyfiles=[self.TomogramsMd])    
        
        # Removing temporary files
        self.insertDeleteTmpDir()
         


    def validate(self):
        errors = []
        md = self.getTomogramMd()
        if md.size() < 1:
            errors.append("There are no tomograms to process in " + self.Input['tomograms'])
        
        return errors

    def summary(self):
        md = self.getTomogramMd()
        size = md.size()
        message = []
        if size > 1:
            message.append("Alignment of <%d> tomograms" % size)
        elif size == 1:
            message.append("Alignment of <%s> tomogram" % removeBasenameExt(md.getValue(xmipp.MDL_TOMOGRAMMD, md.firstObject())))

        return message
    
    def visualize(self):
        resultMd = self.workingDirPath('tomograms.xmd')
        if os.path.exists(resultMd):
            from protlib_utils import runShowJ
            runShowJ(resultMd, extraParams='--mode gallery')
        pass
    
    
# def tiltAlign(log, fnRootIn, fnRootOut, fnRootTmp):
#              
#     tiltAlignParams = {}
#     tiltAlignParams['-tiltfile'] = '%(fnRootIn)s.tlt' 
#     tiltAlignParams['-ModelFile'] = '%(fnRootTmp)s.fid'
#     tiltAlignParams['-ImageFile'] = '%(fnRootTmp)s.preali'
#     tiltAlignParams['-OutputTransformFile'] = '%(fnRootTmp)s.tltxf' 
#     tiltAlignParams['-OutputLocalFile'] = '%(fnRootTmp)s_local.xf'
#     tiltAlignParams['-OutputTiltFile'] = '%(fnRootOut)s.tlt' 
#     tiltAlignParams['-RotationAngle'] = 0.0 
#     tiltAlignParams['-AngleOffset'] = 0.0 
#     tiltAlignParams['-RotOption'] = -1
#     tiltAlignParams['-RotDefaultGrouping'] = 5 
#     tiltAlignParams['-TiltOption'] = 0 
#     tiltAlignParams['-MagReferenceView'] = 1 
#     tiltAlignParams['-MagOption'] = 0 
#     tiltAlignParams['-MagDefaultGrouping'] = 4 
#     tiltAlignParams['-XStretchOption'] = 0 
#     tiltAlignParams['-XStretchDefaultGrouping'] = 7 
#     tiltAlignParams['-SkewOption'] = 0 
#     tiltAlignParams['-SkewDefaultGrouping'] = 11 
#     tiltAlignParams['-ResidualReportCriterion'] = 3.0  
#     tiltAlignParams['-SurfacesToAnalyze'] = 1 
#     tiltAlignParams['-MetroFactor'] = 0.25 
#     tiltAlignParams['-MaximumCycles'] = 1000 
#     tiltAlignParams['-AxisZShift'] = 0.0 
#     tiltAlignParams['-LocalAlignments'] = 0 
#     tiltAlignParams['-MinFidsTotalAndEachSurface'] = '8,3' 
#     tiltAlignParams['-LocalOutputOptions'] = '0,0,0' 
#     tiltAlignParams['-LocalRotOption'] = 3 
#     tiltAlignParams['-LocalRotDefaultGrouping'] = 6 
#     tiltAlignParams['-LocalTiltOption'] = 5 
#     tiltAlignParams['-LocalTiltDefaultGrouping'] = 6 
#     tiltAlignParams['-LocalMagReferenceView'] = 1 
#     tiltAlignParams['-LocalMagOption'] = 3 
#     tiltAlignParams['-LocalMagDefaultGrouping'] = 7 
#     tiltAlignParams['-LocalXStretchOption'] = 0 
#     tiltAlignParams['-LocalXStretchDefaultGrouping'] = 7 
#     tiltAlignParams['-LocalSkewOption'] = 0
#     tiltAlignParams['-LocalSkewDefaultGrouping'] = 11 
#     tiltAlignParams['-BeamTiltOption'] = 0 
#      
#      
#     params = ''
#     for k, v in tiltAlignParams.iteritems():
#         params += " %s %s " % (k, str(v))
#      
#     command = 'tiltalign ' + params % locals()
#     
#     from protlib_xmipp import greenStr
#     
#     printLog("Running command: %s" % greenStr(command),log)
#     
#     try:
#         retcode = call(command, shell=True, stdout=sys.stdout, stderr=sys.stderr)
#         if log:
#             printLog("Process returned with code %d" % retcode,log)
#             if retcode != 0:
#                 raise Exception("Process returned with code %d, command: %s" % (retcode,command))
#     except OSError, e:
#         raise Exception("Execution failed %s, command: %s" % (e, command))
# 
# 
# #     os.system('tiltalign ' + params % locals())
    
    
def createtMd(log, tomoOut, tomoRoot, fnTiltIn, fnTiltOut):
    md = xmipp.MetaData()
    md.read(tomoOut)
    md.addPlain(fnTiltOut, 'angleTilt')
    md.addPlain(fnTiltIn, 'angleTilt2')
    md.write("tomo@%s.xmd" % tomoRoot)  
        

def createResultMd(log, TomogramList, resultMd):
#     md = xmipp.MetaData()
    mdOut = xmipp.MetaData()
     
    for tomogram in TomogramList:
#         tomoRoot = removeFilenameExt(tomogram)
#         tomoBaseName = basename(tomoRoot)
#         md.read(tomoRoot + ".xmd")
#         md.write('tomo_%s@%s' % (tomoBaseName, resultMd), xmipp.MD_APPEND)
        mdOut.setValue(xmipp.MDL_TOMOGRAMMD, tomogram+'.xmd', mdOut.addObject())
     
    mdOut.write('tomograms@%s' % (resultMd), xmipp.MD_APPEND)
