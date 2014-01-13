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
            
            tiltXcorrParams = {}
            tiltXcorrParams['-input'] = '%(fnRootIn)s.mrc'
            tiltXcorrParams['-output'] = '%(fnRootTmp)s.prexf'
            tiltXcorrParams['-tiltfile'] = '%(fnRootIn)s.tlt'
            tiltXcorrParams['-rotation'] = self.rotation
            tiltXcorrParams['-sigma1'] = self.sigma1
            tiltXcorrParams['-radius2'] = self.radius2 
            tiltXcorrParams['-sigma2'] = self.sigma2
            params = ''
            for k, v in tiltXcorrParams.iteritems():
                params += " %s %s " % (k, str(v))
            
#             params = '-input %(fnRootIn)s.mrc -output %(fnRootTmp)s.prexf  -tiltfile %(fnRootIn)s.tlt '\
#                      '-rotation 0.0 -sigma1 0.03 -radius2 0.25 -sigma2 0.05'
            self.insertRunJobStep('tiltxcorr', params=params % locals(), verifyFiles=[fnRootTmp + '.prexf'])

            
            params = '-input %(fnRootTmp)s.prexf -nfit 0 -goutput %(fnRootTmp)s.prexg' 
            self.insertRunJobStep('xftoxg', params=params % locals(), verifyFiles=[fnRootTmp + '.prexg'])
            
            params = '-input %(fnRootIn)s.mrc -output %(fnRootTmp)s.preali -mode 0 -xform %(fnRootTmp)s.prexg -float 2' 
            self.insertRunJobStep('newstack', params=params % locals(), verifyFiles=[fnRootTmp + '.preali'])
            
            
            tiltXcorrParams['-input'] = '%(fnRootTmp)s.preali'
            tiltXcorrParams['-output'] = '%(fnRootTmp)s.fid '
            tiltXcorrParams['-tiltfile'] = '%(fnRootIn)s.tlt'
            tiltXcorrParams['-prexf'] = '%(fnRootTmp)s.prexg'
            tiltXcorrParams['-border'] = self.border
            tiltXcorrParams['-size'] = self.size
            tiltXcorrParams['-overlap'] = self.overlap
            params = ''
            for k, v in tiltXcorrParams.iteritems():
                params += " %s %s " % (k, str(v))
            
#             params = '-input %(fnRootTmp)s.preali -output %(fnRootTmp)s.fid -tiltfile %(fnRootIn)s.tlt '\
#             '-prexf %(fnRootTmp)s.prexg -rotation 0.0 -sigma1 0.03 -radius2 0.25 -sigma2 0.05 -border 49,49 '\
#             '-size 300,300 -overlap 0.33,0.33'
            self.insertRunJobStep('tiltxcorr', params=params % locals(), verifyFiles=[fnRootTmp + '.fid'])
   
   
            tiltAlignParams = {}
            tiltAlignParams['-tiltfile'] = '%(fnRootIn)s.tlt' 
            tiltAlignParams['-ModelFile'] = '%(fnRootTmp)s.fid'
            tiltAlignParams['-ImageFile'] = '%(fnRootTmp)s.preali'
            tiltAlignParams['-OutputTransformFile'] = '%(fnRootTmp)s.tltxf' 
            tiltAlignParams['-OutputLocalFile'] = '%(fnRootTmp)s_local.xf'
            tiltAlignParams['-OutputTiltFile'] = '%(fnRootOut)s.tlt' 
            tiltAlignParams['-RotationAngle'] = self.RotationAngle 
            tiltAlignParams['-AngleOffset'] = self.AngleOffset
            tiltAlignParams['-RotOption'] = self.RotOption
            tiltAlignParams['-RotDefaultGrouping'] = self.RotDefaultGrouping
            tiltAlignParams['-TiltOption'] = self.TiltOption
            tiltAlignParams['-MagReferenceView'] = self.MagReferenceView
            tiltAlignParams['-MagOption'] = self.MagOption
            tiltAlignParams['-MagDefaultGrouping'] = self.MagDefaultGrouping
            tiltAlignParams['-XStretchOption'] = self.XStretchOption
            tiltAlignParams['-XStretchDefaultGrouping'] = self.XStretchDefaultGrouping
            tiltAlignParams['-SkewOption'] = self.SkewOption
            tiltAlignParams['-SkewDefaultGrouping'] = self.SkewDefaultGrouping
            tiltAlignParams['-ResidualReportCriterion'] = self.ResidualReportCriterion
            tiltAlignParams['-SurfacesToAnalyze'] = self.SurfacesToAnalyze
            tiltAlignParams['-MetroFactor'] = self.MetroFactor
            tiltAlignParams['-MaximumCycles'] = self.MaximumCycles
            tiltAlignParams['-AxisZShift'] = self.AxisZShift
            tiltAlignParams['-LocalAlignments'] = self.LocalAlignments
            tiltAlignParams['-MinFidsTotalAndEachSurface'] = self.MinFidsTotalAndEachSurface
            tiltAlignParams['-LocalOutputOptions'] = self.LocalOutputOptions
            tiltAlignParams['-LocalRotOption'] = self.LocalRotOption
            tiltAlignParams['-LocalRotDefaultGrouping'] = self.LocalRotDefaultGrouping
            tiltAlignParams['-LocalTiltOption'] = self.LocalTiltOption
            tiltAlignParams['-LocalTiltDefaultGrouping'] = self.LocalTiltDefaultGrouping
            tiltAlignParams['-LocalMagReferenceView'] = self.LocalMagReferenceView
            tiltAlignParams['-LocalMagOption'] = self.LocalMagOption
            tiltAlignParams['-LocalMagDefaultGrouping'] = self.LocalMagDefaultGrouping
            tiltAlignParams['-LocalXStretchOption'] = self.LocalXStretchOption
            tiltAlignParams['-LocalXStretchDefaultGrouping'] = self.LocalXStretchDefaultGrouping
            tiltAlignParams['-LocalSkewOption'] = self.LocalSkewOption
            tiltAlignParams['-LocalSkewDefaultGrouping'] = self.LocalSkewDefaultGrouping
            tiltAlignParams['-BeamTiltOption'] = self.BeamTiltOption
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
            self.insertStep('createtMd', tomoOut=fnOut, tomoRoot=fnRootOut, fnTiltIn=fnRootIn + '.tlt', fnTiltOut=fnRootOut + '.tlt', verifyfiles=[fnRootOut + '.xmd'])
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
    
    def papers(self):
        papers = []
        papers.append('IMOD: Kremer, JSB (1997) [http://http://www.sciencedirect.com/science/article/pii/S1047847796900131]')
        return papers
    

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
        mdOut.setValue(xmipp.MDL_TOMOGRAMMD, tomogram + '.xmd', mdOut.addObject())
     
    mdOut.write('tomograms@%s' % (resultMd), xmipp.MD_APPEND)
