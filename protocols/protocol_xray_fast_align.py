#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# General script for Xmipp-based pre-processing of micrographs
# Author: Joaquin Oton, Oct 2013

from glob import glob
from protlib_base import *
import xmipp
from protlib_filesystem import  join, basename, createDir, createLink, removeFilenamePrefix, removeBasenameExt, removeFilenameExt, copyFile
from protlib_utils import runJob
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
            fnRootIn = removeFilenameExt(md.getValue(xmipp.MDL_IMAGE, objId))
            fnBaseName = basename(fnRootIn)
            tomoDir = self.workingDirPath(fnBaseName)
            fnRootOut = join(tomoDir, fnBaseName)
            
            self.insertStep('createDir', path=tomoDir)
            
            params = '-input %(fnRootIn)s.mrc -output %(fnRootOut)s.prexf  -tiltfile %(fnRootIn)s.tlt '\
                     '-rotation 0.0 -sigma1 0.03 -radius2 0.25 -sigma2 0.05'
            self.insertRunJobStep('tiltxcorr', params % locals())

            
            params = '-input %(fnRootOut)s.prexf -nfit 0 -goutput %(fnRootOut)s.prexg' 
            self.insertRunJobStep('xftoxg', params % locals())
            
            params = '-input %(fnRootIn)s.mrc -output %(fnRootOut)s.preali -mode 0 -xform %(fnRootOut)s.prexg -float 2' 
            self.insertRunJobStep('newstack', params % locals())
            
            params = '-input %(fnRootOut)s.preali -output %(fnRootOut)s.fid -tiltfile %(fnRootIn)s.tlt '\
            '-prexf %(fnRootOut)s.prexg -rotation 0.0 -sigma1 0.03 -radius2 0.25 -sigma2 0.05 -border 49,49 '\
            '-size 300,300 -overlap 0.33,0.33'
            self.insertRunJobStep('tiltxcorr', params % locals())
   
            self.insertStep('tiltAlign', fnRootIn=fnRootIn, fnRootOut=fnRootOut)
   
            params = '-in1 %(fnRootOut)s.prexg -in2 %(fnRootOut)s.tltxf -output %(fnRootOut)s_fid.xf' 
            self.insertRunJobStep('xfproduct', params % locals())
            
            self.insertStep('copyFile', source=fnRootOut + '_fid.xf', dest=fnRootOut + '.xf')

            params = '-input %(fnRootIn)s.mrc -output %(fnRootOut)s_ali.mrc -offset 0,0 -xform %(fnRootOut)s_fid.xf -origin -taper 1,0' 
            self.insertRunJobStep('newstack', params % locals())
            
            # Create metadata 
            self.insertStep('createtMd', tomoRoot=fnRootOut)
            # Add aligned tomo to list
            self.tomoAlignedList.append(fnRootOut)
            
        self.insertStep('createResultMd', TomogramList=self.tomoAlignedList, resultMd=self.TomogramsMd)     


    def validate(self):
        errors = []
        
        return errors

    def summary(self):
        size = self.getTomogramMd().size()
        message = "Alignment of <%d> tomogram" % size
        if size > 1:
            message += 's' # plural

        return [message]
    
    def visualize(self):
        resultMd = self.workingDirPath('tomograms.xmd')
        if os.path.exists(resultMd):
            from protlib_utils import runShowJ
            runShowJ(resultMd, extraParams='--mode gallery')
        pass
    
    
def tiltAlign(log, fnRootIn, fnRootOut):
             
    tiltAlignParams = {}
    tiltAlignParams['-tiltfile'] = '%(fnRootIn)s.tlt' 
    tiltAlignParams['-ModelFile'] = '%(fnRootOut)s.fid'
    tiltAlignParams['-ImageFile'] = '%(fnRootOut)s.preali'
    tiltAlignParams['-ModelFile'] = '%(fnRootOut)s.fid' 
    tiltAlignParams['-ImageFile'] = '%(fnRootOut)s.preali' 
    tiltAlignParams['-OutputTiltFile'] = '%(fnRootOut)s_ali.tlt' 
    tiltAlignParams['-OutputTransformFile'] = '%(fnRootOut)s.tltxf' 
    tiltAlignParams['-OutputLocalFile'] = '%(fnRootOut)s_local.xf'
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
     
    print params
    print params % locals() 
    
    os.system('tiltalign ' + params % locals())
    
    
def createtMd(log, tomoRoot):
    md = xmipp.MetaData()
    md.read(tomoRoot + '_ali.mrc')
    md.write(tomoRoot + '.xmd')  
        

def createResultMd(log, TomogramList, resultMd):
    md = xmipp.MetaData()
    mdOut = xmipp.MetaData()
     
    for tomogram in TomogramList:
        tomoBaseName = basename(tomogram)
        mdOut.setValue(xmipp.MDL_IMAGE, tomogram + '_ali.mrc', mdOut.addObject())
        md.read(tomogram + ".xmd")
        md.write('tomo_%s@%s' % (tomoBaseName, resultMd), xmipp.MD_APPEND)
     
    mdOut.write('tomograms@%s' % (resultMd), xmipp.MD_APPEND)
