#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Wrapper for relion 1.2
#
#   Author: Roberto Marabini
#

from os.path import join, exists
from os import remove, rename
from protlib_base import protocolMain
from config_protocols import protDict
from xmipp import *
from xmipp import MetaDataInfo
from protlib_utils import runShowJ, getListFromVector, getListFromRangeString, \
                          runJob, runChimera, which, runChimeraClient
from protlib_parser import ProtocolParser
from protlib_xmipp import redStr, cyanStr
from protlib_gui_ext import showWarning, showTable, showError
from protlib_filesystem import xmippExists, findAcquisitionInfo, moveFile, \
                               replaceBasenameExt
from protlib_filesystem import createLink

from protocol_relion_base import ProtRelionBase, runNormalizeRelion, convertImagesMd, renameOutput, \
                                 convertRelionBinaryData, convertRelionMetadata

class ProtRelionClassifierContinue(ProtRelionBase):
    def __init__(self, scriptname, project):
        ProtRelionBase.__init__(self, protDict.relion_classify_continue.name, scriptname, project)
        self.setPreviousRun(self.ImportRun)
        self.Import = 'from protocol_relion_classify_continue import *'
        self.relionType='classify'
        self.NumberOfClasses      = self.PrevRun.NumberOfClasses
        self.SamplingRate         = self.PrevRun.SamplingRate
        self.MaskDiameterA        = self.PrevRun.MaskDiameterA
        self.RegularisationParamT = self.PrevRun.RegularisationParamT

    def summary(self):
        lines         = ProtRelionBase.summary(self)
        lastIteration = self.lastIter()
        lines += ['Continuation from run: <%s>, iter: <%d>' % (self.PrevRunName, self.ContinueFromIteration)]
        if (lastIteration - self.ContinueFromIteration) < 0 :
            performedIteration=0
        else:
            performedIteration=lastIteration - self.ContinueFromIteration
        lines += ['Performed <%d> iterations (number estimated from the files in working directory)' % performedIteration ]
        #lines += ['Input fileName = <%s>' %self.getFilename('optimiserRe', iter=int(self.ContinueFromIteration))]
        #lines += ['test = <%s>'%self.getFilename('optimiserRe',iter=3)]
        #lines += ['test = <%s>'%self.PrevRun.getFilename('optimiserRe',iter=3)]
        #lines += ['WorkingDir = <%s>'%self.WorkingDir]
        #lines += ['WorkingDir2 = <%s>'%self.PrevRun.WorkingDir]
        #lines += ['lastIter = <%s>'%self.lastIter()]
        #lines += ['lastIter2 = <%s>'%self.PrevRun.lastIter()]
        return lines
    
    def validate(self):
        #errors = ProtRelionBase.validate(self)
        lastIterationPrecRun=self.PrevRun.lastIter()
        errors=[]
        if lastIterationPrecRun < self.ContinueFromIteration:
            errors +=['protocol <%s> last iteration is <%d> so I cannot continue from iteration <%d> '%
                      (self.PrevRunName,lastIterationPrecRun,self.ContinueFromIteration)]
        return errors 
    
    def defineSteps(self): 
        self.insertStep('createDir', verifyfiles=[self.ExtraDir], path=self.ExtraDir)
        # launch relion program
        self.insertRelionClassifyContinue()
        #self.ImgStar = self.extraPath(replaceBasenameExt(tmpFileNameXMD, '.star'))
        lastIteration = self.NumberOfIterations
        NumberOfClasses=self.NumberOfClasses
        firstIteration = self.ContinueFromIteration
        ProtRelionBase.defineSteps2(self, firstIteration
                                        , lastIteration 
                                        , NumberOfClasses)

    def createFilenameTemplates(self):
        myDict=ProtRelionBase.createFilenameTemplates(self)        
        #myDict['volume']=self.extraIter + "class%(ref3d)03d.spi"
        #myDict['volumeMRC']=self.extraIter + "class%(ref3d)03d.mrc:mrc"
        myDict['volume']=self.extraIter + "class%(ref3d)03d.mrc:mrc"
        
        self.relionFiles += ['model']
        #relionFiles=['data','model','optimiser','sampling']
        for v in self.relionFiles:
            myDict[v+'Re']=self.extraIter + v +'.star'
            myDict[v+'Xm']=self.extraIter + v +'.xmd'
        myDict['imagesAssignedToClass']='imagesAssignedToClass@'+self.ExtraDir+'/dataForVisualize.xmd'
        extraIter = join('%(workingDir)s', 'extra', 'relion_it%(iter)03d_optimiser.star')
        myDict['optimiserCont'] = extraIter
        return myDict

    def insertRelionClassifyContinue(self):
        import sys
        _workingDir=join(self.PrevRun.projectDir,self.PrevRun.WorkingDir)
        inputFileName = self.getFilename('optimiserCont', iter=int(self.ContinueFromIteration),workingDir=_workingDir)
        args = {
                '--o': '%s/relion' % self.ExtraDir,
                '--continue': inputFileName,
                
                '--iter': self.NumberOfIterations,
                
                '--tau2_fudge': self.RegularisationParamT,# should not be changed 
                '--flatten_solvent': '',# use always
                '--zero_mask': '',# use always. This is confussing since Sjors gui creates the command line with this option
                                  # but then the program complains about it. 
                '--oversampling': '1',
                '--norm': '',
                '--scale': '',
                }
        iover = 1 #TODO: check this DROP THIS
        index = ['30','15','7.5','3.7','1.8',
                 '0.9','0.5','0.2','0.1'].index(self.AngularSamplingDeg)
        args['--healpix_order'] = float(index + 1 - iover)
        
        if self.PerformLocalAngularSearch:
            args['--sigma_ang'] = self.LocalAngularSearchRange / 3.
            
        args['--offset_range'] = self.OffsetSearchRangePix
        args['--offset_step']  = self.OffsetSearchStepPix * pow(2, iover)

        args['--j'] = self.NumberOfThreads
        
        # Join in a single line all key, value pairs of the args dict    
        params = ' '.join(['%s %s' % (k, str(v)) for k, v in args.iteritems()])
        params += self.AdditionalArguments
        verifyFiles=[]
        #relionFiles=['data','model','optimiser','sampling']
        for v in self.relionFiles:
             verifyFiles += [self.getFilename(v+'Re', iter=self.NumberOfIterations, workingDir=self.WorkingDir )]
        self.insertRunJobStep(self.program, params,verifyFiles)
        #############self.insertRunJobStep('echo shortcut', params,verifyFiles)
