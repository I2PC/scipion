#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using maximum-likelihood principles
#
#   Author:  Sjors Scheres, January 2008
#  Updated:  J. M. de la Rosa Trevin July 2011
#

from xmipp import MetaData, MetaDataInfo, MDL_ITER, MDL_LL, MDL_REF, MDValueEQ, getBlocksInMetaDataFile, \
MDL_PMAX, MDL_SIGNALCHANGE, AGGR_MAX, MDL_MAX, MDL_MIRRORFRAC, MDL_WEIGHT, MDL_CLASS_COUNT, MD_APPEND,\
getImageSize, MDL_IMAGE, MDL_SAMPLINGRATE
from protlib_base import XmippProtocol, protocolMain
from config_protocols import protDict
import os
from os.path import exists, join
from protlib_utils import runShowJ
from protlib_gui_ext import showWarning
from protlib_xmipp import greenStr, redStr
from protlib_filesystem import deleteFile, xmippExists, renameFile, findAcquisitionInfo

def lastIteration(self, key='iter_logs'):
    ''' Find the last iteration number '''
    iter = 0        
    while True:
        if not xmippExists(self.getFilename(key, iter=iter+1)):
            break
        iter = iter + 1
    return iter

class ProtML2D(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.ml2d.name, scriptname, project)
        self.Import = 'from protocol_ml2d import *'
        self.acquisionInfo = self.findAcquisitionInfo(self.ImgMd)
        if not self.acquisionInfo is None: 
            md = MetaData(self.acquisionInfo)
            self.SamplingRate = md.getValue(MDL_SAMPLINGRATE, md.firstObject())
        
    def createFilenameTemplates(self):
        extra = self.workingDirPath('extra')
        extraIter = join(extra, 'iter%(iter)03d')
        extraRefs = join(extraIter, "iter_classes.xmd")
        return {
                'iter_logs': "info@" + extraRefs,
                'iter_refs': "classes@" + extraRefs
                }
        #self.fnIterLogs = self.workingDirPath('ml2d_iter_logs.xmd')
        #self.fnIterRefs = self.workingDirPath('ml2d_iter_refs.xmd')
        
            
    def summary(self):
        _, _, _, _, size = MetaDataInfo(self.ImgMd)
        lines = ["Input images:  [%s] (<%u>)" % (self.ImgMd, size)]

        if self.DoMlf:
            if self.DoCorrectAmplitudes:
                suffix = "with CTF correction "
            else:
                suffix = "ignoring CTF effects "
            lines.append("Using a ML in <Fourier-space> " + suffix)
        
        if self.DoGenerateReferences:
            lines.append("Number of references: <%d>" % self.NumberOfReferences)
        else:
            lines.append("Reference image(s): [%s]" % self.RefMd)
        
        
        #logs = self.getFilename('iter_logs')   
        lastIter = lastIteration(self) 
        if lastIter > 0:#exists(logs):
            logs = self.getFilename('iter_logs', iter=lastIter)
            md = MetaData(logs)
            id = md.lastObject()
            #iteration = md.getValue(MDL_ITER, id)
            lines.append("Last iteration:  <%d>" % lastIter)
            LL = md.getValue(MDL_LL, id)
            lines.append("LogLikelihood:  %f" % LL)
            mdRefs = self.getFilename('iter_refs', iter=lastIter)
            lines.append("Last classes: [%s]" % mdRefs)

        return lines
    
    def validate(self):
        errors = []
        
        md = MetaData(self.ImgMd)
        if md.containsLabel(MDL_IMAGE):
            # If using reference check that have same size as images:
            if not self.DoGenerateReferences:
                from protlib_xmipp import validateInputSize
                mdRef = MetaData(self.RefMd)
                references = mdRef.getColumnValues(MDL_IMAGE)
                validateInputSize(references, self.ImgMd, errors)
        else:
            errors.append("Input metadata <%s> doesn't contains image label" % self.ImgMd)

        if self.acquisionInfo is None:
            errors.append("Acquisition info was not found relative to image location")
            
        return errors 
    
    def getId(self):
        progId = "ml"
        if (self.DoMlf):
            progId += "f" 
        return progId
        
    def defineSteps(self):
        self.insertStep("linkAcquisitionInfo",InputFile=self.ImgMd,dirDest=self.WorkingDir)
        progId = self.getId()
        
        program = "xmipp_%s_align2d" % progId

        restart = False
        if (restart):
            pass 
        else: 
            # Dictionary with boolean options and the cmd options
            booleanDict = {'DoMirror': '--mirror', 'DoNorm': '--norm'}
            
            prefix = '%s2d' % progId
            oroot = self.workingDirPath(prefix)
            params = ' -i %s --oroot %s/' % (self.ImgMd, self.WorkingDir)
            # Number of references will be ignored if -ref is passed as expert option
            if self.DoGenerateReferences:
                params += ' --nref %d' % self.NumberOfReferences
            else:
                params += ' --ref %s' % self.RefMd
            
            if (self.DoFast and not self.DoMlf):
                params += ' --fast'
            if (self.NumberOfThreads > 1  and not self.DoMlf):
                params += ' --thr %i' % self.NumberOfThreads
            if (self.DoMlf):
                if not self.DoCorrectAmplitudes:
                    params += ' --no_ctf'                    
                if (not self.ImagesArePhaseFlipped):
                    params += ' --not_phase_flipped'
                if (self.HighResLimit > 0):
                    params += ' --limit_resolution 0 %f' % self.HighResLimit
                params += ' --sampling_rate %f' % self.SamplingRate
            if self.MaxIters != 100:
                params += " --iter %d" % self.MaxIters
            #Add all boolean options if true
            for k, v in booleanDict.iteritems():
                if getattr(self, k):
                    params += " " + v
            
            #Add extra options
            #params += ' ' + self.ExtraParams
        self.insertRunJobStep(program, params, 
                              [self.getFilename(k) for k in ['images', 'classes']])

#        self.insertStep('collectResults', WorkingDir=self.WorkingDir, Prefix=prefix,
#                           verifyfiles=[self.getFilename(k) for k in ['images', 'classes']])

    def visualize(self):
        plots = [k for k in ['DoShowLL', 'DoShowPmax', 'DoShowSignalChange', 'DoShowMirror'] if self.ParamsDict[k]]
        if self.DoShowReferences:
            self.visualizeVar('DoShowReferences')
        if len(plots):
            self.launchPlots(plots)
         
    def visualizeReferences(self):
        refs = self.getFilename('iter_refs', iter=lastIteration(self))
        if xmippExists(refs):
            try:
                runShowJ(refs, extraParams="--mode metadata --render first")
            except Exception, e:
                from protlib_gui_ext import showError
                showError("Error launching java app", str(e))
               
    def visualizeVar(self, varName):
        if varName == 'DoShowReferences':
            self.visualizeReferences()
        else:
            self.launchPlots([varName])
            
    def launchPlots(self, selectedPlots):
        xplotter = launchML2DPlots(self, selectedPlots)
        xplotter.show()
        
def launchML2DPlots(protML, selectedPlots):
    ''' Launch some plot for an ML2D protocol run '''
    #import matplotlib
    import numpy as np
    from protlib_gui_figure import XmippPlotter

    protML._plot_count = 0
    lastIter = lastIteration(protML)
    if lastIter == 0:
        return
    refs = protML.getFilename('iter_refs', iter=lastIter)
#    if not exists(refs):
#        return 
#    blocks = getBlocksInMetaDataFile(refs)
#    lastBlock = blocks[-1]
    
    def doPlot(plotName):
        return plotName in selectedPlots
    
    # Remove 'mirror' from list if DoMirror is false
    if doPlot('DoShowMirror') and not protML.DoMirror:
        selectedPlots.remove('DoShowMirror')
        
    n = len(selectedPlots)
    if n == 0:
        showWarning("ML2D plots", "Nothing to plot")
        return
    elif n == 1:
        gridsize = [1, 1]
    elif n == 2:
        gridsize = [2, 1]
    else:
        gridsize = [2, 2]
        
    xplotter = XmippPlotter(*gridsize)
        
    # Create data to plot
    iters = range(1, lastIter+1)
    ll = []
    pmax = []
    for iter in iters:
#    for id in md:
#        iter = md.getValue(MDL_ITER, id)
#        if iter > 0:
        logs = protML.getFilename('iter_logs', iter=iter)
        md = MetaData(logs)
        id = md.firstObject()
        ll.append(md.getValue(MDL_LL, id))
        pmax.append(md.getValue(MDL_PMAX, id))
            
    if doPlot('DoShowLL'):
        a = xplotter.createSubPlot('Log-likelihood (should increase)', 'iterations', 'LL', yformat=True)
        a.plot(iters, ll)

    #Create plot of mirror for last iteration
    if doPlot('DoShowMirror'):
        from numpy import arange
        from matplotlib.ticker import FormatStrFormatter
        md = MetaData(refs)
        mirrors = [md.getValue(MDL_MIRRORFRAC, id) for id in md]
        nrefs = len(mirrors)
        ind = arange(1, nrefs + 1)
        width = 0.85
        a = xplotter.createSubPlot('Mirror fractions on last iteration', 'references', 'mirror fraction')
        a.set_xticks(ind + 0.45)
        a.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))
        a.bar(ind, mirrors, width, color='b')
        a.set_ylim([0, 1.])
        a.set_xlim([0.8, nrefs + 1])
        
    if doPlot('DoShowPmax'):
        a = xplotter.createSubPlot('Probabilities distribution', 'iterations', 'Pmax/Psum') 
        a.plot(iters, pmax, color='green')
    
    if doPlot('DoShowSignalChange'):
        md = MetaData()
        for iter in iters:
            fn = protML.getFilename('iter_refs', iter=iter)
            md2 = MetaData(fn)
            md2.fillConstant(MDL_ITER, str(iter))
            md.unionAll(md2)
        # 'iter(.*[1-9].*)@2D/ML2D/run_004/ml2d_iter_refs.xmd')
        #a = plt.subplot(gs[1, 1])
        #print "md:", md
        md2 = MetaData()    
        md2.aggregate(md, AGGR_MAX, MDL_ITER, MDL_SIGNALCHANGE, MDL_MAX)
        signal_change = [md2.getValue(MDL_MAX, id) for id in md2]
        xplotter.createSubPlot('Maximum signal change', 'iterations', 'signal change')
        xplotter.plot(iters, signal_change, color='green')
    
    return xplotter
    
             
def collectResults(log, WorkingDir, Prefix):
    for k in ['images', 'classes']:
        src = join(WorkingDir, Prefix  + '_result_%s.xmd' % k)
        dst = join(WorkingDir, 'result_%s.xmd' % k)
        renameFile(log, src, dst)
#        
#    imgs = join(WorkingDir, Prefix  + 'result_images.xmd')
#    refs = join(WorkingDir, Prefix  + 'result_classes.xmd')
#    
#    
#    c
#    for iter in range(LastIter):
#        fnIter = join(orootRoot, "iter%03d" % (iter+1))
#    
#    mdImgs = MetaData(oroot + '_final_images.xmd')
#    outImages = os.path.join(WorkingDir, 'result_images.xmd')
#    # change the sign of the angle, since is the expected one
#    # to align the images agains the reference
#    mdImgs.operate("anglePsi=360-anglePsi")
#    mdImgs.write('images@' + outImages)
#    mdRefs = MetaData(oroot + '_final_refs.xmd')
#    outRefs = os.path.join(WorkingDir, 'result_classes.xmd')
#    mdGroup = MetaData()
#    deleteFile(log, outRefs)
#    for idx in mdRefs:
#        ref = mdRefs.getValue(MDL_REF, idx)
#        mdRefs.setValue(MDL_CLASS_COUNT, long(round(mdRefs.getValue(MDL_WEIGHT, idx))), idx)
#        mdGroup.importObjects( mdImgs, MDValueEQ(MDL_REF, ref))
#        mdGroup.write('class%(ref)06d_images@%(outRefs)s' % locals(),MD_APPEND)
#    mdRefs.write('classes@'+outRefs,MD_APPEND)
    
