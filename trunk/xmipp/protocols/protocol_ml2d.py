#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using maximum-likelihood principles
#
#   Author:  Sjors Scheres, January 2008
#  Updated:  J. M. de la Rosa Trevin July 2011
#

from xmipp import MetaData, MDL_ITER, MDL_LL, MDL_REF, MDValueEQ, getBlocksInMetaDataFile, \
MDL_PMAX, MDL_SIGNALCHANGE, AGGR_MAX, MDL_MAX, MDL_MIRRORFRAC
from protlib_base import XmippProtocol, protocolMain
from config_protocols import protDict
import os
from os.path import exists
from protlib_utils import runImageJPlugin
from protlib_gui_ext import showWarning

class ProtML2D(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.ml2d.name, scriptname, project)
        self.Import = 'from protocol_ml2d import *'
        
    def createFilenameTemplates(self):
        return {
                'iter_logs': self.workingDirPath('ml2d_iter_logs.xmd'),
                'iter_refs': self.workingDirPath('ml2d_iter_refs.xmd')
                }
        #self.fnIterLogs = self.workingDirPath('ml2d_iter_logs.xmd')
        #self.fnIterRefs = self.workingDirPath('ml2d_iter_refs.xmd')
        
    def validate(self):
        errors = []
        if self.DoCorrectAmplitudes and not exists(self.InCtfDatFile):
            errors.append("Missing '%s' file for correcting amplitudes" % self.InCtfDatFile)
        return errors
    
    def summary(self):
        md = MetaData(self.ImgMd)
        lines = [('Input images            ', "%s (%u)" % (self.ImgMd, md.size())),
                 ('Reference image', self.RefMd)]
        
        logs = self.getFilename('iter_logs')    
        if exists(logs):
            md = MetaData(logs)
            id = md.lastObject()
            iteration = md.getValue(MDL_ITER, id)
            lines.append(('Iteration                   ', str(iteration)))
            LL = md.getValue(MDL_LL, id)
            lines.append(('LogLikelihood          ', str(LL)))
        
        output = ["%s : %s" % (k.ljust(20),  v) for k, v in lines]
        return output
    
    def defineSteps(self):
        progId = "ml"
        if (self.DoMlf):
            progId += "f"  
        
        program = "xmipp_%s_align2d" % progId

        restart = False
        if (restart):
            pass 
            #Not yet implemented
            #params= ' --restart ' + utils_xmipp.composeFileName('ml2d_it', RestartIter,'log')
        else: 
            # Dictionary with boolean options and the cmd options
            booleanDict = {'DoMirror': '--mirror', 'DoNorm': '--norm', 'ZeroOffsets': '--zero_offsets',
                           'FixSigmaNoise': '--fix_sigma_noise', 'FixSigmaOffset': '--fix_sigma_offset',
                           'FixFractions': '--fix_fractions'}
            
            prefix = '%s2d' % progId
            oroot = self.workingDirPath(prefix)
            params = ' -i %s --oroot %s' % (self.ImgMd, oroot)
            # Number of references will be ignored if -ref is passed as expert option
            if self.DoGenerateReferences:
                params += ' --nref %d' % self.NumberOfReferences
            
            if (self.DoFast and not self.DoMlf):
                params += ' --fast'
            if (self.NumberOfThreads > 1  and not self.DoMlf):
                params += ' --thr %i' % self.NumberOfThreads
            if (self.DoMlf):
                if (self.DoCorrectAmplitudes):
                    params += ' --ctfdat %s' % self.InCtfDatFile
                else:
                    params += ' --no_ctf --pixel_size %f' % self.PixelSize
                if (not self.ImagesArePhaseFlipped):
                    params += ' --not_phase_flipped'
                if (self.HighResLimit > 0):
                    params += ' --high %f' % self.HighResLimit
            if self.MaxIters != 100:
                params += " --iter %d" % self.MaxIters
            #Add all boolean options if true
            for k, v in booleanDict.iteritems():
                if self.__dict__[k]:
                    params += " " + v
            
            #Add extra options
            #params += ' ' + self.ExtraParams
        self.insertRunJobStep(program, params)

        self.Db.insertStep('collectResults', WorkingDir=self.WorkingDir, Prefix=prefix,
                           verifyfiles=[self.workingDirPath(f) for f in ['result_images.xmd', 'result_classes.xmd']])

    def visualize(self):
        plots = [k for k in ['DoShowLL', 'DoShowPmax', 'DoShowSignalChange', 'DoShowMirror'] if self.ParamsDict[k]]
        if len(plots):
            self.launchML2DPlots(plots)
        if self.DoShowReferences:
            self.visualizeVar('DoShowReferences')
         
    def visualizeReferences(self):
        refs = self.getFilename('iter_refs')
        if exists(refs):
            blocks = getBlocksInMetaDataFile(refs)
            lastBlock = blocks[-1]
            try:
                runImageJPlugin("512m", "XmippBrowser.txt", "-i %(lastBlock)s@%(refs)s --mode metadata --render" 
                            % locals(), batchMode=True)
            except Exception, e:
                from protlib_gui_ext import showError
                showError("Error launching java app", str(e))
               
    def visualizeVar(self, varName):
        if varName == 'DoShowReferences':
            self.visualizeReferences()
        else:
            self.launchML2DPlots([varName])
        
    def launchML2DPlots(self, selectedPlots):
        ''' Launch some plot for an ML2D protocol run '''
        #import matplotlib
        import numpy as np
        #matplotlib.use('TkAgg') # do this before importing pylab
        import matplotlib.ticker as ticker
        import matplotlib.gridspec as gridspec
        import matplotlib.pyplot as plt
    
        self._plot_count = 0
        refs = self.getFilename('iter_refs')
        if not exists(refs):
            return 
        blocks = getBlocksInMetaDataFile(refs)
        lastBlock = blocks[-1]
        
        def doPlot(plotName):
            return plotName in selectedPlots
        
        # Remove 'mirror' from list if DoMirror is false
        if doPlot('DoShowMirror') and not self.DoMirror:
            selectedPlots.remove('DoShowMirror')
            
        n = len(selectedPlots)
        if n == 0:
            showWarning("ML2D plots", "Nothing to plot")
            return
        elif n == 1:
            dims = (1, 1, 6, 5)
        elif n == 2:
            dims = (2, 1, 4, 6)
        else:
            dims = (2, 2, 8, 6)
            
        xg, yg, xs, ys = dims
    
        #f = plt.figure()
        def createSubPlot(title, xlabel, ylabel, yformat=None):
            self._plot_count += 1
            a = fig.add_subplot(xg, yg, self._plot_count)
            #a.get_label().set_fontsize(12)
            a.set_title(title)
            a.set_xlabel(xlabel)
            a.set_ylabel(ylabel)
            
            if yformat:
                formatter = ticker.FormatStrFormatter('%1.2e')
                a.yaxis.set_major_formatter(formatter)
            a.xaxis.get_label().set_fontsize(10)
            a.yaxis.get_label().set_fontsize(10)
            labels = a.xaxis.get_ticklabels() + a.yaxis.get_ticklabels()
            for label in labels:
                label.set_fontsize(8) # Set fontsize
                label.set_text('aa')
                #print label.
                #label.set_visible(False)
            return a
        
        # Create figure and prepare grid    
        gs = gridspec.GridSpec(xg, yg)#, height_ratios=[7,4])
        gs.update(left=0.15, right=0.95, hspace=0.25, wspace=0.4)#, top=0.8, bottom=0.2)    
        fig = plt.figure(figsize=(xs, ys), dpi=100)
        
        # Create data to plot
        logs = self.getFilename('iter_logs')
        md = MetaData(logs)
        iters = []
        ll = []
        pmax = []
        for id in md:
            iter = md.getValue(MDL_ITER, id)
            if iter > 0:
                iters.append(iter)
                ll.append(md.getValue(MDL_LL, id))
                pmax.append(md.getValue(MDL_PMAX, id))
                
        if doPlot('DoShowLL'):
            a = createSubPlot('Log-likelihood (should increase)', 'iterations', 'LL', yformat='%1.2e')
            a.plot(iters, ll)
    
        #Create plot of mirror for last iteration
        if doPlot('DoShowMirror'):
            md = MetaData('%(lastBlock)s@%(refs)s' % locals())
            mirrors = [md.getValue(MDL_MIRRORFRAC, id) for id in md]
            nrefs = len(mirrors)
            ind = np.arange(nrefs)
            width = 0.85#min(0.15, 5/float(nrefs))       # the width of the bars: can also be len(x) sequence
            a = createSubPlot('Mirror fractions on last iteration', 'references', 'mirror fraction')
            a.bar(ind, mirrors, width, color='b')
            a.set_ylim([0, 1.])
            #a.set_xticks(ind + width / 2., [str(r+1) for r in ind])
            #a.set_yticks(np.arange(0, 1, 10))
            
        if doPlot('DoShowPmax'):
            a = createSubPlot('Probabilities distribution', 'iterations', 'Pmax/Psum') 
            a.plot(iters, pmax, color='green')
        
        if doPlot('DoShowSignalChange'):
            # Read all iteration block into one metadata (iter 0 is avoided)
            md = MetaData('iter(.*[1-9].*)@%s' % refs)
            # 'iter(.*[1-9].*)@2D/ML2D/run_004/ml2d_iter_refs.xmd')
            #a = plt.subplot(gs[1, 1])
            md2 = MetaData()    
            md2.aggregate(md, AGGR_MAX, MDL_ITER, MDL_SIGNALCHANGE, MDL_MAX)
            signal_change = [md2.getValue(MDL_MAX, id) for id in md2]
            a = createSubPlot('Maximum signal change', 'iterations', 'signal change')
            a.plot(iters, signal_change, color='green')    
        
        plt.tight_layout()
        plt.show()   
             
def collectResults(log, WorkingDir, Prefix):
    oroot = os.path.join(WorkingDir, Prefix)
    mdImgs = MetaData(oroot + '_final_images.xmd')
    outImages = os.path.join(WorkingDir, 'result_images.xmd')
    mdImgs.write('images@' + outImages)
    mdRefs = MetaData(oroot + '_final_refs.xmd')
    outRefs = os.path.join(WorkingDir, 'result_classes.xmd')
    mdRefs.write('classes@' + outRefs)
    mdGroup = MetaData()
    for idx in mdRefs:
        ref = mdRefs.getValue(MDL_REF, idx)
        mdGroup.importObjects( mdImgs, MDValueEQ(MDL_REF, ref))
        mdGroup.writeBlock(outRefs, 'class%06d_images' % ref)


    
