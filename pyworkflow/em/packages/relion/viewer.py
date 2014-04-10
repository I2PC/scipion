# **************************************************************************
# *
# * Authors:     Jose Gutierrez (jose.gutierrez@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
from pyworkflow.utils.path import cleanPath
"""
This module implement the wrappers around xmipp_showj
visualization program.
"""
import os

from pyworkflow.viewer import Viewer, ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from protocol_classify2d import ProtRelionClassify2D
from protocol_classify3d import ProtRelionClassify3D
from protocol_refine3d import ProtRelionRefine3D
from pyworkflow.protocol.params import *
from convert import addRelionLabels, restoreXmippLabels
import xmipp
from pyworkflow.em.packages.xmipp3.plotter import XmippPlotter

ITER_LAST = 0
ITER_SELECTION = 1

ANGDIST_2DPLOT = 0
ANGDIST_CHIMERA = 1

VOLUME_SLICES = 0
VOLUME_CHIMERA = 1

CLASSES_ALL = 0
CLASSES_SEL = 1

    
class RelionViewer(ProtocolViewer):
    """ Class to visualize Relion protocols """
    _targets = [ProtRelionClassify2D, ProtRelionClassify3D, ProtRelionRefine3D]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    _label = 'viewer relion'
    
    # The following "tricks" with changing _defineParams
    # is to postpone the definition of parameters after
    # the protocol is set, to allow a more dynamic definition
    # if this use become common, we need to move it to base class.
    def setProtocol(self, protocol):
        ProtocolViewer.setProtocol(self, protocol)
        self.__defineParams(self._form)
        self._createVarsFromDefinition()
        
    def _defineParams(self, form):
        self._form = form
        
    def __defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('viewIter', EnumParam, choices=['last', 'selection'], default=ITER_LAST, 
                      display=EnumParam.DISPLAY_LIST,
                      label="Iteration to visualize", 
                      help="""
*last*: only the last iteration will be visualized.
*selection*: you may specify a range of iterations.
Examples:
"1,5-8,10" -> [1,5,6,7,8,10]
"2,6,9-11" -> [2,6,9,10,11]
"2 5, 6-8" -> [2,5,6,7,8]                      
                           """)
        form.addParam('iterSelection', NumericRangeParam, 
                      condition='viewIter==%d' % ITER_SELECTION, 
                      label="Iterations list", 
                      help="Write the iteration list to visualize.")

        changesLabel = 'Changes in Offset and Angles'
        if self.protocol.IS_CLASSIFY:
            form.addParam('showImagesInClasses', BooleanParam, default=True,
                          label='Images assigned to each Class',
                          help='Display the classes and the images associated.')
            changesLabel = 'Changes in Offset, Angles and Classes'
        
        form.addParam('showLL', BooleanParam, label="Show maximum model probability?", default=False, 
                      help='Max likelihood per image may be used to delete images with smaller value.'
                           'The higher, the better. Consider remove particles with low values.')      
        form.addParam('showPMax', BooleanParam, default=True, 
                      label="Show average PMax?", 
                      help='Average (per class) of the maximum value\n of normalized probability function')      
        form.addParam('showChanges', BooleanParam, default=True,
                      label=changesLabel,
                      help='Visualize changes in orientation, offset and\n number images assigned to each class')

         
        if self.protocol.IS_3D:
            group = form.addGroup('3D analysis')
            
            if self.protocol.IS_CLASSIFY:
                group.addParam('showClasses3D', EnumParam, choices=['all', 'selection'], default=CLASSES_ALL,
                              display=EnumParam.DISPLAY_LIST,
                              label='CLASS 3D to visualize',
                              help='')
                group.addParam('class3DSelection', NumericRangeParam, default='1',
                              condition='showClasses3D == %d' % CLASSES_SEL,
                              label='Classes list',
                              help='')
            else:
                group.addParam('showHalves', EnumParam, choices=['half1', 'half2', 'both'], default=0,
                              label='Half to visualize',
                              help='Select which half do you want to visualize.')
            
            group.addParam('displayVol', EnumParam, choices=['slices', 'chimera'], 
                          display=EnumParam.DISPLAY_LIST, default=VOLUME_SLICES,
                          label='Display volume with',
                          help='*slices*: display volumes as 2D slices along z axis.\n'
                               '*chimera*: display volumes as surface with Chimera.')
            group.addParam('displayAngDist', EnumParam, choices=['2D plot', 'chimera'], 
                          display=EnumParam.DISPLAY_LIST, default=ANGDIST_2DPLOT,
                          label='Display angular distribution',
                          help='*2D plot*: display angular distribution as interative 2D in matplotlib.\n'
                               '*chimera*: display angular distribution using Chimera with red spheres.') 
            group.addParam('spheresScale', IntParam, default=-1, 
                          expertLevel=LEVEL_ADVANCED,
                          label='',
                          help='')
            group.addParam('resolutionPlotsSSNR', BooleanParam, default=True,
                          label='Display SSNR plots?',
                          help='Display signal to noise ratio plots (SSNR) ')
            group.addParam('resolutionPlotsFSC', BooleanParam, default=True,
                          label='Display resolution plots (FSC) ?',
                          help='')
            group.addParam('resolutionThresholdFSC', FloatParam, default=0.5, 
                          expertLevel=LEVEL_ADVANCED,
                          label='Threshold in resolution plots',
                          help='')                                      
        
    def _getVisualizeDict(self):
        return {'showImagesInClasses': self._showImagesInClasses,
                'showLL': self._showLL,
                'showPMax': self._showPMax,
                'showChanges': self._showChanges,
                'displayVol': self._showVolumes,
                'displayAngDist': self._showAngularDistribution,
                'resolutionPlotsSSNR': self._showSSNR,
                'resolutionPlotsFSC': self._showFSC
                }
        
    def _viewAll(self, *args):
        pass
        
    def display2D(self, filename, extraParams=''):
        from pyworkflow.em.packages.xmipp3.viewer import runShowJ
        runShowJ(filename, extraParams=extraParams)
        
    def displayScipion(self, filename, extraParams=''):
        inputParticlesId = self.protocol.inputParticles.get().strId()
        from pyworkflow.em.packages.xmipp3.viewer import runScipionShowJ        
        runScipionShowJ(filename, "Particles", self._project.getName(), self.protocol.strId(), inputParticlesId)

    def _load(self):
        """ Load selected iterations and classes 3D for visualization mode. """
        self._refsList = [1] 
        if self.protocol.IS_3D and self.protocol.IS_CLASSIFY:
            if self.showClasses3D == CLASSES_ALL:
                self._refsList = range(1, self.protocol.numberOfClasses.get()+1)
            else:
                self._refsList = self._getListFromRangeString(self.class3DSelection.get())
        self.protocol._initialize() # Load filename templates
        self.firstIter = self.protocol._firstIter()
        self.lastIter = self.protocol._lastIter()
        
        if self.viewIter.get() == ITER_LAST:
            self._iterations = [self.lastIter]
        else:
            self._iterations = self._getListFromRangeString(self.iterSelection.get())
            
        from matplotlib.ticker import FuncFormatter
        self._plotFormatter = FuncFormatter(self._formatFreq) 
        
    def _formatFreq(self, value, pos):
        """ Format function for Matplotlib formatter. """
        inv = 999
        if value:
            inv = int(1/value)
        return "1/%d" % inv

    def _getGridSize(self, n=None):
        """ Figure out the layout of the plots given the number of references. """
        if n is None:
            n = len(self._refsList)
        
        if n == 1:
            gridsize = [1, 1]
        elif n == 2:
            gridsize = [2, 1]
        else:
            gridsize = [(n+1)/2, 2]
            
        return gridsize
    
    def _getPrefixes(self):
        prefixes = self.protocol.PREFIXES
        halves = getattr(self, 'showHalves', None)
        if halves:
            if halves == 0:
                prefixes = ['half1_']
            elif halves == 1:
                prefixes = ['half2_']
        return prefixes
                

#===============================================================================
# showImagesInClasses     
#===============================================================================
    
    def _createImagesInClasses(self, paramName=None):    
        """ Read Relion _data.star images file and 
        generate a new metadata with the Xmipp classification standard:
        a 'classes' block and a 'class00000?_images' block per class.
        If the new metadata was already written, it is just shown.
        """
        self._load()
        
        data_classes = []
        
        for it in self._iterations:
            data_classes.append(self.protocol._getIterClasses(it))
        
        return data_classes
            
    def _showImagesInClasses(self, paramName=None):
            data_classes = self._createImagesInClasses()
            for data in data_classes:
                self.displayScipion(data, extraParams='--mode metadata --render first')
          
#=====================================================================
# showLLRelion
#=====================================================================
          
    def _createLL(self):
        self._load()
        plotters = []
        files = []
        for it in self._iterations:
            # FILES
            fn = self._getFilesPerIterationLL(it)
            files.append(fn)
            # PLOTTERS
            xplotter = self._createPlotPerIterationLL(it, fn)
            plotters.append(xplotter)
            
        return plotters, files
    
    def _getFilesPerIterationLL(self, it):
        fn = self.protocol._getIterSortedData(it)
        addRelionLabels()
        return fn
    
    def _createPlotPerIterationLL(self, it, fn):
        md = xmipp.MetaData(fn)
        xplotter = XmippPlotter(windowTitle="max Likelihood particles sorting Iter_%d" % it)
        xplotter.createSubPlot("Particle sorting: Iter_%d" % it, "Particle number", "maxLL")
        xplotter.plotMd(md, False, mdLabelY=xmipp.MDL_LL)
        
        return xplotter
     
    def _showLL(self, paramName=None):
        plotters, files = self._createLL()
        for xplotter, fn in zip(plotters, files):
            self.display2D(fn)
            xplotter.show()
            
#===============================================================================
# ShowPMax
#===============================================================================
        
    def _createPMax(self):
        self._load()
        
        # FILE
        fn = self._createFilePMax()
        # PLOT
        xplotter = self._createPlotPMax(fn)
        
        return xplotter, fn
    
    def _createFilePMax(self):
        labels = [xmipp.MDL_AVGPMAX, xmipp.MDL_PMAX]
        addRelionLabels(extended=True)  
        
        mdIters = xmipp.MetaData()
        iterations = range(self.firstIter, self.lastIter+1)
        
        for it in iterations: # range (firstIter,self._visualizeLastIteration+1): #alwaya list all iteration
            objId = mdIters.addObject()
            mdIters.setValue(xmipp.MDL_ITER, it, objId)
            for i, prefix in enumerate(self.protocol.PREFIXES):
                fn = 'model_general@'+ self.protocol._getFileName(prefix + 'model', iter=it)
                md = xmipp.MetaData(fn)
                pmax = md.getValue(xmipp.MDL_AVGPMAX, md.firstObject())
                mdIters.setValue(labels[i], pmax, objId)
        fn = self.protocol._getFileName('all_avgPmax_xmipp')
        mdIters.write(fn)
            
        return fn
    
    def _createPlotPMax(self, fn):
        labels = [xmipp.MDL_AVGPMAX, xmipp.MDL_PMAX]
        colors = ['g', 'b']

        md = xmipp.MetaData(fn)
        
        xplotter = XmippPlotter()
        xplotter.createSubPlot("Avg PMax per Iterations", "Iterations", "Avg PMax")
        
        for label, color in zip(labels, colors):
            xplotter.plotMd(md, xmipp.MDL_ITER, label, color)
        
        if len(self.protocol.PREFIXES) > 1:
            xplotter.showLegend(self.protocol.PREFIXES)
        
        return xplotter
        
    def _showPMax(self, paramName=None):
        xplotter, fn = self._createPMax()
        self.display2D(fn)                    
        xplotter.show()
        
    
#===============================================================================
# ShowChanges    
#===============================================================================    

    def _createChanges(self):
        self._load()
        addRelionLabels(extended=True)  
        
        mdIters = xmipp.MetaData()
        iterations = range(self.firstIter, self.lastIter+1)
        
        print " Computing average changes in offset, angles, and class membership"
        for it in iterations:
            print "Computing data for iteration; %03d" % it
            objId = mdIters.addObject()
            mdIters.setValue(xmipp.MDL_ITER, it, objId)
            #agregar por ref3D
            fn = self.protocol._getFileName('optimiser', iter=it )
            md = xmipp.MetaData(fn)
            firstId = md.firstObject()
            for label in self.protocol.CHANGE_LABELS:
                mdIters.setValue(label, md.getValue(label, firstId), objId)
        fn = self.protocol._getFileName('all_changes_xmipp')
        mdIters.write(fn)
        return fn

    def _showChanges(self, paramName=None):
        fn = self._createChanges()
        self.display2D(fn)
        
#===============================================================================
# ShowVolumes
#===============================================================================
        
    def _createVolumesMd(self):
        """ Write a metadata with all volumes selected for visualization. """
        import xmipp
        self._load()
        prefixes = self._getPrefixes()
        
        mdPath = self.protocol._getExtraPath('relion_viewer_volumes.xmd')
        cleanPath(mdPath)
        md = xmipp.MetaData()
        
        for it in self._iterations:
            md.clear()
            for ref3d in self._refsList:
                for prefix in prefixes:
                    volFn = self.protocol._getFileName(prefix + 'volume', iter=it, ref3d=ref3d)
                    md.setValue(xmipp.MDL_IMAGE, volFn, md.addObject())
            md.write('iter%03d@%s' % (it, mdPath), xmipp.MD_APPEND)
        return mdPath
    
    def _showVolumesChimera(self):
        """ Create a chimera script to visualize selected volumes. """
        self._load()
        prefixes = self._getPrefixes()
        volumes = []
        
        for it in self._iterations:
            for ref3d in self._refsList:
                for prefix in prefixes:
                    volFn = self.protocol._getFileName(prefix + 'volume', iter=it, ref3d=ref3d)
                    volumes.append(volFn)
                    
        if len(volumes) > 1:
            cmdFile = self.protocol._getExtraPath('chimera_volumes.cmd')
            f = open(cmdFile, 'w+')
            for volFn in volumes:
                # We assume that the chimera script will be generated
                # at the same folder than relion volumes
                localVol = os.path.basename(volFn.replace(':mrc', ''))
                f.write("open %s\n" % localVol)
            f.write('tile\n')
            f.close()
            os.system('chimera %s &' % cmdFile)                    
        else:
            os.system('xmipp_chimera_client --input "%s" --mode projector 256 &' % volumes[0])
            
    def _showVolumes(self, paramName=None):
        if self.displayVol == VOLUME_CHIMERA:
            self._showVolumesChimera()
        
        elif self.displayVol == VOLUME_SLICES:
            mdPath = self._createVolumesMd()
            self.display2D(mdPath)
            
#===============================================================================
# showAngularDistribution
#===============================================================================
                            
    def _createAngularDistribution(self):
        self._load()
        
        arguments = []
        plotters = []
        
        if self.displayAngDist == ANGDIST_CHIMERA:
            for it in self._iterations:
                arguments = self._createAngDistChimera(it)
                        
        elif self.displayAngDist == ANGDIST_2DPLOT:
            for it in self._iterations:
                xplotter = self._createAngDist2D(it)
                plotters.append(xplotter)
                
        return plotters, arguments
    
    def _createAngDistChimera(self, it):
        arguments = []
        # FIXME
        #outerRadius = int(float(self.MaskRadiusA)/self.SamplingRate)
        outerRadius = 30
        radius = float(outerRadius) * 1.1
        # Common variables to use
        sphere = self.spheresScale.get()
        prefixes = self._getPrefixes()

        data_angularDist = self.protocol._getIterAngularDist(it)
        for ref3d in self._refsList:
            for prefix in prefixes:
                volFn = self.protocol._getFileName(prefix + 'volume', iter=it, ref3d=ref3d)
                args = "--input '%s' --mode projector 256 -a %sclass%06d_angularDist@%s red %f " % (volFn, prefix, ref3d, data_angularDist, radius)
                if sphere > 0:
                    args += ' %f ' % sphere
                arguments.append(args)
                    
        return arguments
        
    
    def _createAngDist2D(self, it):
        # Common variables to use
        prefixes = self._getPrefixes()
        nrefs = len(self._refsList)
        n = nrefs * len(prefixes)
        gridsize = self._getGridSize(n)
        
        data_angularDist = self.protocol._getIterAngularDist(it)
        xplotter = XmippPlotter(*gridsize, mainTitle='Iteration %d' % it, windowTitle="Angular Distribution")
        for ref3d in self._refsList:
            for prefix in prefixes:
                md = xmipp.MetaData("class%06d_angularDist@%s" % (ref3d, data_angularDist))
                plot_title = '%sclass %d' % (prefix, ref3d)
                xplotter.plotMdAngularDistribution(plot_title, md)
        
        return xplotter
    
    def _showAngularDistribution(self, paramName=None):
        plotters, arguments = self._createAngularDistribution()
        
        if arguments:
            for args in arguments:
                os.system('xmipp_chimera_client %s &' % args)
        
        if plotters:
            for xplotter in plotters:
                xplotter.show()
                
#===============================================================================
# plotSSNR              
#===============================================================================
               
    def _plotSSNR(self, a, fn):
        mdOut = xmipp.MetaData(fn)
        md = xmipp.MetaData()
        # only cross by 1 is important
        md.importObjects(mdOut, xmipp.MDValueGT(xmipp.MDL_RESOLUTION_SSNR, 0.9))
        md.operate("resolutionSSNR=log(resolutionSSNR)")
        resolution_inv = [md.getValue(xmipp.MDL_RESOLUTION_FREQ, id) for id in md]
        frc = [md.getValue(xmipp.MDL_RESOLUTION_SSNR, id) for id in md]
        a.plot(resolution_inv, frc)
        a.xaxis.set_major_formatter(self._plotFormatter)               
 
    def _createSSNR(self, paramName=None):
        self._load()
        prefixes = self._getPrefixes()        
        nrefs = len(self._refsList)
        n = nrefs * len(prefixes)
        gridsize = self._getGridSize(n)
        addRelionLabels()
        xmipp.activateMathExtensions()
        xplotter = XmippPlotter(*gridsize)
        
        for prefix in prefixes:
            for ref3d in self._refsList:
                plot_title = 'Resolution SSNR %s, for Class %s' % (prefix, ref3d)
                a = xplotter.createSubPlot(plot_title, 'Armstrongs^-1', 'log(SSNR)', yformat=False)
                blockName = 'model_class_%d@' % ref3d
                legendName = []
                for it in self._iterations:
                    fn = self.protocol._getFileName(prefix + 'model', iter=it)
                    if os.path.exists(fn):
                        self._plotSSNR(a, blockName+fn)
                    legendName.append('iter %d' % it)
                xplotter.showLegend(legendName)
                a.grid(True)
        
        return xplotter
        
    def _showSSNR(self, paramName=None):
        xplotter = self._createSSNR()
        xplotter.show()
            
            
#===============================================================================
# plotFSC            
#===============================================================================

    def _plotFSC(self, a, model_star):
        md = xmipp.MetaData(model_star)
        resolution_inv = [md.getValue(xmipp.MDL_RESOLUTION_FREQ, id) for id in md]
        frc = [md.getValue(xmipp.MDL_RESOLUTION_FRC, id) for id in md]
        self.maxFrc = max(frc)
        self.minInv = min(resolution_inv)
        self.maxInv = max(resolution_inv)
        a.plot(resolution_inv, frc)
        a.xaxis.set_major_formatter(self._plotFormatter)
        a.set_ylim([-0.1, 1.1])
            
    def _createFSC(self):
        self._load()
        threshold = self.resolutionThresholdFSC.get()
        prefixes = self._getPrefixes()        
        nrefs = len(self._refsList)
        n = nrefs * len(prefixes)
        gridsize = self._getGridSize(n)
        
        xmipp.activateMathExtensions()
        addRelionLabels()
        
        xplotter = XmippPlotter(*gridsize, windowTitle='Resolution FSC')

        for prefix in prefixes:
            for ref3d in self._refsList:
                plot_title = prefix + 'class %s' % ref3d
                a = xplotter.createSubPlot(plot_title, 'Armstrongs^-1', 'FSC', yformat=False)
                legends = []
                blockName = 'model_class_%d@' % ref3d
                for it in self._iterations:
                    model_star = self.protocol._getFileName(prefix + 'model', iter=it)
                    if os.path.exists(model_star):
                        self._plotFSC(a, blockName + model_star)
                        legends.append('iter %d' % it)
                xplotter.showLegend(legends)
                if threshold < self.maxFrc:
                    a.plot([self.minInv, self.maxInv],[threshold, threshold], color='black', linestyle='--')
                a.grid(True)
            
        return xplotter
        
    def _showFSC(self, paramName=None):
        plots = [self._createFSC()]
        for xplotter in plots:
            xplotter.show()
            
### WEB UTILS ################################################# 

    def getVisualizeDictWeb(self):
        return {'showImagesInClasses': 'doShowImagesInClasses',
                'showLL': 'doShowLLRelion',
                'showPMax': 'doShowPMax',
                'showChanges': 'doShowChanges',
                'displayVol': 'doShowVolumes',
                'displayAngDist': 'doShowAngularDistributionRelion',
                'resolutionPlotsSSNR': 'doPlotsSSNR',
                'resolutionPlotsFSC': 'doPlotsFSC'
                }

    @classmethod
    def getView(cls):
        """ This function will notify the web viewer for this protocol"""
        return "viewerForm"
    
    @classmethod
    def getViewFunction(cls):
        """ This will return the name of the function to view
        in web one (or all) params of the protocol"""
        return "viewerRelion"

        
