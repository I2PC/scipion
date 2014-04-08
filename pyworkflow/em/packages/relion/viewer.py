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


ITER_LAST = 0
ITER_SELECTION = 1

    
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
                group.addParam('showClasses3D', EnumParam, choices=['selection', 'all'], default=0,
                              display=EnumParam.DISPLAY_LIST,
                              label='CLASS 3D to visualize',
                              help='')
                group.addParam('class3DSelection', NumericRangeParam, default='1',
                              condition='showClasses3D == 0',
                              label='Classes list',
                              help='')
            else:
                group.addParam('showHalves', EnumParam, choices=['half1', 'half2', 'both'], default=0,
                              label='Half to visualize',
                              help='Select which half do you want to visualize.')
            
            group.addParam('displayVol', EnumParam, choices=['slices', 'chimera'], 
                          display=EnumParam.DISPLAY_LIST, default=0,
                          label='Display volume with',
                          help='*slices*: display volumes as 2D slices along z axis.\n'
                               '*chimera*: display volumes as surface with Chimera.')
            group.addParam('displayAngDist', EnumParam, choices=['2D plot', 'chimera'], 
                          display=EnumParam.DISPLAY_LIST, default=0,
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
                'displayAngDist': self._showAngularDistribution
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
#        DisplayRef3DNo = self.parser.getTkValue('DisplayRef3DNo')
#        VisualizeIter = self.parser.getTkValue('VisualizeIter')
#        
#        if DisplayRef3DNo == 'all':
#            self._visualizeRef3Ds = range(1, self.NumberOfClasses + 1)
#        else:
#            self._visualizeRef3Ds = getListFromRangeString(self.parser.getTkValue('SelectedRef3DNo'))
#        
#        self._visualizeNrefs = len(self._visualizeRef3Ds)

        self._refsList = [1] #FIXME: Change this for classify 2d and 3d
        self.protocol._initialize() # Load filename templates
        self.firstIter = self.protocol._firstIter()
        self.lastIter = self.protocol._lastIter()
        
        if self.viewIter == ITER_LAST:
            self._iterations = [self.lastIter]
        else:
            self._iterations = self._getListFromRangeString(self.iterSelection.get())

#        self._visualizeLastIteration = self._iterations[-1]
#        self._visualizeLastRef3D = self._visualizeRef3Ds[-1]
#        
#        self._visualizeVolumesMode = self.parser.getTkValue('DisplayReconstruction')
#        
#        from matplotlib.ticker import FuncFormatter
#        self._plotFormatter = FuncFormatter(self._formatFreq) 

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
                
    def _showImagesInClasses(self, paramName):
        """ Read Relion _data.star images file and 
        generate a new metadata with the Xmipp classification standard:
        a 'classes' block and a 'class00000?_images' block per class.
        If the new metadata was already written, it is just shown.
        """
        self._load()
        
        for it in self._iterations:
            data_classes = self.protocol._getIterClasses(it)
            self.displayScipion(data_classes, extraParams='--mode metadata --render first')
          
    def _createLL(self, paramName):
        self._load()
        import xmipp
        from pyworkflow.em.packages.xmipp3.plotter import XmippPlotter
        from convert import addRelionLabels
        plotters = []
        files = []
        for it in self._iterations:
            fn = self.protocol._getIterSortedData(it)
            addRelionLabels()
            md = xmipp.MetaData(fn)
            files.append(fn)
            xplotter = XmippPlotter(windowTitle="max Likelihood particles sorting Iter_%d" % it)
            xplotter.createSubPlot("Particle sorting: Iter_%d" % it, "Particle number", "maxLL")
            xplotter.plotMd(md, False, mdLabelY=xmipp.MDL_LL)
            plotters.append(xplotter)
            
        return plotters, files
     
    def _showLL(self, paramName):
        plotters, files = self._createLL(paramName)
        for xplotter, fn in zip(plotters, files):
            self.display2D(fn)
            xplotter.show()
        
    def _createPMax(self, paramName):
        self._load()
        import xmipp
        from pyworkflow.em.packages.xmipp3.plotter import XmippPlotter
        from convert import addRelionLabels
        addRelionLabels(extended=True)  
                  
        mdIters = xmipp.MetaData()
        iterations = range(self.firstIter, self.lastIter+1)
        labels = [xmipp.MDL_AVGPMAX, xmipp.MDL_PMAX]
        colors = ['g', 'b']
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
        #self.display2D(fn, extraParams)
        xplotter = XmippPlotter()
        xplotter.createSubPlot("Avg PMax per Iterations", "Iterations", "Avg PMax")
        for label, color in zip(labels, colors):
            xplotter.plotMd(mdIters, xmipp.MDL_ITER, label, color)
        
        if len(self.protocol.PREFIXES) > 1:
            xplotter.showLegend(self.protocol.PREFIXES)
        
        return xplotter, fn
        
    def _showPMax(self, paramName):
        xplotter, fn = self._createPMax(paramName)
        self.display2D(fn)                    
        xplotter.show()
        
    def _createChanges(self, paramName):
        self._load()
        import xmipp
        from pyworkflow.em.packages.xmipp3.plotter import XmippPlotter
        from convert import addRelionLabels
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

    def _showChanges(self, paramName):
        fn = self._createChanges(paramName)
        self.display2D(fn)
        
    def _createVolumes(self, paramName):
        files = []
        self._load()
        prefixes = self._getPrefixes()
        for it in self._iterations:
            for ref3d in self._refsList:
                for prefix in prefixes:
                    volFn = self.protocol._getFileName(prefix + 'volume', iter=it, ref3d=ref3d)
                    files.append(volFn)
        return files
        
    def _showVolumes(self, paramName=None):
        files = self._createVolumes(paramName)
        for volFn in files:
            os.system('xmipp_chimera_client --input "%s" --mode projector 256 &' % volFn)
                            
    def _createAngularDistribution(self, paramName):
        self._load()
        import xmipp
        from pyworkflow.em.packages.xmipp3.plotter import XmippPlotter
        sphere = self.spheresScale.get()
        prefixes = self._getPrefixes()
        arguments = []
        plotters = []
        
        if self.displayAngDist == 1:
            # FIXME
            #outerRadius = int(float(self.MaskRadiusA)/self.SamplingRate)
            outerRadius = 30
            radius = float(outerRadius) * 1.1

            for it in self._iterations:
                data_angularDist = self.protocol._getIterAngularDist(it)
                for ref3d in self._refsList:
                    for prefix in prefixes:
                        volFn = self.protocol._getFileName(prefix + 'volume', iter=it, ref3d=ref3d)
                        args = "--input '%s' --mode projector 256 -a %sclass%06d_angularDist@%s red %f " % (volFn, prefix, ref3d, data_angularDist, radius)
                        if sphere > 0:
                            args += ' %f ' % sphere
                        arguments.append(args)
                        
        elif self.displayAngDist == 0:
            nrefs = len(self._refsList)
            n = nrefs * len(prefixes)
            gridsize = self._getGridSize(n)
            
            for it in self._iterations:
                data_angularDist = self.protocol._getIterAngularDist(it)
                xplotter = XmippPlotter(*gridsize, mainTitle='Iteration %d' % it, windowTitle="Angular Distribution")
                for ref3d in self._refsList:
                    for prefix in prefixes:
                        md = xmipp.MetaData("class%06d_angularDist@%s" % (ref3d, data_angularDist))
                        plot_title = '%sclass %d' % (prefix, ref3d)
                        xplotter.plotMdAngularDistribution(plot_title, md)
                plotters.append(xplotter)
                
        return plotters, arguments
    
    def _showAngularDistribution(self, paramName):
        plotters, arguments = self._createAngularDistribution(paramName)
        
        if arguments:
            for args in arguments:
                os.system('xmipp_chimera_client %s &' % args)
        
        if plotters:
            for xplotter in plotters:
                xplotter.show()
            
            
    def getVisualizeDictWeb(self):
        return {'showImagesInClasses': 'doShowImagesInClasses',
                'showLL': 'doShowLL',
                'showPMax': 'doShowPMax',
                'showChanges': 'doShowChanges',
                'displayVol': 'doShowVolumes',
                'displayAngDist': 'doShowAngularDistribution'
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

        
