# **************************************************************************
# *
# * Authors:     L. del Cano (ldelcano@cnb.csic.es)
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
This module implement the wrappers aroung Xmipp ML3D protocol
visualization program.
"""
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em import *
from pyworkflow.gui.form import FormWindow
from protocol_ml3d import XmippProtML3D
from viewer import runShowJ
from plotter import XmippPlotter
from xmipp import MetaData, AGGR_SUM, MDL_REF3D, MDL_WEIGHT, MDL_WEIGHT, MDValueEQ

import numpy as np

LAST_ITER = 0
ALL_ITER = 1
SELECTED_ITERS = 2

class XmippML3DViewer(ProtocolViewer):
    """ Wrapper to visualize different type of data objects
    with the Xmipp program xmipp_showj
    """
    _targets = [XmippProtML3D]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    _label = 'Xmipp Viewer ML3D'

    
    def _defineParams(self, form):
        form.addSection(label='Preparation')
        form.addParam('doShowGreyScaleRefVol', BooleanParam, label="Visualize the grey-scale corrected reference volume?", default=True)
        form.addParam('doShowFilterRefVol', BooleanParam, label="Visualize the low-pass filtered reference volume?", default=True)
        form.addParam('doShowSeedsVols', BooleanParam, label="Visualize the generated seeds volumes?", default=True)
        
        form.addSection(label='Results per Iteration')
        form.addParam('iterToShow', EnumParam, label="Which iteration do you want to visualize?", default=0, 
                      choices=['last','all','selection'], display=EnumParam.DISPLAY_LIST)
        form.addParam('selectedIters', StringParam, default='',
              label='Selected iterations', condition='iterToShow==2',
              help='Which iteration do you want to visualize.')     
        form.addParam('doShow2DAvgs', BooleanParam, label="Visualize weighted 2D-averages?.", default=True)
        form.addParam('doShow3DRefsVolumes', BooleanParam, label="Visualize the 3D-references volumes?", default=True)
        form.addParam('doShowAngDist', BooleanParam, label="Plot angular distribution of 3D-references?", default=True)
        form.addParam('doShowDataDist', BooleanParam, label="Plot data distribution over 3d-references?", default=True)     
        
        form.addSection(label='Overall Results')
        form.addParam('doShowStatistics', BooleanParam, label="Plot overall convergence statistics?", default=True)
    
    def _getVisualizeDict(self):
        return {'doShowGreyScaleRefVol': self._viewCorrectedVols,
                'doShowFilterRefVol': self._viewFilteredVols,
                'doShowSeedsVols': self._viewGeneratedVols,
                'doShow2DAvgs': self._view2DAvgs,
                'doShow3DRefsVolumes': self._view3DRefsVolumes,
                'doShowAngDist': self._plotAngularDistribution,
                'doShowDataDist': self._plotClassDistribution,
                'doShowStatistics': self._plotStatistics
                }
        
    def _viewAll(self, *args):
        if self.doShowGreyScaleRefVol:
            self._viewCorrectedVols()
        if self.doShowFilterRefVol:
            self._viewFilteredVols()
        if self.doShowSeedsVols:
            self._viewGeneratedVols()
        if self.doShow2DAvgs:
            self._view2DAvgs()
        if self.doShow3DRefsVolumes:
            self._view3DRefsVolumes()
        if self.doShowAngDist:
            self._plotAngularDistribution()
        if self.doShowDataDist:
            self._plotClassDistribution()
        if self.doShowStatistics: 
            self._plotStatistics()
        
    def _viewCorrectedVols(self, e=None):
        runShowJ(self.protocol._getExtraPath("corrected_volumes.stk"))
        
    def _viewFilteredVols(self, e=None):
        runShowJ(self.protocol._getExtraPath("filtered_volumes.stk"))
        
    def _viewGeneratedVols(self, e=None):
        runShowJ(self.protocol._getExtraPath("generated_volumes.stk"))
        
    def _view2DAvgs(self, e=None):
        self._viewIterationFile("ml2dextra/iter%03d/iter_classes.xmd")
        
    def _view3DRefsVolumes(self, e=None):
        self._viewIterationFile("extra/iter%03d/vol000001.vol")
    
    def _plotAngularDistribution(self, e=None):
        self.setVisualizeIterations()
        for iter in self.visualizeIters:
            extraPath = self.protocol._getPath("ml2dextra/iter%03d/iter_classes.xmd" % iter)
            print "extraPath=%s" % extraPath
            if not os.path.exists(extraPath):
                self.formWindow.showError('Iteration %s does not exist.' % iter)        
            else:
                md = MetaData("classes@%s" % extraPath)
                md2 = MetaData()
                md2.aggregate(md, AGGR_SUM, MDL_REF3D, MDL_WEIGHT, MDL_WEIGHT)
                nrefs = md2.size()
                figsize = None
                if nrefs == 1:
                    gridsize = [1, 1]
                    figsize = (4, 4)
                elif nrefs == 2:
                    gridsize = [1, 2]
                    figsize = (8, 4)
                else:
                    gridsize = [(nrefs+1)/2, 2]
                    figsize = (8, 12)
        
                xplotter = XmippPlotter(*gridsize, figsize=figsize, 
                                        windowTitle="Angular distribution - iteration %d" % iter)
                
                for r in range(1, nrefs+1):
                    md2.importObjects(md, MDValueEQ(MDL_REF3D, r))  
                    plot_title = 'ref %d' % r
                    xplotter.plotAngularDistribution(plot_title, md2)
                if xplotter is not None:
                    xplotter.show()
#                    return self._showOrReturn(xplotter)
                    

    def _plotClassDistribution(self, e=None):
        from numpy import arange
        from matplotlib.ticker import FormatStrFormatter
        self.setVisualizeIterations()
        for iter in self.visualizeIters:
            extraPath = self.protocol._getPath("ml2dextra/iter%03d/iter_classes.xmd" % iter)
            print "extraPath=%s" % extraPath
            if not os.path.exists(extraPath):
                self.formWindow.showError('Iteration %s does not exist.' % iter)        
            else:
                xplotter = XmippPlotter(1, 1, figsize=(4,4),
                                        windowTitle="Images distribution - iteration %d" % iter)
                md = MetaData("classes@%s" % extraPath)
                md2 = MetaData()    
                md2.aggregate(md, AGGR_SUM, MDL_REF3D, MDL_WEIGHT, MDL_WEIGHT)
                weights = [md2.getValue(MDL_WEIGHT, objId) for objId in md2]
                nrefs = len(weights)
                refs3d = arange(1, nrefs + 1)
                width = 0.85
                a = xplotter.createSubPlot('3D references weights on last iteration', 'references', 'weight')
                a.set_xticks(refs3d + 0.45)
                a.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))
                a.set_xlim([0.8, nrefs + 1])
                a.bar(refs3d, weights, width, color='b')
                if xplotter is not None:
                    return self._showOrReturn(xplotter)
                    
    def _plotStatistics(self, e=None):
        from viewer_ml2d import createPlots
        xplotter = createPlots(self.protocol, ['doShowLL', 'doShowPmax'])
        if xplotter is not None:
            return self._showOrReturn(xplotter)
                
    def _viewIterationFile(self, filePath):
        self.setVisualizeIterations()
        for iter in self.visualizeIters:
            extraPath = self.protocol._getPath(filePath % iter)
            print "extraPath=%s" % extraPath
            if os.path.exists(extraPath):
                runShowJ(extraPath)        
            else:
                self.formWindow.showError('Iteration %s does not exist.' % iter)        
        
    def setVisualizeIterations(self):
        '''Validate and set the set of iterations to visualize.
        If not set is selected, only use the last iteration'''
        self.lastIter = self.protocol._lastIteration()
        self.visualizeIters = []
        
        if self.iterToShow.get() == LAST_ITER:
            self.visualizeIters = [self.lastIter]
            
        elif self.iterToShow.get() == ALL_ITER:
            self.visualizeIters = range(1, self.lastIter + 1)
        elif self.iterToShow.get() == SELECTED_ITERS:
            if self.selectedIters.empty():
                self.formWindow.showError('Please select the iterations that you want to visualize.')
            else:
                try:
                    self.visualizeIters = self._getListFromRangeString(self.selectedIters.get())
                except Exception, ex:
                    self.formWindow.showError('Invalid iterations range.')
        

    def getVisualizeDictWeb(self):
        return {'doShowGreyScaleRefVol': 'viewCorrectedVols',
                'doShowFilterRefVol': 'viewFilteredVolsProtocol',
                'doShowSeedsVols': 'viewGeneratedVols',
                'doShow2DAvgs': 'view2DAvgs',
                'doShow3DRefsVolumes': 'view3DRefsVolumes',
                'doShowAngDist': 'doPlotAngularDistribution',
                'doShowDataDist': 'doPlotClassDistribution',
                'doShowStatistics': 'doPlotStatistics'
                }

    @classmethod
    def getView(cls):
        """ This function will notify the web viewer for this protocol"""
        return "viewerForm"
    
    @classmethod
    def getViewFunction(cls):
        """ This will return the name of the function to view
        in web one (or all) params of the protocol"""
        return "viewerML3D"
    
